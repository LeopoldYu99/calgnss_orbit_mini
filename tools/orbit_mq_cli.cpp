#include "orbit_mq_protocol.h"

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mqueue.h>
#include <stdexcept>
#include <string>
#include <thread>
#include <time.h>
#include <unistd.h>
#include <utility>
#include <vector>

namespace {

using orbit_prediction::MqFrame;
using orbit_prediction::MqRecvType;
using orbit_prediction::MqSendType;

std::string EnvironmentString(const char *name, const char *fallback)
{
    const char *value = std::getenv(name);
    return value && *value ? value : fallback;
}

long long ParseInteger(const char *text, const char *label, long long minimum,
                       long long maximum)
{
    if (!text || !*text) {
        throw std::invalid_argument(std::string("missing ") + label);
    }
    char *end = nullptr;
    errno = 0;
    const long long value = std::strtoll(text, &end, 10);
    if (errno != 0 || !end || *end != '\0' || value < minimum || value > maximum) {
        throw std::invalid_argument(std::string("invalid ") + label + ": " + text);
    }
    return value;
}

timespec DeadlineAfter(std::chrono::milliseconds delay)
{
    timespec deadline{};
    if (clock_gettime(CLOCK_REALTIME, &deadline) != 0) {
        throw std::runtime_error(std::string("clock_gettime failed: ") +
                                 std::strerror(errno));
    }
    const auto seconds = std::chrono::duration_cast<std::chrono::seconds>(delay);
    const auto nanoseconds =
        std::chrono::duration_cast<std::chrono::nanoseconds>(delay - seconds);
    deadline.tv_sec += seconds.count();
    deadline.tv_nsec += nanoseconds.count();
    if (deadline.tv_nsec >= 1000000000L) {
        ++deadline.tv_sec;
        deadline.tv_nsec -= 1000000000L;
    }
    return deadline;
}

std::string ReadFile(const std::string &path)
{
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        throw std::runtime_error("cannot open input file: " + path);
    }
    return std::string(std::istreambuf_iterator<char>(input),
                       std::istreambuf_iterator<char>());
}

std::uint32_t NewSequence()
{
    static std::atomic<std::uint32_t> next([]() {
        const auto ticks = static_cast<std::uint32_t>(
            std::chrono::steady_clock::now().time_since_epoch().count());
        const auto pid = static_cast<std::uint32_t>(getpid());
        return ticks ^ (pid << 16U);
    }());
    return next.fetch_add(1, std::memory_order_relaxed);
}

class MqClient final {
public:
    MqClient(std::string request_name, std::string response_name)
        : request_name_(std::move(request_name)), response_name_(std::move(response_name))
    {
        request_queue_ = OpenWithRetry(request_name_, O_WRONLY);
        response_queue_ = OpenWithRetry(response_name_, O_RDONLY);
        mq_attr request_attributes{};
        mq_attr response_attributes{};
        if (mq_getattr(request_queue_, &request_attributes) != 0 ||
            mq_getattr(response_queue_, &response_attributes) != 0) {
            throw std::runtime_error(std::string("mq_getattr failed: ") +
                                     std::strerror(errno));
        }
        request_size_ = static_cast<std::size_t>(request_attributes.mq_msgsize);
        response_size_ = static_cast<std::size_t>(response_attributes.mq_msgsize);
    }

    ~MqClient()
    {
        if (request_queue_ != static_cast<mqd_t>(-1)) {
            mq_close(request_queue_);
        }
        if (response_queue_ != static_cast<mqd_t>(-1)) {
            mq_close(response_queue_);
        }
    }

    std::size_t MaximumPayloadSize() const
    {
        return request_size_ > orbit_prediction::kMqFrameHeaderSize
            ? request_size_ - orbit_prediction::kMqFrameHeaderSize
            : 0;
    }

    void Send(const MqFrame &frame)
    {
        std::vector<std::uint8_t> bytes;
        std::string error;
        if (!orbit_prediction::BuildMqFrame(frame, request_size_, &bytes, &error)) {
            throw std::runtime_error(error);
        }
        const timespec deadline = DeadlineAfter(std::chrono::seconds(5));
        if (mq_timedsend(request_queue_,
                         reinterpret_cast<const char *>(bytes.data()),
                         bytes.size(), 0, &deadline) != 0) {
            throw std::runtime_error(std::string("mq_timedsend failed: ") +
                                     std::strerror(errno));
        }
    }

    bool Receive(std::chrono::milliseconds timeout, MqFrame *frame)
    {
        std::vector<std::uint8_t> bytes(response_size_);
        const timespec deadline = DeadlineAfter(timeout);
        const ssize_t length = mq_timedreceive(
            response_queue_, reinterpret_cast<char *>(bytes.data()), bytes.size(),
            nullptr, &deadline);
        if (length < 0) {
            if (errno == EINTR || errno == EAGAIN || errno == ETIMEDOUT) {
                return false;
            }
            throw std::runtime_error(std::string("mq_timedreceive failed: ") +
                                     std::strerror(errno));
        }
        std::string error;
        if (!orbit_prediction::ParseMqFrame(
                bytes.data(), static_cast<std::size_t>(length), frame, &error)) {
            throw std::runtime_error("invalid received frame: " + error);
        }
        return true;
    }

private:
    static mqd_t OpenWithRetry(const std::string &name, int flags)
    {
        for (int attempt = 0; attempt < 50; ++attempt) {
            const mqd_t queue = mq_open(name.c_str(), flags);
            if (queue != static_cast<mqd_t>(-1)) {
                return queue;
            }
            if (errno != ENOENT && errno != EINTR) {
                break;
            }
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
        throw std::runtime_error("cannot open MQ " + name + ": " +
                                 std::strerror(errno));
    }

    std::string request_name_;
    std::string response_name_;
    mqd_t request_queue_{static_cast<mqd_t>(-1)};
    mqd_t response_queue_{static_cast<mqd_t>(-1)};
    std::size_t request_size_{};
    std::size_t response_size_{};
};

MqFrame NewRequest(MqSendType type)
{
    MqFrame frame;
    frame.message_type = static_cast<std::uint8_t>(type);
    frame.message_seq = NewSequence();
    frame.stream_end = orbit_prediction::kMqStreamEnd;
    return frame;
}

std::chrono::milliseconds TimeoutSeconds(const char *text, long long fallback)
{
    const long long seconds = text
        ? ParseInteger(text, "timeout_s", 1, 86400)
        : fallback;
    return std::chrono::milliseconds(seconds * 1000);
}

int SendOnly(MqClient &client, const MqFrame &frame)
{
    client.Send(frame);
    std::cout << "RESULT queued=1 type="
              << static_cast<unsigned>(frame.message_type)
              << " sequence=" << frame.message_seq << '\n';
    return 0;
}

int RunStatus(MqClient &client, int argc, char **argv)
{
    const MqFrame request = NewRequest(MqSendType::GetStatus);
    client.Send(request);
    MqFrame response;
    if (!client.Receive(TimeoutSeconds(argc > 2 ? argv[2] : nullptr, 5), &response)) {
        std::cout << "RESULT reply_received=0 success=unknown message=\"response timeout\"\n";
        return 0;
    }
    if (response.message_type != static_cast<std::uint8_t>(MqRecvType::Heartbeat) ||
        response.message_seq != request.message_seq) {
        throw std::runtime_error("received frame is not the requested status response");
    }
    orbit_prediction::StatusReply status;
    std::string error;
    if (!orbit_prediction::DecodeStatusReply(response.payload, &status, &error)) {
        throw std::runtime_error(error);
    }
    std::cout << "STATUS value=" << static_cast<unsigned>(status.status)
              << " timestamp_ms=" << status.timestamp_ms
              << " message=" << std::quoted(status.message) << '\n'
              << "RESULT reply_received=1 success=1\n";
    return 0;
}

int RunRtcm(MqClient &client, int argc, char **argv)
{
    if (argc < 3) {
        throw std::invalid_argument("rtcm requires FILE");
    }
    const std::string bytes = ReadFile(argv[2]);
    if (bytes.empty()) {
        throw std::invalid_argument("RTCM file is empty");
    }
    const std::size_t chunk_size = client.MaximumPayloadSize();
    if (chunk_size == 0) {
        throw std::runtime_error("request MQ is smaller than the CCSM frame header");
    }
    std::size_t chunks = 0;
    for (std::size_t offset = 0; offset < bytes.size(); offset += chunk_size) {
        const std::size_t length = std::min(chunk_size, bytes.size() - offset);
        MqFrame frame;
        frame.message_type = static_cast<std::uint8_t>(MqSendType::ReceiveRtcmData);
        frame.message_seq = static_cast<std::uint32_t>(chunks);
        frame.stream_end = offset + length == bytes.size()
            ? orbit_prediction::kMqStreamEnd
            : orbit_prediction::kMqStreamMore;
        frame.payload.assign(bytes.begin() + static_cast<std::ptrdiff_t>(offset),
                             bytes.begin() + static_cast<std::ptrdiff_t>(offset + length));
        client.Send(frame);
        ++chunks;
    }
    std::cout << "RESULT chunks=" << chunks << " queued=" << chunks
              << " success=1\n";
    return 0;
}

int RunPredict(MqClient &client, int argc, char **argv)
{
    if (argc < 5) {
        throw std::invalid_argument(
            "predict requires START_TIME_S DURATION_S OUTPUT.csv");
    }
    orbit_prediction::OrbitPredictionRequest prediction;
    prediction.start_time_s = ParseInteger(
        argv[2], "start_time_s", 0, 253402300799LL);
    prediction.duration_s = static_cast<std::uint32_t>(
        ParseInteger(argv[3], "duration_s", 1, 3600));
    const auto timeout = TimeoutSeconds(argc > 5 ? argv[5] : nullptr, 600);

    std::ofstream csv(argv[4], std::ios::binary | std::ios::trunc);
    if (!csv) {
        throw std::runtime_error(std::string("cannot create CSV: ") + argv[4]);
    }
    csv << "timestamp_ms,x,y,z,vx,vy,vz\n" << std::setprecision(17);

    MqFrame request = NewRequest(MqSendType::PredictOrbit);
    if (!orbit_prediction::EncodePredictionRequest(prediction, &request.payload)) {
        throw std::runtime_error("cannot encode prediction request");
    }
    const auto started = std::chrono::steady_clock::now();
    client.Send(request);

    std::size_t points = 0;
    std::size_t frames = 0;
    std::uint32_t expected_sequence = 0;
    bool terminal_received = false;
    for (;;) {
        MqFrame response;
        if (!client.Receive(timeout, &response)) {
            break;
        }
        if (response.message_type !=
            static_cast<std::uint8_t>(MqRecvType::ReceiveTelemetryData)) {
            std::cerr << "warning: ignored non-prediction end-to-main frame\n";
            continue;
        }
        if (response.message_seq != expected_sequence) {
            throw std::runtime_error("prediction frame sequence gap");
        }
        ++expected_sequence;
        ++frames;
        std::string error;
        if (!response.payload.empty()) {
            orbit_prediction::OrbitData data;
            if (!orbit_prediction::DecodeOrbitData(response.payload, &data, &error)) {
                throw std::runtime_error(error);
            }
            for (const auto &point : data.points) {
                csv << point.timestamp_ms << ',' << point.x << ',' << point.y << ','
                    << point.z << ',' << point.vx << ',' << point.vy << ','
                    << point.vz << '\n';
                ++points;
            }
        }
        if (orbit_prediction::IsMqStreamEnd(response.stream_end)) {
            terminal_received = true;
            break;
        }
    }
    csv.close();
    const auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                                std::chrono::steady_clock::now() - started)
                                .count();
    std::cout << "RESULT reply_received=" << (frames > 0 ? 1 : 0)
              << " terminal_received=" << (terminal_received ? 1 : 0)
              << " success=" << (terminal_received ? "1" : "unknown")
              << " points=" << points << " frames=" << frames
              << " elapsed_ms=" << elapsed_ms << " output=" << std::quoted(argv[4])
              << '\n';
    return 0;
}

void PrintUsage(const char *program)
{
    std::cerr
        << "Usage:\n"
        << "  " << program << " status [timeout_s]\n"
        << "  " << program << " time-sync TIMESTAMP_MS\n"
        << "  " << program << " uplink FILE [port]\n"
        << "  " << program << " rtcm FILE\n"
        << "  " << program << " predict START_TIME_S DURATION_S OUTPUT.csv [timeout_s]\n"
        << "  " << program << " stop\n"
        << "  " << program << " reset\n";
}

}  // namespace

int main(int argc, char **argv)
{
    try {
        if (argc < 2) {
            PrintUsage(argv[0]);
            return 1;
        }
        MqClient client(
            EnvironmentString("ORBIT_MQ_REQUEST_QUEUE", "/csm_main_to_end0"),
            EnvironmentString("ORBIT_MQ_RESPONSE_QUEUE", "/csm_end0_to_main"));
        const std::string command = argv[1];
        if (command == "status") {
            return RunStatus(client, argc, argv);
        }
        if (command == "time-sync") {
            if (argc < 3) {
                throw std::invalid_argument("time-sync requires TIMESTAMP_MS");
            }
            const auto timestamp_ms = ParseInteger(
                argv[2], "timestamp_ms", 0, std::numeric_limits<long long>::max());
            orbit_prediction::TimeSyncData time_sync;
            time_sync.timestamp_s = timestamp_ms / 1000;
            MqFrame request = NewRequest(MqSendType::ReceiveTimeSync);
            orbit_prediction::EncodeTimeSyncData(time_sync, &request.payload);
            return SendOnly(client, request);
        }
        if (command == "uplink") {
            if (argc < 3) {
                throw std::invalid_argument("uplink requires FILE");
            }
            orbit_prediction::UplinkPacket packet;
            packet.port = static_cast<std::uint32_t>(
                ParseInteger(argc > 3 ? argv[3] : "1", "port", 0, 65535));
            packet.data = ReadFile(argv[2]);
            MqFrame request = NewRequest(MqSendType::ReceiveUplinkData);
            std::string error;
            if (!orbit_prediction::EncodeUplinkPacket(packet, &request.payload, &error)) {
                throw std::runtime_error(error);
            }
            return SendOnly(client, request);
        }
        if (command == "rtcm") {
            return RunRtcm(client, argc, argv);
        }
        if (command == "predict") {
            return RunPredict(client, argc, argv);
        }
        if (command == "stop" || command == "reset") {
            return SendOnly(client, NewRequest(
                command == "stop" ? MqSendType::Stop : MqSendType::Reset));
        }
        PrintUsage(argv[0]);
        throw std::invalid_argument("unknown command: " + command);
    } catch (const std::exception &error) {
        std::cerr << "orbit_mq_cli: " << error.what() << '\n';
        return 1;
    }
}

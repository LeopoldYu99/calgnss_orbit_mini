#include "orbit_mq_protocol.h"

#include <cerrno>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <iomanip>
#include <iostream>
#include <mqueue.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <time.h>
#include <utility>
#include <vector>

namespace {

void Require(bool condition, const std::string &message)
{
    if (!condition) throw std::runtime_error(message);
}

std::string EnvironmentString(const char *name, const char *fallback)
{
    const char *value = std::getenv(name);
    return value && *value ? value : fallback;
}

timespec DeadlineAfter(std::chrono::seconds delay)
{
    timespec deadline{};
    Require(clock_gettime(CLOCK_REALTIME, &deadline) == 0, "clock_gettime failed");
    deadline.tv_sec += delay.count();
    return deadline;
}

mqd_t OpenQueueWithRetry(const std::string &name, int flags)
{
    for (int attempt = 0; attempt < 50; ++attempt) {
        const mqd_t queue = mq_open(name.c_str(), flags);
        if (queue != static_cast<mqd_t>(-1)) return queue;
        if (errno != ENOENT) break;
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
    throw std::runtime_error("cannot open MQ " + name + ": " + std::strerror(errno));
}

std::string TestNmeaArc()
{
    std::ostringstream data;
    data << "$GNRMC,040403.00,A,1724.2594222,N,15304.2955483,E,0.0,0.0,230626,,,A\n";
    for (int index = 0; index < 12; ++index) {
        const int second = 3 + index;
        const double longitude = 15300.0 + 4.2955483 + index * 4.2;
        data << "$GNGGA,0404" << std::setw(2) << std::setfill('0') << second
             << ".00,1724.2594222,N," << std::fixed << std::setprecision(7)
             << longitude << ",E,1,14,0.69,271441.347,M,0.989,M,,\n";
    }
    return data.str();
}

std::string TestRtcmFrame()
{
    const std::vector<std::uint8_t> frame{
        0xD3, 0x00, 0x04, 0x43, 0x50, 0x00, 0x00, 0x44, 0xFE, 0x2E};
    return std::string(reinterpret_cast<const char *>(frame.data()), frame.size());
}

class Client final {
public:
    Client(std::string request_name, std::string response_name)
        : request_name_(std::move(request_name)), response_name_(std::move(response_name))
    {
        request_queue_ = OpenQueueWithRetry(request_name_, O_WRONLY);
        response_queue_ = OpenQueueWithRetry(response_name_, O_RDONLY);
        mq_attr attributes{};
        Require(mq_getattr(request_queue_, &attributes) == 0, "request mq_getattr failed");
        request_size_ = static_cast<std::size_t>(attributes.mq_msgsize);
        Require(mq_getattr(response_queue_, &attributes) == 0, "response mq_getattr failed");
        response_size_ = static_cast<std::size_t>(attributes.mq_msgsize);
    }

    ~Client()
    {
        if (request_queue_ != static_cast<mqd_t>(-1)) mq_close(request_queue_);
        if (response_queue_ != static_cast<mqd_t>(-1)) mq_close(response_queue_);
        if (EnvironmentString("ORBIT_MQ_TEST_CLEANUP", "0") == "1") {
            mq_unlink(request_name_.c_str());
            mq_unlink(response_name_.c_str());
        }
    }

    void Send(const orbit_prediction::MqFrame &frame)
    {
        std::vector<std::uint8_t> bytes;
        std::string error;
        Require(orbit_prediction::BuildMqFrame(frame, request_size_, &bytes, &error), error);
        const timespec deadline = DeadlineAfter(std::chrono::seconds(5));
        Require(mq_timedsend(request_queue_,
                             reinterpret_cast<const char *>(bytes.data()),
                             bytes.size(), 0, &deadline) == 0,
                std::string("mq_timedsend failed: ") + std::strerror(errno));
    }

    orbit_prediction::MqFrame Receive()
    {
        std::vector<std::uint8_t> bytes(response_size_);
        const timespec deadline = DeadlineAfter(std::chrono::seconds(15));
        const ssize_t length = mq_timedreceive(
            response_queue_, reinterpret_cast<char *>(bytes.data()), bytes.size(),
            nullptr, &deadline);
        Require(length >= 0, std::string("mq_timedreceive failed: ") +
                                 std::strerror(errno));
        orbit_prediction::MqFrame frame;
        std::string error;
        Require(orbit_prediction::ParseMqFrame(
                    bytes.data(), static_cast<std::size_t>(length), &frame, &error), error);
        return frame;
    }

private:
    std::string request_name_;
    std::string response_name_;
    mqd_t request_queue_{static_cast<mqd_t>(-1)};
    mqd_t response_queue_{static_cast<mqd_t>(-1)};
    std::size_t request_size_{};
    std::size_t response_size_{};
};

orbit_prediction::MqFrame NewFrame(orbit_prediction::MqSendType type,
                                   std::uint32_t sequence)
{
    orbit_prediction::MqFrame frame;
    frame.message_type = static_cast<std::uint8_t>(type);
    frame.message_seq = sequence;
    frame.stream_end = orbit_prediction::kMqStreamEnd;
    return frame;
}

}  // namespace

int main()
{
    try {
        Client client(
            EnvironmentString("ORBIT_MQ_REQUEST_QUEUE", "/orbit_prediction_test_requests"),
            EnvironmentString("ORBIT_MQ_RESPONSE_QUEUE", "/orbit_prediction_test_responses"));

        orbit_prediction::MqFrame status =
            NewFrame(orbit_prediction::MqSendType::GetStatus, 100);
        client.Send(status);
        auto status_response = client.Receive();
        Require(status_response.message_type == static_cast<std::uint8_t>(
                    orbit_prediction::MqRecvType::Heartbeat) &&
                    status_response.message_seq == 100,
                "status response header mismatch");
        orbit_prediction::StatusReply decoded_status;
        std::string error;
        Require(orbit_prediction::DecodeStatusReply(
                    status_response.payload, &decoded_status, &error), error);
        Require(decoded_status.status == orbit_prediction::ServiceStatus::Idle,
                "service must initially be idle");

        auto time_frame = NewFrame(orbit_prediction::MqSendType::ReceiveTimeSync, 1);
        // ReceiveTimeSync is deliberately ignored, including malformed/empty
        // payloads, and must not affect RTCM or prediction time handling.
        client.Send(time_frame);

        orbit_prediction::UplinkPacket uplink;
        uplink.port = 1;
        uplink.data = TestNmeaArc();
        auto uplink_frame = NewFrame(orbit_prediction::MqSendType::ReceiveUplinkData, 2);
        Require(orbit_prediction::EncodeUplinkPacket(uplink, &uplink_frame.payload, &error),
                error);
        client.Send(uplink_frame);

        auto rtcm_frame = NewFrame(orbit_prediction::MqSendType::ReceiveRtcmData, 3);
        const std::string rtcm = TestRtcmFrame();
        rtcm_frame.payload.assign(rtcm.begin(), rtcm.end());
        client.Send(rtcm_frame);

        orbit_prediction::OrbitPredictionRequest prediction{1782187454LL, 2};
        auto predict_frame = NewFrame(orbit_prediction::MqSendType::PredictOrbit, 4);
        Require(orbit_prediction::EncodePredictionRequest(
                    prediction, &predict_frame.payload),
                "cannot encode prediction request");
        client.Send(predict_frame);

        std::size_t point_count = 0;
        std::uint32_t expected_sequence = 0;
        for (;;) {
            const auto response = client.Receive();
            Require(response.message_type == static_cast<std::uint8_t>(
                        orbit_prediction::MqRecvType::ReceiveTelemetryData),
                    "prediction telemetry response type mismatch");
            Require(response.message_seq == expected_sequence++,
                    "prediction sequence mismatch");
            if (!response.payload.empty()) {
                orbit_prediction::OrbitData data;
                Require(orbit_prediction::DecodeOrbitData(response.payload, &data, &error),
                        error);
                point_count += data.points.size();
            }
            if (response.stream_end == orbit_prediction::kMqStreamEnd) break;
        }
        Require(point_count == 0,
                "prediction without an RTCM position must only return a terminal packet");

        client.Send(NewFrame(orbit_prediction::MqSendType::Stop, 5));
        client.Send(NewFrame(orbit_prediction::MqSendType::Reset, 6));

        std::cout << "CCSM fixed-binary MQ smoke test passed: 7 inbound methods, "
                  << point_count << " orbit points (RTCM start-time override active)\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "CCSM fixed-binary MQ smoke test failed: " << error.what() << '\n';
        return 1;
    }
}

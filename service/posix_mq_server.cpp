#include "posix_mq_server.h"

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cstring>
#include <fcntl.h>
#include <iostream>
#include <stdexcept>
#include <sys/stat.h>
#include <time.h>
#include <utility>

namespace orbit_prediction {
namespace {

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

void ValidateQueueName(const std::string &name, const char *label)
{
    if (name.size() < 2 || name.front() != '/' ||
        name.find('/', 1) != std::string::npos) {
        throw std::invalid_argument(std::string(label) +
                                    " must begin with '/' and contain no other '/'");
    }
}

OrbitPoint CopyPoint(const PredictionPoint &source)
{
    return {source.timestamp_ms, source.x, source.y, source.z,
            source.vx, source.vy, source.vz};
}

}  // namespace

PosixMqServer::PosixMqServer(OrbitPredictionHandler &handler, PosixMqOptions options)
    : handler_(handler), options_(std::move(options))
{
    ValidateQueueName(options_.request_queue, "request queue name");
    ValidateQueueName(options_.response_queue, "response queue name");
    if (options_.request_queue == options_.response_queue) {
        throw std::invalid_argument("request and response queues must be different");
    }
    if (options_.max_messages < 1 ||
        options_.message_size < static_cast<long>(kMqFrameHeaderSize)) {
        throw std::invalid_argument("invalid POSIX MQ size options");
    }

    mq_attr attributes{};
    attributes.mq_maxmsg = options_.max_messages;
    attributes.mq_msgsize = options_.message_size;

    const int receive_flags = O_RDONLY | (options_.create_queues ? O_CREAT : 0);
    request_queue_ = options_.create_queues
        ? mq_open(options_.request_queue.c_str(), receive_flags, 0666, &attributes)
        : mq_open(options_.request_queue.c_str(), receive_flags);
    if (request_queue_ == static_cast<mqd_t>(-1)) {
        throw std::runtime_error("cannot open request queue " + options_.request_queue +
                                 ": " + std::strerror(errno));
    }

    const int send_flags = O_WRONLY | (options_.create_queues ? O_CREAT : 0);
    response_queue_ = options_.create_queues
        ? mq_open(options_.response_queue.c_str(), send_flags, 0666, &attributes)
        : mq_open(options_.response_queue.c_str(), send_flags);
    if (response_queue_ == static_cast<mqd_t>(-1)) {
        const std::string error = std::strerror(errno);
        mq_close(request_queue_);
        request_queue_ = static_cast<mqd_t>(-1);
        throw std::runtime_error("cannot open response queue " + options_.response_queue +
                                 ": " + error);
    }

    mq_attr request_attributes{};
    mq_attr response_attributes{};
    if (mq_getattr(request_queue_, &request_attributes) != 0 ||
        mq_getattr(response_queue_, &response_attributes) != 0) {
        const std::string error = std::strerror(errno);
        mq_close(request_queue_);
        mq_close(response_queue_);
        request_queue_ = static_cast<mqd_t>(-1);
        response_queue_ = static_cast<mqd_t>(-1);
        throw std::runtime_error("mq_getattr failed: " + error);
    }
    request_message_size_ = static_cast<std::size_t>(request_attributes.mq_msgsize);
    response_message_size_ = static_cast<std::size_t>(response_attributes.mq_msgsize);
}

PosixMqServer::~PosixMqServer()
{
    StopAndJoinPrediction();
    if (request_queue_ != static_cast<mqd_t>(-1)) {
        mq_close(request_queue_);
    }
    if (response_queue_ != static_cast<mqd_t>(-1)) {
        mq_close(response_queue_);
    }
}

void PosixMqServer::Run(const std::atomic<bool> &shutdown_requested)
{
    std::vector<std::uint8_t> buffer(request_message_size_);
    while (!shutdown_requested.load(std::memory_order_acquire)) {
        const timespec deadline = DeadlineAfter(std::chrono::milliseconds(250));
        const ssize_t received = mq_timedreceive(
            request_queue_, reinterpret_cast<char *>(buffer.data()), buffer.size(),
            nullptr, &deadline);
        if (received < 0) {
            if (errno == EINTR || errno == EAGAIN || errno == ETIMEDOUT) {
                continue;
            }
            throw std::runtime_error(std::string("mq_timedreceive failed: ") +
                                     std::strerror(errno));
        }

        MqFrame request;
        std::string error;
        if (!ParseMqFrame(buffer.data(), static_cast<std::size_t>(received),
                          &request, &error)) {
            std::cerr << "invalid CCSM MQ frame ignored: " << error << '\n';
            continue;
        }
        if (options_.console_logging) {
            const auto type = static_cast<MqSendType>(request.message_type);
            std::cout << "Received MqSendType: type="
                      << static_cast<unsigned>(request.message_type)
                      << '(' << MqSendTypeName(type) << ')'
                      << " seq=" << request.message_seq
                      << " stream_end=" << static_cast<unsigned>(request.stream_end)
                      << " payload_bytes=" << request.payload.size() << '\n';
        }
        Dispatch(request);
    }
    StopAndJoinPrediction();
}

void PosixMqServer::Dispatch(const MqFrame &request)
{
    std::string error;
    CommonReply result;
    switch (static_cast<MqSendType>(request.message_type)) {
    case MqSendType::ReceiveUplinkData: {
        UplinkPacket packet;
        if (!DecodeUplinkPacket(request.payload, &packet, &error)) {
            LogRejected(request, error);
            return;
        }
        result = handler_.ReceiveUplinkData(packet);
        break;
    }
    case MqSendType::ReceiveTimeSync:
        // Intentionally ignored. RTCM processing and prediction must not use
        // MQ time synchronization or the host system clock.
        return;
    case MqSendType::ReceiveRtcmData: {
        RtcmData rtcm;
        if (!request.payload.empty()) {
            rtcm.data.assign(reinterpret_cast<const char *>(request.payload.data()),
                             request.payload.size());
        }
        result = handler_.ReceiveRTCMData(rtcm);
        break;
    }
    case MqSendType::GetStatus:
        if (!request.payload.empty()) {
            LogRejected(request, "GetStatus payload must be empty");
            return;
        }
        SendStatus(request.message_seq);
        return;
    case MqSendType::Stop:
        if (!request.payload.empty()) {
            LogRejected(request, "Stop payload must be empty");
            return;
        }
        result = handler_.Stop();
        break;
    case MqSendType::Reset:
        if (!request.payload.empty()) {
            LogRejected(request, "Reset payload must be empty");
            return;
        }
        handler_.Stop();
        {
            std::lock_guard<std::mutex> lock(prediction_thread_mutex_);
            if (prediction_thread_.joinable()) {
                prediction_thread_.join();
            }
        }
        result = handler_.Reset();
        break;
    case MqSendType::PredictOrbit: {
        OrbitPredictionRequest prediction_request;
        if (!DecodePredictionRequest(request.payload, &prediction_request, &error)) {
            LogRejected(request, error);
            return;
        }
        StartPrediction(request, prediction_request);
        return;
    }
    case MqSendType::Unspecified:
    default:
        LogRejected(request, "unsupported main-to-orbit message type");
        return;
    }

    if (!result.success) {
        LogRejected(request, result.message);
    }
}

void PosixMqServer::StartPrediction(
    const MqFrame &request, const OrbitPredictionRequest &prediction_request)
{
    if (!handler_.BeginPrediction()) {
        LogRejected(request, "another orbit prediction is already running");
        return;
    }

    std::lock_guard<std::mutex> lock(prediction_thread_mutex_);
    if (prediction_thread_.joinable()) {
        prediction_thread_.join();
    }
    prediction_thread_ = std::thread([this, prediction_request]() {
        std::uint32_t sequence = 0;
        CommonReply terminal;
        try {
            terminal = handler_.RunPrediction(
                prediction_request,
                [this, &sequence](const std::vector<PredictionPoint> &batch) {
                    SendPredictionBatch(&sequence, batch);
                });
        } catch (const std::exception &error) {
            terminal.success = false;
            terminal.message = std::string("prediction worker failed: ") + error.what();
        } catch (...) {
            terminal.success = false;
            terminal.message = "prediction worker failed with an unknown exception";
        }
        if (!terminal.success) {
            std::cerr << terminal.message << '\n';
        }
        MqFrame final_frame;
        final_frame.message_type =
            static_cast<std::uint8_t>(MqRecvType::ReceiveTelemetryData);
        final_frame.message_seq = sequence;
        final_frame.stream_end = kMqStreamEnd;
        SendFrameBestEffort(final_frame);
    });
}

void PosixMqServer::StopAndJoinPrediction()
{
    handler_.Stop();
    std::lock_guard<std::mutex> lock(prediction_thread_mutex_);
    if (prediction_thread_.joinable()) {
        prediction_thread_.join();
    }
}

bool PosixMqServer::SendFrameBestEffort(const MqFrame &frame)
{
    std::vector<std::uint8_t> bytes;
    std::string error;
    if (!BuildMqFrame(frame, response_message_size_, &bytes, &error)) {
        std::cerr << "cannot build end-to-main MQ frame: " << error << '\n';
        return false;
    }
    std::lock_guard<std::mutex> lock(send_mutex_);
    const timespec deadline = DeadlineAfter(std::chrono::milliseconds(100));
    if (mq_timedsend(response_queue_, reinterpret_cast<const char *>(bytes.data()),
                     bytes.size(), 0, &deadline) != 0) {
        std::cerr << "best-effort end-to-main MQ frame dropped: "
                  << std::strerror(errno) << '\n';
        return false;
    }
    if (options_.console_logging) {
        const auto type = static_cast<MqRecvType>(frame.message_type);
        std::cout << "Sent MqRecvType: type="
                  << static_cast<unsigned>(frame.message_type)
                  << '(' << MqRecvTypeName(type) << ')'
                  << " seq=" << frame.message_seq
                  << " stream_end=" << static_cast<unsigned>(frame.stream_end)
                  << " payload_bytes=" << frame.payload.size() << '\n';
    }
    return true;
}

void PosixMqServer::SendStatus(std::uint32_t sequence)
{
    MqFrame frame;
    frame.message_type = static_cast<std::uint8_t>(MqRecvType::Heartbeat);
    frame.message_seq = sequence;
    frame.stream_end = kMqStreamEnd;
    std::string error;
    if (!EncodeStatusReply(handler_.GetStatus(), &frame.payload, &error)) {
        std::cerr << "cannot encode status response: " << error << '\n';
        return;
    }
    SendFrameBestEffort(frame);
}

void PosixMqServer::SendPredictionBatch(
    std::uint32_t *sequence, const std::vector<PredictionPoint> &batch)
{
    if (!sequence || response_message_size_ < kMqFrameHeaderSize + 56U) {
        return;
    }
    const std::size_t maximum_points =
        (response_message_size_ - kMqFrameHeaderSize) / 56U;
    for (std::size_t offset = 0; offset < batch.size(); offset += maximum_points) {
        const std::size_t count = std::min(maximum_points, batch.size() - offset);
        OrbitData data;
        data.points.reserve(count);
        for (std::size_t index = 0; index < count; ++index) {
            data.points.push_back(CopyPoint(batch[offset + index]));
        }
        MqFrame frame;
        frame.message_type =
            static_cast<std::uint8_t>(MqRecvType::ReceiveTelemetryData);
        frame.message_seq = (*sequence)++;
        frame.stream_end = kMqStreamMore;
        std::string error;
        if (!EncodeOrbitData(data, &frame.payload, &error)) {
            std::cerr << "cannot encode raw OrbitPoint payload: " << error << '\n';
            continue;
        }
        SendFrameBestEffort(frame);
    }
}

void PosixMqServer::LogRejected(const MqFrame &request,
                                const std::string &message) const
{
    std::cerr << "main-to-orbit MQ message type="
              << static_cast<unsigned>(request.message_type)
              << " seq=" << request.message_seq << " rejected: " << message << '\n';
}

}  // namespace orbit_prediction

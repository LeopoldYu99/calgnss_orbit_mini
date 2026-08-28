#include "orbit_prediction_service.h"

#include <chrono>
#include <iostream>
#include <sstream>
#include <utility>

namespace orbit_prediction {
namespace {

constexpr double kPredictionStepSeconds = 1.0;

std::int64_t CurrentTimestampMs()
{
    return std::chrono::duration_cast<std::chrono::milliseconds>(
               std::chrono::system_clock::now().time_since_epoch())
        .count();
}

}  // namespace

OrbitPredictionHandler::OrbitPredictionHandler(EngineOptions options)
    : engine_(std::move(options))
{
}

CommonReply OrbitPredictionHandler::Reply(bool success, const std::string &message)
{
    CommonReply reply;
    reply.timestamp_ms = CurrentTimestampMs();
    reply.success = success;
    reply.message = message;
    return reply;
}

CommonReply OrbitPredictionHandler::ReceiveUplinkData(const UplinkPacket &request)
{
    std::string error;
    const bool success = engine_.ReceiveUplinkData(
        request.port, request.packet_header, request.data, &error);
    if (!success) {
        std::cerr << "ReceiveUplinkData rejected: " << error << '\n';
    }
    return Reply(success, success ? "uplink NMEA data accepted" : error);
}

CommonReply OrbitPredictionHandler::ReceiveRTCMData(const RtcmData &request)
{
    std::string message;
    const bool success = engine_.ReceiveRTCMData(request.data, &message);
    if (!success) {
        std::cerr << "ReceiveRTCMData rejected: " << message << '\n';
    }
    return Reply(success, message);
}

StatusReply OrbitPredictionHandler::GetStatus() const
{
    StatusReply reply;
    reply.timestamp_ms = CurrentTimestampMs();
    std::ostringstream message;
    if (prediction_running_.load(std::memory_order_acquire)) {
        reply.status = ServiceStatus::Running;
        message << "orbit prediction is running";
    } else {
        reply.status = ServiceStatus::Idle;
        message << (engine_.IsStopped() ? "orbit prediction is stopped" :
                                             "orbit prediction service is idle");
    }
    message << "; observations=" << engine_.ObservationCount()
            << "; rtcm_frames=" << engine_.RTCMFrameCount()
            << "; rtcm_positions=" << engine_.RTCMPositionCount();
    reply.message = message.str();
    return reply;
}

CommonReply OrbitPredictionHandler::Stop()
{
    engine_.Stop();
    return Reply(true, "orbit prediction stopped");
}

CommonReply OrbitPredictionHandler::Reset()
{
    std::string error;
    const bool success = engine_.Reset(&error);
    if (!success) {
        std::cerr << "Reset failed: " << error << '\n';
    }
    return Reply(success, success ? "orbit prediction service reset" : error);
}

bool OrbitPredictionHandler::BeginPrediction()
{
    bool expected = false;
    return prediction_running_.compare_exchange_strong(
        expected, true, std::memory_order_acq_rel, std::memory_order_acquire);
}

CommonReply OrbitPredictionHandler::RunPrediction(
    const OrbitPredictionRequest &request, const BatchCallback &on_batch)
{
    if (!prediction_running_.load(std::memory_order_acquire)) {
        return Reply(false, "prediction was not started");
    }
    struct PredictionCompletion {
        std::atomic<bool> &running;
        ~PredictionCompletion()
        {
            running.store(false, std::memory_order_release);
        }
    } completion{prediction_running_};

    std::int64_t latest_rtcm_timestamp_ms = 0;
    if (!engine_.LatestRTCMObservationTimestamp(&latest_rtcm_timestamp_ms)) {
        return Reply(false,
                     "temporary RTCM test override requires at least one RTCM position");
    }

    // TEMPORARY TEST OVERRIDE: ignore MQ_METHOD_PREDICT_ORBIT.start_time_s and
    // start at the timestamp of the most recently solved RTCM ECEF position.
    // Remove this override after the RTCM/MQ field test is complete.
    const std::int64_t rtcm_start_time_s = latest_rtcm_timestamp_ms / 1000;
    std::cout << "TEMPORARY RTCM prediction-time override: requested_start_time_s="
              << request.start_time_s << " rtcm_start_time_s=" << rtcm_start_time_s
              << '\n';

    std::string error;
    const bool success = engine_.PredictOrbit(
        rtcm_start_time_s, request.duration_s, kPredictionStepSeconds,
        [&on_batch](const std::vector<PredictionPoint> &batch) {
            if (on_batch) {
                on_batch(batch);
            }
            // MQ responses are best-effort. Missing responses do not make the
            // underlying prediction calculation fail.
            return true;
        },
        [this]() { return engine_.IsStopped(); }, &error);

    if (!success) {
        std::cerr << "PredictOrbit failed: " << error << '\n';
        return Reply(false, error);
    }
    return Reply(true, engine_.IsStopped() ? "orbit prediction stopped" :
                                               "orbit prediction completed");
}

bool OrbitPredictionHandler::IsPredictionRunning() const
{
    return prediction_running_.load(std::memory_order_acquire);
}

}  // namespace orbit_prediction

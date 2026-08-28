#ifndef ORBIT_PREDICTION_ENGINE_H
#define ORBIT_PREDICTION_ENGINE_H

#include "calgnss_orbit_mini.h"

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <functional>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

namespace orbit_prediction {

#ifdef ORBIT_ENABLE_RTCM_POSITIONING
class RtcmPositionSolver;
#endif

struct PredictionPoint {
    std::int64_t timestamp_ms{};
    double x{};
    double y{};
    double z{};
    double vx{};
    double vy{};
    double vz{};
};

struct EngineOptions {
    std::size_t observation_capacity{600};
    int fit_degree{10};
    std::size_t stream_batch_size{100};
    std::size_t rtcm_cache_bytes{4 * 1024 * 1024};
};

class OrbitPredictionEngine {
public:
    explicit OrbitPredictionEngine(EngineOptions options = {});
    ~OrbitPredictionEngine();

    OrbitPredictionEngine(const OrbitPredictionEngine &) = delete;
    OrbitPredictionEngine &operator=(const OrbitPredictionEngine &) = delete;

    bool ReceiveUplinkData(
        std::uint32_t port,
        const std::string &packet_header,
        const std::string &data,
        std::string *error);

    bool ReceiveTimeSync(std::int64_t timestamp_ms, std::string *error);
    bool ReceiveRTCMData(const std::string &data, std::string *message);
    void Stop();
    bool Reset(std::string *error);

    bool PredictOrbit(
        std::int64_t start_time_s,
        std::uint32_t duration_s,
        double step_s,
        const std::function<bool(const std::vector<PredictionPoint> &)> &on_batch,
        const std::function<bool()> &is_cancelled,
        std::string *error);

    std::size_t ObservationCount() const;
    std::uint64_t RTCMFrameCount() const;
    std::uint64_t RTCMPositionCount() const;
    bool LatestRTCMObservationTimestamp(std::int64_t *timestamp_ms) const;
    bool IsStopped() const;

private:
    bool RecreateContext(std::string *error);

    EngineOptions options_;
    std::vector<cg_observation_t> observation_buffer_;
    cg_context_t *context_{nullptr};
    mutable std::mutex mutex_;
    std::atomic<bool> stopped_{false};
    bool has_time_sync_{false};
    std::int64_t time_sync_ms_{0};
    bool has_nmea_date_{false};
    int nmea_year_{0};
    int nmea_month_{0};
    int nmea_day_{0};
    bool has_latest_observation_{false};
    std::int64_t latest_observation_ms_{0};
    std::vector<std::uint8_t> rtcm_input_buffer_;
    std::deque<std::vector<std::uint8_t>> rtcm_frames_;
    std::size_t rtcm_frame_bytes_{0};
    std::uint64_t rtcm_frames_received_{0};
    std::uint64_t rtcm_solutions_received_{0};
    bool has_latest_rtcm_observation_{false};
    std::int64_t latest_rtcm_observation_ms_{0};
#ifdef ORBIT_ENABLE_RTCM_POSITIONING
    std::unique_ptr<RtcmPositionSolver> rtcm_position_solver_;
#endif
};

}  // namespace orbit_prediction

#endif

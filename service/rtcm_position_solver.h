#ifndef ORBIT_PREDICTION_RTCM_POSITION_SOLVER_H
#define ORBIT_PREDICTION_RTCM_POSITION_SOLVER_H

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace orbit_prediction {

struct EcefPositionObservation {
    std::int64_t timestamp_ms{};
    double x_m{};
    double y_m{};
    double z_m{};
    int quality{};
    int satellites{};
};

// Resolve RTCM 1019's 10-bit GPS week in the mission's current rollover era
// without consulting MQ time synchronization or the host system clock.
int ResolveGpsWeek10(std::uint16_t week10);

class RtcmPositionSolver final {
public:
    RtcmPositionSolver();
    ~RtcmPositionSolver();

    RtcmPositionSolver(const RtcmPositionSolver &) = delete;
    RtcmPositionSolver &operator=(const RtcmPositionSolver &) = delete;

    void Reset();
    bool FeedFrame(const std::vector<std::uint8_t> &frame,
                   std::vector<EcefPositionObservation> *positions,
                   std::string *diagnostic);

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace orbit_prediction

#endif

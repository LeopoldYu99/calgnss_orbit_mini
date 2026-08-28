#include "rtcm_position_solver.h"

extern "C" {
#include "rtklib.h"
}

#include <cmath>
#include <cstring>
#include <stdexcept>

namespace orbit_prediction {
namespace {

constexpr int kGpsWeekRolloverOffset = 2048;

std::uint16_t MessageType(const std::vector<std::uint8_t> &frame)
{
    if (frame.size() < 5) return 0;
    return static_cast<std::uint16_t>(
        (static_cast<std::uint16_t>(frame[3]) << 4U) | (frame[4] >> 4U));
}

std::uint32_t UnsignedBits(const std::vector<std::uint8_t> &frame,
                           std::size_t position, std::size_t length)
{
    std::uint32_t value = 0;
    for (std::size_t index = 0; index < length; ++index) {
        value = (value << 1U) |
                ((frame[(position + index) / 8U] >>
                  (7U - (position + index) % 8U)) & 1U);
    }
    return value;
}

std::int64_t GpsTimeToUnixMilliseconds(gtime_t gps_time)
{
    const gtime_t utc = gpst2utc(gps_time);
    return static_cast<std::int64_t>(utc.time) * 1000 +
           static_cast<std::int64_t>(std::llround(utc.sec * 1000.0));
}

}  // namespace

int ResolveGpsWeek10(std::uint16_t week10)
{
    // RTCM 1019 only carries 10 bits. The deployed mission is in the third
    // GPS rollover era (weeks 2048..3071), so resolve the week entirely from
    // the RTCM field and never from the host clock or MQ time-sync messages.
    return static_cast<int>(week10 & 0x03FFU) + kGpsWeekRolloverOffset;
}

class RtcmPositionSolver::Impl final {
public:
    Impl() { Initialize(); }
    ~Impl() { free_rtcm(&rtcm_); }

    void Initialize()
    {
        std::memset(&rtcm_, 0, sizeof(rtcm_));
        if (!init_rtcm(&rtcm_)) {
            throw std::runtime_error("RTKLIB init_rtcm failed");
        }
        options_ = prcopt_default;
        options_.mode = PMODE_SINGLE;
        options_.navsys = SYS_GPS | SYS_CMP;
        options_.nf = 1;
        options_.elmin = 10.0 * D2R;
        options_.ionoopt = IONOOPT_BRDC;
        options_.tropopt = TROPOPT_SAAS;
        options_.sateph = EPHOPT_BRDC;
    }

    void Reset()
    {
        free_rtcm(&rtcm_);
        Initialize();
    }

    bool FeedFrame(const std::vector<std::uint8_t> &frame,
                   std::vector<EcefPositionObservation> *positions,
                   std::string *diagnostic)
    {
        if (!positions) {
            if (diagnostic) *diagnostic = "RTCM position output is required";
            return false;
        }
        const std::uint16_t message_type = MessageType(frame);
        PrepareRtcmTimeFromEphemeris(frame, message_type);
        if (IsObservationMessage(message_type) && !has_rtcm_time_) {
            if (diagnostic) {
                *diagnostic = "waiting for RTCM ephemeris time reference";
            }
            return true;
        }

        bool ephemeris_decoded = false;
        for (const std::uint8_t byte : frame) {
            const int decode_result = input_rtcm3(&rtcm_, byte);
            if (decode_result == -1) {
                if (diagnostic) *diagnostic = "RTKLIB rejected an RTCM3 frame";
                return false;
            }
            if (decode_result == 2) ephemeris_decoded = true;
            if (decode_result != 1) continue;

            sol_t solution{};
            ssat_t satellite_status[MAXSAT]{};
            double azimuth_elevation[MAXOBS * 2]{};
            char message[256]{};
            if (!pntpos(rtcm_.obs.data, rtcm_.obs.n, &rtcm_.nav, &options_,
                        &solution, azimuth_elevation, satellite_status, message)) {
                if (diagnostic && message[0] != '\0') *diagnostic = message;
                continue;
            }
            if (solution.stat == SOLQ_NONE || solution.ns < 4 ||
                !std::isfinite(solution.rr[0]) || !std::isfinite(solution.rr[1]) ||
                !std::isfinite(solution.rr[2])) {
                continue;
            }
            positions->push_back({GpsTimeToUnixMilliseconds(solution.time),
                                  solution.rr[0], solution.rr[1], solution.rr[2],
                                  static_cast<int>(solution.stat),
                                  static_cast<int>(solution.ns)});
            has_rtcm_time_ = true;
            if (diagnostic) diagnostic->clear();
        }
        if (ephemeris_decoded) {
            CorrectDecodedEphemeris(frame, message_type);
        }
        return true;
    }

private:
    static bool IsObservationMessage(std::uint16_t message_type)
    {
        return message_type >= 1071 && message_type <= 1137;
    }

    void PrepareRtcmTimeFromEphemeris(const std::vector<std::uint8_t> &frame,
                                      std::uint16_t message_type)
    {
        if (message_type == 1019 && frame.size() * 8U >= 328U) {
            const int week = ResolveGpsWeek10(
                static_cast<std::uint16_t>(UnsignedBits(frame, 42, 10)));
            const double toe = static_cast<double>(UnsignedBits(frame, 312, 16)) * 16.0;
            rtcm_.time = gpst2time(week, toe);
            has_rtcm_time_ = true;
        } else if (message_type == 1042 && frame.size() * 8U >= 340U) {
            const int week = static_cast<int>(UnsignedBits(frame, 42, 13));
            const double toe = static_cast<double>(UnsignedBits(frame, 323, 17)) * 8.0;
            rtcm_.time = bdt2gpst(bdt2time(week, toe));
            has_rtcm_time_ = true;
        }
    }

    void CorrectDecodedEphemeris(const std::vector<std::uint8_t> &frame,
                                 std::uint16_t message_type)
    {
        if (rtcm_.ephsat <= 0 || rtcm_.ephsat > MAXSAT) return;
        eph_t &ephemeris = rtcm_.nav.eph[rtcm_.ephsat - 1];
        if (message_type == 1019 && frame.size() * 8U >= 328U) {
            const int week = ResolveGpsWeek10(
                static_cast<std::uint16_t>(UnsignedBits(frame, 42, 10)));
            const double toc = static_cast<double>(UnsignedBits(frame, 80, 16)) * 16.0;
            ephemeris.week = week;
            ephemeris.toe = gpst2time(week, ephemeris.toes);
            ephemeris.toc = gpst2time(week, toc);
            ephemeris.ttr = rtcm_.time;
        } else if (message_type == 1042 && frame.size() * 8U >= 340U) {
            const int week = static_cast<int>(UnsignedBits(frame, 42, 13));
            const double toc = static_cast<double>(UnsignedBits(frame, 78, 17)) * 8.0;
            ephemeris.week = week;
            ephemeris.toe = bdt2gpst(bdt2time(week, ephemeris.toes));
            ephemeris.toc = bdt2gpst(bdt2time(week, toc));
            ephemeris.ttr = rtcm_.time;
        }
    }

    rtcm_t rtcm_{};
    prcopt_t options_{};
    bool has_rtcm_time_{false};
};

RtcmPositionSolver::RtcmPositionSolver() : impl_(std::make_unique<Impl>()) {}
RtcmPositionSolver::~RtcmPositionSolver() = default;
void RtcmPositionSolver::Reset() { impl_->Reset(); }
bool RtcmPositionSolver::FeedFrame(const std::vector<std::uint8_t> &frame,
                                   std::vector<EcefPositionObservation> *positions,
                                   std::string *diagnostic)
{
    return impl_->FeedFrame(frame, positions, diagnostic);
}

}  // namespace orbit_prediction

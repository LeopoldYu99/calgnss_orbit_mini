#include "orbit_prediction_engine.h"

#include <cstdlib>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

void Require(bool condition, const std::string &message)
{
    if (!condition) {
        std::cerr << "FAILED: " << message << '\n';
        std::exit(1);
    }
}

std::string TestNmeaArc()
{
    std::ostringstream data;
    data << "$GNRMC,040403.00,A,1724.2594222,N,15304.2955483,E,0.0,0.0,230626,,,A\n";
    for (int index = 0; index < 12; ++index) {
        const int second = 3 + index;
        const double longitude = 15300.0 + 4.2955483 + index * 4.2;
        data << "$GNGGA,0404" << std::setw(2) << std::setfill('0') << second
             << ".00,1724.2594222,N," << std::fixed << std::setprecision(7) << longitude
             << ",E,1,14,0.69,271441.347,M,0.989,M,,\n";
    }
    return data.str();
}

std::string TestRTCMFrame()
{
    // RTCM3 type 1077 test frame. CRC24Q for d3000443500000 is 44fe2e.
    const std::vector<std::uint8_t> frame{
        0xD3, 0x00, 0x04, 0x43, 0x50, 0x00, 0x00, 0x44, 0xFE, 0x2E};
    return std::string(reinterpret_cast<const char *>(frame.data()), frame.size());
}

}  // namespace

int main()
{
    orbit_prediction::EngineOptions options;
    options.observation_capacity = 32;
    options.fit_degree = 3;
    options.stream_batch_size = 1;
    orbit_prediction::OrbitPredictionEngine engine(options);
    std::string error;

    Require(engine.ReceiveTimeSync(1782187443000LL, &error), "time sync: " + error);
    Require(engine.ReceiveUplinkData(1, "", TestNmeaArc(), &error), "NMEA uplink: " + error);
    Require(engine.ObservationCount() == 12, "all GGA observations must be stored");

    const std::string rtcm_frame = TestRTCMFrame();
    Require(engine.ReceiveRTCMData(rtcm_frame.substr(0, 5), &error),
            "partial RTCM frame: " + error);
    Require(engine.RTCMFrameCount() == 0, "partial RTCM frame must not be counted yet");
    Require(engine.ReceiveRTCMData(rtcm_frame.substr(5), &error),
            "completed RTCM frame: " + error);
    Require(engine.RTCMFrameCount() == 1, "one RTCM frame must be accepted");

    std::string corrupt_rtcm = rtcm_frame;
    corrupt_rtcm.back() = static_cast<char>(corrupt_rtcm.back() ^ 0x01);
    Require(!engine.ReceiveRTCMData(corrupt_rtcm, &error),
            "RTCM frame with invalid CRC must be rejected");

    std::vector<orbit_prediction::PredictionPoint> predictions;
    Require(engine.PredictOrbit(1782187454LL, 2, 0.5,
                [&predictions](const auto &batch) {
                    predictions.insert(predictions.end(), batch.begin(), batch.end());
                    return true;
                },
                []() { return false; }, &error),
            "orbit prediction: " + error);
    Require(predictions.size() == 5, "five prediction points must be streamed");
    Require(predictions[0].timestamp_ms == 1782187454000LL,
            "first point must match start_time_s");
    Require(predictions[1].timestamp_ms == 1782187454500LL,
            "fractional prediction steps must use millisecond timestamps");

    engine.Stop();
    Require(engine.IsStopped(), "Stop must set stopped state");
    Require(!engine.PredictOrbit(1782187454LL, 2, 1.0, [](const auto &) { return true; },
                                 []() { return false; }, &error),
            "PredictOrbit must reject calls after Stop");

    Require(engine.Reset(&error), "reset: " + error);
    Require(!engine.IsStopped(), "Reset must clear stopped state");
    Require(engine.ObservationCount() == 0, "Reset must clear observations");
    Require(engine.RTCMFrameCount() == 0, "Reset must clear RTCM frames");
    Require(!engine.ReceiveUplinkData(1, "", "not NMEA", &error),
            "invalid uplink data must be rejected");

    std::cout << "orbit_prediction_engine_test passed\n";
    return 0;
}

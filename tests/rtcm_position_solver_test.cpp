#include "rtcm_position_solver.h"

#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <vector>

int main(int argc, char **argv)
{
    try {
        if (orbit_prediction::ResolveGpsWeek10(385) != 2433) {
            throw std::runtime_error("GPS week 385 must resolve to full week 2433");
        }
        orbit_prediction::RtcmPositionSolver solver;
        if (argc > 1) {
            std::ifstream input(argv[1], std::ios::binary);
            if (!input) throw std::runtime_error("cannot open RTCM capture");
            const std::vector<std::uint8_t> bytes{
                std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>()};
            std::vector<orbit_prediction::EcefPositionObservation> positions;
            std::size_t cursor = 0;
            while (cursor + 6 <= bytes.size()) {
                if (bytes[cursor] != 0xD3U || (bytes[cursor + 1] & 0xFCU) != 0) {
                    ++cursor;
                    continue;
                }
                const std::size_t payload_size =
                    (static_cast<std::size_t>(bytes[cursor + 1] & 0x03U) << 8) |
                    bytes[cursor + 2];
                const std::size_t end = cursor + payload_size + 6;
                if (end > bytes.size()) break;
                const std::vector<std::uint8_t> frame(bytes.begin() + cursor,
                                                      bytes.begin() + end);
                std::string diagnostic;
                if (!solver.FeedFrame(frame, &positions, &diagnostic)) {
                    throw std::runtime_error("RTCM decode failed: " + diagnostic);
                }
                cursor = end;
            }
            if (positions.size() < 1000) {
                throw std::runtime_error("expected at least 1000 ECEF positions");
            }
            const auto &first = positions.front();
            const auto &last = positions.back();
            std::cout << std::fixed << std::setprecision(4)
                      << "RTCM SPP positions=" << positions.size()
                      << " first=" << first.timestamp_ms << ',' << first.x_m << ','
                      << first.y_m << ',' << first.z_m
                      << " last=" << last.timestamp_ms << ',' << last.x_m << ','
                      << last.y_m << ',' << last.z_m << '\n';
        }
        solver.Reset();
        std::cout << "RTCM position solver lifecycle test passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "RTCM position solver lifecycle test failed: " << error.what() << '\n';
        return 1;
    }
}

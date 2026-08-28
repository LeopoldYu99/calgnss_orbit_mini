#include "rtcm_position_solver.h"

#include <cstdint>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::uint16_t MessageType(const std::vector<std::uint8_t> &frame)
{
    if (frame.size() < 5) return 0;
    return static_cast<std::uint16_t>(
        (static_cast<std::uint16_t>(frame[3]) << 4U) | (frame[4] >> 4U));
}

}  // namespace

int main(int argc, char **argv)
{
    try {
        if (argc != 2) {
            std::cerr << "usage: " << argv[0] << " RTCM_FILE\n";
            return 2;
        }
        std::ifstream input(argv[1], std::ios::binary);
        if (!input) throw std::runtime_error("cannot open RTCM file");
        const std::vector<std::uint8_t> bytes{
            std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>()};

        orbit_prediction::RtcmPositionSolver solver;
        std::size_t offset = 0;
        std::size_t frames = 0;
        std::size_t positions = 0;
        while (offset + 6 <= bytes.size()) {
            if (bytes[offset] != 0xD3U) {
                ++offset;
                continue;
            }
            const std::size_t payload_length =
                (static_cast<std::size_t>(bytes[offset + 1] & 0x03U) << 8U) |
                bytes[offset + 2];
            const std::size_t frame_length = payload_length + 6U;
            if (offset + frame_length > bytes.size()) break;
            std::vector<std::uint8_t> frame(bytes.begin() + offset,
                                            bytes.begin() + offset + frame_length);
            std::vector<orbit_prediction::EcefPositionObservation> decoded;
            std::string diagnostic;
            solver.FeedFrame(frame, &decoded, &diagnostic);
            ++frames;
            if (!diagnostic.empty()) {
                std::cout << "frame=" << frames << " type=" << MessageType(frame)
                          << " diagnostic=\"" << diagnostic << "\"\n";
            }
            for (const auto &position : decoded) {
                ++positions;
                std::cout << "position=" << positions
                          << " timestamp_ms=" << position.timestamp_ms
                          << " x=" << position.x_m << " y=" << position.y_m
                          << " z=" << position.z_m << " ns=" << position.satellites
                          << " quality=" << position.quality << '\n';
            }
            offset += frame_length;
        }
        std::cout << "summary frames=" << frames << " positions=" << positions << '\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "rtcm_position_debug: " << error.what() << '\n';
        return 1;
    }
}

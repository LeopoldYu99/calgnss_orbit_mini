#include "orbit_mq_protocol.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void Require(bool condition, const std::string &message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

}  // namespace

int main()
{
    try {
        orbit_prediction::MqFrame source;
        source.message_type = 3;
        source.message_seq = 0x12345678U;
        source.stream_end = orbit_prediction::kMqStreamEnd;
        source.payload = {0xD3, 0x00, 0x04};
        std::vector<std::uint8_t> bytes;
        std::string error;
        Require(orbit_prediction::BuildMqFrame(source, 8192, &bytes, &error), error);
        Require(bytes.size() == 13, "frame length must be header plus payload");
        Require(bytes[0] == 3 && bytes[1] == 0x78 && bytes[2] == 0x56 &&
                    bytes[3] == 0x34 && bytes[4] == 0x12 && bytes[5] == 0xFF,
                "CCSM frame prefix mismatch");
        Require(bytes[6] == 3 && bytes[7] == 0 && bytes[8] == 0 && bytes[9] == 0,
                "CCSM uint32 payload length mismatch");

        orbit_prediction::MqFrame decoded;
        Require(orbit_prediction::ParseMqFrame(bytes.data(), bytes.size(), &decoded, &error),
                error);
        Require(decoded.message_type == source.message_type &&
                    decoded.message_seq == source.message_seq &&
                    decoded.stream_end == source.stream_end &&
                    decoded.payload == source.payload,
                "frame round trip mismatch");

        bytes[5] = orbit_prediction::kMqStreamEndLegacy;
        Require(orbit_prediction::ParseMqFrame(bytes.data(), bytes.size(), &decoded, &error) &&
                    orbit_prediction::IsMqStreamEnd(decoded.stream_end),
                "live CCSM 0x01 stream-end marker was rejected");
        bytes[5] = orbit_prediction::kMqStreamEnd;
        bytes[6] = 4;
        Require(!orbit_prediction::ParseMqFrame(bytes.data(), bytes.size(), &decoded, &error),
                "mismatched payload length must be rejected");

        orbit_prediction::StatusReply status;
        status.timestamp_ms = 1787283546000LL;
        status.status = orbit_prediction::ServiceStatus::Running;
        status.message = "running";
        std::vector<std::uint8_t> payload;
        Require(orbit_prediction::EncodeStatusReply(status, &payload, &error), error);
        orbit_prediction::StatusReply decoded_status;
        Require(orbit_prediction::DecodeStatusReply(payload, &decoded_status, &error), error);
        Require(decoded_status.timestamp_ms == status.timestamp_ms &&
                    decoded_status.status == status.status &&
                    decoded_status.message == status.message,
                "StatusReply round trip mismatch");

        orbit_prediction::OrbitData orbit;
        orbit.points.push_back({1787283546000LL, 1, 2, 3, 4, 5, 6});
        Require(orbit_prediction::EncodeOrbitData(orbit, &payload, &error), error);
        Require(payload.size() == 56, "one raw OrbitPoint payload must be 56 bytes");
        orbit_prediction::OrbitData decoded_orbit;
        Require(orbit_prediction::DecodeOrbitData(payload, &decoded_orbit, &error), error);
        Require(decoded_orbit.points.size() == 1 &&
                    std::fabs(decoded_orbit.points[0].vz - 6.0) < 1.0e-12,
                "OrbitData round trip mismatch");

        std::cout << "CCSM fixed-binary MQ protocol test passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "CCSM fixed-binary MQ protocol test failed: " << error.what() << '\n';
        return 1;
    }
}

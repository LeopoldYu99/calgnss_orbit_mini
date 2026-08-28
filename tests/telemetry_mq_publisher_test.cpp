#include "telemetry_mq_publisher.h"

#include <fcntl.h>
#include <iostream>
#include <mqueue.h>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr const char *kQueueName = "/orbit_telemetry_publisher_test";
void Require(bool condition, const std::string &message)
{
    if (!condition) throw std::runtime_error(message);
}
}  // namespace

int main()
{
    mq_unlink(kQueueName);
    mq_attr attributes{};
    attributes.mq_maxmsg = 10;
    attributes.mq_msgsize = 8192;
    const mqd_t queue = mq_open(kQueueName, O_CREAT | O_RDONLY, 0600, &attributes);
    if (queue == static_cast<mqd_t>(-1)) return 1;
    int result = 1;
    try {
        orbit_prediction::TelemetryPacket packet;
        packet.port = 7;
        packet.packet_header = "header";
        packet.data = "telemetry";
        orbit_prediction::TelemetryMqPublisher publisher(kQueueName, 8192);
        Require(publisher.Publish(packet, 42), "telemetry frame was not queued");
        std::vector<std::uint8_t> bytes(8192);
        const ssize_t length = mq_receive(
            queue, reinterpret_cast<char *>(bytes.data()), bytes.size(), nullptr);
        Require(length > 0, "telemetry frame was not received");
        orbit_prediction::MqFrame frame;
        std::string error;
        Require(orbit_prediction::ParseMqFrame(
                    bytes.data(), static_cast<std::size_t>(length), &frame, &error), error);
        Require(frame.message_type ==
                    static_cast<std::uint8_t>(orbit_prediction::MqRecvType::TelemetryData) &&
                    frame.message_seq == 42,
                "telemetry frame header mismatch");
        orbit_prediction::TelemetryPacket decoded;
        Require(orbit_prediction::DecodeTelemetryPacket(frame.payload, &decoded, &error), error);
        Require(decoded.port == packet.port && decoded.packet_header == packet.packet_header &&
                    decoded.data == packet.data,
                "telemetry payload mismatch");
        std::cout << "CCSM telemetry MQ publisher test passed\n";
        result = 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
    }
    mq_close(queue);
    mq_unlink(kQueueName);
    return result;
}

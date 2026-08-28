#include "heartbeat_mq_publisher.h"

#include <fcntl.h>
#include <iostream>
#include <mqueue.h>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr const char *kQueueName = "/orbit_heartbeat_publisher_test";
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
        orbit_prediction::StatusReply status;
        status.timestamp_ms = 1787283546000LL;
        status.status = orbit_prediction::ServiceStatus::Running;
        status.message = "orbit prediction is running";
        orbit_prediction::HeartbeatMqPublisher publisher(kQueueName, 8192);
        std::string error;
        Require(publisher.Publish(status, 99, &error), error);
        std::vector<std::uint8_t> bytes(8192);
        const ssize_t length = mq_receive(
            queue, reinterpret_cast<char *>(bytes.data()), bytes.size(), nullptr);
        Require(length > 0, "heartbeat frame was not received");
        orbit_prediction::MqFrame frame;
        Require(orbit_prediction::ParseMqFrame(
                    bytes.data(), static_cast<std::size_t>(length), &frame, &error), error);
        Require(frame.message_type ==
                    static_cast<std::uint8_t>(orbit_prediction::MqRecvType::Heartbeat) &&
                    frame.message_seq == 99 &&
                    frame.stream_end == orbit_prediction::kMqStreamEnd,
                "heartbeat frame header mismatch");
        orbit_prediction::StatusReply decoded;
        Require(orbit_prediction::DecodeStatusReply(frame.payload, &decoded, &error), error);
        Require(decoded.timestamp_ms == status.timestamp_ms &&
                    decoded.status == status.status && decoded.message == status.message,
                "heartbeat payload mismatch");
        std::cout << "CCSM heartbeat MQ publisher test passed\n";
        result = 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
    }
    mq_close(queue);
    mq_unlink(kQueueName);
    return result;
}

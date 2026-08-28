#include "telemetry_mq_publisher.h"

#include <cerrno>
#include <cstring>
#include <fcntl.h>
#include <iostream>
#include <mqueue.h>
#include <utility>

namespace orbit_prediction {
namespace {

bool IsValidQueueName(const std::string &name)
{
    return name.size() >= 2 && name.front() == '/' &&
           name.find('/', 1) == std::string::npos;
}

}  // namespace

TelemetryMqPublisher::TelemetryMqPublisher(std::string lifecycle_request_queue,
                                           std::size_t maximum_message_size)
    : lifecycle_request_queue_(std::move(lifecycle_request_queue)),
      maximum_message_size_(maximum_message_size)
{
}

bool TelemetryMqPublisher::Publish(const TelemetryPacket &packet,
                                   std::uint32_t message_seq) const
{
    if (!IsValidQueueName(lifecycle_request_queue_)) {
        std::cerr << "telemetry MQ publish rejected: invalid queue name\n";
        return false;
    }

    MqFrame frame;
    frame.message_type = static_cast<std::uint8_t>(MqRecvType::TelemetryData);
    frame.message_seq = message_seq;
    frame.stream_end = kMqStreamEnd;
    std::string error;
    if (!EncodeTelemetryPacket(packet, &frame.payload, &error)) {
        std::cerr << "cannot encode telemetry payload: " << error << '\n';
        return false;
    }
    std::vector<std::uint8_t> bytes;
    if (!BuildMqFrame(frame, maximum_message_size_, &bytes, &error)) {
        std::cerr << "cannot build telemetry MQ frame: " << error << '\n';
        return false;
    }

    const mqd_t queue = mq_open(lifecycle_request_queue_.c_str(), O_WRONLY | O_NONBLOCK);
    if (queue == static_cast<mqd_t>(-1)) {
        std::cerr << "telemetry MQ is unavailable; sequence " << message_seq
                  << " was not queued: " << std::strerror(errno) << '\n';
        return false;
    }
    const int status = mq_send(
        queue, reinterpret_cast<const char *>(bytes.data()), bytes.size(), 0);
    const int send_error = errno;
    mq_close(queue);
    if (status != 0) {
        std::cerr << "telemetry MQ sequence " << message_seq
                  << " was not queued: " << std::strerror(send_error) << '\n';
        return false;
    }
    return true;
}

}  // namespace orbit_prediction

#include "heartbeat_mq_publisher.h"

#include <cerrno>
#include <cstring>
#include <fcntl.h>
#include <mqueue.h>
#include <utility>

namespace orbit_prediction {
namespace {

bool IsValidQueueName(const std::string &name)
{
    return name.size() >= 2 && name.front() == '/' &&
           name.find('/', 1) == std::string::npos;
}

void SetError(std::string *error, const std::string &message)
{
    if (error) {
        *error = message;
    }
}

}  // namespace

HeartbeatMqPublisher::HeartbeatMqPublisher(std::string lifecycle_request_queue,
                                           std::size_t maximum_message_size)
    : lifecycle_request_queue_(std::move(lifecycle_request_queue)),
      maximum_message_size_(maximum_message_size)
{
}

bool HeartbeatMqPublisher::Publish(const StatusReply &status,
                                   std::uint32_t message_seq,
                                   std::string *error) const
{
    if (!IsValidQueueName(lifecycle_request_queue_)) {
        SetError(error, "invalid lifecycle queue name");
        return false;
    }

    MqFrame frame;
    frame.message_type = static_cast<std::uint8_t>(MqRecvType::Heartbeat);
    frame.message_seq = message_seq;
    frame.stream_end = kMqStreamEnd;
    if (!EncodeStatusReply(status, &frame.payload, error)) {
        return false;
    }
    std::vector<std::uint8_t> bytes;
    if (!BuildMqFrame(frame, maximum_message_size_, &bytes, error)) {
        return false;
    }

    const mqd_t queue = mq_open(lifecycle_request_queue_.c_str(), O_WRONLY | O_NONBLOCK);
    if (queue == static_cast<mqd_t>(-1)) {
        SetError(error, std::string("cannot open lifecycle MQ: ") + std::strerror(errno));
        return false;
    }
    const int send_status = mq_send(
        queue, reinterpret_cast<const char *>(bytes.data()), bytes.size(), 0);
    const int send_error = errno;
    mq_close(queue);
    if (send_status != 0) {
        SetError(error, std::string("cannot queue heartbeat: ") +
                            std::strerror(send_error));
        return false;
    }
    if (error) {
        error->clear();
    }
    return true;
}

}  // namespace orbit_prediction

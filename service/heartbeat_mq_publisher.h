#ifndef ORBIT_PREDICTION_HEARTBEAT_MQ_PUBLISHER_H
#define ORBIT_PREDICTION_HEARTBEAT_MQ_PUBLISHER_H

#include "orbit_mq_protocol.h"

#include <cstddef>
#include <cstdint>
#include <string>

namespace orbit_prediction {

// Best-effort, fire-and-forget heartbeat publisher to the lifecycle service.
class HeartbeatMqPublisher final {
public:
    explicit HeartbeatMqPublisher(
        std::string lifecycle_request_queue = "/csm_end0_to_main",
        std::size_t maximum_message_size = 2048);

    bool Publish(const StatusReply &status, std::uint32_t message_seq,
                 std::string *error = nullptr) const;

private:
    std::string lifecycle_request_queue_;
    std::size_t maximum_message_size_;
};

}  // namespace orbit_prediction

#endif

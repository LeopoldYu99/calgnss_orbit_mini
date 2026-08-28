#ifndef ORBIT_PREDICTION_TELEMETRY_MQ_PUBLISHER_H
#define ORBIT_PREDICTION_TELEMETRY_MQ_PUBLISHER_H

#include "orbit_mq_protocol.h"

#include <cstddef>
#include <cstdint>
#include <string>

namespace orbit_prediction {

// Best-effort, fire-and-forget publisher for LifecycleService.ReceiveTelemetryData.
// It only confirms whether the request entered the lifecycle MQ; it never waits for a reply.
class TelemetryMqPublisher final {
public:
    explicit TelemetryMqPublisher(
        std::string lifecycle_request_queue = "/csm_end0_to_main",
        std::size_t maximum_message_size = 2048);

    bool Publish(const TelemetryPacket &packet, std::uint32_t message_seq) const;

private:
    std::string lifecycle_request_queue_;
    std::size_t maximum_message_size_;
};

}  // namespace orbit_prediction

#endif

#ifndef ORBIT_PREDICTION_SERVICE_H
#define ORBIT_PREDICTION_SERVICE_H

#include "orbit_prediction_engine.h"
#include "orbit_mq_protocol.h"

#include <atomic>
#include <functional>

namespace orbit_prediction {

class OrbitPredictionHandler final {
public:
    using BatchCallback = std::function<void(const std::vector<PredictionPoint> &)>;

    explicit OrbitPredictionHandler(EngineOptions options = {});

    CommonReply ReceiveUplinkData(const UplinkPacket &request);
    CommonReply ReceiveRTCMData(const RtcmData &request);
    StatusReply GetStatus() const;
    CommonReply Stop();
    CommonReply Reset();

    bool BeginPrediction();
    CommonReply RunPrediction(const OrbitPredictionRequest &request,
                              const BatchCallback &on_batch);
    bool IsPredictionRunning() const;

private:
    static CommonReply Reply(bool success, const std::string &message);

    OrbitPredictionEngine engine_;
    std::atomic<bool> prediction_running_{false};
};

}  // namespace orbit_prediction

#endif

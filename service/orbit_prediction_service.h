#ifndef ORBIT_PREDICTION_SERVICE_H
#define ORBIT_PREDICTION_SERVICE_H

#include "orbit_prediction.grpc.pb.h"
#include "orbit_prediction_engine.h"

namespace orbit_prediction {

class OrbitPredictionServiceImpl final : public OrbitPredictionService::Service {
public:
    explicit OrbitPredictionServiceImpl(EngineOptions options = {});

    grpc::Status ReceiveUplinkData(
        grpc::ServerContext *context,
        const UplinkPacket *request,
        CommonReply *reply) override;

    grpc::Status ReceiveTimeSync(
        grpc::ServerContext *context,
        const TimeSyncData *request,
        CommonReply *reply) override;

    grpc::Status ReceiveRTCMData(
        grpc::ServerContext *context,
        const RTCMData *request,
        CommonReply *reply) override;

    grpc::Status Stop(
        grpc::ServerContext *context,
        const google::protobuf::Empty *request,
        CommonReply *reply) override;

    grpc::Status Reset(
        grpc::ServerContext *context,
        const google::protobuf::Empty *request,
        CommonReply *reply) override;

    grpc::Status PredictOrbit(
        grpc::ServerContext *context,
        const OrbitPredictionRequest *request,
        grpc::ServerWriter<OrbitData> *writer) override;

private:
    OrbitPredictionEngine engine_;
};

}  // namespace orbit_prediction

#endif

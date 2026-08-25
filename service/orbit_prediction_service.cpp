#include "orbit_prediction_service.h"

#include <iostream>
#include <utility>

namespace orbit_prediction {

OrbitPredictionServiceImpl::OrbitPredictionServiceImpl(EngineOptions options)
    : engine_(std::move(options))
{
}

grpc::Status OrbitPredictionServiceImpl::ReceiveUplinkData(
    grpc::ServerContext *context, const UplinkPacket *request, CommonReply *reply)
{
    (void)context;
    std::string error;
    const bool success = request && engine_.ReceiveUplinkData(
                                        request->port(), request->packet_header(), request->data(),
                                        &error);
    reply->set_success(success);
    reply->set_message(success ? "uplink NMEA data accepted" :
                                 (request ? error : "missing request"));
    if (!success) {
        std::cerr << "ReceiveUplinkData rejected: "
                  << (request ? error : "missing request") << '\n';
    }
    return grpc::Status::OK;
}

grpc::Status OrbitPredictionServiceImpl::ReceiveTimeSync(
    grpc::ServerContext *context, const TimeSyncData *request, CommonReply *reply)
{
    (void)context;
    std::string error;
    const bool success = request && engine_.ReceiveTimeSync(request->timestamp_ms(), &error);
    reply->set_success(success);
    reply->set_message(success ? "time synchronization accepted" :
                                 (request ? error : "missing request"));
    if (!success) {
        std::cerr << "ReceiveTimeSync rejected: "
                  << (request ? error : "missing request") << '\n';
    }
    return grpc::Status::OK;
}

grpc::Status OrbitPredictionServiceImpl::ReceiveRTCMData(
    grpc::ServerContext *context, const RTCMData *request, CommonReply *reply)
{
    (void)context;
    std::string message;
    const bool success = request && engine_.ReceiveRTCMData(request->data(), &message);
    reply->set_success(success);
    reply->set_message(request ? message : "missing request");
    if (!success) {
        std::cerr << "ReceiveRTCMData rejected: " << reply->message() << '\n';
    }
    return grpc::Status::OK;
}

grpc::Status OrbitPredictionServiceImpl::Stop(
    grpc::ServerContext *context, const google::protobuf::Empty *request, CommonReply *reply)
{
    (void)context;
    (void)request;
    engine_.Stop();
    reply->set_success(true);
    reply->set_message("orbit prediction stopped");
    return grpc::Status::OK;
}

grpc::Status OrbitPredictionServiceImpl::Reset(
    grpc::ServerContext *context, const google::protobuf::Empty *request, CommonReply *reply)
{
    (void)context;
    (void)request;
    std::string error;
    const bool success = engine_.Reset(&error);
    reply->set_success(success);
    reply->set_message(success ? "orbit prediction service reset" : error);
    if (!success) {
        std::cerr << "Reset failed: " << error << '\n';
    }
    return grpc::Status::OK;
}

grpc::Status OrbitPredictionServiceImpl::PredictOrbit(
    grpc::ServerContext *context, const OrbitPredictionRequest *request,
    grpc::ServerWriter<OrbitData> *writer)
{
    if (!request) {
        return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, "missing request");
    }
    std::string error;
    const bool success = engine_.PredictOrbit(
        request->start_time_s(), request->duration_s(), request->step_s(),
        [writer](const std::vector<PredictionPoint> &batch) {
            OrbitData response;
            for (const auto &value : batch) {
                OrbitPoint *point = response.add_points();
                point->set_timestamp_ms(value.timestamp_ms);
                point->set_x(value.x);
                point->set_y(value.y);
                point->set_z(value.z);
                point->set_vx(value.vx);
                point->set_vy(value.vy);
                point->set_vz(value.vz);
            }
            return writer->Write(response);
        },
        [context]() { return context->IsCancelled(); }, &error);

    if (!success) {
        std::cerr << "PredictOrbit failed: " << error << '\n';
        return grpc::Status(grpc::StatusCode::FAILED_PRECONDITION, error);
    }
    return grpc::Status::OK;
}

}  // namespace orbit_prediction

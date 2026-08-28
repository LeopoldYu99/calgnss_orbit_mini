#ifndef ORBIT_PREDICTION_MQ_PROTOCOL_H
#define ORBIT_PREDICTION_MQ_PROTOCOL_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace orbit_prediction {

// CCSM frame layout (little-endian multi-byte fields):
//   u8 message_type | u32 message_seq | u8 stream_end | u32 payload_length
constexpr std::size_t kMqFrameHeaderSize = 10;
constexpr std::uint8_t kMqStreamMore = 0x00;
// The live CCSM producer currently uses 0x01, while the written interface
// document specifies 0xFF. Receivers accept both; this service sends 0xFF.
constexpr std::uint8_t kMqStreamEndLegacy = 0x01;
constexpr std::uint8_t kMqStreamEnd = 0xFF;

constexpr bool IsMqStreamEnd(std::uint8_t value)
{
    return value == kMqStreamEndLegacy || value == kMqStreamEnd;
}

enum class MqSendType : std::uint8_t {
    Unspecified = 0,
    ReceiveUplinkData = 1,
    ReceiveTimeSync = 2,
    ReceiveRtcmData = 3,
    GetStatus = 4,
    Stop = 5,
    Reset = 6,
    PredictOrbit = 7,
};

enum class MqRecvType : std::uint8_t {
    Unspecified = 0,
    Heartbeat = 1,
    TelemetryData = 2,
    TelemetryStatus = 3,
    ReceiveTelemetryData = 4,
};

constexpr const char *MqSendTypeName(MqSendType type)
{
    switch (type) {
    case MqSendType::ReceiveUplinkData: return "ReceiveUplinkData";
    case MqSendType::ReceiveTimeSync: return "ReceiveTimeSync";
    case MqSendType::ReceiveRtcmData: return "ReceiveRtcmData";
    case MqSendType::GetStatus: return "GetStatus";
    case MqSendType::Stop: return "Stop";
    case MqSendType::Reset: return "Reset";
    case MqSendType::PredictOrbit: return "PredictOrbit";
    case MqSendType::Unspecified: return "Unspecified";
    }
    return "Unknown";
}

constexpr const char *MqRecvTypeName(MqRecvType type)
{
    switch (type) {
    case MqRecvType::Heartbeat: return "Heartbeat";
    case MqRecvType::TelemetryData: return "TelemetryData";
    case MqRecvType::TelemetryStatus: return "TelemetryStatus";
    case MqRecvType::ReceiveTelemetryData: return "ReceiveTelemetryData";
    case MqRecvType::Unspecified: return "Unspecified";
    }
    return "Unknown";
}

enum class ServiceStatus : std::uint8_t {
    Unspecified = 0,
    Idle = 1,
    Running = 2,
    Error = 3,
};

struct MqFrame {
    std::uint8_t message_type{};
    std::uint32_t message_seq{};
    std::uint8_t stream_end{kMqStreamEnd};
    std::vector<std::uint8_t> payload;
};

struct UplinkPacket {
    std::uint32_t port{};
    std::string packet_header;
    std::string data;
};

// CCSM specification uses Unix seconds for the time-sync payload.
struct TimeSyncData {
    std::int64_t timestamp_s{};
};

struct RtcmData {
    std::string data;
};

struct OrbitPredictionRequest {
    std::int64_t start_time_s{};
    std::uint32_t duration_s{};
};

struct OrbitPoint {
    std::int64_t timestamp_ms{};
    double x{};
    double y{};
    double z{};
    double vx{};
    double vy{};
    double vz{};
};

struct OrbitData {
    std::vector<OrbitPoint> points;
};

struct TelemetryPacket {
    std::uint32_t port{};
    std::string packet_header;
    std::string data;
};

struct StatusReply {
    std::int64_t timestamp_ms{};
    ServiceStatus status{ServiceStatus::Unspecified};
    std::string message;
};

struct CommonReply {
    std::int64_t timestamp_ms{};
    bool success{};
    std::string message;
};

bool BuildMqFrame(const MqFrame &frame, std::size_t maximum_size,
                  std::vector<std::uint8_t> *bytes, std::string *error = nullptr);
bool ParseMqFrame(const void *bytes, std::size_t length, MqFrame *frame,
                  std::string *error = nullptr);

bool EncodeUplinkPacket(const UplinkPacket &packet, std::vector<std::uint8_t> *payload,
                        std::string *error = nullptr);
bool DecodeUplinkPacket(const std::vector<std::uint8_t> &payload, UplinkPacket *packet,
                        std::string *error = nullptr);
bool EncodeTimeSyncData(const TimeSyncData &time_sync,
                        std::vector<std::uint8_t> *payload);
bool DecodeTimeSyncData(const std::vector<std::uint8_t> &payload, TimeSyncData *time_sync,
                        std::string *error = nullptr);
bool EncodePredictionRequest(const OrbitPredictionRequest &request,
                             std::vector<std::uint8_t> *payload);
bool DecodePredictionRequest(const std::vector<std::uint8_t> &payload,
                             OrbitPredictionRequest *request,
                             std::string *error = nullptr);
bool EncodeStatusReply(const StatusReply &status, std::vector<std::uint8_t> *payload,
                       std::string *error = nullptr);
bool DecodeStatusReply(const std::vector<std::uint8_t> &payload, StatusReply *status,
                       std::string *error = nullptr);
bool EncodeTelemetryPacket(const TelemetryPacket &packet,
                           std::vector<std::uint8_t> *payload,
                           std::string *error = nullptr);
bool DecodeTelemetryPacket(const std::vector<std::uint8_t> &payload,
                           TelemetryPacket *packet, std::string *error = nullptr);
bool EncodeOrbitData(const OrbitData &data, std::vector<std::uint8_t> *payload,
                     std::string *error = nullptr);
bool DecodeOrbitData(const std::vector<std::uint8_t> &payload, OrbitData *data,
                     std::string *error = nullptr);

}  // namespace orbit_prediction

#endif

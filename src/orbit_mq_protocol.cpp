#include "orbit_mq_protocol.h"

#include <cstring>
#include <limits>
#include <type_traits>

namespace orbit_prediction {
namespace {

void SetError(std::string *error, const std::string &message)
{
    if (error) {
        *error = message;
    }
}

template <typename T>
void AppendUnsigned(std::vector<std::uint8_t> *output, T value)
{
    static_assert(std::is_unsigned<T>::value, "unsigned integer required");
    for (std::size_t index = 0; index < sizeof(T); ++index) {
        output->push_back(static_cast<std::uint8_t>(value & 0xFFU));
        value >>= 8U;
    }
}

void AppendInt64(std::vector<std::uint8_t> *output, std::int64_t value)
{
    std::uint64_t bits{};
    std::memcpy(&bits, &value, sizeof(bits));
    AppendUnsigned(output, bits);
}

void AppendDouble(std::vector<std::uint8_t> *output, double value)
{
    static_assert(sizeof(double) == sizeof(std::uint64_t), "64-bit double required");
    std::uint64_t bits{};
    std::memcpy(&bits, &value, sizeof(bits));
    AppendUnsigned(output, bits);
}

class Reader final {
public:
    Reader(const std::uint8_t *data, std::size_t size) : data_(data), size_(size) {}

    template <typename T>
    bool ReadUnsigned(T *value)
    {
        static_assert(std::is_unsigned<T>::value, "unsigned integer required");
        if (!value || remaining() < sizeof(T)) {
            return false;
        }
        T decoded{};
        for (std::size_t index = 0; index < sizeof(T); ++index) {
            decoded |= static_cast<T>(data_[offset_ + index]) << (index * 8U);
        }
        offset_ += sizeof(T);
        *value = decoded;
        return true;
    }

    bool ReadInt64(std::int64_t *value)
    {
        std::uint64_t bits{};
        if (!ReadUnsigned(&bits) || !value) {
            return false;
        }
        std::memcpy(value, &bits, sizeof(bits));
        return true;
    }

    bool ReadDouble(double *value)
    {
        std::uint64_t bits{};
        if (!ReadUnsigned(&bits) || !value) {
            return false;
        }
        std::memcpy(value, &bits, sizeof(bits));
        return true;
    }

    bool ReadString(std::size_t length, std::string *value)
    {
        if (!value || remaining() < length) {
            return false;
        }
        value->assign(reinterpret_cast<const char *>(data_ + offset_), length);
        offset_ += length;
        return true;
    }

    std::size_t remaining() const { return size_ - offset_; }
    bool done() const { return offset_ == size_; }

private:
    const std::uint8_t *data_{};
    std::size_t size_{};
    std::size_t offset_{};
};

bool AppendString(std::vector<std::uint8_t> *output, const std::string &value,
                  std::string *error)
{
    if (value.size() > std::numeric_limits<std::uint32_t>::max()) {
        SetError(error, "binary string exceeds uint32 length");
        return false;
    }
    AppendUnsigned(output, static_cast<std::uint32_t>(value.size()));
    output->insert(output->end(), value.begin(), value.end());
    return true;
}

template <typename Packet>
bool EncodePortPacket(const Packet &packet, std::vector<std::uint8_t> *payload,
                      std::string *error)
{
    if (!payload) {
        SetError(error, "output payload is null");
        return false;
    }
    payload->clear();
    AppendUnsigned(payload, packet.port);
    return AppendString(payload, packet.packet_header, error) &&
           AppendString(payload, packet.data, error);
}

template <typename Packet>
bool DecodePortPacket(const std::vector<std::uint8_t> &payload, Packet *packet,
                      std::string *error)
{
    if (!packet) {
        SetError(error, "output packet is null");
        return false;
    }
    Reader reader(payload.data(), payload.size());
    std::uint32_t header_length{};
    std::uint32_t data_length{};
    if (!reader.ReadUnsigned(&packet->port) || !reader.ReadUnsigned(&header_length) ||
        !reader.ReadString(header_length, &packet->packet_header) ||
        !reader.ReadUnsigned(&data_length) ||
        !reader.ReadString(data_length, &packet->data) || !reader.done()) {
        SetError(error, "invalid port/header/data payload");
        return false;
    }
    return true;
}

}  // namespace

bool BuildMqFrame(const MqFrame &frame, std::size_t maximum_size,
                  std::vector<std::uint8_t> *bytes, std::string *error)
{
    if (!bytes) {
        SetError(error, "output frame is null");
        return false;
    }
    if (frame.stream_end != kMqStreamMore && !IsMqStreamEnd(frame.stream_end)) {
        SetError(error, "stream_end must be 0x00, 0x01, or 0xFF");
        return false;
    }
    if (frame.payload.size() > std::numeric_limits<std::uint32_t>::max() ||
        frame.payload.size() + kMqFrameHeaderSize > maximum_size) {
        SetError(error, "MQ frame exceeds message size limit");
        return false;
    }
    bytes->clear();
    bytes->reserve(kMqFrameHeaderSize + frame.payload.size());
    bytes->push_back(frame.message_type);
    AppendUnsigned(bytes, frame.message_seq);
    bytes->push_back(frame.stream_end);
    AppendUnsigned(bytes, static_cast<std::uint32_t>(frame.payload.size()));
    bytes->insert(bytes->end(), frame.payload.begin(), frame.payload.end());
    return true;
}

bool ParseMqFrame(const void *bytes, std::size_t length, MqFrame *frame,
                  std::string *error)
{
    if (!bytes || !frame || length < kMqFrameHeaderSize) {
        SetError(error, "MQ frame is shorter than the 10-byte header");
        return false;
    }
    const auto *raw = static_cast<const std::uint8_t *>(bytes);
    Reader reader(raw, length);
    std::uint32_t payload_length{};
    frame->message_type = raw[0];
    reader = Reader(raw + 1, length - 1);
    if (!reader.ReadUnsigned(&frame->message_seq) || reader.remaining() < 1) {
        SetError(error, "invalid MQ frame header");
        return false;
    }
    frame->stream_end = raw[5];
    Reader length_reader(raw + 6, length - 6);
    if (!length_reader.ReadUnsigned(&payload_length) ||
        payload_length != length - kMqFrameHeaderSize) {
        SetError(error, "MQ frame payload length does not match mq_receive length");
        return false;
    }
    if (frame->stream_end != kMqStreamMore && !IsMqStreamEnd(frame->stream_end)) {
        SetError(error, "invalid MQ stream_end marker");
        return false;
    }
    frame->payload.assign(raw + kMqFrameHeaderSize, raw + length);
    return true;
}

bool EncodeUplinkPacket(const UplinkPacket &packet, std::vector<std::uint8_t> *payload,
                        std::string *error)
{
    return EncodePortPacket(packet, payload, error);
}

bool DecodeUplinkPacket(const std::vector<std::uint8_t> &payload, UplinkPacket *packet,
                        std::string *error)
{
    return DecodePortPacket(payload, packet, error);
}

bool EncodeTimeSyncData(const TimeSyncData &time_sync,
                        std::vector<std::uint8_t> *payload)
{
    if (!payload) {
        return false;
    }
    payload->clear();
    AppendInt64(payload, time_sync.timestamp_s);
    return true;
}

bool DecodeTimeSyncData(const std::vector<std::uint8_t> &payload, TimeSyncData *time_sync,
                        std::string *error)
{
    Reader reader(payload.data(), payload.size());
    if (!time_sync || !reader.ReadInt64(&time_sync->timestamp_s) || !reader.done()) {
        SetError(error, "time-sync payload must contain one int64 Unix timestamp in seconds");
        return false;
    }
    return true;
}

bool EncodePredictionRequest(const OrbitPredictionRequest &request,
                             std::vector<std::uint8_t> *payload)
{
    if (!payload) {
        return false;
    }
    payload->clear();
    AppendInt64(payload, request.start_time_s);
    AppendUnsigned(payload, request.duration_s);
    return true;
}

bool DecodePredictionRequest(const std::vector<std::uint8_t> &payload,
                             OrbitPredictionRequest *request, std::string *error)
{
    Reader reader(payload.data(), payload.size());
    if (!request || !reader.ReadInt64(&request->start_time_s) ||
        !reader.ReadUnsigned(&request->duration_s) || !reader.done()) {
        SetError(error, "prediction payload must be int64 start_time_s + uint32 duration_s");
        return false;
    }
    return true;
}

bool EncodeStatusReply(const StatusReply &status, std::vector<std::uint8_t> *payload,
                       std::string *error)
{
    if (!payload) {
        SetError(error, "output payload is null");
        return false;
    }
    payload->clear();
    AppendInt64(payload, status.timestamp_ms);
    payload->push_back(static_cast<std::uint8_t>(status.status));
    return AppendString(payload, status.message, error);
}

bool DecodeStatusReply(const std::vector<std::uint8_t> &payload, StatusReply *status,
                       std::string *error)
{
    Reader reader(payload.data(), payload.size());
    std::uint8_t status_value{};
    std::uint32_t message_length{};
    if (!status || !reader.ReadInt64(&status->timestamp_ms) ||
        !reader.ReadUnsigned(&status_value) ||
        status_value > static_cast<std::uint8_t>(ServiceStatus::Error) ||
        !reader.ReadUnsigned(&message_length) ||
        !reader.ReadString(message_length, &status->message) || !reader.done()) {
        SetError(error, "invalid StatusReply payload");
        return false;
    }
    status->status = static_cast<ServiceStatus>(status_value);
    return true;
}

bool EncodeTelemetryPacket(const TelemetryPacket &packet,
                           std::vector<std::uint8_t> *payload, std::string *error)
{
    return EncodePortPacket(packet, payload, error);
}

bool DecodeTelemetryPacket(const std::vector<std::uint8_t> &payload,
                           TelemetryPacket *packet, std::string *error)
{
    return DecodePortPacket(payload, packet, error);
}

bool EncodeOrbitData(const OrbitData &data, std::vector<std::uint8_t> *payload,
                     std::string *error)
{
    if (!payload || data.points.size() >
                        std::numeric_limits<std::size_t>::max() / 56U) {
        SetError(error, "invalid OrbitData output or point count");
        return false;
    }
    payload->clear();
    payload->reserve(data.points.size() * 56U);
    for (const OrbitPoint &point : data.points) {
        AppendInt64(payload, point.timestamp_ms);
        AppendDouble(payload, point.x);
        AppendDouble(payload, point.y);
        AppendDouble(payload, point.z);
        AppendDouble(payload, point.vx);
        AppendDouble(payload, point.vy);
        AppendDouble(payload, point.vz);
    }
    return true;
}

bool DecodeOrbitData(const std::vector<std::uint8_t> &payload, OrbitData *data,
                     std::string *error)
{
    if (!data || payload.size() % 56U != 0) {
        SetError(error, "raw OrbitPoint payload length must be a multiple of 56 bytes");
        return false;
    }
    Reader reader(payload.data(), payload.size());
    const std::size_t point_count = payload.size() / 56U;
    data->points.clear();
    data->points.reserve(point_count);
    for (std::size_t index = 0; index < point_count; ++index) {
        OrbitPoint point;
        if (!reader.ReadInt64(&point.timestamp_ms) || !reader.ReadDouble(&point.x) ||
            !reader.ReadDouble(&point.y) || !reader.ReadDouble(&point.z) ||
            !reader.ReadDouble(&point.vx) || !reader.ReadDouble(&point.vy) ||
            !reader.ReadDouble(&point.vz)) {
            SetError(error, "truncated OrbitPoint payload");
            return false;
        }
        data->points.push_back(point);
    }
    if (!reader.done()) {
        SetError(error, "unexpected bytes after OrbitData payload");
        return false;
    }
    return true;
}

}  // namespace orbit_prediction

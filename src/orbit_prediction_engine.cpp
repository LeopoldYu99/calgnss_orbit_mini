#include "orbit_prediction_engine.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace orbit_prediction {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kWgs84A = 6378137.0;
constexpr double kWgs84F = 1.0 / 298.257223563;
constexpr double kWgs84E2 = kWgs84F * (2.0 - kWgs84F);
constexpr std::int64_t kMillisecondsPerDay = 86400000;
constexpr std::size_t kMaximumUplinkBytes = 1024 * 1024;
constexpr std::size_t kMaximumRtcmChunkBytes = 1024 * 1024;
constexpr std::size_t kMaximumPredictionPoints = 100000;
constexpr std::uint32_t kCrc24QPolynomial = 0x1864CFB;

struct CivilDate {
    int year{};
    int month{};
    int day{};
};

std::int64_t DaysFromCivil(int year, unsigned month, unsigned day)
{
    year -= month <= 2;
    const int era = (year >= 0 ? year : year - 399) / 400;
    const unsigned year_of_era = static_cast<unsigned>(year - era * 400);
    const unsigned adjusted_month = month > 2 ? month - 3 : month + 9;
    const unsigned day_of_year = (153 * adjusted_month + 2) / 5 + day - 1;
    const unsigned day_of_era = year_of_era * 365 + year_of_era / 4 - year_of_era / 100 +
                                day_of_year;
    return static_cast<std::int64_t>(era) * 146097 + static_cast<std::int64_t>(day_of_era) -
           719468;
}

CivilDate CivilFromDays(std::int64_t days)
{
    days += 719468;
    const std::int64_t era = (days >= 0 ? days : days - 146096) / 146097;
    const unsigned day_of_era = static_cast<unsigned>(days - era * 146097);
    const unsigned year_of_era =
        (day_of_era - day_of_era / 1460 + day_of_era / 36524 - day_of_era / 146096) / 365;
    int year = static_cast<int>(year_of_era) + static_cast<int>(era) * 400;
    const unsigned day_of_year = day_of_era -
                                 (365 * year_of_era + year_of_era / 4 - year_of_era / 100);
    const unsigned month_prime = (5 * day_of_year + 2) / 153;
    const unsigned day = day_of_year - (153 * month_prime + 2) / 5 + 1;
    const unsigned month = month_prime < 10 ? month_prime + 3 : month_prime - 9;
    year += month <= 2;
    return {year, static_cast<int>(month), static_cast<int>(day)};
}

std::int64_t FloorDiv(std::int64_t numerator, std::int64_t denominator)
{
    std::int64_t quotient = numerator / denominator;
    const std::int64_t remainder = numerator % denominator;
    if (remainder != 0 && ((remainder < 0) != (denominator < 0))) {
        --quotient;
    }
    return quotient;
}

cg_time_t MakeTime(std::int64_t timestamp_ms)
{
    const std::int64_t days = FloorDiv(timestamp_ms, kMillisecondsPerDay);
    const std::int64_t milliseconds_of_day = timestamp_ms - days * kMillisecondsPerDay;
    const CivilDate date = CivilFromDays(days);
    const int hour = static_cast<int>(milliseconds_of_day / 3600000);
    const int minute = static_cast<int>((milliseconds_of_day % 3600000) / 60000);
    const double second = static_cast<double>(milliseconds_of_day % 60000) / 1000.0;
    const double unix_seconds = static_cast<double>(timestamp_ms) / 1000.0;

    cg_time_t result{};
    result.year = date.year;
    result.month = date.month;
    result.day = date.day;
    result.hour = hour;
    result.minute = minute;
    result.second = second;
    result.unix_seconds = unix_seconds;
    result.jd_utc = unix_seconds / 86400.0 + 2440587.5;
    return result;
}

std::vector<std::string> Split(const std::string &text, char delimiter)
{
    std::vector<std::string> fields;
    std::size_t begin = 0;
    for (;;) {
        const std::size_t end = text.find(delimiter, begin);
        fields.emplace_back(text.substr(begin, end == std::string::npos ? end : end - begin));
        if (end == std::string::npos) {
            break;
        }
        begin = end + 1;
    }
    return fields;
}

bool EndsWith(const std::string &text, const char *suffix)
{
    const std::size_t length = std::char_traits<char>::length(suffix);
    return text.size() >= length && text.compare(text.size() - length, length, suffix) == 0;
}

bool ParseInt(const std::string &text, int *value)
{
    if (!value || text.empty()) {
        return false;
    }
    char *end = nullptr;
    const long parsed = std::strtol(text.c_str(), &end, 10);
    if (!end || *end != '\0' || parsed < std::numeric_limits<int>::min() ||
        parsed > std::numeric_limits<int>::max()) {
        return false;
    }
    *value = static_cast<int>(parsed);
    return true;
}

bool ParseDouble(const std::string &text, double *value)
{
    if (!value || text.empty()) {
        return false;
    }
    char *end = nullptr;
    const double parsed = std::strtod(text.c_str(), &end);
    if (!end || *end != '\0' || !std::isfinite(parsed)) {
        return false;
    }
    *value = parsed;
    return true;
}

int HexValue(char value)
{
    if (value >= '0' && value <= '9') {
        return value - '0';
    }
    value = static_cast<char>(std::toupper(static_cast<unsigned char>(value)));
    if (value >= 'A' && value <= 'F') {
        return value - 'A' + 10;
    }
    return -1;
}

bool NormalizeNmeaSentence(const std::string &line, std::string *sentence)
{
    const std::size_t dollar = line.find('$');
    if (dollar == std::string::npos || !sentence) {
        return false;
    }

    std::string value = line.substr(dollar);
    while (!value.empty() && (value.back() == '\r' || value.back() == '\n' ||
                              std::isspace(static_cast<unsigned char>(value.back())))) {
        value.pop_back();
    }
    const std::size_t checksum_marker = value.find('*');
    if (checksum_marker != std::string::npos) {
        if (checksum_marker + 2 >= value.size()) {
            return false;
        }
        const int high = HexValue(value[checksum_marker + 1]);
        const int low = HexValue(value[checksum_marker + 2]);
        if (high < 0 || low < 0) {
            return false;
        }
        unsigned char checksum = 0;
        for (std::size_t index = 1; index < checksum_marker; ++index) {
            checksum ^= static_cast<unsigned char>(value[index]);
        }
        if (checksum != static_cast<unsigned char>((high << 4) | low)) {
            return false;
        }
        value.resize(checksum_marker);
    }
    *sentence = std::move(value);
    return true;
}

bool ParseNmeaDate(const std::string &text, CivilDate *date)
{
    if (!date || text.size() != 6) {
        return false;
    }
    int day = 0;
    int month = 0;
    int year = 0;
    if (!ParseInt(text.substr(0, 2), &day) || !ParseInt(text.substr(2, 2), &month) ||
        !ParseInt(text.substr(4, 2), &year)) {
        return false;
    }
    year += year >= 80 ? 1900 : 2000;
    if (month < 1 || month > 12 || day < 1 || day > 31) {
        return false;
    }
    const std::int64_t days = DaysFromCivil(year, static_cast<unsigned>(month),
                                            static_cast<unsigned>(day));
    const CivilDate round_trip = CivilFromDays(days);
    if (round_trip.year != year || round_trip.month != month || round_trip.day != day) {
        return false;
    }
    *date = round_trip;
    return true;
}

bool ParseNmeaTimeOfDay(const std::string &text, std::int64_t *milliseconds)
{
    if (!milliseconds || text.size() < 6) {
        return false;
    }
    int hour = 0;
    int minute = 0;
    double second = 0.0;
    if (!ParseInt(text.substr(0, 2), &hour) || !ParseInt(text.substr(2, 2), &minute) ||
        !ParseDouble(text.substr(4), &second)) {
        return false;
    }
    if (hour < 0 || hour > 23 || minute < 0 || minute > 59 || second < 0.0 ||
        second >= 60.0) {
        return false;
    }
    *milliseconds = static_cast<std::int64_t>(hour) * 3600000 +
                    static_cast<std::int64_t>(minute) * 60000 +
                    static_cast<std::int64_t>(std::llround(second * 1000.0));
    return *milliseconds < kMillisecondsPerDay;
}

bool ParseNmeaCoordinate(const std::string &text, const std::string &hemisphere,
                         double *degrees)
{
    double raw = 0.0;
    if (!degrees || !ParseDouble(text, &raw) || hemisphere.size() != 1 || raw < 0.0) {
        return false;
    }
    const double whole_degrees = std::floor(raw / 100.0);
    const double minutes = raw - whole_degrees * 100.0;
    if (minutes < 0.0 || minutes >= 60.0) {
        return false;
    }
    double result = whole_degrees + minutes / 60.0;
    const char direction = static_cast<char>(std::toupper(
        static_cast<unsigned char>(hemisphere.front())));
    if (direction == 'S' || direction == 'W') {
        result = -result;
    } else if (direction != 'N' && direction != 'E') {
        return false;
    }
    if ((direction == 'N' || direction == 'S') && std::abs(result) > 90.0) {
        return false;
    }
    if ((direction == 'E' || direction == 'W') && std::abs(result) > 180.0) {
        return false;
    }
    *degrees = result;
    return true;
}

cg_vec3_t LlaToEcef(double latitude_degrees, double longitude_degrees, double height_m)
{
    const double latitude = latitude_degrees * kPi / 180.0;
    const double longitude = longitude_degrees * kPi / 180.0;
    const double sin_latitude = std::sin(latitude);
    const double cos_latitude = std::cos(latitude);
    const double radius = kWgs84A / std::sqrt(1.0 - kWgs84E2 * sin_latitude * sin_latitude);

    return {
        (radius + height_m) * cos_latitude * std::cos(longitude),
        (radius + height_m) * cos_latitude * std::sin(longitude),
        (radius * (1.0 - kWgs84E2) + height_m) * sin_latitude,
    };
}

std::int64_t NearestTimestampForTimeOfDay(
    const CivilDate &date, std::int64_t time_of_day_ms, std::int64_t reference_ms)
{
    const std::int64_t base =
        DaysFromCivil(date.year, static_cast<unsigned>(date.month), static_cast<unsigned>(date.day)) *
            kMillisecondsPerDay +
        time_of_day_ms;
    std::array<std::int64_t, 3> candidates{
        base - kMillisecondsPerDay,
        base,
        base + kMillisecondsPerDay,
    };
    return *std::min_element(candidates.begin(), candidates.end(), [reference_ms](auto left, auto right) {
        return std::llabs(left - reference_ms) < std::llabs(right - reference_ms);
    });
}

void SetError(std::string *error, const std::string &message)
{
    if (error) {
        *error = message;
    }
}

std::uint32_t CRC24Q(const std::uint8_t *data, std::size_t size)
{
    std::uint32_t crc = 0;
    for (std::size_t index = 0; index < size; ++index) {
        crc ^= static_cast<std::uint32_t>(data[index]) << 16;
        for (int bit = 0; bit < 8; ++bit) {
            crc <<= 1;
            if ((crc & 0x1000000U) != 0) {
                crc ^= kCrc24QPolynomial;
            }
        }
    }
    return crc & 0xFFFFFFU;
}

}  // namespace

OrbitPredictionEngine::OrbitPredictionEngine(EngineOptions options)
    : options_(options), observation_buffer_(options.observation_capacity)
{
    if (options_.observation_capacity < 2 || options_.fit_degree < 1 ||
        options_.fit_degree > CG_MAX_DEGREE || options_.stream_batch_size < 1 ||
        options_.rtcm_cache_bytes < 1029) {
        throw std::invalid_argument("invalid orbit prediction engine options");
    }
    std::string error;
    if (!RecreateContext(&error)) {
        throw std::runtime_error(error);
    }
}

OrbitPredictionEngine::~OrbitPredictionEngine()
{
    cg_context_destroy(context_);
}

bool OrbitPredictionEngine::RecreateContext(std::string *error)
{
    cg_context_destroy(context_);
    context_ = nullptr;
    cg_options_t core_options = cg_default_options();
    core_options.degree = options_.fit_degree;
    const int status = cg_context_create(&context_, observation_buffer_.data(),
                                         observation_buffer_.size(), &core_options);
    if (status != CG_OK) {
        SetError(error, std::string("cannot create orbit context: ") + cg_status_string(status));
        return false;
    }
    return true;
}

bool OrbitPredictionEngine::ReceiveUplinkData(
    std::uint32_t port, const std::string &packet_header, const std::string &data,
    std::string *error)
{
    (void)port;
    (void)packet_header;
    if (data.empty() || data.size() > kMaximumUplinkBytes) {
        SetError(error, "uplink data must contain 1 to 1048576 bytes");
        return false;
    }

    std::lock_guard<std::mutex> lock(mutex_);
    CivilDate current_date{nmea_year_, nmea_month_, nmea_day_};
    bool has_current_date = has_nmea_date_;
    bool date_is_from_rmc = has_nmea_date_;
    if (!has_current_date && has_time_sync_) {
        current_date = CivilFromDays(FloorDiv(time_sync_ms_, kMillisecondsPerDay));
        has_current_date = true;
    }

    std::vector<cg_observation_t> parsed_observations;
    std::istringstream input(data);
    std::string line;
    bool saw_supported_sentence = false;
    bool accepted_sentence = false;
    while (std::getline(input, line)) {
        std::string sentence;
        if (!NormalizeNmeaSentence(line, &sentence)) {
            continue;
        }
        const std::vector<std::string> fields = Split(sentence, ',');
        if (fields.empty()) {
            continue;
        }

        if (EndsWith(fields[0], "RMC")) {
            saw_supported_sentence = true;
            if (fields.size() < 10 || fields[2] != "A") {
                continue;
            }
            CivilDate date;
            if (!ParseNmeaDate(fields[9], &date)) {
                SetError(error, "invalid RMC date");
                return false;
            }
            current_date = date;
            has_current_date = true;
            date_is_from_rmc = true;
            accepted_sentence = true;
            continue;
        }

        if (!EndsWith(fields[0], "GGA")) {
            continue;
        }
        saw_supported_sentence = true;
        if (fields.size() < 12 || !has_current_date) {
            SetError(error, fields.size() < 12 ? "invalid GGA sentence" :
                                                "GGA requires RMC date or ReceiveTimeSync");
            return false;
        }

        int fix_quality = 0;
        std::int64_t time_of_day_ms = 0;
        double latitude = 0.0;
        double longitude = 0.0;
        double altitude_msl = 0.0;
        double geoid_separation = 0.0;
        if (!ParseInt(fields[6], &fix_quality) || fix_quality <= 0 ||
            !ParseNmeaTimeOfDay(fields[1], &time_of_day_ms) ||
            !ParseNmeaCoordinate(fields[2], fields[3], &latitude) ||
            !ParseNmeaCoordinate(fields[4], fields[5], &longitude) ||
            !ParseDouble(fields[9], &altitude_msl) ||
            (!fields[11].empty() && !ParseDouble(fields[11], &geoid_separation))) {
            SetError(error, "invalid GGA fix, time, coordinate, or height");
            return false;
        }

        std::int64_t timestamp_ms =
            DaysFromCivil(current_date.year, static_cast<unsigned>(current_date.month),
                          static_cast<unsigned>(current_date.day)) *
                kMillisecondsPerDay +
            time_of_day_ms;
        if (!date_is_from_rmc && has_time_sync_) {
            timestamp_ms = NearestTimestampForTimeOfDay(current_date, time_of_day_ms, time_sync_ms_);
            current_date = CivilFromDays(FloorDiv(timestamp_ms, kMillisecondsPerDay));
        }

        cg_observation_t observation{};
        observation.time_utc = MakeTime(timestamp_ms);
        observation.r_ecef_m = LlaToEcef(latitude, longitude, altitude_msl + geoid_separation);
        parsed_observations.push_back(observation);
        accepted_sentence = true;
    }

    if (!saw_supported_sentence || !accepted_sentence) {
        SetError(error, "uplink data contains no valid NMEA RMC or GGA data");
        return false;
    }

    double previous_time = has_latest_observation_
                               ? static_cast<double>(latest_observation_ms_) / 1000.0
                               : -std::numeric_limits<double>::infinity();
    for (const auto &observation : parsed_observations) {
        if (observation.time_utc.unix_seconds <= previous_time) {
            SetError(error, "NMEA observations must be strictly increasing in time");
            return false;
        }
        previous_time = observation.time_utc.unix_seconds;
    }

    for (const auto &observation : parsed_observations) {
        const int status = cg_context_push(context_, &observation);
        if (status != CG_OK) {
            SetError(error, std::string("cannot store observation: ") + cg_status_string(status));
            return false;
        }
    }

    has_nmea_date_ = date_is_from_rmc;
    nmea_year_ = current_date.year;
    nmea_month_ = current_date.month;
    nmea_day_ = current_date.day;
    if (!parsed_observations.empty()) {
        has_latest_observation_ = true;
        latest_observation_ms_ = static_cast<std::int64_t>(
            std::llround(parsed_observations.back().time_utc.unix_seconds * 1000.0));
    }
    if (error) {
        error->clear();
    }
    return true;
}

bool OrbitPredictionEngine::ReceiveTimeSync(std::int64_t timestamp_ms, std::string *error)
{
    if (timestamp_ms < 0 || timestamp_ms > 253402300799999LL) {
        SetError(error, "timestamp_ms is outside the supported UTC range");
        return false;
    }
    const cg_time_t time = MakeTime(timestamp_ms);
    if (time.year < 1970 || time.year > 9999) {
        SetError(error, "timestamp_ms is outside the supported UTC range");
        return false;
    }
    std::lock_guard<std::mutex> lock(mutex_);
    has_time_sync_ = true;
    time_sync_ms_ = timestamp_ms;
    if (error) {
        error->clear();
    }
    return true;
}

bool OrbitPredictionEngine::ReceiveRTCMData(const std::string &data, std::string *message)
{
    if (data.empty() || data.size() > kMaximumRtcmChunkBytes) {
        SetError(message, "RTCM data must contain 1 to 1048576 raw binary bytes");
        return false;
    }

    std::lock_guard<std::mutex> lock(mutex_);
    rtcm_input_buffer_.insert(rtcm_input_buffer_.end(), data.begin(), data.end());

    std::size_t cursor = 0;
    std::size_t accepted_frames = 0;
    std::size_t bad_crc_frames = 0;
    std::size_t discarded_bytes = 0;
    std::map<std::uint16_t, std::size_t> accepted_types;
    while (cursor < rtcm_input_buffer_.size()) {
        const auto preamble = std::find(rtcm_input_buffer_.begin() +
                                            static_cast<std::ptrdiff_t>(cursor),
                                        rtcm_input_buffer_.end(), 0xD3U);
        if (preamble == rtcm_input_buffer_.end()) {
            discarded_bytes += rtcm_input_buffer_.size() - cursor;
            cursor = rtcm_input_buffer_.size();
            break;
        }
        const std::size_t frame_begin =
            static_cast<std::size_t>(std::distance(rtcm_input_buffer_.begin(), preamble));
        discarded_bytes += frame_begin - cursor;
        cursor = frame_begin;
        if (rtcm_input_buffer_.size() - frame_begin < 3) {
            break;
        }
        if ((rtcm_input_buffer_[frame_begin + 1] & 0xFCU) != 0) {
            ++discarded_bytes;
            ++cursor;
            continue;
        }

        const std::size_t payload_size =
            (static_cast<std::size_t>(rtcm_input_buffer_[frame_begin + 1] & 0x03U) << 8) |
            rtcm_input_buffer_[frame_begin + 2];
        const std::size_t frame_size = payload_size + 6;
        if (rtcm_input_buffer_.size() - frame_begin < frame_size) {
            break;
        }

        const std::size_t crc_offset = frame_begin + frame_size - 3;
        const std::uint32_t expected_crc =
            (static_cast<std::uint32_t>(rtcm_input_buffer_[crc_offset]) << 16) |
            (static_cast<std::uint32_t>(rtcm_input_buffer_[crc_offset + 1]) << 8) |
            rtcm_input_buffer_[crc_offset + 2];
        const std::uint32_t actual_crc =
            CRC24Q(rtcm_input_buffer_.data() + frame_begin, frame_size - 3);
        if (actual_crc != expected_crc) {
            ++bad_crc_frames;
            ++discarded_bytes;
            ++cursor;
            continue;
        }

        std::vector<std::uint8_t> frame(
            rtcm_input_buffer_.begin() + static_cast<std::ptrdiff_t>(frame_begin),
            rtcm_input_buffer_.begin() + static_cast<std::ptrdiff_t>(frame_begin + frame_size));
        if (payload_size >= 2) {
            const std::uint16_t message_type = static_cast<std::uint16_t>(
                (static_cast<std::uint16_t>(frame[3]) << 4) | (frame[4] >> 4));
            ++accepted_types[message_type];
        }
        rtcm_frame_bytes_ += frame.size();
        rtcm_frames_.push_back(std::move(frame));
        while (rtcm_frame_bytes_ > options_.rtcm_cache_bytes && !rtcm_frames_.empty()) {
            rtcm_frame_bytes_ -= rtcm_frames_.front().size();
            rtcm_frames_.pop_front();
        }
        ++accepted_frames;
        ++rtcm_frames_received_;
        cursor = frame_begin + frame_size;
    }

    if (cursor > 0) {
        rtcm_input_buffer_.erase(
            rtcm_input_buffer_.begin(),
            rtcm_input_buffer_.begin() + static_cast<std::ptrdiff_t>(cursor));
    }

    std::ostringstream result;
    if (accepted_frames > 0) {
        result << "accepted " << accepted_frames << " RTCM3 frame(s)";
        if (!accepted_types.empty()) {
            result << "; message types=";
            bool first = true;
            for (const auto &entry : accepted_types) {
                if (!first) {
                    result << ',';
                }
                first = false;
                result << entry.first << 'x' << entry.second;
            }
        }
    } else if (!rtcm_input_buffer_.empty() && bad_crc_frames == 0) {
        result << "buffered " << rtcm_input_buffer_.size() << " byte(s) of a partial RTCM3 frame";
    } else {
        result << "no valid RTCM3 frame found";
    }
    if (bad_crc_frames > 0) {
        result << "; CRC failures=" << bad_crc_frames;
    }
    if (discarded_bytes > 0) {
        result << "; discarded bytes=" << discarded_bytes;
    }
    SetError(message, result.str());
    return accepted_frames > 0 || (!rtcm_input_buffer_.empty() && bad_crc_frames == 0);
}

void OrbitPredictionEngine::Stop()
{
    stopped_.store(true, std::memory_order_release);
}

bool OrbitPredictionEngine::Reset(std::string *error)
{
    std::lock_guard<std::mutex> lock(mutex_);
    cg_context_reset(context_);
    has_time_sync_ = false;
    time_sync_ms_ = 0;
    has_nmea_date_ = false;
    nmea_year_ = 0;
    nmea_month_ = 0;
    nmea_day_ = 0;
    has_latest_observation_ = false;
    latest_observation_ms_ = 0;
    rtcm_input_buffer_.clear();
    rtcm_frames_.clear();
    rtcm_frame_bytes_ = 0;
    rtcm_frames_received_ = 0;
    stopped_.store(false, std::memory_order_release);
    if (error) {
        error->clear();
    }
    return true;
}

bool OrbitPredictionEngine::PredictOrbit(
    std::int64_t start_time_s, std::uint32_t duration_s, double step_s,
    const std::function<bool(const std::vector<PredictionPoint> &)> &on_batch,
    const std::function<bool()> &is_cancelled, std::string *error)
{
    if (!on_batch) {
        SetError(error, "prediction callback is required");
        return false;
    }
    if (start_time_s < 0 || start_time_s > 253402300799LL || duration_s == 0 ||
        duration_s > 3600 || !std::isfinite(step_s) || step_s < 0.001 ||
        step_s > static_cast<double>(duration_s)) {
        SetError(error,
                 "request requires a valid start_time_s, duration_s in [1,3600], and step_s "
                 "in [0.001,duration_s]");
        return false;
    }
    const std::size_t point_count = static_cast<std::size_t>(
                                        std::floor(static_cast<double>(duration_s) / step_s + 1.0e-12)) +
                                    1;
    if (point_count > kMaximumPredictionPoints) {
        SetError(error, "prediction request exceeds the 100000 point limit");
        return false;
    }
    if (stopped_.load(std::memory_order_acquire)) {
        SetError(error, "prediction is stopped; call Reset before PredictOrbit");
        return false;
    }

    {
        std::lock_guard<std::mutex> lock(mutex_);
        if (!context_ || cg_context_count(context_) < 2 || !has_latest_observation_) {
            SetError(error, "at least two observations are required");
            return false;
        }
    }

    std::vector<PredictionPoint> batch;
    batch.reserve(options_.stream_batch_size);
    const std::int64_t prediction_epoch_ms = start_time_s * 1000;
    for (std::size_t point_index = 0; point_index < point_count; ++point_index) {
        if (stopped_.load(std::memory_order_acquire) || (is_cancelled && is_cancelled())) {
            break;
        }

        const std::int64_t query_timestamp_ms =
            prediction_epoch_ms + static_cast<std::int64_t>(
                                      std::llround(static_cast<double>(point_index) * step_s * 1000.0));
        const cg_time_t query_time = MakeTime(query_timestamp_ms);
        cg_state_t state{};
        int status = CG_OK;
        {
            std::lock_guard<std::mutex> lock(mutex_);
            status = cg_context_query_state(context_, &query_time, &state);
        }
        if (status != CG_OK) {
            SetError(error, std::string("orbit prediction failed: ") + cg_status_string(status));
            return false;
        }

        batch.push_back({
            query_timestamp_ms,
            state.r_j2000_m.x,
            state.r_j2000_m.y,
            state.r_j2000_m.z,
            state.v_j2000_mps.x,
            state.v_j2000_mps.y,
            state.v_j2000_mps.z,
        });
        if (batch.size() >= options_.stream_batch_size) {
            if (!on_batch(batch)) {
                return true;
            }
            batch.clear();
        }
    }

    if (!batch.empty() && !stopped_.load(std::memory_order_acquire) &&
        !(is_cancelled && is_cancelled())) {
        if (!on_batch(batch)) {
            return true;
        }
    }
    if (error) {
        error->clear();
    }
    return true;
}

std::size_t OrbitPredictionEngine::ObservationCount() const
{
    std::lock_guard<std::mutex> lock(mutex_);
    return context_ ? cg_context_count(context_) : 0;
}

std::uint64_t OrbitPredictionEngine::RTCMFrameCount() const
{
    std::lock_guard<std::mutex> lock(mutex_);
    return rtcm_frames_received_;
}

bool OrbitPredictionEngine::IsStopped() const
{
    return stopped_.load(std::memory_order_acquire);
}

}  // namespace orbit_prediction

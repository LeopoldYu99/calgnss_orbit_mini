#include "heartbeat_mq_publisher.h"
#include "orbit_prediction_service.h"
#include "posix_mq_server.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <csignal>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>

namespace {

std::atomic<bool> shutdown_requested{false};

void HandleSignal(int signal)
{
    if (signal == SIGINT || signal == SIGTERM) {
        shutdown_requested.store(true, std::memory_order_release);
    }
}

std::string EnvironmentString(const char *name, const char *fallback)
{
    const char *value = std::getenv(name);
    return value && *value ? value : fallback;
}

long long EnvironmentInteger(const char *name, long long fallback, long long minimum,
                             long long maximum)
{
    const char *value = std::getenv(name);
    if (!value || !*value) {
        return fallback;
    }
    char *end = nullptr;
    const long long parsed = std::strtoll(value, &end, 10);
    if (!end || *end != '\0' || parsed < minimum || parsed > maximum) {
        throw std::invalid_argument(std::string("invalid environment variable ") + name);
    }
    return parsed;
}

orbit_prediction::EngineOptions LoadOptions()
{
    orbit_prediction::EngineOptions options;
    options.observation_capacity = static_cast<std::size_t>(EnvironmentInteger(
        "ORBIT_OBSERVATION_CAPACITY", 600, 2, 100000));
    options.fit_degree = static_cast<int>(EnvironmentInteger(
        "ORBIT_FIT_DEGREE", 10, 1, CG_MAX_DEGREE));
    options.stream_batch_size = static_cast<std::size_t>(EnvironmentInteger(
        "ORBIT_STREAM_BATCH_SIZE", 100, 1, 10000));
    options.rtcm_cache_bytes = static_cast<std::size_t>(EnvironmentInteger(
        "ORBIT_RTCM_CACHE_BYTES", 4 * 1024 * 1024, 1029, 256 * 1024 * 1024));
    return options;
}

bool EnvironmentBoolean(const char *name, bool fallback)
{
    const std::string value = EnvironmentString(name, fallback ? "1" : "0");
    if (value == "1" || value == "true" || value == "TRUE" || value == "yes") {
        return true;
    }
    if (value == "0" || value == "false" || value == "FALSE" || value == "no") {
        return false;
    }
    throw std::invalid_argument(std::string("invalid environment variable ") + name);
}

orbit_prediction::PosixMqOptions LoadMqOptions()
{
    orbit_prediction::PosixMqOptions options;
    options.request_queue =
        EnvironmentString("ORBIT_MQ_REQUEST_QUEUE", "/csm_main_to_end0");
    options.response_queue =
        EnvironmentString("ORBIT_MQ_RESPONSE_QUEUE", "/csm_end0_to_main");
    options.max_messages = static_cast<long>(
        EnvironmentInteger("ORBIT_MQ_MAX_MESSAGES", 10, 1, 1024));
    options.message_size = static_cast<long>(
        EnvironmentInteger("ORBIT_MQ_MESSAGE_SIZE", 2048, 256, 1024 * 1024));
    options.create_queues = EnvironmentBoolean("ORBIT_MQ_CREATE", false);
    options.console_logging = EnvironmentBoolean("ORBIT_CONSOLE_LOG", false);
    return options;
}

std::uint32_t InitialHeartbeatSequence()
{
    const auto ticks = static_cast<std::uint32_t>(
        std::chrono::steady_clock::now().time_since_epoch().count());
    return ticks;
}

void RunHeartbeat(orbit_prediction::OrbitPredictionHandler &handler,
                  const orbit_prediction::HeartbeatMqPublisher &publisher,
                  std::chrono::milliseconds interval, bool console_logging)
{
    std::uint32_t message_seq = InitialHeartbeatSequence();
    bool first_attempt = true;
    bool previously_available = false;
    while (!shutdown_requested.load(std::memory_order_acquire)) {
        std::string error;
        const orbit_prediction::StatusReply status = handler.GetStatus();
        const std::uint32_t current_sequence = message_seq++;
        const bool available = publisher.Publish(status, current_sequence, &error);
        if (available && console_logging) {
            std::cout << "Sent MqRecvType: type="
                      << static_cast<unsigned>(orbit_prediction::MqRecvType::Heartbeat)
                      << '(' << orbit_prediction::MqRecvTypeName(
                                      orbit_prediction::MqRecvType::Heartbeat) << ')'
                      << " seq=" << current_sequence
                      << " timestamp_ms=" << status.timestamp_ms
                      << " status=" << static_cast<unsigned>(status.status)
                      << " message=\"" << status.message << "\"\n";
        } else if (available && (first_attempt || !previously_available)) {
            std::cout << "Lifecycle heartbeat MQ is available\n";
        } else if (!available && (first_attempt || previously_available)) {
            std::cerr << "Lifecycle heartbeat is not being queued: " << error << '\n';
        }
        first_attempt = false;
        previously_available = available;

        const auto wake_at = std::chrono::steady_clock::now() + interval;
        while (!shutdown_requested.load(std::memory_order_acquire) &&
               std::chrono::steady_clock::now() < wake_at) {
            const auto remaining = wake_at - std::chrono::steady_clock::now();
            std::this_thread::sleep_for(
                std::min(remaining, std::chrono::steady_clock::duration(
                                        std::chrono::milliseconds(100))));
        }
    }
}

}  // namespace

int main()
{
    try {
        std::cout << std::unitbuf;
        std::cerr << std::unitbuf;
        std::signal(SIGINT, HandleSignal);
        std::signal(SIGTERM, HandleSignal);

        const orbit_prediction::PosixMqOptions mq_options = LoadMqOptions();
        orbit_prediction::OrbitPredictionHandler handler(LoadOptions());
        orbit_prediction::PosixMqServer server(handler, mq_options);
        const bool heartbeat_enabled =
            EnvironmentBoolean("ORBIT_HEARTBEAT_ENABLED", true);
        const auto heartbeat_interval = std::chrono::milliseconds(EnvironmentInteger(
            "ORBIT_HEARTBEAT_INTERVAL_MS", 1000, 100, 60000));
        const orbit_prediction::HeartbeatMqPublisher heartbeat_publisher(
            EnvironmentString("ORBIT_LIFECYCLE_REQUEST_QUEUE", "/csm_end0_to_main"),
            static_cast<std::size_t>(mq_options.message_size));
        std::thread heartbeat_thread;
        if (heartbeat_enabled) {
            heartbeat_thread = std::thread(RunHeartbeat, std::ref(handler),
                                           std::cref(heartbeat_publisher),
                                           heartbeat_interval,
                                           mq_options.console_logging);
        }
        std::cout << "Orbit prediction MQ service receiving on " << mq_options.request_queue
                  << " and sending best-effort responses on " << mq_options.response_queue
                  << '\n';
        try {
            server.Run(shutdown_requested);
        } catch (...) {
            shutdown_requested.store(true, std::memory_order_release);
            if (heartbeat_thread.joinable()) {
                heartbeat_thread.join();
            }
            throw;
        }
        shutdown_requested.store(true, std::memory_order_release);
        if (heartbeat_thread.joinable()) {
            heartbeat_thread.join();
        }
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "orbit_prediction_mq_server: " << error.what() << '\n';
        return 1;
    }
}

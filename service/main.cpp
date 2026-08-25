#include "orbit_prediction_service.h"

#include <grpcpp/grpcpp.h>
#include <grpcpp/health_check_service_interface.h>

#include <atomic>
#include <chrono>
#include <csignal>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
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
        "ORBIT_STREAM_BATCH_SIZE", 10, 1, 10000));
    options.rtcm_cache_bytes = static_cast<std::size_t>(EnvironmentInteger(
        "ORBIT_RTCM_CACHE_BYTES", 4 * 1024 * 1024, 1029, 256 * 1024 * 1024));
    return options;
}

}  // namespace

int main()
{
    try {
        const std::string address =
            EnvironmentString("ORBIT_GRPC_ADDRESS", "192.168.104.100:50051");
        orbit_prediction::OrbitPredictionServiceImpl service(LoadOptions());

        grpc::EnableDefaultHealthCheckService(true);
        grpc::ServerBuilder builder;
        int selected_port = 0;
        builder.SetMaxReceiveMessageSize(2 * 1024 * 1024);
        builder.AddListeningPort(address, grpc::InsecureServerCredentials(), &selected_port);
        builder.RegisterService(&service);
        std::unique_ptr<grpc::Server> server = builder.BuildAndStart();
        if (!server || selected_port == 0) {
            std::cerr << "Cannot listen on " << address
                      << ". Ensure this Linux host owns the configured IP address.\n";
            return 1;
        }

        std::signal(SIGINT, HandleSignal);
        std::signal(SIGTERM, HandleSignal);
        std::thread shutdown_thread([&server]() {
            while (!shutdown_requested.load(std::memory_order_acquire)) {
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
            }
            server->Shutdown(std::chrono::system_clock::now() + std::chrono::seconds(5));
        });

        std::cout << "orbit_prediction.OrbitPredictionService listening on " << address << '\n';
        server->Wait();
        shutdown_requested.store(true, std::memory_order_release);
        shutdown_thread.join();
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "orbit_prediction_server: " << error.what() << '\n';
        return 1;
    }
}

#ifndef ORBIT_PREDICTION_POSIX_MQ_SERVER_H
#define ORBIT_PREDICTION_POSIX_MQ_SERVER_H

#include "orbit_prediction_service.h"

#include <atomic>
#include <cstddef>
#include <mqueue.h>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace orbit_prediction {

struct PosixMqOptions {
    std::string request_queue{"/csm_main_to_end0"};
    std::string response_queue{"/csm_end0_to_main"};
    long max_messages{10};
    long message_size{2048};
    bool create_queues{false};
    bool console_logging{false};
};

class PosixMqServer final {
public:
    PosixMqServer(OrbitPredictionHandler &handler, PosixMqOptions options = {});
    ~PosixMqServer();

    PosixMqServer(const PosixMqServer &) = delete;
    PosixMqServer &operator=(const PosixMqServer &) = delete;

    void Run(const std::atomic<bool> &shutdown_requested);

private:
    void Dispatch(const MqFrame &request);
    void StartPrediction(const MqFrame &request,
                         const OrbitPredictionRequest &prediction_request);
    void StopAndJoinPrediction();
    bool SendFrameBestEffort(const MqFrame &frame);
    void SendStatus(std::uint32_t sequence);
    void SendPredictionBatch(std::uint32_t *sequence,
                             const std::vector<PredictionPoint> &batch);
    void LogRejected(const MqFrame &request, const std::string &message) const;

    OrbitPredictionHandler &handler_;
    PosixMqOptions options_;
    mqd_t request_queue_{static_cast<mqd_t>(-1)};
    mqd_t response_queue_{static_cast<mqd_t>(-1)};
    std::size_t request_message_size_{0};
    std::size_t response_message_size_{0};
    std::mutex send_mutex_;
    std::mutex prediction_thread_mutex_;
    std::thread prediction_thread_;
};

}  // namespace orbit_prediction

#endif

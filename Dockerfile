FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        ninja-build \
        pkg-config \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /src
COPY . .
RUN cmake -S . -B build-container -G Ninja \
        -DCMAKE_BUILD_TYPE=Release \
        -DCALGNSS_BUILD_MQ=ON \
        -DCALGNSS_BUILD_TESTS=ON \
    && cmake --build build-container \
    && ctest --test-dir build-container --output-on-failure

ENV ORBIT_MQ_REQUEST_QUEUE=/orbit_prediction_requests \
    ORBIT_MQ_RESPONSE_QUEUE=/orbit_prediction_responses \
    ORBIT_MQ_MESSAGE_SIZE=8192
ENTRYPOINT ["/src/build-container/orbit_prediction_server"]

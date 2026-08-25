FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        ninja-build \
        pkg-config \
        libprotobuf-dev \
        protobuf-compiler \
        libgrpc++-dev \
        protobuf-compiler-grpc \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /src
COPY . .
RUN cmake -S . -B build-container -G Ninja \
        -DCMAKE_BUILD_TYPE=Release \
        -DCALGNSS_BUILD_GRPC=ON \
        -DCALGNSS_BUILD_TESTS=ON \
    && cmake --build build-container \
    && ctest --test-dir build-container --output-on-failure

ENV ORBIT_GRPC_ADDRESS=0.0.0.0:50051
EXPOSE 50051
ENTRYPOINT ["/src/build-container/orbit_prediction_server"]

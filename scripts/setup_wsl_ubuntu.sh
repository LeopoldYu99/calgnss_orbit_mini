#!/usr/bin/env bash
set -euo pipefail

sudo apt-get update
sudo apt-get install -y \
    build-essential \
    cmake \
    ninja-build \
    pkg-config \
    libprotobuf-dev \
    protobuf-compiler \
    libgrpc++-dev \
    protobuf-compiler-grpc

cmake -S . -B build-wsl -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCALGNSS_BUILD_GRPC=ON \
    -DCALGNSS_BUILD_TESTS=ON
cmake --build build-wsl
ctest --test-dir build-wsl --output-on-failure

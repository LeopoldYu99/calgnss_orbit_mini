# Standalone consumer: compiles only the C caller and links installed libod.a.
cmake_minimum_required(VERSION 3.16)
project(od_prebuilt_consumer LANGUAGES C)
set(CMAKE_C_STANDARD 99)
find_package(od CONFIG REQUIRED)
add_executable(od_prebuilt_replay ../od_replay.c)
target_link_libraries(od_prebuilt_replay PRIVATE od::orbit_determination)

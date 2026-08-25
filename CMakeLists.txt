cmake_minimum_required(VERSION 3.16)

project(calgnss_orbit_mini LANGUAGES C CXX)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(CALGNSS_BUILD_GRPC "Build the OrbitPredictionService gRPC microservice" ON)
option(CALGNSS_BUILD_TESTS "Build project tests" ON)

set(CALGNSS_ORBIT_MINI_SOURCES
    src/calgnss_orbit_mini.c
)

add_library(calgnss_orbit_mini SHARED
    ${CALGNSS_ORBIT_MINI_SOURCES}
)

add_library(calgnss_orbit_mini_static STATIC
    ${CALGNSS_ORBIT_MINI_SOURCES}
)

set_target_properties(calgnss_orbit_mini PROPERTIES
    WINDOWS_EXPORT_ALL_SYMBOLS ON
)

set_target_properties(calgnss_orbit_mini_static PROPERTIES
    OUTPUT_NAME calgnss_orbit_mini_static
    POSITION_INDEPENDENT_CODE ON
)

target_include_directories(calgnss_orbit_mini PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_include_directories(calgnss_orbit_mini_static PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

if (UNIX)
    target_link_libraries(calgnss_orbit_mini PUBLIC m)
    target_link_libraries(calgnss_orbit_mini_static PUBLIC m)
endif()

find_package(Threads REQUIRED)

add_library(orbit_prediction_engine STATIC
    src/orbit_prediction_engine.cpp
)
target_include_directories(orbit_prediction_engine PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(orbit_prediction_engine PUBLIC calgnss_orbit_mini_static Threads::Threads)
if (MSVC)
    target_compile_options(orbit_prediction_engine PRIVATE /W4 /permissive-)
else()
    target_compile_options(orbit_prediction_engine PRIVATE -Wall -Wextra -Wpedantic)
endif()

if (CALGNSS_BUILD_TESTS)
    enable_testing()
    add_executable(orbit_prediction_engine_test tests/orbit_prediction_engine_test.cpp)
    target_link_libraries(orbit_prediction_engine_test PRIVATE orbit_prediction_engine)
    add_test(NAME orbit_prediction_engine_test COMMAND orbit_prediction_engine_test)
endif()

if (CALGNSS_BUILD_GRPC)
    find_package(Protobuf REQUIRED)
    find_package(gRPC CONFIG REQUIRED)

    set(ORBIT_PROTO_FILE ${CMAKE_CURRENT_SOURCE_DIR}/proto/orbit_prediction.proto)
    set(ORBIT_GENERATED_DIR ${CMAKE_CURRENT_BINARY_DIR}/generated)
    set(ORBIT_PROTO_PATH_ARGS --proto_path=${CMAKE_CURRENT_SOURCE_DIR}/proto)
    foreach(PROTOBUF_INCLUDE_DIR IN LISTS Protobuf_INCLUDE_DIRS)
        list(APPEND ORBIT_PROTO_PATH_ARGS --proto_path=${PROTOBUF_INCLUDE_DIR})
    endforeach()
    file(MAKE_DIRECTORY ${ORBIT_GENERATED_DIR})

    if (CMAKE_CROSSCOMPILING)
        find_program(ORBIT_PROTOC_EXECUTABLE protoc REQUIRED NO_CMAKE_FIND_ROOT_PATH)
        find_program(ORBIT_GRPC_CPP_PLUGIN_EXECUTABLE grpc_cpp_plugin REQUIRED NO_CMAKE_FIND_ROOT_PATH)
        set(ORBIT_PROTO_TOOL_DEPENDENCIES)
    else()
        set(ORBIT_PROTOC_EXECUTABLE protobuf::protoc)
        set(ORBIT_GRPC_CPP_PLUGIN_EXECUTABLE $<TARGET_FILE:gRPC::grpc_cpp_plugin>)
        set(ORBIT_PROTO_TOOL_DEPENDENCIES protobuf::protoc gRPC::grpc_cpp_plugin)
    endif()

    add_custom_command(
        OUTPUT
            ${ORBIT_GENERATED_DIR}/orbit_prediction.pb.cc
            ${ORBIT_GENERATED_DIR}/orbit_prediction.pb.h
            ${ORBIT_GENERATED_DIR}/orbit_prediction.grpc.pb.cc
            ${ORBIT_GENERATED_DIR}/orbit_prediction.grpc.pb.h
        COMMAND ${ORBIT_PROTOC_EXECUTABLE}
        ARGS
            ${ORBIT_PROTO_PATH_ARGS}
            --cpp_out=${ORBIT_GENERATED_DIR}
            --grpc_out=${ORBIT_GENERATED_DIR}
            --plugin=protoc-gen-grpc=${ORBIT_GRPC_CPP_PLUGIN_EXECUTABLE}
            ${ORBIT_PROTO_FILE}
        DEPENDS
            ${ORBIT_PROTO_FILE}
            ${ORBIT_PROTO_TOOL_DEPENDENCIES}
        VERBATIM
    )

    add_library(orbit_prediction_proto STATIC
        ${ORBIT_GENERATED_DIR}/orbit_prediction.pb.cc
        ${ORBIT_GENERATED_DIR}/orbit_prediction.grpc.pb.cc
    )
    target_include_directories(orbit_prediction_proto PUBLIC ${ORBIT_GENERATED_DIR})
    target_link_libraries(orbit_prediction_proto PUBLIC protobuf::libprotobuf gRPC::grpc++)

    add_library(orbit_prediction_grpc_service STATIC
        service/orbit_prediction_service.cpp
    )
    target_include_directories(orbit_prediction_grpc_service PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/service
        ${ORBIT_GENERATED_DIR}
    )
    target_link_libraries(orbit_prediction_grpc_service PUBLIC
        orbit_prediction_engine
        orbit_prediction_proto
        gRPC::grpc++
    )
    if (MSVC)
        target_compile_options(orbit_prediction_grpc_service PRIVATE /W4 /permissive-)
    else()
        target_compile_options(orbit_prediction_grpc_service PRIVATE -Wall -Wextra -Wpedantic)
    endif()

    add_executable(orbit_prediction_server service/main.cpp)
    target_link_libraries(orbit_prediction_server PRIVATE orbit_prediction_grpc_service)
    if (MSVC)
        target_compile_options(orbit_prediction_server PRIVATE /W4 /permissive-)
    else()
        target_compile_options(orbit_prediction_server PRIVATE -Wall -Wextra -Wpedantic)
    endif()

    install(TARGETS orbit_prediction_server RUNTIME DESTINATION bin)

endif()

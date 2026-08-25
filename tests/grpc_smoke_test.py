#!/usr/bin/env python3
"""End-to-end smoke test for OrbitPredictionService."""

import argparse
import sys

import grpc
from google.protobuf import empty_pb2

import orbit_prediction_pb2 as messages
import orbit_prediction_pb2_grpc as services


def require_reply(reply, operation):
    if not reply.success:
        raise RuntimeError(f"{operation} failed: {reply.message}")


def test_nmea_arc():
    lines = ["$GNRMC,040403.00,A,1724.2594222,N,15304.2955483,E,0.0,0.0,230626,,,A"]
    for index in range(12):
        second = 3 + index
        longitude = 15304.2955483 + index * 4.2
        lines.append(
            f"$GNGGA,0404{second:02d}.00,1724.2594222,N,{longitude:.7f},"
            "E,1,14,0.69,271441.347,M,0.989,M,,"
        )
    return ("\n".join(lines) + "\n").encode("ascii")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--address", default="127.0.0.1:50051")
    args = parser.parse_args()

    channel = grpc.insecure_channel(args.address)
    grpc.channel_ready_future(channel).result(timeout=20)
    stub = services.OrbitPredictionServiceStub(channel)

    require_reply(stub.Reset(empty_pb2.Empty(), timeout=10), "Reset")
    require_reply(
        stub.ReceiveTimeSync(messages.TimeSyncData(timestamp_ms=1782187443000), timeout=10),
        "ReceiveTimeSync",
    )
    require_reply(
        stub.ReceiveUplinkData(messages.UplinkPacket(port=1, data=test_nmea_arc()), timeout=10),
        "ReceiveUplinkData",
    )

    rtcm = bytes.fromhex("d300044350000044fe2e")
    require_reply(stub.ReceiveRTCMData(messages.RTCMData(data=rtcm[:5]), timeout=10),
                  "ReceiveRTCMData(partial)")
    require_reply(stub.ReceiveRTCMData(messages.RTCMData(data=rtcm[5:]), timeout=10),
                  "ReceiveRTCMData(complete)")

    batches = stub.PredictOrbit(
        messages.OrbitPredictionRequest(
            start_time_s=1782187454, duration_s=2, step_s=0.5
        ),
        timeout=20,
    )
    points = [point for batch in batches for point in batch.points]
    if len(points) != 5:
        raise RuntimeError(f"PredictOrbit returned {len(points)} points instead of 5")
    if [points[0].timestamp_ms, points[1].timestamp_ms] != [1782187454000, 1782187454500]:
        raise RuntimeError("PredictOrbit returned unexpected timestamps")

    require_reply(stub.Stop(empty_pb2.Empty(), timeout=10), "Stop")
    channel.close()
    print(f"gRPC smoke test passed against {args.address}: all 6 RPCs, 5 orbit points")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:  # Keep CI output concise and actionable.
        print(f"gRPC smoke test failed: {error}", file=sys.stderr)
        sys.exit(1)

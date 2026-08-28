"""Capture raw orbit-prediction MQ frames and decode type-4 OrbitPoint payloads.

The queue is opened for non-blocking reads. Every received message is removed
from the queue, preserved byte-for-byte, and described in frames.csv. Type-4
payloads are decoded as a packed sequence of ``int64_t + 6 * double`` records.
"""

from __future__ import annotations

import argparse
import csv
import ctypes
import errno
import os
from pathlib import Path
import struct
import time


FRAME_HEADER = struct.Struct("<BIBI")
ORBIT_POINT = struct.Struct("<q6d")
TYPE_RECEIVE_TELEMETRY_DATA = 4


class MqAttr(ctypes.Structure):
    _fields_ = [
        ("mq_flags", ctypes.c_long),
        ("mq_maxmsg", ctypes.c_long),
        ("mq_msgsize", ctypes.c_long),
        ("mq_curmsgs", ctypes.c_long),
    ]


class MqReader:
    def __init__(self, queue_name: str) -> None:
        self.queue_name = queue_name
        self.libc = ctypes.CDLL(None, use_errno=True)
        self.libc.mq_open.argtypes = [ctypes.c_char_p, ctypes.c_int]
        self.libc.mq_open.restype = ctypes.c_int
        self.libc.mq_getattr.argtypes = [ctypes.c_int, ctypes.POINTER(MqAttr)]
        self.libc.mq_getattr.restype = ctypes.c_int
        self.libc.mq_receive.argtypes = [
            ctypes.c_int,
            ctypes.c_void_p,
            ctypes.c_size_t,
            ctypes.POINTER(ctypes.c_uint),
        ]
        self.libc.mq_receive.restype = ctypes.c_ssize_t
        self.libc.mq_close.argtypes = [ctypes.c_int]
        self.libc.mq_close.restype = ctypes.c_int

        self.descriptor = self.libc.mq_open(
            queue_name.encode(), os.O_RDONLY | os.O_NONBLOCK
        )
        if self.descriptor < 0:
            code = ctypes.get_errno()
            raise OSError(code, os.strerror(code), queue_name)

        attributes = MqAttr()
        if self.libc.mq_getattr(self.descriptor, ctypes.byref(attributes)) != 0:
            code = ctypes.get_errno()
            self.close()
            raise OSError(code, os.strerror(code), queue_name)
        self.message_size = attributes.mq_msgsize

    def close(self) -> None:
        if self.descriptor >= 0:
            self.libc.mq_close(self.descriptor)
            self.descriptor = -1

    def receive(self) -> tuple[bytes, int] | None:
        buffer = ctypes.create_string_buffer(self.message_size)
        priority = ctypes.c_uint()
        length = self.libc.mq_receive(
            self.descriptor, buffer, len(buffer), ctypes.byref(priority)
        )
        if length >= 0:
            return buffer.raw[:length], priority.value
        code = ctypes.get_errno()
        if code == errno.EAGAIN:
            return None
        raise OSError(code, os.strerror(code), self.queue_name)


def decode_frame(frame: bytes) -> tuple[int, int, int, bytes]:
    if len(frame) < FRAME_HEADER.size:
        raise ValueError(f"frame is only {len(frame)} bytes")
    message_type, sequence, stream_end, payload_length = FRAME_HEADER.unpack_from(frame)
    payload = frame[FRAME_HEADER.size :]
    if payload_length != len(payload):
        raise ValueError(
            f"declared payload {payload_length} != actual payload {len(payload)}"
        )
    return message_type, sequence, stream_end, payload


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("queue", help="POSIX MQ name, for example /csm_end0_to_main")
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--max-messages", type=int, default=100)
    parser.add_argument(
        "--idle-seconds",
        type=float,
        default=2.0,
        help="stop after the queue remains empty for this long",
    )
    args = parser.parse_args()
    if args.max_messages <= 0 or args.idle_seconds < 0:
        parser.error("--max-messages must be positive and --idle-seconds non-negative")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    frames_path = args.output_dir / "frames.csv"
    points_path = args.output_dir / "orbit_points.csv"
    frame_count = 0
    point_count = 0
    invalid_count = 0
    idle_started = time.monotonic()

    reader = MqReader(args.queue)
    try:
        with frames_path.open("w", newline="", encoding="utf-8") as frames_file, \
                points_path.open("w", newline="", encoding="utf-8") as points_file:
            frames = csv.writer(frames_file)
            points = csv.writer(points_file)
            frames.writerow(
                ["capture_index", "priority", "frame_bytes", "type", "sequence",
                 "stream_end", "payload_bytes", "valid", "error"]
            )
            points.writerow(["frame_index", "point_index", "timestamp_ms", "x", "y", "z",
                             "vx", "vy", "vz"])

            while frame_count < args.max_messages:
                received = reader.receive()
                if received is None:
                    if time.monotonic() - idle_started >= args.idle_seconds:
                        break
                    time.sleep(0.05)
                    continue
                idle_started = time.monotonic()
                frame, priority = received
                frame_path = args.output_dir / f"frame_{frame_count:04d}.bin"
                frame_path.write_bytes(frame)
                try:
                    message_type, sequence, stream_end, payload = decode_frame(frame)
                    valid = True
                    error = ""
                    if message_type == TYPE_RECEIVE_TELEMETRY_DATA:
                        if len(payload) % ORBIT_POINT.size:
                            raise ValueError(
                                f"type-4 payload {len(payload)} is not divisible by 56"
                            )
                        for point_index, values in enumerate(
                            ORBIT_POINT.iter_unpack(payload)
                        ):
                            points.writerow([frame_count, point_index, *values])
                            point_count += 1
                except ValueError as exc:
                    message_type = sequence = stream_end = -1
                    payload = b""
                    valid = False
                    error = str(exc)
                    invalid_count += 1
                frames.writerow(
                    [frame_count, priority, len(frame), message_type, sequence,
                     stream_end, len(payload), int(valid), error]
                )
                frame_count += 1
    finally:
        reader.close()

    print(
        f"CAPTURE queue={args.queue} frames={frame_count} points={point_count} "
        f"invalid={invalid_count} output={args.output_dir}"
    )
    return 0 if invalid_count == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())

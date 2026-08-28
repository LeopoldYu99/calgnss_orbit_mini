"""Send one CCSM PredictOrbit request without consuming the response queue."""

from __future__ import annotations

import argparse
import ctypes
import os
import struct


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("queue", help="request POSIX MQ name")
    parser.add_argument("duration_s", type=int)
    parser.add_argument("--start-time-s", type=int, default=0)
    parser.add_argument("--sequence", type=int, default=0)
    args = parser.parse_args()
    if not 0 <= args.duration_s <= 0xFFFFFFFF:
        parser.error("duration_s must fit uint32_t")
    if not 0 <= args.sequence <= 0xFFFFFFFF:
        parser.error("sequence must fit uint32_t")

    payload = struct.pack("<qI", args.start_time_s, args.duration_s)
    frame = struct.pack("<BIBI", 7, args.sequence, 0xFF, len(payload)) + payload

    libc = ctypes.CDLL(None, use_errno=True)
    libc.mq_open.argtypes = [ctypes.c_char_p, ctypes.c_int]
    libc.mq_open.restype = ctypes.c_int
    libc.mq_send.argtypes = [ctypes.c_int, ctypes.c_void_p, ctypes.c_size_t, ctypes.c_uint]
    libc.mq_send.restype = ctypes.c_int
    libc.mq_close.argtypes = [ctypes.c_int]
    libc.mq_close.restype = ctypes.c_int

    descriptor = libc.mq_open(args.queue.encode(), os.O_WRONLY)
    if descriptor < 0:
        code = ctypes.get_errno()
        raise OSError(code, os.strerror(code), args.queue)
    try:
        buffer = ctypes.create_string_buffer(frame)
        if libc.mq_send(descriptor, buffer, len(frame), 0) != 0:
            code = ctypes.get_errno()
            raise OSError(code, os.strerror(code), args.queue)
    finally:
        libc.mq_close(descriptor)

    print(
        f"SENT queue={args.queue} type=7 sequence={args.sequence} "
        f"start_time_s={args.start_time_s} duration_s={args.duration_s} bytes={len(frame)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

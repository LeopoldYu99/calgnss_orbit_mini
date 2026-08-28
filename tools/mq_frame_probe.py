"""Consume and inspect one CCSM fixed-binary POSIX MQ frame.

This tool removes one message from the queue. Use it only on a queue whose
message is intentionally being handed to this process.
"""

from __future__ import annotations

import argparse
import ctypes
import errno
import os
from pathlib import Path
import struct


class MqAttr(ctypes.Structure):
    _fields_ = [
        ("mq_flags", ctypes.c_long),
        ("mq_maxmsg", ctypes.c_long),
        ("mq_msgsize", ctypes.c_long),
        ("mq_curmsgs", ctypes.c_long),
    ]


def crc24q(data: bytes) -> int:
    crc = 0
    for byte in data:
        crc ^= byte << 16
        for _ in range(8):
            crc = ((crc << 1) ^ 0x1864CFB) & 0xFFFFFF if crc & 0x800000 else (crc << 1) & 0xFFFFFF
    return crc


def count_rtcm3(data: bytes) -> tuple[int, int, int]:
    offset = 0
    valid = 0
    invalid = 0
    skipped = 0
    while offset + 6 <= len(data):
        if data[offset] != 0xD3:
            offset += 1
            skipped += 1
            continue
        body_length = ((data[offset + 1] & 0x03) << 8) | data[offset + 2]
        frame_length = body_length + 6
        if offset + frame_length > len(data):
            break
        frame = data[offset : offset + frame_length]
        expected = int.from_bytes(frame[-3:], "big")
        if crc24q(frame[:-3]) == expected:
            valid += 1
        else:
            invalid += 1
        offset += frame_length
    return valid, invalid, skipped


def receive_one(queue_name: str) -> bytes:
    libc = ctypes.CDLL(None, use_errno=True)
    libc.mq_open.argtypes = [ctypes.c_char_p, ctypes.c_int]
    libc.mq_open.restype = ctypes.c_int
    libc.mq_getattr.argtypes = [ctypes.c_int, ctypes.POINTER(MqAttr)]
    libc.mq_getattr.restype = ctypes.c_int
    libc.mq_receive.argtypes = [ctypes.c_int, ctypes.c_void_p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_uint)]
    libc.mq_receive.restype = ctypes.c_ssize_t
    libc.mq_close.argtypes = [ctypes.c_int]
    libc.mq_close.restype = ctypes.c_int

    descriptor = libc.mq_open(queue_name.encode(), os.O_RDONLY | os.O_NONBLOCK)
    if descriptor < 0:
        code = ctypes.get_errno()
        raise OSError(code, os.strerror(code), queue_name)
    try:
        attributes = MqAttr()
        if libc.mq_getattr(descriptor, ctypes.byref(attributes)) != 0:
            code = ctypes.get_errno()
            raise OSError(code, os.strerror(code), queue_name)
        buffer = ctypes.create_string_buffer(attributes.mq_msgsize)
        priority = ctypes.c_uint()
        length = libc.mq_receive(descriptor, buffer, len(buffer), ctypes.byref(priority))
        if length < 0:
            code = ctypes.get_errno()
            if code == errno.EAGAIN:
                raise RuntimeError(f"queue {queue_name} is empty")
            raise OSError(code, os.strerror(code), queue_name)
        print(
            f"QUEUE name={queue_name} maxmsg={attributes.mq_maxmsg} "
            f"msgsize={attributes.mq_msgsize} remaining_before_receive={attributes.mq_curmsgs} "
            f"priority={priority.value}"
        )
        return buffer.raw[:length]
    finally:
        libc.mq_close(descriptor)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("queue")
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    frame = receive_one(args.queue)
    args.output.write_bytes(frame)
    print(f"MESSAGE bytes={len(frame)} backup={args.output}")
    print(f"HEX first64={frame[:64].hex(' ')}")
    if len(frame) < 10:
        print("FRAME invalid=shorter-than-10-byte-header")
        return 2

    message_type, message_seq, stream_end, payload_length = struct.unpack_from("<BIBI", frame)
    payload = frame[10:]
    print(
        f"FRAME type={message_type} sequence={message_seq} stream_end=0x{stream_end:02X} "
        f"payload_length={payload_length} actual_payload_length={len(payload)} "
        f"length_match={int(payload_length == len(payload))}"
    )
    if message_type == 3:
        valid, invalid, skipped = count_rtcm3(payload)
        print(
            f"RTCM valid_frames={valid} invalid_frames={invalid} "
            f"skipped_bytes={skipped} payload_starts_d3={int(payload.startswith(bytes([0xD3])))}"
        )
    elif message_type == 1 and len(payload) >= 13:
        timestamp_ms, service_status, message_length = struct.unpack_from("<qBI", payload)
        message_bytes = payload[13:]
        if message_length == len(message_bytes):
            message = message_bytes.decode("utf-8", "replace")
            print(
                f"STATUS timestamp_ms={timestamp_ms} service_status={service_status} "
                f"message={message!r}"
            )
        else:
            print(
                f"STATUS invalid_length={message_length} "
                f"actual_message_length={len(message_bytes)}"
            )
    return 0 if payload_length == len(payload) else 3


if __name__ == "__main__":
    raise SystemExit(main())

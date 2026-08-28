"""Extract CRC-valid RTCM3 frames from raw bytes or a MobaXterm hex dump."""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass, replace
from datetime import datetime, timezone
from pathlib import Path
import re


HEX_BYTE = re.compile(r"(?<![0-9A-Fa-f])[0-9A-Fa-f]{2}(?![0-9A-Fa-f])")


@dataclass(frozen=True)
class CaptureStats:
    source_bytes: int
    decoded_bytes: int
    valid_frames: int
    valid_bytes: int
    invalid_candidates: int
    message_types: dict[int, int]
    first_gps_tow_ms: int | None = None
    last_gps_tow_ms: int | None = None
    suggested_time_sync_ms: int | None = None
    suggested_last_utc_ms: int | None = None

    def summary(self) -> str:
        types = ", ".join(f"{key}×{value}" for key, value in self.message_types.items())
        time_hint = (
            f"，建议UTC授时={self.suggested_time_sync_ms}"
            if self.suggested_time_sync_ms is not None
            else ""
        )
        return (
            f"RTCM3有效帧={self.valid_frames}，有效字节={self.valid_bytes}，"
            f"CRC失败候选={self.invalid_candidates}，类型={types or '无'}{time_hint}"
        )


def crc24q(data: bytes) -> int:
    crc = 0
    for value in data:
        crc ^= value << 16
        for _ in range(8):
            if crc & 0x800000:
                crc = ((crc << 1) ^ 0x1864CFB) & 0xFFFFFF
            else:
                crc = (crc << 1) & 0xFFFFFF
    return crc


def _decode_source(source: bytes) -> bytes:
    # A raw stream normally contains a binary 0xD3 preamble. MobaXterm's
    # "terminal output" capture instead stores printable hexadecimal bytes.
    if b"\xd3" in source[:4096]:
        return source
    text = source.decode("utf-8", "replace")
    decoded = bytearray()
    for line in text.splitlines():
        # CPU/fEuler diagnostics were written in the middle of some hex lines.
        # Only the hexadecimal prefix belongs to the RTCM byte stream.
        cut = len(line)
        for marker in ("[", "{"):
            position = line.find(marker)
            if position >= 0:
                cut = min(cut, position)
        decoded.extend(int(token, 16) for token in HEX_BYTE.findall(line[:cut]))
    return bytes(decoded)


def _get_unsigned_bits(data: bytes, position: int, length: int) -> int:
    value = 0
    for index in range(length):
        value = (value << 1) | (
            (data[(position + index) // 8] >> (7 - (position + index) % 8)) & 1
        )
    return value


def _resolve_gps_tow_utc_ms(tow_ms: int, reference_utc: float) -> int:
    # GPS time has been 18 seconds ahead of UTC since 2017. Pick the GPS week
    # nearest the capture file timestamp; RTCM MSM carries TOW but not the week.
    gps_epoch_unix = datetime(1980, 1, 6, tzinfo=timezone.utc).timestamp()
    leap_seconds = 18
    reference_gps_seconds = reference_utc - gps_epoch_unix + leap_seconds
    week = round((reference_gps_seconds - tow_ms / 1000.0) / 604800.0)
    utc_seconds = gps_epoch_unix + week * 604800.0 + tow_ms / 1000.0 - leap_seconds
    return round(utc_seconds * 1000.0)


def _capture_reference_timestamp(path: Path) -> float:
    date_match = re.search(r"(20\d{6})", path.name)
    if date_match:
        capture_date = datetime.strptime(date_match.group(1), "%Y%m%d").replace(
            hour=12, tzinfo=timezone.utc
        )
        return capture_date.timestamp()
    return path.stat().st_mtime


def extract_frames(source: bytes) -> tuple[bytes, CaptureStats]:
    decoded = _decode_source(source)
    output = bytearray()
    message_types: Counter[int] = Counter()
    gps_tows: list[int] = []
    invalid = 0
    cursor = 0
    while cursor + 6 <= len(decoded):
        if decoded[cursor] != 0xD3 or decoded[cursor + 1] & 0xFC:
            cursor += 1
            continue
        payload_size = ((decoded[cursor + 1] & 0x03) << 8) | decoded[cursor + 2]
        end = cursor + payload_size + 6
        if end > len(decoded):
            break
        frame = decoded[cursor:end]
        expected_crc = int.from_bytes(frame[-3:], "big")
        if crc24q(frame[:-3]) != expected_crc:
            invalid += 1
            cursor += 1
            continue
        output.extend(frame)
        if payload_size >= 2:
            message_type = (frame[3] << 4) | (frame[4] >> 4)
            message_types[message_type] += 1
            if 1071 <= message_type <= 1077 and payload_size * 8 >= 54:
                gps_tows.append(_get_unsigned_bits(frame[3:-3], 24, 30))
        cursor = end
    stats = CaptureStats(
        source_bytes=len(source),
        decoded_bytes=len(decoded),
        valid_frames=sum(message_types.values()),
        valid_bytes=len(output),
        invalid_candidates=invalid,
        message_types=dict(sorted(message_types.items())),
        first_gps_tow_ms=min(gps_tows) if gps_tows else None,
        last_gps_tow_ms=max(gps_tows) if gps_tows else None,
    )
    return bytes(output), stats


def convert_file(source_path: Path, destination_path: Path) -> CaptureStats:
    frames, stats = extract_frames(source_path.read_bytes())
    if not frames:
        raise ValueError("文件中没有找到通过 CRC24Q 校验的 RTCM3 帧")
    destination_path.write_bytes(frames)
    if stats.first_gps_tow_ms is not None:
        stats = replace(
            stats,
            suggested_time_sync_ms=_resolve_gps_tow_utc_ms(
                stats.first_gps_tow_ms, _capture_reference_timestamp(source_path)
            ),
            suggested_last_utc_ms=_resolve_gps_tow_utc_ms(
                stats.last_gps_tow_ms, _capture_reference_timestamp(source_path)
            ) if stats.last_gps_tow_ms is not None else None,
        )
    return stats


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    arguments = parser.parse_args()
    stats = convert_file(arguments.input, arguments.output)
    print(stats.summary())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

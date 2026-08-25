import csv
import math
from datetime import datetime


# ============================================================
# 配置
# ============================================================

input_log = "270km.log"
output_csv = "270km_NMEA_to_ECEF_position_only.csv"

# 坐标输出单位：
# True  = km，适合 STK / SimGEN 常见 Fixed_Position 文件
# False = m
OUTPUT_KM = True

# 质量筛选条件
MIN_SATS = 8          # 最少参与解算卫星数
MAX_HDOP = 1.5        # 最大 HDOP，越小越好

# 是否严格要求 RMC 状态为 A
REQUIRE_RMC_VALID = True


# ============================================================
# WGS84 椭球参数
# ============================================================

WGS84_A = 6378137.0                  # 长半轴，m
WGS84_F = 1.0 / 298.257223563
WGS84_E2 = WGS84_F * (2.0 - WGS84_F)


# ============================================================
# 工具函数
# ============================================================

def strip_checksum(nmea_line: str) -> str:
    """
    去掉 NMEA 校验和部分。
    $GNGGA,...*5B -> $GNGGA,...
    """
    if "*" in nmea_line:
        return nmea_line.split("*", 1)[0]
    return nmea_line


def extract_nmea_sentence(line: str) -> str | None:
    """
    从 log 行中提取 NMEA 句子。
    兼容这种格式：
    [2026-06-30 ...]# RECV ASCII ... <<<
    $GNGGA,...
    """
    if "$" not in line:
        return None

    sentence = line[line.find("$"):].strip()
    sentence = strip_checksum(sentence)
    return sentence


def nmea_latlon_to_deg(value: str, hemi: str) -> float | None:
    """
    NMEA 经纬度格式转十进制度。

    纬度格式：ddmm.mmmm
    经度格式：dddmm.mmmm

    例如：
    1724.2594222,N -> 17 + 24.2594222 / 60
    15304.2955483,W -> -(153 + 04.2955483 / 60)
    """
    if not value or not hemi:
        return None

    raw = float(value)
    degrees = int(raw // 100)
    minutes = raw - degrees * 100
    decimal_deg = degrees + minutes / 60.0

    if hemi.upper() in ["S", "W"]:
        decimal_deg = -decimal_deg

    return decimal_deg


def lla_to_ecef(lat_deg: float, lon_deg: float, h_m: float):
    """
    WGS84 经纬高转 ECEF。

    输入：
        lat_deg: 纬度，单位 deg
        lon_deg: 经度，单位 deg
        h_m: 椭球高，单位 m

    输出：
        x, y, z，单位 m
    """
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)

    sin_lat = math.sin(lat)
    cos_lat = math.cos(lat)

    sin_lon = math.sin(lon)
    cos_lon = math.cos(lon)

    N = WGS84_A / math.sqrt(1.0 - WGS84_E2 * sin_lat * sin_lat)

    x = (N + h_m) * cos_lat * cos_lon
    y = (N + h_m) * cos_lat * sin_lon
    z = (N * (1.0 - WGS84_E2) + h_m) * sin_lat

    return x, y, z


def parse_rmc(fields):
    """
    解析 RMC 句子，主要取日期和有效状态。

    RMC 示例：
    $GNRMC,040403.00,A,1724.2594222,N,15304.2955483,W,15184.663,349.74,230626,,,A,V

    字段含义：
    fields[1] = UTC 时间
    fields[2] = 状态，A 有效，V 无效
    fields[9] = 日期 ddmmyy
    """
    if len(fields) < 10:
        return None

    status = fields[2].strip().upper()
    date_str = fields[9].strip()

    if REQUIRE_RMC_VALID and status != "A":
        return None

    if not date_str or len(date_str) != 6:
        return None

    day = int(date_str[0:2])
    month = int(date_str[2:4])
    year = int(date_str[4:6]) + 2000

    return {
        "date": (year, month, day),
        "valid": status == "A"
    }


def parse_gga(fields, current_date):
    """
    解析 GGA 句子，并做质量筛选。

    GGA 示例：
    $GNGGA,040403.00,1724.2594222,N,15304.2955483,W,1,14,0.69,271441.347,M,0.989,M,,

    字段含义：
    fields[1]  = UTC 时间
    fields[2]  = 纬度
    fields[3]  = N/S
    fields[4]  = 经度
    fields[5]  = E/W
    fields[6]  = Fix Quality，0 表示无效
    fields[7]  = 参与解算卫星数
    fields[8]  = HDOP
    fields[9]  = 海拔高度，MSL，高于平均海平面，m
    fields[11] = geoid separation，大地水准面相对 WGS84 椭球的差值，m
    """
    if len(fields) < 12:
        return None, "bad_gga_length"

    if current_date is None:
        return None, "no_date"

    utc_time = fields[1].strip()
    lat_raw = fields[2].strip()
    lat_hemi = fields[3].strip()
    lon_raw = fields[4].strip()
    lon_hemi = fields[5].strip()
    fix_quality = fields[6].strip()
    num_sats = fields[7].strip()
    hdop = fields[8].strip()
    alt_msl = fields[9].strip()
    geoid_sep = fields[11].strip()

    # 1. 定位质量筛选
    # fix_quality = 0 表示无效定位
    if not fix_quality or fix_quality == "0":
        return None, "invalid_fix"

    # 2. 时间字段检查
    if not utc_time:
        return None, "missing_time"

    # 3. 卫星数和 HDOP 筛选
    try:
        num_sats_int = int(num_sats)
        hdop_float = float(hdop)
    except ValueError:
        return None, "bad_quality_fields"

    if num_sats_int < MIN_SATS:
        return None, "too_few_sats"

    if hdop_float > MAX_HDOP:
        return None, "hdop_too_large"

    # 4. 经纬度解析
    try:
        lat_deg = nmea_latlon_to_deg(lat_raw, lat_hemi)
        lon_deg = nmea_latlon_to_deg(lon_raw, lon_hemi)
    except ValueError:
        return None, "bad_latlon"

    if lat_deg is None or lon_deg is None:
        return None, "missing_latlon"

    # 5. 高度解析
    # GGA 的 altitude 是 MSL 高，也就是高于平均海平面的高度。
    # ECEF 需要 WGS84 椭球高。
    #
    # 椭球高 h = MSL 高度 + geoid separation
    try:
        h_m = float(alt_msl)
        if geoid_sep:
            h_m += float(geoid_sep)
    except ValueError:
        return None, "bad_height"

    # 6. 时间解析
    try:
        hour = int(utc_time[0:2])
        minute = int(utc_time[2:4])

        second_float = float(utc_time[4:])
        second = int(second_float)
        microsecond = int(round((second_float - second) * 1_000_000))

        # 防止 59.999999 四舍五入成 1000000 微秒
        if microsecond >= 1_000_000:
            second += 1
            microsecond -= 1_000_000

        year, month, day = current_date
        dt = datetime(year, month, day, hour, minute, second, microsecond)

    except Exception:
        return None, "bad_time"

    # 7. LLA -> ECEF
    x_m, y_m, z_m = lla_to_ecef(lat_deg, lon_deg, h_m)

    if OUTPUT_KM:
        x = x_m / 1000.0
        y = y_m / 1000.0
        z = z_m / 1000.0
    else:
        x = x_m
        y = y_m
        z = z_m

    return {
        "time": dt,
        "x": x,
        "y": y,
        "z": z,
        "lat_deg": lat_deg,
        "lon_deg": lon_deg,
        "h_m": h_m,
        "num_sats": num_sats_int,
        "hdop": hdop_float,
        "fix_quality": fix_quality,
    }, None


def format_stk_utcg(dt: datetime) -> str:
    """
    输出 STK 常见 UTCG 时间格式：
    23 Jun 2026 04:04:03.000
    """
    return dt.strftime("%d %b %Y %H:%M:%S.") + f"{int(dt.microsecond / 1000):03d}"


# ============================================================
# 主程序
# ============================================================

def main():
    rows = []
    current_date = None

    stats = {
        "total_lines": 0,
        "nmea_lines": 0,
        "rmc_lines": 0,
        "gga_lines": 0,
        "output_rows": 0,
        "bad_gga_length": 0,
        "no_date": 0,
        "invalid_fix": 0,
        "missing_time": 0,
        "bad_quality_fields": 0,
        "too_few_sats": 0,
        "hdop_too_large": 0,
        "bad_latlon": 0,
        "missing_latlon": 0,
        "bad_height": 0,
        "bad_time": 0,
    }

    with open(input_log, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            stats["total_lines"] += 1

            sentence = extract_nmea_sentence(line)
            if sentence is None:
                continue

            stats["nmea_lines"] += 1

            fields = sentence.split(",")
            sentence_type = fields[0]

            # RMC 提供日期
            if sentence_type.endswith("RMC"):
                stats["rmc_lines"] += 1

                rmc_info = parse_rmc(fields)
                if rmc_info is not None:
                    current_date = rmc_info["date"]

            # GGA 提供经纬高、卫星数、HDOP
            elif sentence_type.endswith("GGA"):
                stats["gga_lines"] += 1

                record, reason = parse_gga(fields, current_date)

                if record is None:
                    if reason in stats:
                        stats[reason] += 1
                    continue

                rows.append(record)
                stats["output_rows"] += 1

    # 写出 CSV
    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)

        # 只输出坐标，不输出速度，也不输出 HDOP/卫星数
        writer.writerow(["Time (UTCG)", "x", "y", "z"])

        for r in rows:
            writer.writerow([
                format_stk_utcg(r["time"]),
                f"{r['x']:.9f}",
                f"{r['y']:.9f}",
                f"{r['z']:.9f}",
            ])

    print("转换完成")
    print(f"输入文件: {input_log}")
    print(f"输出文件: {output_csv}")
    print(f"输出坐标单位: {'km' if OUTPUT_KM else 'm'}")
    print(f"最少卫星数 MIN_SATS = {MIN_SATS}")
    print(f"最大 HDOP MAX_HDOP = {MAX_HDOP}") 

    print("\n统计信息:")
    for k, v in stats.items():
        print(f"{k}: {v}")


if __name__ == "__main__":
    main()
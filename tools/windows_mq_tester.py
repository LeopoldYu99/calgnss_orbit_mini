"""Windows GUI for testing OrbitPredictionService through an SSH-hosted MQ CLI."""

from __future__ import annotations

import csv
import queue
import shlex
import threading
import time
import tempfile
from dataclasses import dataclass
from pathlib import Path
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

from rtcm_capture import convert_file

try:
    import paramiko
except ImportError:  # Display a useful GUI error instead of failing at import time.
    paramiko = None


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SERVER = ROOT / "build-riscv64-static" / "orbit_prediction_server_static"
DEFAULT_CLI = ROOT / "build-riscv64-static" / "orbit_mq_cli_static"
REMOTE_SERVER = "/tmp/orbit_gui_server"
REMOTE_CLI = "/tmp/orbit_gui_cli"
REMOTE_PID = "/tmp/orbit_gui_server.pid"
REMOTE_LOG = "/tmp/orbit_gui_server.log"
REQUEST_QUEUE = "/orbit_gui_requests"
RESPONSE_QUEUE = "/orbit_gui_responses"


@dataclass
class OrbitPoint:
    timestamp_ms: int
    x: float
    y: float
    z: float
    vx: float
    vy: float
    vz: float


class RemoteSession:
    def __init__(self, host: str, username: str, password: str) -> None:
        if paramiko is None:
            raise RuntimeError(
                "缺少 Paramiko。请先运行：py -3 -m pip install -r tools/requirements-windows.txt"
            )
        self.client = paramiko.SSHClient()
        self.client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        self.client.connect(
            hostname=host,
            username=username,
            password=password,
            timeout=10,
            banner_timeout=15,
            auth_timeout=15,
            look_for_keys=True,
            allow_agent=True,
        )

    def close(self) -> None:
        self.client.close()

    def run(self, command: str, timeout: int = 30) -> tuple[int, str, str]:
        stdin, stdout, stderr = self.client.exec_command(command, timeout=timeout)
        stdin.close()
        code = stdout.channel.recv_exit_status()
        return code, stdout.read().decode("utf-8", "replace"), stderr.read().decode(
            "utf-8", "replace"
        )

    def upload(self, local_path: Path, remote_path: str) -> None:
        with self.client.open_sftp() as sftp:
            sftp.put(str(local_path), remote_path)

    def download(self, remote_path: str, local_path: Path) -> None:
        with self.client.open_sftp() as sftp:
            sftp.get(remote_path, str(local_path))


class OrbitMqTester(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.title("轨道预报 MQ 测试台")
        self.geometry("1280x820")
        self.minsize(980, 680)
        self.events: queue.Queue[tuple[str, object]] = queue.Queue()
        self.points: list[OrbitPoint] = []
        self._build_style()
        self._build_ui()
        self.after(100, self._drain_events)

    def _build_style(self) -> None:
        style = ttk.Style(self)
        if "vista" in style.theme_names():
            style.theme_use("vista")
        style.configure("Title.TLabel", font=("Microsoft YaHei UI", 17, "bold"))
        style.configure("Sub.TLabel", foreground="#52606d")
        style.configure("Metric.TLabel", font=("Consolas", 11, "bold"))
        style.configure("Accent.TButton", font=("Microsoft YaHei UI", 10, "bold"))

    def _build_ui(self) -> None:
        header = ttk.Frame(self, padding=(16, 12, 16, 6))
        header.pack(fill="x")
        ttk.Label(header, text="轨道预报 MQ 测试台", style="Title.TLabel").pack(side="left")
        ttk.Label(
            header,
            text="Windows → SSH → RISC-V MQ CLI → OrbitPredictionService",
            style="Sub.TLabel",
        ).pack(side="left", padx=18, pady=(6, 0))

        connection = ttk.LabelFrame(self, text="目标机与隔离测试服务", padding=10)
        connection.pack(fill="x", padx=16, pady=6)
        self.host = tk.StringVar(value="192.168.104.100")
        self.user = tk.StringVar(value="root")
        self.password = tk.StringVar(value="root")
        self.server_binary = tk.StringVar(value=str(DEFAULT_SERVER))
        self.cli_binary = tk.StringVar(value=str(DEFAULT_CLI))
        fields = [
            ("IP", self.host, 16, False),
            ("用户", self.user, 9, False),
            ("密码", self.password, 9, True),
            ("服务程序", self.server_binary, 42, False),
            ("MQ CLI", self.cli_binary, 42, False),
        ]
        for column, (label, variable, width, secret) in enumerate(fields):
            ttk.Label(connection, text=label).grid(row=0, column=column * 2, padx=(2, 4))
            ttk.Entry(connection, textvariable=variable, width=width, show="*" if secret else "").grid(
                row=0, column=column * 2 + 1, padx=(0, 8), sticky="ew"
            )
        connection.columnconfigure(7, weight=1)
        connection.columnconfigure(9, weight=1)
        ttk.Button(connection, text="连接测试", command=self.check_connection).grid(
            row=1, column=0, columnspan=2, padx=2, pady=(9, 0), sticky="ew"
        )
        ttk.Button(
            connection, text="部署并启动隔离服务", style="Accent.TButton", command=self.deploy
        ).grid(row=1, column=2, columnspan=3, padx=4, pady=(9, 0), sticky="ew")
        ttk.Button(connection, text="查询状态", command=self.status).grid(
            row=1, column=5, columnspan=2, padx=4, pady=(9, 0), sticky="ew"
        )
        ttk.Button(connection, text="复位", command=self.reset_service).grid(
            row=1, column=7, padx=4, pady=(9, 0), sticky="ew"
        )
        ttk.Button(connection, text="停止测试服务", command=self.stop_service).grid(
            row=1, column=8, columnspan=2, padx=4, pady=(9, 0), sticky="ew"
        )

        body = ttk.Panedwindow(self, orient="horizontal")
        body.pack(fill="both", expand=True, padx=16, pady=6)
        controls = ttk.Frame(body, padding=(0, 0, 10, 0))
        output = ttk.Frame(body)
        body.add(controls, weight=1)
        body.add(output, weight=3)

        self._build_inputs(controls)
        self._build_output(output)

    def _build_inputs(self, parent: ttk.Frame) -> None:
        time_box = ttk.LabelFrame(parent, text="1. 授时", padding=10)
        time_box.pack(fill="x", pady=(0, 8))
        self.timestamp_ms = tk.StringVar(value="1787283546000")
        ttk.Entry(time_box, textvariable=self.timestamp_ms).pack(fill="x")
        ttk.Button(time_box, text="发送授时", command=self.send_time).pack(fill="x", pady=(7, 0))

        input_box = ttk.LabelFrame(parent, text="2. RTCM 输入", padding=10)
        input_box.pack(fill="x", pady=8)
        self.rtcm_file = tk.StringVar()
        ttk.Label(input_box, text="RTCM3 二进制文件").pack(anchor="w")
        self._file_row(input_box, self.rtcm_file, [("发送 RTCM", self.send_rtcm)])
        ttk.Label(
            input_box,
            text="RTCM 将由服务端解码/解算为带时间戳的 ECEF 观测，再用于轨道拟合。",
            style="Sub.TLabel",
            wraplength=310,
        ).pack(anchor="w", pady=(8, 0))

        predict_box = ttk.LabelFrame(parent, text="3. 轨道位置/速度输出（J2000）", padding=10)
        predict_box.pack(fill="x", pady=8)
        self.start_time_s = tk.StringVar(value="1787284669")
        self.duration_s = tk.StringVar(value="3600")
        grid = ttk.Frame(predict_box)
        grid.pack(fill="x")
        ttk.Label(grid, text="起始 Unix 秒").grid(row=0, column=0, sticky="w")
        ttk.Entry(grid, textvariable=self.start_time_s, width=19).grid(
            row=1, column=0, padx=(0, 8), sticky="ew"
        )
        ttk.Label(grid, text="时长（秒）").grid(row=0, column=1, sticky="w")
        ttk.Entry(grid, textvariable=self.duration_s, width=9).grid(row=1, column=1, sticky="ew")
        grid.columnconfigure(0, weight=1)
        ttk.Button(
            predict_box, text="开始预报并下载 CSV", style="Accent.TButton", command=self.predict
        ).pack(fill="x", pady=(9, 0))

        log_box = ttk.LabelFrame(parent, text="运行日志", padding=6)
        log_box.pack(fill="both", expand=True, pady=(8, 0))
        self.log_text = tk.Text(log_box, height=12, wrap="word", font=("Consolas", 9))
        self.log_text.pack(fill="both", expand=True)

    def _file_row(self, parent: ttk.Frame, variable: tk.StringVar, actions: list[tuple[str, object]]) -> None:
        row = ttk.Frame(parent)
        row.pack(fill="x", pady=(3, 0))
        ttk.Entry(row, textvariable=variable).pack(side="left", fill="x", expand=True)
        ttk.Button(
            row,
            text="…",
            width=3,
            command=lambda: self.choose_file(variable),
        ).pack(side="left", padx=4)
        for text, command in actions:
            ttk.Button(row, text=text, command=command).pack(side="left")

    def _build_output(self, parent: ttk.Frame) -> None:
        metrics = ttk.Frame(parent)
        metrics.pack(fill="x", pady=(0, 7))
        self.point_metric = tk.StringVar(value="点数：—")
        self.elapsed_metric = tk.StringVar(value="端到端耗时：—")
        self.range_metric = tk.StringVar(value="时间范围：—")
        for variable in (self.point_metric, self.elapsed_metric, self.range_metric):
            ttk.Label(metrics, textvariable=variable, style="Metric.TLabel").pack(
                side="left", padx=(0, 24)
            )

        chart_box = ttk.LabelFrame(parent, text="轨迹投影", padding=6)
        chart_box.pack(fill="both", expand=True)
        toolbar = ttk.Frame(chart_box)
        toolbar.pack(fill="x")
        ttk.Label(toolbar, text="投影").pack(side="left")
        self.projection = tk.StringVar(value="X-Y")
        projection = ttk.Combobox(
            toolbar,
            textvariable=self.projection,
            values=("X-Y", "X-Z", "Y-Z"),
            width=7,
            state="readonly",
        )
        projection.pack(side="left", padx=6)
        projection.bind("<<ComboboxSelected>>", lambda _event: self.draw_plot())
        ttk.Label(toolbar, text="单位：m；速度输出见下表", style="Sub.TLabel").pack(side="left")
        self.canvas = tk.Canvas(chart_box, background="#0c1724", highlightthickness=0, height=360)
        self.canvas.pack(fill="both", expand=True, pady=(5, 0))
        self.canvas.bind("<Configure>", lambda _event: self.draw_plot())

        table_box = ttk.LabelFrame(parent, text="定位/轨道数据采样", padding=6)
        table_box.pack(fill="both", expand=False, pady=(8, 0))
        columns = ("time", "x", "y", "z", "vx", "vy", "vz")
        self.table = ttk.Treeview(table_box, columns=columns, show="headings", height=8)
        widths = (125, 105, 105, 105, 90, 90, 90)
        labels = ("timestamp_ms", "X(m)", "Y(m)", "Z(m)", "Vx(m/s)", "Vy(m/s)", "Vz(m/s)")
        for column, width, label in zip(columns, widths, labels):
            self.table.heading(column, text=label)
            self.table.column(column, width=width, anchor="e")
        scrollbar = ttk.Scrollbar(table_box, orient="vertical", command=self.table.yview)
        self.table.configure(yscrollcommand=scrollbar.set)
        self.table.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

    def choose_file(self, variable: tk.StringVar) -> None:
        path = filedialog.askopenfilename()
        if path:
            variable.set(path)

    def _credentials(self) -> tuple[str, str, str]:
        return self.host.get().strip(), self.user.get().strip(), self.password.get()

    def _with_session(self, action):
        session = RemoteSession(*self._credentials())
        try:
            return action(session)
        finally:
            session.close()

    def _run_async(self, label: str, action, done=None) -> None:
        self._log(f"▶ {label}")

        def worker() -> None:
            try:
                result = action()
                self.events.put(("done", (label, result, done)))
            except Exception as exc:  # Show operational failures in the GUI.
                self.events.put(("error", (label, str(exc))))

        threading.Thread(target=worker, daemon=True).start()

    def _drain_events(self) -> None:
        try:
            while True:
                kind, payload = self.events.get_nowait()
                if kind == "done":
                    label, result, callback = payload
                    self._log(f"✓ {label}")
                    # Prediction returns (CLI text, 3601 parsed points).  Logging
                    # the tuple would stringify every point and stall Tk for a
                    # large response; let its completion callback summarize it.
                    if result and (isinstance(result, str) or callback is None):
                        self._log(str(result).strip())
                    if callback:
                        callback(result)
                elif kind == "error":
                    label, error = payload
                    self._log(f"✗ {label}: {error}")
                    messagebox.showerror(label, error)
                elif kind == "rtcm_times":
                    sync_ms, last_ms = payload
                    self.timestamp_ms.set(str(sync_ms))
                    self.start_time_s.set(str(last_ms // 1000))
                    self._log(
                        f"已按RTCM历元设置授时={sync_ms}、预报起点={last_ms // 1000}"
                    )
        except queue.Empty:
            pass
        self.after(100, self._drain_events)

    def _log(self, message: str) -> None:
        stamp = time.strftime("%H:%M:%S")
        self.log_text.insert("end", f"[{stamp}] {message}\n")
        self.log_text.see("end")

    @staticmethod
    def _checked(code: int, stdout: str, stderr: str) -> str:
        text = "\n".join(part.strip() for part in (stdout, stderr) if part.strip())
        if code != 0:
            raise RuntimeError(text or f"远程命令退出码 {code}")
        return text

    @staticmethod
    def _mq_command(arguments: str) -> str:
        return (
            f"ORBIT_MQ_REQUEST_QUEUE={REQUEST_QUEUE} "
            f"ORBIT_MQ_RESPONSE_QUEUE={RESPONSE_QUEUE} {REMOTE_CLI} {arguments}"
        )

    def check_connection(self) -> None:
        self._run_async(
            "连接测试",
            lambda: self._with_session(
                lambda session: self._checked(*session.run("uname -a; cat /proc/sys/fs/mqueue/msgsize_max"))
            ),
        )

    def deploy(self) -> None:
        server = Path(self.server_binary.get())
        cli = Path(self.cli_binary.get())

        def action() -> str:
            if not server.is_file() or not cli.is_file():
                raise RuntimeError("静态服务程序或 MQ CLI 不存在，请先执行静态构建")

            def remote(session: RemoteSession) -> str:
                session.upload(server, REMOTE_SERVER)
                session.upload(cli, REMOTE_CLI)
                command = (
                    f"set -e; if read oldpid <{REMOTE_PID} 2>/dev/null; then "
                    f"kill \"$oldpid\" 2>/dev/null || true; fi; "
                    f"rm -f /dev/mqueue/{REQUEST_QUEUE[1:]} /dev/mqueue/{RESPONSE_QUEUE[1:]} "
                    f"{REMOTE_LOG}; chmod 755 {REMOTE_SERVER} {REMOTE_CLI}; "
                    f"ORBIT_MQ_REQUEST_QUEUE={REQUEST_QUEUE} "
                    f"ORBIT_MQ_RESPONSE_QUEUE={RESPONSE_QUEUE} "
                    f"nohup {REMOTE_SERVER} >{REMOTE_LOG} 2>&1 </dev/null & echo $! >{REMOTE_PID}; "
                    "sleep 1; read pid <" + REMOTE_PID + "; kill -0 \"$pid\"; "
                    + self._mq_command("status")
                )
                return self._checked(*session.run(command, timeout=40))

            return self._with_session(remote)

        self._run_async("部署并启动隔离服务", action)

    def status(self) -> None:
        self._simple_remote("查询状态", "status")

    def reset_service(self) -> None:
        self._simple_remote("复位服务", "reset")

    def send_time(self) -> None:
        timestamp = self.timestamp_ms.get().strip()
        if not timestamp.isdigit():
            messagebox.showerror("授时", "timestamp_ms 必须是非负整数")
            return
        self._simple_remote("发送授时", f"time-sync {timestamp}")

    def _simple_remote(self, label: str, arguments: str) -> None:
        self._run_async(
            label,
            lambda: self._with_session(
                lambda session: self._checked(*session.run(self._mq_command(arguments), timeout=30))
            ),
        )

    def _upload_and_run(self, label: str, local_text: str, command: str) -> None:
        local = Path(local_text)

        def action() -> str:
            if not local.is_file():
                raise RuntimeError(f"文件不存在：{local}")
            remote_file = f"/tmp/orbit_gui_input_{int(time.time() * 1000)}{local.suffix}"

            def remote(session: RemoteSession) -> str:
                session.upload(local, remote_file)
                try:
                    full = self._mq_command(command.format(file=shlex.quote(remote_file)))
                    return self._checked(*session.run(full, timeout=120))
                finally:
                    session.run(f"rm -f {shlex.quote(remote_file)}")

            return self._with_session(remote)

        self._run_async(label, action)

    def send_rtcm(self) -> None:
        local = Path(self.rtcm_file.get())

        def action() -> str:
            if not local.is_file():
                raise RuntimeError(f"文件不存在：{local}")
            temporary = tempfile.NamedTemporaryFile(suffix=".rtcm3", delete=False)
            temporary.close()
            cleaned = Path(temporary.name)
            try:
                stats = convert_file(local, cleaned)
                if (
                    stats.suggested_time_sync_ms is not None
                    and stats.suggested_last_utc_ms is not None
                ):
                    self.events.put(
                        (
                            "rtcm_times",
                            (stats.suggested_time_sync_ms, stats.suggested_last_utc_ms),
                        )
                    )
                remote_file = f"/tmp/orbit_gui_input_{int(time.time() * 1000)}.rtcm3"

                def remote(session: RemoteSession) -> str:
                    messages = [stats.summary()]
                    if stats.suggested_time_sync_ms is not None:
                        messages.append(
                            self._checked(
                                *session.run(
                                    self._mq_command(
                                        f"time-sync {stats.suggested_time_sync_ms}"
                                    ),
                                    timeout=30,
                                )
                            )
                        )
                    session.upload(cleaned, remote_file)
                    try:
                        result = self._checked(
                            *session.run(
                                self._mq_command(f"rtcm {shlex.quote(remote_file)}"),
                                timeout=300,
                            )
                        )
                        messages.append(result)
                        return "\n".join(messages)
                    finally:
                        session.run(f"rm -f {shlex.quote(remote_file)}")

                return self._with_session(remote)
            finally:
                cleaned.unlink(missing_ok=True)

        self._run_async("清洗并发送 RTCM", action)

    def stop_service(self) -> None:
        def action() -> str:
            command = (
                f"if read pid <{REMOTE_PID} 2>/dev/null; then kill \"$pid\" 2>/dev/null || true; fi; "
                f"rm -f /dev/mqueue/{REQUEST_QUEUE[1:]} /dev/mqueue/{RESPONSE_QUEUE[1:]} {REMOTE_PID}; "
                "echo isolated_test_service_stopped"
            )
            return self._with_session(
                lambda session: self._checked(*session.run(command, timeout=20))
            )

        self._run_async("停止测试服务", action)

    def predict(self) -> None:
        start = self.start_time_s.get().strip()
        duration = self.duration_s.get().strip()
        if not start.isdigit() or not duration.isdigit() or not 1 <= int(duration) <= 3600:
            messagebox.showerror("预报", "起始时间必须为非负整数，时长必须为 1～3600 秒")
            return
        output = filedialog.asksaveasfilename(
            title="保存轨道预报 CSV",
            defaultextension=".csv",
            initialfile=f"orbit_prediction_{duration}s.csv",
            filetypes=[("CSV", "*.csv")],
        )
        if not output:
            return
        local_output = Path(output)
        remote_output = f"/tmp/orbit_gui_prediction_{int(time.time() * 1000)}.csv"

        def action():
            def remote(session: RemoteSession):
                command = self._mq_command(
                    f"predict {start} {duration} {shlex.quote(remote_output)} 600"
                )
                code, stdout, stderr = session.run(command, timeout=660)
                text = self._checked(code, stdout, stderr)
                session.download(remote_output, local_output)
                session.run(f"rm -f {shlex.quote(remote_output)}")
                return text, load_csv(local_output)

            return self._with_session(remote)

        self._run_async("轨道预报", action, self._prediction_done)

    def _prediction_done(self, result) -> None:
        text, points = result
        if text:
            self._log(text.strip())
        self.points = points
        elapsed = "—"
        for token in text.replace("\n", " ").split():
            if token.startswith("elapsed_ms="):
                elapsed = f"{int(token.split('=', 1)[1]) / 1000:.3f} s"
        self.point_metric.set(f"点数：{len(points)}")
        self.elapsed_metric.set(f"端到端耗时：{elapsed}")
        if points:
            self.range_metric.set(
                f"时间范围：{points[0].timestamp_ms} ～ {points[-1].timestamp_ms}"
            )
        self._fill_table()
        self.draw_plot()

    def _fill_table(self) -> None:
        self.table.delete(*self.table.get_children())
        if not self.points:
            return
        stride = max(1, len(self.points) // 200)
        selected = self.points[::stride]
        if selected[-1] is not self.points[-1]:
            selected.append(self.points[-1])
        for point in selected:
            self.table.insert(
                "",
                "end",
                values=(
                    point.timestamp_ms,
                    f"{point.x:.3f}",
                    f"{point.y:.3f}",
                    f"{point.z:.3f}",
                    f"{point.vx:.6f}",
                    f"{point.vy:.6f}",
                    f"{point.vz:.6f}",
                ),
            )

    def draw_plot(self) -> None:
        canvas = self.canvas
        canvas.delete("all")
        width = max(canvas.winfo_width(), 200)
        height = max(canvas.winfo_height(), 150)
        margin = 48
        canvas.create_text(
            width / 2,
            height / 2,
            text="运行预报后显示轨迹",
            fill="#708399",
            font=("Microsoft YaHei UI", 13),
            tags="placeholder",
        )
        if len(self.points) < 2:
            return
        canvas.delete("placeholder")
        projection = self.projection.get()
        axes = {"X-Y": ("x", "y"), "X-Z": ("x", "z"), "Y-Z": ("y", "z")}
        horizontal, vertical = axes[projection]
        values = [(getattr(p, horizontal), getattr(p, vertical)) for p in self.points]
        min_x, max_x = min(v[0] for v in values), max(v[0] for v in values)
        min_y, max_y = min(v[1] for v in values), max(v[1] for v in values)
        span_x = max(max_x - min_x, 1.0)
        span_y = max(max_y - min_y, 1.0)
        for index in range(6):
            x = margin + (width - 2 * margin) * index / 5
            y = margin + (height - 2 * margin) * index / 5
            canvas.create_line(x, margin, x, height - margin, fill="#203246")
            canvas.create_line(margin, y, width - margin, y, fill="#203246")
        coordinates: list[float] = []
        for x_value, y_value in values:
            coordinates.extend(
                [
                    margin + (x_value - min_x) / span_x * (width - 2 * margin),
                    height - margin - (y_value - min_y) / span_y * (height - 2 * margin),
                ]
            )
        canvas.create_line(*coordinates, fill="#4fd1c5", width=2, smooth=False)
        canvas.create_oval(
            coordinates[0] - 4,
            coordinates[1] - 4,
            coordinates[0] + 4,
            coordinates[1] + 4,
            fill="#68d391",
            outline="",
        )
        canvas.create_oval(
            coordinates[-2] - 4,
            coordinates[-1] - 4,
            coordinates[-2] + 4,
            coordinates[-1] + 4,
            fill="#fc8181",
            outline="",
        )
        canvas.create_text(margin, 18, text=f"{projection} 绿色=起点 红色=终点", fill="#d9e2ec", anchor="w")
        canvas.create_text(margin, height - 18, text=f"{horizontal.upper()}: {min_x:.3f} … {max_x:.3f} m", fill="#9fb3c8", anchor="w")
        canvas.create_text(width - margin, 18, text=f"{vertical.upper()}: {min_y:.3f} … {max_y:.3f} m", fill="#9fb3c8", anchor="e")


def load_csv(path: Path) -> list[OrbitPoint]:
    points: list[OrbitPoint] = []
    with path.open("r", newline="", encoding="utf-8") as source:
        for row in csv.DictReader(source):
            points.append(
                OrbitPoint(
                    timestamp_ms=int(row["timestamp_ms"]),
                    x=float(row["x"]),
                    y=float(row["y"]),
                    z=float(row["z"]),
                    vx=float(row["vx"]),
                    vy=float(row["vy"]),
                    vz=float(row["vz"]),
                )
            )
    return points


if __name__ == "__main__":
    app = OrbitMqTester()
    if paramiko is None:
        app.after(
            200,
            lambda: messagebox.showwarning(
                "缺少依赖",
                "请先执行：\npy -3 -m pip install -r tools/requirements-windows.txt",
            ),
        )
    app.mainloop()

"""Cross-platform runtime and peak-memory measurement for EasyHap.

Runtime is wall-clock time measured with :func:`time.perf_counter` over the
actual analysis window. Peak memory is the maximum *simultaneous resident
memory* (RSS / Windows Working Set) of the EasyHap process tree observed during
that same window. The process tree includes child processes such as a bundled
``bcftools.exe`` backend.

The preferred implementation uses psutil when it is available.  On Windows and
Linux a standard-library/native fallback is provided so a PyInstaller build does
not silently report zero memory if psutil was omitted.  Windows native API calls
have explicit ctypes signatures, which is important for 64-bit executables.
"""
from __future__ import annotations

from dataclasses import dataclass
import os
import platform
import threading
import time
from typing import Optional


_GIB = 1024.0 ** 3

try:  # Preferred cross-platform process-tree implementation.
    import psutil as _psutil  # type: ignore
except Exception:  # pragma: no cover - exercised on minimal installations.
    _psutil = None


# ---------------------------------------------------------------------------
# Windows native memory helpers
# ---------------------------------------------------------------------------

def _windows_memory_info(pid: Optional[int] = None) -> tuple[int, int]:
    """Return (current_working_set, lifetime_peak_working_set) for *pid*.

    Returns (0, 0) if the process cannot be queried.  K32GetProcessMemoryInfo is
    used first (Windows 7+), with psapi.GetProcessMemoryInfo as a compatibility
    fallback.  Function signatures are declared explicitly to avoid pointer /
    HANDLE truncation in 64-bit frozen Python applications.
    """
    if os.name != "nt":
        return (0, 0)

    try:
        import ctypes
        from ctypes import wintypes

        SIZE_T = ctypes.c_size_t

        class PROCESS_MEMORY_COUNTERS_EX(ctypes.Structure):
            _fields_ = [
                ("cb", wintypes.DWORD),
                ("PageFaultCount", wintypes.DWORD),
                ("PeakWorkingSetSize", SIZE_T),
                ("WorkingSetSize", SIZE_T),
                ("QuotaPeakPagedPoolUsage", SIZE_T),
                ("QuotaPagedPoolUsage", SIZE_T),
                ("QuotaPeakNonPagedPoolUsage", SIZE_T),
                ("QuotaNonPagedPoolUsage", SIZE_T),
                ("PagefileUsage", SIZE_T),
                ("PeakPagefileUsage", SIZE_T),
                ("PrivateUsage", SIZE_T),
            ]

        kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
        kernel32.GetCurrentProcess.argtypes = []
        kernel32.GetCurrentProcess.restype = wintypes.HANDLE
        kernel32.OpenProcess.argtypes = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD]
        kernel32.OpenProcess.restype = wintypes.HANDLE
        kernel32.CloseHandle.argtypes = [wintypes.HANDLE]
        kernel32.CloseHandle.restype = wintypes.BOOL

        # QUERY_INFORMATION + VM_READ gives broad compatibility.  If that is
        # denied, retry with QUERY_LIMITED_INFORMATION on modern Windows.
        PROCESS_QUERY_INFORMATION = 0x0400
        PROCESS_VM_READ = 0x0010
        PROCESS_QUERY_LIMITED_INFORMATION = 0x1000

        current_pid = os.getpid()
        need_close = False
        if pid is None or int(pid) == current_pid:
            handle = kernel32.GetCurrentProcess()
        else:
            pid = int(pid)
            handle = kernel32.OpenProcess(
                PROCESS_QUERY_INFORMATION | PROCESS_VM_READ, False, pid
            )
            if not handle:
                handle = kernel32.OpenProcess(
                    PROCESS_QUERY_LIMITED_INFORMATION | PROCESS_VM_READ, False, pid
                )
            if not handle:
                return (0, 0)
            need_close = True

        counters = PROCESS_MEMORY_COUNTERS_EX()
        counters.cb = ctypes.sizeof(counters)

        ok = False
        try:
            k32 = getattr(kernel32, "K32GetProcessMemoryInfo", None)
            if k32 is not None:
                k32.argtypes = [
                    wintypes.HANDLE,
                    ctypes.POINTER(PROCESS_MEMORY_COUNTERS_EX),
                    wintypes.DWORD,
                ]
                k32.restype = wintypes.BOOL
                ok = bool(k32(handle, ctypes.byref(counters), counters.cb))

            if not ok:
                psapi = ctypes.WinDLL("psapi", use_last_error=True)
                get_mem = psapi.GetProcessMemoryInfo
                get_mem.argtypes = [
                    wintypes.HANDLE,
                    ctypes.POINTER(PROCESS_MEMORY_COUNTERS_EX),
                    wintypes.DWORD,
                ]
                get_mem.restype = wintypes.BOOL
                ok = bool(get_mem(handle, ctypes.byref(counters), counters.cb))
        finally:
            if need_close:
                kernel32.CloseHandle(handle)

        if ok:
            return (int(counters.WorkingSetSize), int(counters.PeakWorkingSetSize))
    except Exception:
        pass

    return (0, 0)


def _windows_descendant_pids(root_pid: int) -> list[int]:
    """Return root PID plus all descendants using Toolhelp32Snapshot."""
    if os.name != "nt":
        return [root_pid]

    try:
        import ctypes
        from ctypes import wintypes

        TH32CS_SNAPPROCESS = 0x00000002
        MAX_PATH = 260
        ULONG_PTR = ctypes.c_size_t

        class PROCESSENTRY32W(ctypes.Structure):
            _fields_ = [
                ("dwSize", wintypes.DWORD),
                ("cntUsage", wintypes.DWORD),
                ("th32ProcessID", wintypes.DWORD),
                ("th32DefaultHeapID", ULONG_PTR),
                ("th32ModuleID", wintypes.DWORD),
                ("cntThreads", wintypes.DWORD),
                ("th32ParentProcessID", wintypes.DWORD),
                ("pcPriClassBase", wintypes.LONG),
                ("dwFlags", wintypes.DWORD),
                ("szExeFile", wintypes.WCHAR * MAX_PATH),
            ]

        kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
        kernel32.CreateToolhelp32Snapshot.argtypes = [wintypes.DWORD, wintypes.DWORD]
        kernel32.CreateToolhelp32Snapshot.restype = wintypes.HANDLE
        kernel32.Process32FirstW.argtypes = [wintypes.HANDLE, ctypes.POINTER(PROCESSENTRY32W)]
        kernel32.Process32FirstW.restype = wintypes.BOOL
        kernel32.Process32NextW.argtypes = [wintypes.HANDLE, ctypes.POINTER(PROCESSENTRY32W)]
        kernel32.Process32NextW.restype = wintypes.BOOL
        kernel32.CloseHandle.argtypes = [wintypes.HANDLE]
        kernel32.CloseHandle.restype = wintypes.BOOL

        snapshot = kernel32.CreateToolhelp32Snapshot(TH32CS_SNAPPROCESS, 0)
        invalid = ctypes.c_void_p(-1).value
        if not snapshot or ctypes.cast(snapshot, ctypes.c_void_p).value == invalid:
            return [root_pid]

        parent_to_children: dict[int, list[int]] = {}
        try:
            entry = PROCESSENTRY32W()
            entry.dwSize = ctypes.sizeof(entry)
            ok = kernel32.Process32FirstW(snapshot, ctypes.byref(entry))
            while ok:
                pid = int(entry.th32ProcessID)
                ppid = int(entry.th32ParentProcessID)
                parent_to_children.setdefault(ppid, []).append(pid)
                entry.dwSize = ctypes.sizeof(entry)
                ok = kernel32.Process32NextW(snapshot, ctypes.byref(entry))
        finally:
            kernel32.CloseHandle(snapshot)

        found: list[int] = []
        stack = [int(root_pid)]
        seen: set[int] = set()
        while stack:
            pid = stack.pop()
            if pid in seen:
                continue
            seen.add(pid)
            found.append(pid)
            stack.extend(parent_to_children.get(pid, ()))
        return found
    except Exception:
        return [root_pid]


# ---------------------------------------------------------------------------
# Linux native process-tree fallback
# ---------------------------------------------------------------------------

def _linux_descendant_pids(root_pid: int) -> list[int]:
    if not (os.name == "posix" and os.path.isdir("/proc")):
        return [root_pid]

    found: list[int] = []
    stack = [int(root_pid)]
    seen: set[int] = set()
    while stack:
        pid = stack.pop()
        if pid in seen:
            continue
        seen.add(pid)
        if not os.path.exists(f"/proc/{pid}"):
            continue
        found.append(pid)
        children_file = f"/proc/{pid}/task/{pid}/children"
        try:
            with open(children_file, "r", encoding="ascii") as fh:
                stack.extend(int(x) for x in fh.read().split())
        except (OSError, ValueError):
            pass
    return found


def _linux_rss_for_pid(pid: int) -> int:
    try:
        with open(f"/proc/{pid}/statm", "r", encoding="ascii") as fh:
            fields = fh.readline().split()
        if len(fields) >= 2:
            return int(fields[1]) * int(os.sysconf("SC_PAGE_SIZE"))
    except (OSError, ValueError, AttributeError):
        pass
    return 0


# ---------------------------------------------------------------------------
# Cross-platform process-tree RSS
# ---------------------------------------------------------------------------

def _process_tree_rss_bytes() -> int:
    """Return simultaneous RSS/working-set bytes for EasyHap + descendants."""
    root_pid = os.getpid()

    # psutil is both robust and fast enough for a 10-ms benchmark sampler.
    if _psutil is not None:
        try:
            root = _psutil.Process(root_pid)
            processes = [root] + root.children(recursive=True)
            total = 0
            for proc in processes:
                try:
                    total += int(proc.memory_info().rss)
                except (_psutil.NoSuchProcess, _psutil.AccessDenied, _psutil.ZombieProcess):
                    continue
            if total > 0:
                return total
        except (_psutil.NoSuchProcess, _psutil.AccessDenied, _psutil.ZombieProcess):
            pass
        except Exception:
            pass

    # Native fallbacks.  Windows includes descendants so a bcftools.exe backend
    # is still counted even when psutil was not bundled into the EXE.
    if os.name == "nt":
        total = 0
        for pid in _windows_descendant_pids(root_pid):
            current, _peak = _windows_memory_info(pid)
            total += current
        if total > 0:
            return total

    if os.name == "posix" and os.path.isdir("/proc"):
        total = sum(_linux_rss_for_pid(pid) for pid in _linux_descendant_pids(root_pid))
        if total > 0:
            return total

    # Last-resort process-only measurement.
    peak = _lifetime_peak_rss_bytes()
    return int(peak or 0)


def _lifetime_peak_rss_bytes() -> Optional[int]:
    """Return this process's OS-maintained RSS/working-set high-water mark."""
    if os.name == "nt":
        _current, peak = _windows_memory_info(os.getpid())
        return peak if peak > 0 else None

    # Linux exposes VmHWM directly in bytes after conversion from KiB.
    if os.name == "posix" and os.path.exists("/proc/self/status"):
        try:
            with open("/proc/self/status", "r", encoding="ascii") as fh:
                for line in fh:
                    if line.startswith("VmHWM:"):
                        parts = line.split()
                        if len(parts) >= 2:
                            return int(parts[1]) * 1024
        except (OSError, ValueError):
            pass

    try:
        import resource

        value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        # Linux reports KiB; macOS/BSD commonly report bytes.
        if platform.system() == "Darwin":
            return int(value)
        return int(value * 1024.0)
    except Exception:
        return None


@dataclass(frozen=True)
class PerformanceMetrics:
    runtime_s: float
    peak_ram_bytes: int

    @property
    def peak_ram_gb(self) -> float:
        return self.peak_ram_bytes / _GIB


class PerformanceMonitor:
    """Measure one EasyHap analysis window.

    Peak RAM is the maximum sampled simultaneous RSS/working set of the EasyHap
    process tree.  The default 10-ms interval balances transient-peak capture
    against benchmark overhead.  The main process's native high-water mark is
    also retained when it establishes a new peak during the analysis window.
    """

    def __init__(self, sample_interval_s: float = 0.01) -> None:
        self.sample_interval_s = max(float(sample_interval_s), 0.001)
        self._stop = threading.Event()
        self._thread: Optional[threading.Thread] = None
        self._start_time = 0.0
        self._start_lifetime_peak: Optional[int] = None
        self._sampled_peak = 0
        self.metrics: Optional[PerformanceMetrics] = None

    def _sample_loop(self) -> None:
        while not self._stop.is_set():
            rss = _process_tree_rss_bytes()
            if rss > self._sampled_peak:
                self._sampled_peak = rss
            self._stop.wait(self.sample_interval_s)
        # Capture once more at the exact end of the analysis window.
        rss = _process_tree_rss_bytes()
        if rss > self._sampled_peak:
            self._sampled_peak = rss

    def start(self) -> "PerformanceMonitor":
        if self._thread is not None:
            raise RuntimeError("PerformanceMonitor has already been started")
        self._stop.clear()
        self._sampled_peak = _process_tree_rss_bytes()
        self._start_lifetime_peak = _lifetime_peak_rss_bytes()
        self._start_time = time.perf_counter()
        self._thread = threading.Thread(
            target=self._sample_loop,
            name="EasyHapPerformanceMonitor",
            daemon=True,
        )
        self._thread.start()
        return self

    def stop(self) -> PerformanceMetrics:
        if self._thread is None:
            raise RuntimeError("PerformanceMonitor has not been started")
        runtime_s = time.perf_counter() - self._start_time
        self._stop.set()
        self._thread.join()

        end_lifetime_peak = _lifetime_peak_rss_bytes()
        peak = self._sampled_peak

        # If the main process itself established a new native high-water mark,
        # preserve it.  This supplements the sampler for very brief spikes.
        if (
            end_lifetime_peak is not None
            and self._start_lifetime_peak is not None
            and end_lifetime_peak > self._start_lifetime_peak
        ):
            peak = max(peak, end_lifetime_peak)
        elif end_lifetime_peak is not None and self._start_lifetime_peak is None:
            peak = max(peak, end_lifetime_peak)

        if peak <= 0:
            # Never silently publish a scientifically misleading 0.000000 GB.
            raise RuntimeError(
                "EasyHap could not read process memory usage on this system. "
                "On Windows, rebuild the EXE with psutil installed "
                "(`pip install psutil`) or ensure Windows process-query APIs "
                "are not blocked by security software."
            )

        self.metrics = PerformanceMetrics(runtime_s=runtime_s, peak_ram_bytes=int(peak))
        return self.metrics

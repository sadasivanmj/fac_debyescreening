#!/usr/bin/env python3
"""Keep at most N sfac processes of the sweep running; pause the rest.

Memory guard for the large paper-active-space CI runs (~2.7 GB each on a
15 GB machine).  sfac processes whose working directory is under runs/ are
paused with SIGSTOP when more than N are running and resumed (oldest first)
with SIGCONT when slots free up.  Pausing does not change any result.
Exits when no sweep sfac process is left for 10 minutes.
"""

import os
import signal
import sys
import time
from pathlib import Path

RUNS = str(Path(__file__).resolve().parent / "runs")
N = int(sys.argv[1]) if len(sys.argv) > 1 else 4


def sweep_procs():
    out = []
    for p in Path("/proc").iterdir():
        if not p.name.isdigit():
            continue
        try:
            if p.joinpath("comm").read_text().strip() != "sfac":
                continue
            if not os.readlink(p / "cwd").startswith(RUNS):
                continue
            state = p.joinpath("stat").read_text().split(")")[1].split()[0]
            start = int(p.joinpath("stat").read_text().split(")")[1].split()[19])
            out.append((start, int(p.name), state))
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            continue
    return sorted(out)


idle_since = None
while True:
    procs = sweep_procs()
    running = [p for p in procs if p[2] != "T"]
    stopped = [p for p in procs if p[2] == "T"]
    if len(running) > N:
        for _, pid, _ in running[N:]:
            os.kill(pid, signal.SIGSTOP)
            print(time.strftime("%H:%M:%S"), "paused", pid, flush=True)
    elif stopped and len(running) < N:
        for _, pid, _ in stopped[: N - len(running)]:
            os.kill(pid, signal.SIGCONT)
            print(time.strftime("%H:%M:%S"), "resumed", pid, flush=True)
    if not procs:
        idle_since = idle_since or time.time()
        if time.time() - idle_since > 600:
            break
    else:
        idle_since = None
    time.sleep(3)

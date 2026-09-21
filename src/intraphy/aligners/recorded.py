"""Argument-array external calls with retained output, timeout and exit status."""
from __future__ import annotations
from datetime import datetime, timezone
import json
import os
import signal
from pathlib import Path
import shutil
import subprocess


def run_recorded(argv: list[str], output_dir: str | Path, label: str, timeout: int = 600) -> Path:
    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    if not label or Path(label).name != label:
        raise ValueError("An external-call label must be a simple filename")
    exe = shutil.which(argv[0])
    if not exe:
        raise FileNotFoundError(f"Required executable not found: {argv[0]}")
    if timeout < 1:
        raise ValueError("External timeout must be positive")
    command = [exe, *map(str, argv[1:])]
    stdout, stderr = directory / f"{label}.stdout", directory / f"{label}.stderr"
    record = {"command": command, "shell": False, "timeout_seconds": timeout,
              "started_at_utc": datetime.now(timezone.utc).isoformat(), "status": "running"}
    try:
        with stdout.open("w") as out, stderr.open("w") as err:
            process = subprocess.Popen(command, stdout=out, stderr=err, start_new_session=(os.name == "posix"))
            try:
                returncode = process.wait(timeout=timeout)
            except BaseException:
                if process.poll() is None:
                    if os.name == "posix":
                        os.killpg(process.pid, signal.SIGKILL)
                    else:
                        process.kill()
                process.wait()
                raise
            result = subprocess.CompletedProcess(command, returncode)
        record.update(returncode=result.returncode, status="completed" if result.returncode == 0 else "failed")
        if result.returncode:
            raise RuntimeError(f"{argv[0]} failed ({result.returncode}); see {stderr}")
    except subprocess.TimeoutExpired:
        record.update(status="timeout")
        raise RuntimeError(f"{argv[0]} timed out; outputs retained in {directory}") from None
    except BaseException as exc:
        record.update(status="failed", error=str(exc))
        raise
    finally:
        record["finished_at_utc"] = datetime.now(timezone.utc).isoformat()
        (directory / f"{label}.command.json").write_text(json.dumps(record, indent=2)+"\n")
    return stdout

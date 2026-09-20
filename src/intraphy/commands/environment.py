"""Execution provenance limited to software actually relevant to IntraPhy."""
from datetime import datetime, timezone
from importlib import metadata
import platform
import shutil
import subprocess
import sys

TOOLS = ("mafft", "minimap2", "miniprot", "lastz")


def inspect_tools():
    rows = []
    for name in TOOLS:
        path = shutil.which(name)
        row = {"tool": name, "path": path, "available": path is not None, "version": None}
        if path:
            try:
                result = subprocess.run([path, "--version"], capture_output=True, text=True,
                                        timeout=10, check=False)
                output = (result.stdout + "\n" + result.stderr).strip()
                row["version"] = output.splitlines()[0] if output else "unreported"
                row["version_exit_code"] = result.returncode
            except (OSError, subprocess.TimeoutExpired) as exc:
                row["version"] = "unavailable"
                row["version_error"] = str(exc)
        rows.append(row)
    return rows


def environment_report():
    from .. import __version__
    dependencies = {}
    for name in ("numpy", "scipy", "biopython", "networkx"):
        try:
            dependencies[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            dependencies[name] = "unavailable"
    return {"intraphy": __version__, "python": platform.python_version(),
            "python_executable": sys.executable, "platform": platform.platform(),
            "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
            "dependencies": dependencies, "alignment_tools": inspect_tools()}

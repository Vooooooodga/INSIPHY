"""Opt-in AGAT normalization. Repaired annotation is not independent evidence."""
from __future__ import annotations

import json
from pathlib import Path
import shutil
import subprocess

from ..preparation.annotation_index import load_annotation_index, clear_annotation_cache
from .resources import GFF_SUFFIXES, expand_files, resource_name

TOOL = "agat_convert_sp_gxf2gxf.pl"


def normalize_annotations(inputs: list[str], output_dir: str,
                          config: str | None = None, timeout: int = 600) -> list[dict]:
    executable = shutil.which(TOOL)
    if executable is None:
        raise RuntimeError(f"{TOOL} is not installed; AGAT is optional. See docs/inputs.md.")
    if timeout < 1:
        raise ValueError("AGAT timeout must be positive")
    root = Path(output_dir).resolve()
    root.mkdir(parents=True, exist_ok=True)
    records = []
    for path in expand_files(inputs, GFF_SUFFIXES):
        destination = root / (resource_name(path, GFF_SUFFIXES) + ".gff3")
        if destination.exists():
            raise ValueError(f"AGAT output already exists: {destination}")
        command = [executable, "--gff", str(path), "--output", str(destination)]
        if config:
            command += ["--config", str(Path(config).resolve())]
        try:
            result = subprocess.run(command, capture_output=True, text=True, timeout=timeout,
                                    cwd=root, check=False)
        except subprocess.TimeoutExpired as exc:
            def decoded(value):
                return value.decode(errors="replace") if isinstance(value, bytes) else (value or "")
            (root / (destination.stem + ".agat.log")).write_text(
                decoded(exc.stdout) + "\n" + decoded(exc.stderr) + f"\nTimed out after {timeout} seconds\n")
            records.append({"source_gff": str(path), "output_gff": str(destination), "command": command,
                            "returncode": None, "status": "timeout", "annotation_modified": True})
            (root / "agat_commands.json").write_text(json.dumps(records, indent=2) + "\n")
            raise RuntimeError(f"AGAT timed out for {path}; inspect {destination.stem}.agat.log") from exc
        (root / (destination.stem + ".agat.log")).write_text(result.stdout + "\n" + result.stderr)
        record = {"source_gff": str(path), "output_gff": str(destination), "command": command,
                  "returncode": result.returncode, "annotation_modified": True,
                  "interpretation": "normalized_annotation_not_independent_biological_validation"}
        records.append(record)
        (root / "agat_commands.json").write_text(json.dumps(records, indent=2) + "\n")
        if result.returncode or not destination.is_file():
            raise RuntimeError(f"AGAT failed for {path}; inspect {destination.stem}.agat.log")
        clear_annotation_cache()
        if not load_annotation_index(destination).rows:
            raise ValueError(f"AGAT produced an empty annotation: {destination}")
    return records

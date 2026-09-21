#!/usr/bin/env python3
"""Exercise the installed V19 wheel from outside the source tree; preserve output."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time
import xml.etree.ElementTree as ET


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    root = Path(args.output_dir).resolve()
    root.mkdir(parents=True, exist_ok=True)
    if any(root.iterdir()):
        raise SystemExit("Smoke output must be empty")
    commands = []
    def invoke(*argv):
        started = time.monotonic()
        command = [sys.executable, "-I", "-m", "intraphy", *map(str, argv)]
        result = subprocess.run(command, cwd=root, capture_output=True, text=True, timeout=600)
        index = len(commands)+1
        (root/f"command-{index:02d}.stdout").write_text(result.stdout)
        (root/f"command-{index:02d}.stderr").write_text(result.stderr)
        commands.append({"argv": command, "returncode": result.returncode, "seconds": time.monotonic()-started})
        (root/"commands.json").write_text(json.dumps(commands, indent=2)+"\n")
        if result.returncode:
            raise RuntimeError(f"Installed command failed: {command}\n{result.stderr}")
    invoke("--version")
    invoke("inspect-aligners")
    invoke("example-exons", "--scenario", "split_insertion", "--output-dir", "raw")
    (root/"raw/truth.json").unlink()
    inputs = ("--fasta", "raw", "--gff", "raw", "--species-tree", "raw/species_tree.nwk")
    invoke("check", *inputs)
    invoke("extract-loci", *inputs, "--flank", 50, "--output-dir", "portable")
    invoke("analyze", *inputs, "--output-dir", "result", "--threads", 2)
    invoke("infer-phylogeny", "--input-dir", "result", "--exon-configurations", "result/exon_configurations.jsonl",
           "--output-dir", "reloaded")
    invoke("exon-rate-template", "--output", "rates.json")
    invoke("infer-phylogeny", "--input-dir", "result", "--exon-configurations", "result/exon_configurations.jsonl",
           "--model", "exon-ctmc", "--exon-rates", "rates.json", "--expected-edits", "--output-dir", "ctmc")
    invoke("visualize", "--input-dir", "result/prepared_inputs", "--result-dir", "result", "--output-dir", "figures")
    invoke("explain", "--output-dir", "guide")
    data = json.loads((root/"result/exon_history.json").read_text())
    assert data["model"] == "exon_configuration_v2"
    edits = [e for u in data["units"] for e in u["views"]["evidence"]["events"] if e["support"] == "required"]
    assert len(edits) == 1 and edits[0]["operation"] == "dna_insertion", edits
    assert "exon_split" in edits[0]["consequences"]
    diag = json.loads((root/"ctmc/model_diagnostics.json").read_text())
    assert all(r.get("probability_status") == "conditional_on_fixed_parameters_and_declared_catalogue" for r in diag["units"]), diag
    svgs = list(root.rglob("*.svg"))
    for path in svgs:
        ET.parse(path)
    assert len(svgs) >= 2
    summary = {"commands": len(commands), "svg_files": len(svgs), "status": "passed",
               "validation_scope": "installed_v19_package_synthetic_raw_input_and_conditional_ctmc"}
    (root/"summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

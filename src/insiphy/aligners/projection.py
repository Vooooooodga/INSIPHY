"""aligners / projection: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from insiphy.aligners.formats import _parse_paf_tags
from insiphy.aligners.formats import _write_temp_fasta
from insiphy.aligners.formats import _write_temp_fasta_records
from insiphy.aligners.types import AlignmentBackendError
from pathlib import Path
from urllib.parse import unquote
import shutil
import subprocess
import tempfile


def _parse_gff_attributes(raw: str) -> dict[str, str]:
    attrs = {}
    for item in raw.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        elif " " in item:
            key, value = item.split(" ", 1)
        else:
            continue
        attrs[unquote(key)] = unquote(value.strip().strip('"'))
    return attrs


def protein_locus_exons(proteins: dict[str, str], locus: str, threads=1) -> list[dict]:
    """Return miniprot CDS exon evidence rows for a batch of proteins against one locus."""
    exe = shutil.which("miniprot")
    if not exe:
        raise AlignmentBackendError("miniprot was requested but is not available on PATH")
    clean_proteins = {str(name): (seq or "").upper() for name, seq in proteins.items() if seq}
    locus = (locus or "").upper()
    if not clean_proteins or not locus:
        return []

    with tempfile.TemporaryDirectory(prefix="insiphy_miniprot_gff_") as tmp:
        tmp = Path(tmp)
        protein_path = tmp / "proteins.faa"
        locus_path = tmp / "locus.fa"
        _write_temp_fasta_records(protein_path, clean_proteins)
        _write_temp_fasta(locus_path, "locus", locus)
        cmd = [exe, "-t", str(max(1, int(threads or 1))), "--gff", str(locus_path), str(protein_path)]
        proc = subprocess.run(cmd, text=True, capture_output=True, check=False)
    if proc.returncode != 0:
        raise AlignmentBackendError(proc.stderr.strip() or "miniprot --gff failed")

    parent_info: dict[str, dict] = {}
    cds_rows = []
    paf_cigars: dict[str, list[str]] = {}
    for line in proc.stdout.splitlines():
        if not line.strip():
            continue
        if line.startswith("##PAF"):
            paf_fields = line.split("\t")[1:]
            if len(paf_fields) >= 12:
                tags = _parse_paf_tags(paf_fields)
                cigar = tags.get("cg")
                if cigar:
                    paf_cigars.setdefault(paf_fields[0], []).append(cigar)
            continue
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 9:
            continue
        seqid, _source, feature, start, end, _score, strand, phase, raw_attrs = fields
        attrs = _parse_gff_attributes(raw_attrs)
        if feature == "mRNA":
            parent_id = attrs.get("ID", "")
            target = attrs.get("Target", "")
            target_parts = target.split()
            protein_id = target_parts[0] if target_parts else attrs.get("ID", "")
            q_start = int(target_parts[1]) if len(target_parts) >= 3 and target_parts[1].isdigit() else 1
            q_end = int(target_parts[2]) if len(target_parts) >= 3 and target_parts[2].isdigit() else 0
            query_coverage = (q_end - q_start + 1) / max(1, len(clean_proteins.get(protein_id, ""))) if q_end else 0.0
            parent_info[parent_id] = {
                "protein_id": protein_id,
                "alignment_protein_id": protein_id,
                "query_start": q_start,
                "query_end": q_end,
                "identity": float(attrs.get("Identity", 0.0) or 0.0),
                "query_coverage": query_coverage,
                "cigar": "NA",
            }
        elif feature == "CDS":
            cds_rows.append((seqid, int(start), int(end), strand, phase, attrs))

    for info in parent_info.values():
        cigars = paf_cigars.get(info["protein_id"], [])
        if len(cigars) == 1:
            info["cigar"] = cigars[0]

    rows = []
    for _seqid, start, end, strand, phase, attrs in cds_rows:
        parent_id = attrs.get("Parent", "")
        info = parent_info.get(parent_id, {})
        cds_target = attrs.get("Target", "").split()
        has_cds_target = len(cds_target) >= 3 and cds_target[1].isdigit() and cds_target[2].isdigit()
        protein_id = cds_target[0] if cds_target else info.get("protein_id", parent_id)
        query_start = int(cds_target[1]) if has_cds_target else "NA"
        query_end = int(cds_target[2]) if has_cds_target else "NA"
        if attrs.get("Identity", "") not in {"", None}:
            identity = float(attrs.get("Identity", 0.0) or 0.0)
            identity_scope = "cds"
        else:
            identity = info.get("identity", 0.0)
            identity_scope = "parent_alignment"
        rows.append({
            "protein_id": protein_id,
            "alignment_protein_id": info.get("alignment_protein_id", protein_id),
            "query_start": query_start,
            "query_end": query_end,
            "parent_query_start": info.get("query_start", "NA"),
            "parent_query_end": info.get("query_end", "NA"),
            "target_start": start,
            "target_end": end,
            "strand": strand,
            "identity": identity,
            "identity_scope": identity_scope,
            "query_coverage": info.get("query_coverage", 0.0),
            "parent_query_coverage": info.get("query_coverage", 0.0),
            "parent_id": parent_id,
            "phase": phase,
            "cigar": info.get("cigar", "NA"),
            "projection_status": "cds_target" if has_cds_target else "missing_cds_target",
        })
    return rows

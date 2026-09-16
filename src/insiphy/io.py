"""TSV and FASTA helpers used by INSIPHY."""

import csv
import gzip
from pathlib import Path


UNKNOWN = {"", "NA", "unknown", "ambiguous", "truncated", "unresolved"}


def open_text(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def read_tsv(path, required=None, optional=False):
    path = Path(path)
    if optional and not path.exists():
        return []
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = [field for field in (required or []) if field not in (reader.fieldnames or [])]
        if missing:
            raise SystemExit(f"{path} missing required fields: {','.join(missing)}")
        return list(reader)


def write_tsv(path, rows, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "NA") for field in fields})


def parse_fasta(path):
    path = Path(path)
    if not path.exists():
        return {}
    records = {}
    name = None
    seq = []
    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(seq).upper()
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
    if name is not None:
        records[name] = "".join(seq).upper()
    return records


def read_fasta_record(path, record_id):
    """Read one FASTA record without retaining the remaining assembly in memory."""
    path = Path(path)
    sequence = []
    collecting = False
    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                if collecting:
                    break
                collecting = current == record_id
            elif collecting:
                sequence.append(line)
    if not sequence:
        raise SystemExit(f"FASTA record not found: {record_id} in {path}")
    return "".join(sequence).upper()


def norm_state(value):
    value = (value or "unknown").strip()
    return "unknown" if value in UNKNOWN else value


def to_float(value, default=0.0):
    try:
        if value in ("", "NA", None):
            return default
        return float(value)
    except ValueError:
        return default


def uniq(values):
    return len({v for v in values if v not in ("", "NA", None)})

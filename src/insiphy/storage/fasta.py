"""storage / fasta: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from functools import lru_cache
from insiphy.storage.tabular import open_text
from pathlib import Path


def _fai_path(path):
    return Path(f"{path}.fai")


@lru_cache(maxsize=64)
def _read_fai(path):
    path = Path(path)
    index_path = _fai_path(path)
    if not index_path.exists():
        return None
    records = {}
    with index_path.open() as handle:
        for raw in handle:
            if not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            name, length, offset, line_bases, line_width = fields[:5]
            records[name] = {
                "length": int(length),
                "offset": int(offset),
                "line_bases": int(line_bases),
                "line_width": int(line_width),
            }
    return records


def _indexed_interval(path, record_id, start, end):
    path = Path(path)
    if path.suffix == ".gz":
        return None
    fai = _read_fai(str(path))
    if not fai or record_id not in fai:
        return None
    rec = fai[record_id]
    if start < 1 or end > rec["length"] or start > end:
        raise SystemExit(f"invalid FASTA interval for {record_id}: {start}-{end}")
    zero_start = start - 1
    zero_end = end
    first_line = zero_start // rec["line_bases"]
    first_col = zero_start % rec["line_bases"]
    byte_start = rec["offset"] + first_line * rec["line_width"] + first_col
    last_base = zero_end - 1
    last_line = last_base // rec["line_bases"]
    last_col = last_base % rec["line_bases"]
    byte_end = rec["offset"] + last_line * rec["line_width"] + last_col + 1
    read_len = byte_end - byte_start
    with path.open("rb") as handle:
        handle.seek(byte_start)
        chunk = handle.read(read_len)
    sequence = chunk.replace(b"\n", b"").replace(b"\r", b"")
    return sequence[: end - start + 1].decode("ascii").upper()


def _stream_fasta_interval(path, record_id, start, end):
    if start < 1 or start > end:
        raise SystemExit(f"invalid FASTA interval for {record_id}: {start}-{end}")
    collecting = False
    position = 0
    pieces = []
    found = False
    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if collecting:
                    break
                current = line[1:].split()[0]
                collecting = current == record_id
                found = found or collecting
                position = 0
                continue
            if not collecting:
                continue
            line_start = position + 1
            line_end = position + len(line)
            if line_end >= start and line_start <= end:
                left = max(start, line_start) - line_start
                right = min(end, line_end) - line_start + 1
                pieces.append(line[left:right])
            position = line_end
            if position >= end:
                break
    if not found:
        raise SystemExit(f"FASTA record not found: {record_id} in {path}")
    sequence = "".join(pieces)
    if len(sequence) != end - start + 1:
        raise SystemExit(f"FASTA interval exceeds record length for {record_id}: {start}-{end}")
    return sequence.upper()


def read_fasta_interval(path, record_id, start, end):
    """Read a 1-based inclusive FASTA interval without indexing side effects."""
    path = Path(path)
    if not path.exists():
        raise SystemExit(f"FASTA not found: {path}")
    indexed = _indexed_interval(path, record_id, int(start), int(end))
    if indexed is not None:
        return indexed
    return _stream_fasta_interval(path, record_id, int(start), int(end))


@lru_cache(maxsize=256)
def fasta_record_length(path, record_id):
    """Return one FASTA record length using an existing .fai when available."""
    path = Path(path)
    fai = _read_fai(str(path))
    if fai and record_id in fai:
        return fai[record_id]["length"]
    found = False
    length = 0
    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if found:
                    break
                found = line[1:].split()[0] == record_id
                continue
            if found:
                length += len(line)
    if not found:
        raise SystemExit(f"FASTA record not found: {record_id} in {path}")
    return length


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
    length = fasta_record_length(str(Path(path)), record_id)
    return read_fasta_interval(path, record_id, 1, length)

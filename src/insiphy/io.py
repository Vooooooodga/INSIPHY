"""TSV and FASTA helpers used by INSIPHY."""

import csv
import gzip
from collections import Counter
from functools import lru_cache
from pathlib import Path


UNKNOWN = {
    "",
    ".",
    "NA",
    "?",
    "unknown",
    "ambiguous",
    "truncated",
    "unavailable",
    "unresolved",
}

STRUCTURAL_SITE_SCHEMA_VERSION = "3"
STRUCTURAL_SITE_SCHEMA_VERSIONS = {STRUCTURAL_SITE_SCHEMA_VERSION, "2", "1-legacy"}
STRUCTURAL_SITE_ANNOTATION_VIEWS = {
    "repertoire",
    "canonical",
    "view_independent",
}
LEGACY_STRUCTURAL_SITE_SCHEMA_VERSION = "1-legacy"
STRUCTURAL_SITE_OBSERVATION_MASKS = {
    "observed",
    "missing",
    "generated_all_zero",
    "explicit_all_zero",
    "masked",
    "excluded",
    "unobserved",
    "unknown",
}
MISSING_OBSERVATION_MASKS = {"missing", "masked", "excluded", "unobserved", "unknown"}
STRUCTURAL_SITE_STATE_LABELS = {
    "exon_presence": ("absent", "present"),
    "exon_role": ("not_exonic", "exonic"),
    "splice_junction": ("absent", "present"),
}
STRUCTURAL_SITE_CONSTANT_FIELDS = (
    "linked_group_id",
    "site_kind",
    "discovery_rule",
    "annotation_view",
)
STRUCTURAL_SITE_FIELDS = [
    "family_id",
    "layer",
    "site_id",
    "species",
    "state",
    "state_0",
    "state_1",
    "evidence",
    "schema_version",
    "annotation_view",
    "applicability",
    "observation_reason",
    "transcript_scope",
    "parent_feature_ids",
    "member_interval_ids",
    "evidence_ids",
    "discovery_rule",
    "discovery_species",
    "observation_mask",
    "linked_group_id",
    "site_kind",
    "observation_source",
    "annotation_completeness",
    "confidence_flag",
    "conclusion_flag",
]


def open_text(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


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


def normalize_structural_site_row(row):
    """Adapt a structural observation to the shared schema without hardening unknowns."""

    normalized = dict(row)
    state = norm_state(normalized.get("state"))
    layer = normalized.get("layer", "NA")
    evidence = normalized.get("evidence", "NA")
    normalized["state"] = state
    normalized.setdefault("schema_version", STRUCTURAL_SITE_SCHEMA_VERSION)
    if "applicability" not in normalized:
        if "not_applicable" in str(evidence):
            normalized["applicability"] = "inapplicable"
        elif state == "unknown":
            normalized["applicability"] = "undetermined"
        else:
            normalized["applicability"] = "applicable"
    normalized.setdefault("observation_reason", evidence)
    if not normalized.get("transcript_scope"):
        view = str(normalized.get("annotation_view", "")).strip().lower()
        normalized["transcript_scope"] = (
            "gene_locus" if layer == "exon_presence"
            else "explicit_canonical_transcript" if view == "canonical"
            else "annotated_transcript_repertoire"
        )
    normalized.setdefault("parent_feature_ids", "NA")
    normalized.setdefault("member_interval_ids", "NA")
    normalized.setdefault("evidence_ids", "NA")
    normalized.setdefault("discovery_rule", "derived_from_element_correspondence")
    normalized.setdefault("discovery_species", "NA")
    if not normalized.get("observation_mask"):
        normalized["observation_mask"] = "missing" if state == "unknown" else "observed"
    normalized["observation_mask"] = normalize_observation_mask(normalized["observation_mask"])
    normalized.setdefault("linked_group_id", "NA")
    linked_group_ids = {
        token
        for token in str(normalized["linked_group_id"]).split(";")
        if token not in {"", "NA"}
    }
    normalized["linked_group_id"] = ";".join(sorted(linked_group_ids)) or "NA"
    normalized.setdefault("site_kind", "unspecified")
    normalized.setdefault("observation_source", "genome_and_supplied_annotation")
    normalized.setdefault("annotation_completeness", "unassessed")
    normalized.setdefault("confidence_flag", "unassessed")
    normalized.setdefault("conclusion_flag", "unassessed")
    return normalized


def _validate_structural_site_rows(rows):
    seen = set()
    labels_by_site = {}
    metadata_by_site = {}
    for row in rows:
        key = tuple(
            row.get(field, "NA")
            for field in ("family_id", "layer", "site_id", "species")
        )
        if key in seen:
            raise ValueError("duplicate structural observation: " + "/".join(key))
        seen.add(key)
        state_0 = row.get("state_0", "NA")
        state_1 = row.get("state_1", "NA")
        if state_0 in UNKNOWN or state_1 in UNKNOWN or state_0 == state_1:
            raise ValueError(
                f"invalid structural-state labels for {key[0]}/{key[1]}/{key[2]}: "
                f"{state_0},{state_1}"
            )
        state = norm_state(row.get("state"))
        if state not in {"unknown", state_0, state_1}:
            raise ValueError(
                f"state {state!r} is outside the declared labels for "
                f"{key[0]}/{key[1]}/{key[2]}/{key[3]}"
            )
        site_key = key[:3]
        labels = (state_0, state_1)
        expected_labels = STRUCTURAL_SITE_STATE_LABELS.get(key[1])
        if expected_labels is not None and labels != expected_labels:
            raise ValueError(
                f"noncanonical structural-state labels for "
                f"{site_key[0]}/{site_key[1]}/{site_key[2]}: "
                f"expected {expected_labels[0]},{expected_labels[1]}; "
                f"observed {state_0},{state_1}"
            )
        previous = labels_by_site.setdefault(site_key, labels)
        if previous != labels:
            raise ValueError(
                f"inconsistent structural-state labels for "
                f"{site_key[0]}/{site_key[1]}/{site_key[2]}"
            )
        metadata = tuple(str(row.get(field, "NA")) for field in STRUCTURAL_SITE_CONSTANT_FIELDS)
        previous_metadata = metadata_by_site.setdefault(site_key, metadata)
        if previous_metadata != metadata:
            differing = [
                field
                for field, old, new in zip(
                    STRUCTURAL_SITE_CONSTANT_FIELDS, previous_metadata, metadata
                )
                if old != new
            ]
            raise ValueError(
                "inconsistent structural-site metadata for "
                f"{site_key[0]}/{site_key[1]}/{site_key[2]}: "
                + ",".join(differing)
            )
        applicability = row.get("applicability")
        if applicability not in {"applicable", "inapplicable", "undetermined"}:
            raise ValueError(
                f"invalid applicability for {key[0]}/{key[1]}/{key[2]}/{key[3]}: "
                f"{applicability!r}"
            )
        observation_mask = str(row.get("observation_mask", "")).strip().lower()
        if observation_mask not in STRUCTURAL_SITE_OBSERVATION_MASKS:
            raise ValueError(
                f"invalid observation_mask for {'/'.join(key)}: {observation_mask!r}"
            )
        if applicability == "inapplicable" and (
            state != "unknown" or observation_mask != "missing"
        ):
            raise ValueError(
                "inapplicable observations must have state=unknown and "
                f"observation_mask=missing for {'/'.join(key)}"
            )
        if observation_mask == "observed" and state not in {state_0, state_1}:
            raise ValueError(
                f"observed structural-site row lacks a declared state for {'/'.join(key)}"
            )
        if state == "unknown" and observation_mask == "observed":
            raise ValueError(
                f"unknown structural-site state cannot be observed for {'/'.join(key)}"
            )
        if observation_mask == "missing":
            reason = row.get("observation_reason")
            if reason is None or str(reason).strip() in UNKNOWN:
                raise ValueError(
                    f"missing structural-site observation requires observation_reason for {'/'.join(key)}"
                )
        if observation_mask in {"generated_all_zero", "explicit_all_zero"} and state != state_0:
            raise ValueError(
                f"all-zero structural-site mask requires state_0 for {'/'.join(key)}"
            )


def _structural_site_layer(row):
    return str(row.get("site_layer", row.get("layer", ""))).strip().lower()


def _scope_annotation_view(row):
    scope = str(row.get("transcript_scope", "")).strip().lower()
    if scope in {"canonical", "canonical_transcript", "explicit_canonical_transcript"}:
        return "canonical"
    if scope in {
        "repertoire",
        "transcript_repertoire",
        "annotated_transcript_repertoire",
        "all_transcripts",
    }:
        return "repertoire"
    return ""


def structural_site_annotation_view(row):
    """Resolve a row's view and reject contradictions between view and transcript scope."""
    layer = _structural_site_layer(row)
    explicit = str(row.get("annotation_view", "")).strip().lower()
    scope_view = _scope_annotation_view(row)
    is_presence = layer in {"presence", "sequence_presence", "element_presence", "exon_presence"}

    if explicit and explicit not in STRUCTURAL_SITE_ANNOTATION_VIEWS:
        raise ValueError(f"invalid annotation_view: {explicit!r}")
    if is_presence:
        if explicit and explicit != "view_independent":
            raise ValueError("Presence observations require annotation_view=view_independent")
        # Presence is defined at the homologous sequence-unit/gene-locus
        # level.  ``gene_locus`` is the canonical scope synthesized for
        # legacy presence rows and carries no transcript-specific view.
        if scope_view and str(row.get("transcript_scope", "")).strip().lower() not in {
            "gene_locus",
            "gene",
            "locus",
        }:
            raise ValueError("Presence observations cannot have a transcript-specific transcript_scope")
        return "view_independent"
    if explicit == "view_independent":
        raise ValueError("Role and junction rows require canonical or repertoire annotation_view")
    if explicit and scope_view and explicit != scope_view:
        raise ValueError(
            f"annotation_view={explicit} conflicts with transcript_scope={row.get('transcript_scope')}"
        )
    return explicit or scope_view


def normalize_observation_mask(value):
    """Normalize retained legacy missing masks and reject undefined mask labels."""
    mask = str(value or "").strip().lower()
    if mask in MISSING_OBSERVATION_MASKS:
        return "missing"
    if mask not in STRUCTURAL_SITE_OBSERVATION_MASKS:
        raise ValueError(f"invalid observation_mask: {mask!r}")
    return mask


def _normalize_structural_site_contract(row, *, upgrade_schema=True):
    schema_version = str(row.get("schema_version", "")).strip()
    if schema_version and schema_version not in STRUCTURAL_SITE_SCHEMA_VERSIONS:
        raise ValueError(f"unsupported structural-site schema_version: {schema_version!r}")
    layer = _structural_site_layer(row)
    is_presence = layer in {"presence", "sequence_presence", "element_presence", "exon_presence"}
    if schema_version in {"2", LEGACY_STRUCTURAL_SITE_SCHEMA_VERSION} and (
        not is_presence and not _scope_annotation_view(row)
    ):
        raise ValueError(
            "legacy role and junction rows require an inferable transcript_scope"
        )
    explicit_view = str(row.get("annotation_view", "")).strip()
    if schema_version == STRUCTURAL_SITE_SCHEMA_VERSION and not explicit_view:
        raise ValueError("schema-v3 structural-site rows require annotation_view")
    annotation_view = structural_site_annotation_view(row)
    if annotation_view not in STRUCTURAL_SITE_ANNOTATION_VIEWS:
        raise ValueError(
            "Structural-site row lacks an unambiguous annotation_view; "
            "legacy role and junction rows require an inferable transcript_scope"
        )

    state = norm_state(row.get("state"))
    row["state"] = state
    raw_mask = str(row.get("observation_mask") or "").strip()
    mask = normalize_observation_mask(raw_mask) if raw_mask else ""
    if mask:
        row["observation_mask"] = mask
        if (
            mask == "missing"
            and raw_mask.lower() in (MISSING_OBSERVATION_MASKS - {"missing"})
            and str(row.get("observation_reason", "")).strip() in UNKNOWN
        ):
            row["observation_reason"] = "legacy_missing_mask"
    applicability = str(row.get("applicability", "")).strip().lower()
    state_0 = str(row.get("state_0", row.get("state_0_label", ""))).strip()
    state_1 = str(row.get("state_1", row.get("state_1_label", ""))).strip()

    if applicability == "inapplicable" and (state != "unknown" or mask != "missing"):
        raise ValueError("inapplicable observations must have state=unknown and observation_mask=missing")
    if mask == "observed" and state not in {state_0, state_1}:
        raise ValueError("observed structural-site rows must contain state_0 or state_1")
    if state == "unknown" and mask == "observed":
        raise ValueError("unknown structural-site states cannot be marked observed")
    if mask in {"generated_all_zero", "explicit_all_zero"} and state != state_0:
        raise ValueError("all-zero structural-site masks require state_0")

    if upgrade_schema:
        row["schema_version"] = STRUCTURAL_SITE_SCHEMA_VERSION
    row["annotation_view"] = annotation_view
    return row


def write_structural_site_matrix(path, rows):
    prepared = [_normalize_structural_site_contract(dict(row)) for row in rows]
    normalized = [normalize_structural_site_row(row) for row in prepared]
    normalized = [
        _normalize_structural_site_contract(row, upgrade_schema=False)
        for row in normalized
    ]
    _validate_structural_site_rows(normalized)
    normalized.sort(
        key=lambda row: tuple(
            str(row.get(field, "NA"))
            for field in ("family_id", "layer", "site_id", "species")
        )
    )
    write_tsv(path, normalized, STRUCTURAL_SITE_FIELDS)


def read_structural_site_matrix(path, optional=False):
    rows = read_tsv(
        path,
        ["family_id", "layer", "site_id", "species", "state", "state_0", "state_1"],
        optional=optional,
    )
    normalized = []
    for row in rows:
        if not str(row.get("schema_version", "")).strip():
            row = {**row, "schema_version": LEGACY_STRUCTURAL_SITE_SCHEMA_VERSION}
        prepared = _normalize_structural_site_contract(dict(row), upgrade_schema=False)
        normalized_row = normalize_structural_site_row(prepared)
        normalized.append(_normalize_structural_site_contract(normalized_row, upgrade_schema=False))
    _validate_structural_site_rows(normalized)
    return normalized


def validate_structural_site_tip_rows(rows, tip_labels):
    """Require an explicit matrix row for every selected tree tip and site."""

    expected = set(tip_labels)
    by_site = {}
    for row in rows:
        key = tuple(row.get(field, "NA") for field in ("family_id", "layer", "site_id"))
        by_site.setdefault(key, Counter())[row.get("species", "NA")] += 1

    violations = []
    for key, counts in sorted(by_site.items()):
        observed = set(counts)
        missing = sorted(expected - observed)
        extra = sorted(observed - expected)
        duplicates = sorted(species for species, count in counts.items() if count != 1)
        if missing:
            violations.append("/".join(key) + ":missing=" + ",".join(missing))
        if extra:
            violations.append("/".join(key) + ":outside_tree=" + ",".join(extra))
        if duplicates:
            violations.append(
                "/".join(key) + ":duplicate="
                + ",".join(f"{species}x{counts[species]}" for species in duplicates)
            )
    if violations:
        preview = "; ".join(violations[:8])
        suffix = "; ..." if len(violations) > 8 else ""
        raise SystemExit(
            "structural-site matrix must contain exactly one explicit row for every "
            "selected tree tip and site: " + preview + suffix
        )


def structural_site_observed_state(row):
    """Return the declared state only when the normalized observation mask consumes it."""

    raw_mask = row.get("observation_mask")
    if raw_mask in {None, ""}:
        raw_mask = "missing" if norm_state(row.get("state")) == "unknown" else "observed"
    mask = normalize_observation_mask(raw_mask)
    if mask == "missing":
        return "unknown"
    return row.get("state", "unknown")


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

"""Project one family protein MSA onto transcript-specific CDS coordinates.

Every annotated transcript remains an independent biological path. Identical
proteins are deduplicated only for the MSA call and retain all transcript
aliases. Protein identity and local coverage are reported as measurements;
neither is a universal homology gate. No protein projection assigns an exonic
role by itself.
"""

import json
from collections import OrderedDict, defaultdict
from dataclasses import dataclass
from urllib.parse import unquote

from Bio.Align import substitution_matrices
from Bio.Data import CodonTable

from . import alignment
from .coordinates import CoordinateBlock, Interval0


_STANDARD_CODE = CodonTable.unambiguous_dna_by_id[1]
_BLOSUM62 = substitution_matrices.load("BLOSUM62")
_EMPTY = {"", ".", "NA", None}
_ANCHOR_WINDOW = 15
_MIN_ANCHOR_PAIRS = 8


@dataclass(frozen=True)
class CodingBaseProjection:
    """One CDS base with transcript, occurrence and genomic provenance."""

    occurrence_id: str
    occurrence_interval: Interval0
    cds_interval: Interval0
    contig: str
    genome_interval: Interval0
    strand: str
    source_path_id: str
    source_feature_ids: tuple


@dataclass(frozen=True)
class CodingResidueProjection:
    """One translated residue placed on the family MSA coordinate system."""

    transcript_key: tuple
    residue_index0: int
    msa_column0: int
    amino_acid: str
    cds_bases: tuple


@dataclass
class CodingTranscript:
    key: tuple
    protein: str
    codons: tuple
    coding_lengths: dict
    unavailable_reason: str = ""
    codon_sources: tuple = ()
    source_ids: tuple = ()


@dataclass(frozen=True)
class FamilyCodingProjection:
    """A family MSA plus the transcript aliases represented by each record."""

    family_id: str
    mode: str
    aligned_records: dict
    aliases_by_record: dict
    record_by_transcript: dict
    residue_columns: dict

    def aligned_for(self, transcript_key):
        return self.aligned_records[self.record_by_transcript[transcript_key]]

    def columns_for(self, transcript_key):
        return self.residue_columns[self.record_by_transcript[transcript_key]]


@dataclass(frozen=True)
class CodingProjectionCandidate:
    """One transcript-pair projection for a requested occurrence pair."""

    query_transcript_key: tuple
    target_transcript_key: tuple
    coordinate_blocks: tuple
    aa_identity: float
    query_cds_coverage: float
    target_cds_coverage: float
    known_aa_pairs: int
    blosum62_score: float
    gap_fraction: float
    left_anchor_pairs: int
    right_anchor_pairs: int
    left_anchor_score: float
    right_anchor_score: float
    left_anchor_supported: bool
    right_anchor_supported: bool
    terminal_side: str
    msa_column_interval: Interval0
    anchor_resolved: bool
    position_monotonic: bool
    competing_occurrences: tuple

    @property
    def query_transcript_id(self):
        return self.query_transcript_key[-1]

    @property
    def target_transcript_id(self):
        return self.target_transcript_key[-1]

    @property
    def candidate_id(self):
        query = ":".join(str(value) for value in self.query_transcript_key)
        target = ":".join(str(value) for value in self.target_transcript_key)
        return f"protein:{query}>{target}"

    @property
    def position_eligible(self):
        return (
            self.anchor_resolved
            and self.position_monotonic
            and not self.competing_occurrences
        )


@dataclass(frozen=True)
class CodingProjectionCandidateSet:
    """All independent transcript-pair projections for one occurrence pair."""

    candidates: tuple

    @property
    def candidate_evidence_available(self):
        return bool(self.candidates)

    @property
    def coordinate_consensus(self):
        return bool(self.candidates) and len({
            candidate.coordinate_blocks for candidate in self.candidates
        }) == 1

    @property
    def competing_occurrences(self):
        return tuple(sorted({
            occurrence
            for candidate in self.candidates
            for occurrence in candidate.competing_occurrences
        }))

    @property
    def position_eligible(self):
        return (
            self.coordinate_consensus
            and not self.competing_occurrences
            and any(candidate.position_eligible for candidate in self.candidates)
        )

    @property
    def membership_eligible(self):
        """Compatibility name for hard, coordinate-resolved membership."""

        return self.position_eligible


def _copy_key(row):
    return row.get("family_id", ""), row.get("species", ""), row.get("gene_copy_id", "")


def _transcript_features(rows, transcript_id):
    by_id = defaultdict(list)
    children = defaultdict(list)
    for row in rows:
        if row.get("id") not in _EMPTY:
            by_id[row["id"]].append(row)
        for parent in str(row.get("parent", "") or "").split(","):
            parent = parent.strip()
            if parent not in _EMPTY:
                children[parent].append(row)
    selected = []
    visited = set()
    pending = [transcript_id]
    while pending:
        feature_id = pending.pop()
        if feature_id in _EMPTY or feature_id in visited:
            continue
        visited.add(feature_id)
        selected.extend(by_id.get(feature_id, []))
        pending.extend(row["id"] for row in children.get(feature_id, []) if row.get("id") not in _EMPTY)
        selected.extend(row for row in children.get(feature_id, []) if row.get("id") in _EMPTY)
    # Ancestor attributes can specify a translation table or a recoding event.
    pending = [row for row in by_id.get(transcript_id, [])]
    while pending:
        row = pending.pop()
        for parent in str(row.get("parent", "") or "").split(","):
            parent = parent.strip()
            if parent not in _EMPTY and parent not in visited:
                visited.add(parent)
                ancestors = by_id.get(parent, [])
                selected.extend(ancestors)
                pending.extend(ancestors)
    return selected


def _translation_reason(features):
    for feature in features:
        attrs = feature.get("attrs", {})
        for key, value in attrs.items():
            key = key.lower()
            value = unquote(str(value)).strip()
            if key in {"transl_table", "translation_table"} and value != "1":
                return "unsupported_translation_table"
            if key in {"transl_except", "translation_exception"} and value not in _EMPTY:
                return "translation_exception"
            if key == "exception" and any(word in value.lower() for word in (
                "ribosomal", "slippage", "frameshift", "recod", "selenocyst",
                "pyrrolys", "editing", "translation",
            )):
                return "translation_exception"
    return ""


def _synthetic_transcript_features(rows, transcript_id, paths, occurrences):
    """Resolve a synthetic path only through owned gene/exon/CDS coordinates."""
    owned = [row for row in rows if row.get("ownership") == "target_gene_descendant"]
    transcript_types = {"mrna", "transcript", "lnc_rna", "ncrna", "rrna", "trna"}
    if any(row.get("type", "").lower() in transcript_types for row in owned):
        return []
    genes = [
        row for row in owned
        if row.get("type", "").lower() == "gene" and row.get("id") not in _EMPTY
        and transcript_id == f"{row['id']}.synthetic_tx"
    ]
    if len(genes) != 1:
        return []
    gene = genes[0]
    expected, exon_intervals = set(), set()
    for path in paths:
        intervals = path.get("cds_intervals")
        if intervals in _EMPTY:
            continue
        occurrence = occurrences[path["occurrence_id"]]
        contig, strand = occurrence.get("contig"), occurrence.get("strand")
        exon_intervals.add((contig, int(occurrence["start"]), int(occurrence["end"]), strand))
        for interval in intervals.split(";"):
            start, end = (int(value) for value in interval.split("-"))
            expected.add((contig, start, end, strand))
    parents = {gene["id"]}
    selected = [gene]
    for row in owned:
        if row.get("type", "").lower() != "exon" or row.get("id") in _EMPTY:
            continue
        location = (row.get("seqid"), int(row["start"]), int(row["end"]), row.get("strand"))
        row_parents = {parent.strip() for parent in str(row.get("parent") or "").split(",")} - _EMPTY
        if row_parents == {gene["id"]} and location in exon_intervals:
            selected.append(row)
            parents.add(row["id"])
    found = set()
    for row in owned:
        if row.get("type", "").lower() != "cds":
            continue
        location = (row.get("seqid"), int(row["start"]), int(row["end"]), row.get("strand"))
        row_parents = {parent.strip() for parent in str(row.get("parent") or "").split(",")} - _EMPTY
        if row_parents and row_parents <= parents and location in expected:
            selected.append(row)
            found.add(location)
    return selected if expected and found == expected else []


def build_coding_transcript(key, paths, occurrences, sequences, raw_features=()):
    """Use path CDS intervals; raw feature attrs, when provided, are parsed dicts."""
    lengths = defaultdict(int)

    def unavailable(reason):
        return CodingTranscript(key, "", (), dict(lengths), reason)

    if raw_features and not any(row.get("id") == key[-1] for row in raw_features):
        features = _synthetic_transcript_features(raw_features, key[-1], paths, occurrences)
        if not features:
            return unavailable("unresolved_raw_transcript_ownership")
    else:
        features = _transcript_features(raw_features, key[-1])
    reason = _translation_reason(features)
    if reason:
        return unavailable(reason)
    raw_phases = defaultdict(set)
    raw_sources = defaultdict(set)
    for feature in features:
        if feature.get("type", "").lower() == "cds":
            location = (feature.get("seqid"), int(feature["start"]), int(feature["end"]), feature.get("strand"))
            raw_phases[location].add(str(feature.get("phase", ".")))
            if feature.get("id") not in _EMPTY:
                raw_sources[location].add(feature["id"])

    segments = []
    for path in paths:
        intervals = path.get("cds_intervals")
        if intervals in _EMPTY:
            if path.get("cds_length") not in _EMPTY and int(path["cds_length"]) > 0:
                return unavailable("missing_CDS_intervals")
            continue
        occurrence = occurrences[path["occurrence_id"]]
        strand = occurrence.get("strand")
        if strand not in {"+", "-"}:
            return unavailable("unknown_CDS_strand")
        coords = [tuple(int(value) for value in token.split("-")) for token in intervals.split(";")]
        coords = sorted(set(coords), reverse=strand == "-")
        for index, (start, end) in enumerate(coords):
            if start < int(occurrence["start"]) or end > int(occurrence["end"]) or start > end:
                return unavailable("CDS_outside_occurrence")
            phases = raw_phases.get((occurrence.get("contig"), start, end, strand), set())
            if len(phases) > 1:
                return unavailable("conflicting_CDS_phase")
            phase = next(iter(phases)) if phases else str(path.get("cds_phase", ".")) if index == 0 else "."
            if phase not in {"0", "1", "2"}:
                return unavailable("unknown_CDS_phase")
            location = (occurrence.get("contig"), start, end, strand)
            segments.append((
                occurrence, start, end, int(phase), str(path.get("path_id", "")),
                tuple(sorted(raw_sources.get(location, ()))),
            ))
            lengths[occurrence["occurrence_id"]] += end - start + 1
    if not segments:
        return unavailable("no_annotated_CDS")
    locations = {
        (row.get("contig"), row.get("strand"))
        for row, _start, _end, _phase, _path_id, _source_ids in segments
    }
    if len(locations) != 1:
        return unavailable("inconsistent_CDS_locus")
    strand = segments[0][0]["strand"]
    segments.sort(key=lambda segment: segment[1], reverse=strand == "-")
    ordered_intervals = sorted((start, end) for _row, start, end, _phase, _path_id, _source_ids in segments)
    if any(left[1] >= right[0] for left, right in zip(ordered_intervals, ordered_intervals[1:])):
        return unavailable("overlapping_CDS_in_transcript")

    bases, positions, base_sources = [], [], []
    initial_phase = segments[0][3]
    for occurrence, start, end, phase, path_id, source_feature_ids in segments:
        if phase != (initial_phase - len(bases)) % 3:
            return unavailable("inconsistent_CDS_phase")
        occurrence_id = occurrence["occurrence_id"]
        sequence = sequences.get(occurrence_id, "").upper()
        if len(sequence) != int(occurrence["end"]) - int(occurrence["start"]) + 1:
            return unavailable("missing_or_incomplete_occurrence_sequence")
        if strand == "-":
            relative_start = int(occurrence["end"]) - end + 1
            relative_end = int(occurrence["end"]) - start + 1
        else:
            relative_start = start - int(occurrence["start"]) + 1
            relative_end = end - int(occurrence["start"]) + 1
        segment_sequence = sequence[relative_start - 1:relative_end]
        bases.extend(segment_sequence)
        for relative_position in range(relative_start, relative_end + 1):
            positions.append((occurrence_id, relative_position))
            genome_position = (
                int(occurrence["start"]) + relative_position - 1
                if strand == "+"
                else int(occurrence["end"]) - relative_position + 1
            )
            cds_position0 = len(base_sources)
            base_sources.append(CodingBaseProjection(
                occurrence_id=occurrence_id,
                occurrence_interval=Interval0(relative_position - 1, relative_position),
                cds_interval=Interval0(cds_position0, cds_position0 + 1),
                contig=str(occurrence.get("contig", "")),
                genome_interval=Interval0(genome_position - 1, genome_position),
                strand=strand,
                source_path_id=path_id,
                source_feature_ids=source_feature_ids,
            ))

    protein, codons, codon_sources = [], [], []
    # Initial partial codons are skipped once. Internal exon phases retain all bases.
    for offset in range(initial_phase, len(bases) - 2, 3):
        codon = "".join(bases[offset:offset + 3])
        amino_acid = "*" if codon in _STANDARD_CODE.stop_codons else _STANDARD_CODE.forward_table.get(codon, "X")
        protein.append(amino_acid)
        codons.append(tuple(positions[offset:offset + 3]) if amino_acid not in {"X", "*"} else None)
        codon_sources.append(tuple(base_sources[offset:offset + 3]))
    if protein and protein[-1] == "*":
        protein.pop()
        codons.pop()
        codon_sources.pop()
    protein = "".join(protein).replace("*", "X")
    if not protein or all(codon is None for codon in codons):
        return unavailable("no_complete_standard_codons")
    source_ids = tuple(sorted({
        source_id
        for codon in codon_sources
        for base in codon
        for source_id in base.source_feature_ids
    }))
    return CodingTranscript(
        key, protein, tuple(codons), dict(lengths), "", tuple(codon_sources), source_ids,
    )


def _blosum62_score(query_aa, target_aa):
    return float(_BLOSUM62[query_aa, target_aa])


def _aligned_occurrence_pairs(query, target, aligned_query, aligned_target):
    """Index aligned codons in 0-based MSA, occurrence and genome coordinates."""
    pairs = {}
    known_columns = []
    query_index = target_index = 0
    for column0, (query_aa, target_aa) in enumerate(zip(aligned_query, aligned_target)):
        query_residue = query_index if query_aa != "-" else None
        target_residue = target_index if target_aa != "-" else None
        query_codon = query.codon_sources[query_residue] if query_residue is not None else None
        target_codon = target.codon_sources[target_residue] if target_residue is not None else None
        query_known = query_residue is not None and query.codons[query_residue] is not None
        target_known = target_residue is not None and target.codons[target_residue] is not None
        query_index += query_aa != "-"
        target_index += target_aa != "-"
        if not query_known or not target_known:
            continue
        known_columns.append((column0, query_aa, target_aa, query_codon, target_codon))
        touched = set()
        for query_base, target_base in zip(query_codon, target_codon):
            key = (query_base.occurrence_id, target_base.occurrence_id)
            record = pairs.setdefault(key, {
                "positions0": [], "aa_pairs": 0, "aa_matches": 0,
                "blosum62_score": 0.0, "msa_columns0": [],
            })
            record["positions0"].append((
                query_base.occurrence_interval.start0,
                target_base.occurrence_interval.start0,
            ))
            touched.add(key)
        for key in touched:
            pairs[key]["aa_pairs"] += 1
            pairs[key]["aa_matches"] += query_aa == target_aa
            pairs[key]["blosum62_score"] += _blosum62_score(query_aa, target_aa)
            pairs[key]["msa_columns0"].append(column0)
    return pairs, tuple(known_columns)


def _coordinate_blocks0(positions0):
    blocks = []
    for query0, target0 in sorted(positions0):
        if (
            blocks
            and query0 == blocks[-1].query.end0
            and target0 == blocks[-1].target.end0
        ):
            blocks[-1] = CoordinateBlock(
                Interval0(blocks[-1].query.start0, query0 + 1),
                Interval0(blocks[-1].target.start0, target0 + 1),
            )
        else:
            blocks.append(CoordinateBlock(
                Interval0(query0, query0 + 1), Interval0(target0, target0 + 1),
            ))
    return tuple(blocks)

def _anchor_metrics(
    record, known_columns, aligned_query, aligned_target,
    query_occurrence, target_occurrence,
):
    columns = sorted(set(record["msa_columns0"]))
    first, last = columns[0], columns[-1]

    def independent_of_focal_occurrences(item):
        query_codon, target_codon = item[3], item[4]
        return (
            all(base.occurrence_id != query_occurrence for base in query_codon)
            and all(base.occurrence_id != target_occurrence for base in target_codon)
        )

    independent_columns = [
        item for item in known_columns if independent_of_focal_occurrences(item)
    ]
    left = [item for item in independent_columns if item[0] < first][-_ANCHOR_WINDOW:]
    right = [item for item in independent_columns if item[0] > last][:_ANCHOR_WINDOW]
    left_score = sum(_blosum62_score(item[1], item[2]) for item in left)
    right_score = sum(_blosum62_score(item[1], item[2]) for item in right)
    left_supported = len(left) >= _MIN_ANCHOR_PAIRS and left_score > 0
    right_supported = len(right) >= _MIN_ANCHOR_PAIRS and right_score > 0
    left_terminal = first == 0
    right_terminal = last + 1 == len(aligned_query) == len(aligned_target)
    span = max(1, last - first + 1)
    gap_columns = sum(
        aligned_query[column0] == "-" or aligned_target[column0] == "-"
        for column0 in range(first, last + 1)
    )
    resolved = (
        (left_supported and right_supported)
        or (left_terminal and right_supported)
        or (right_terminal and left_supported)
    )
    terminal_side = (
        "both" if left_terminal and right_terminal
        else "left" if left_terminal
        else "right" if right_terminal
        else "none"
    )
    return {
        "left_pairs": len(left), "right_pairs": len(right),
        "left_score": left_score, "right_score": right_score,
        "left_supported": left_supported, "right_supported": right_supported,
        "terminal_side": terminal_side, "gap_fraction": gap_columns / span,
        "resolved": resolved,
        "column_interval": Interval0(first, last + 1),
    }


def _candidate_json(candidate):
    """Return one stable, scalar-only record for table serialization."""
    interval = candidate.msa_column_interval
    return {
        "aa_identity": candidate.aa_identity,
        "anchor_resolved": candidate.anchor_resolved,
        "blosum62_score": candidate.blosum62_score,
        "candidate_id": candidate.candidate_id,
        "competing_occurrences": list(candidate.competing_occurrences),
        "known_aa_pairs": candidate.known_aa_pairs,
        "left_anchor_pairs": candidate.left_anchor_pairs,
        "left_anchor_score": candidate.left_anchor_score,
        "left_anchor_supported": candidate.left_anchor_supported,
        "msa_column_end0": interval.end0,
        "msa_column_start0": interval.start0,
        "position_eligible": candidate.position_eligible,
        "position_monotonic": candidate.position_monotonic,
        "projected_blocks": [
            {
                "query_end0": block.query.end0,
                "query_start0": block.query.start0,
                "target_end0": block.target.end0,
                "target_start0": block.target.start0,
            }
            for block in candidate.coordinate_blocks
        ],
        "query_cds_coverage": candidate.query_cds_coverage,
        "query_transcript_id": candidate.query_transcript_id,
        "right_anchor_pairs": candidate.right_anchor_pairs,
        "right_anchor_score": candidate.right_anchor_score,
        "right_anchor_supported": candidate.right_anchor_supported,
        "score_scheme": "BLOSUM62",
        "target_cds_coverage": candidate.target_cds_coverage,
        "target_transcript_id": candidate.target_transcript_id,
        "terminal_side": candidate.terminal_side,
    }


class CodingProjectionIndex:
    """Family MSA index with transcript-specific CDS and genome projections."""

    MAX_CACHED_PAIRS = 32
    MAX_CACHED_BASE_PAIRS = 1_000_000

    def __init__(
        self, occurrences, sequences, transcript_paths, raw_features=(), threads=1,
        msa_mode="linsi",
    ):
        self.threads = threads
        self.msa_mode = msa_mode
        self.cache = OrderedDict()
        self.cached_bases = 0
        self.transcripts = {}
        self.by_occurrence = defaultdict(list)
        self.family_alignments = {}
        self.occurrences = {row["occurrence_id"]: row for row in occurrences}
        self.occurrences_by_copy = defaultdict(list)
        for row in occurrences:
            self.occurrences_by_copy[_copy_key(row)].append(row["occurrence_id"])
        by_copy = defaultdict(list)
        for row in raw_features:
            if row.get("ownership") == "target_gene_descendant":
                by_copy[_copy_key(row)].append(row)
        by_transcript = defaultdict(list)
        for row in transcript_paths:
            by_transcript[(*_copy_key(row), row["transcript_id"])].append(row)
        by_id = self.occurrences
        for key, paths in sorted(by_transcript.items()):
            transcript = build_coding_transcript(key, paths, by_id, sequences, by_copy[key[:-1]])
            self.transcripts[key] = transcript
            for occurrence_id in {row["occurrence_id"] for row in paths}:
                self.by_occurrence[occurrence_id].append(key)

    def _family_alignment(self, family_id):
        if family_id in self.family_alignments:
            return self.family_alignments[family_id]
        family_transcripts = [
            transcript for key, transcript in sorted(self.transcripts.items())
            if key[0] == family_id and not transcript.unavailable_reason and transcript.protein
        ]
        aliases_by_protein = defaultdict(list)
        for transcript in family_transcripts:
            aliases_by_protein[transcript.protein].append(transcript.key)
        records = []
        aliases_by_record = {}
        record_by_transcript = {}
        for index, protein in enumerate(sorted(aliases_by_protein), 1):
            record_id = f"coding_{index:06d}"
            aliases = tuple(sorted(aliases_by_protein[protein]))
            records.append((record_id, protein))
            aliases_by_record[record_id] = aliases
            for transcript_key in aliases:
                record_by_transcript[transcript_key] = record_id
        aligned_records = dict(alignment.protein_multiple_alignment(
            tuple(records), mode=self.msa_mode, threads=self.threads,
        )) if records else {}
        expected_ids = {record_id for record_id, _protein in records}
        if set(aligned_records) != expected_ids:
            raise alignment.AlignmentBackendError("family protein MSA returned a different record set")
        lengths = {len(sequence) for sequence in aligned_records.values()}
        if len(lengths) > 1:
            raise alignment.AlignmentBackendError("family protein MSA returned unequal alignment lengths")
        residue_columns = {}
        for record_id, protein in records:
            aligned = aligned_records[record_id]
            if aligned.replace("-", "") != protein:
                raise alignment.AlignmentBackendError("family protein MSA did not preserve input residues")
            residue_columns[record_id] = tuple(
                column0 for column0, amino_acid in enumerate(aligned) if amino_acid != "-"
            )
        projection = FamilyCodingProjection(
            family_id=family_id, mode=self.msa_mode,
            aligned_records=aligned_records,
            aliases_by_record=aliases_by_record,
            record_by_transcript=record_by_transcript,
            residue_columns=residue_columns,
        )
        self.family_alignments[family_id] = projection
        return projection

    def _pair(self, query_key, target_key):
        key = tuple(sorted((query_key, target_key)))
        if key in self.cache:
            self.cache.move_to_end(key)
            pairs, known_columns, aligned_left, aligned_right, _size = self.cache[key]
        else:
            query, target = (self.transcripts[item] for item in key)
            if query.key[0] != target.key[0]:
                return {}, (), "", "", False
            family = self._family_alignment(query.key[0])
            aligned_left = family.aligned_for(query.key)
            aligned_right = family.aligned_for(target.key)
            pairs, known_columns = _aligned_occurrence_pairs(
                query, target, aligned_left, aligned_right,
            )
            size = sum(len(record["positions0"]) for record in pairs.values())
            if size <= self.MAX_CACHED_BASE_PAIRS:
                while self.cache and (len(self.cache) >= self.MAX_CACHED_PAIRS or self.cached_bases + size > self.MAX_CACHED_BASE_PAIRS):
                    _old_key, (_old_pairs, _old_columns, _left, _right, old_size) = self.cache.popitem(last=False)
                    self.cached_bases -= old_size
                self.cache[key] = (pairs, known_columns, aligned_left, aligned_right, size)
                self.cached_bases += size
        inverted = query_key != key[0]
        if not inverted:
            return pairs, known_columns, aligned_left, aligned_right, False
        reversed_pairs = {}
        for (left_occurrence, right_occurrence), record in pairs.items():
            reversed_pairs[(right_occurrence, left_occurrence)] = {
                **record,
                "positions0": [(right0, left0) for left0, right0 in record["positions0"]],
            }
        reversed_columns = tuple(
            (column0, right_aa, left_aa, right_codon, left_codon)
            for column0, left_aa, right_aa, left_codon, right_codon in known_columns
        )
        return reversed_pairs, reversed_columns, aligned_right, aligned_left, True

    def project(self, occurrence_id):
        """Return residue-to-MSA-to-CDS/genome projections for one occurrence."""
        projections = []
        for transcript_key in sorted(self.by_occurrence.get(occurrence_id, ())):
            transcript = self.transcripts[transcript_key]
            if transcript.unavailable_reason:
                continue
            family = self._family_alignment(transcript_key[0])
            columns = family.columns_for(transcript_key)
            for residue_index0, codon in enumerate(transcript.codon_sources):
                if not any(base.occurrence_id == occurrence_id for base in codon):
                    continue
                projections.append(CodingResidueProjection(
                    transcript_key=transcript_key,
                    residue_index0=residue_index0,
                    msa_column0=columns[residue_index0],
                    amino_acid=transcript.protein[residue_index0],
                    cds_bases=codon,
                ))
        return tuple(projections)

    def family_projection(self, family_id):
        """Expose the run-local family MSA and its complete alias mapping."""
        return self._family_alignment(family_id)

    def candidates(self, query_occurrence, target_copy_id):
        """Return local coding evidence for every occurrence in a target copy."""
        query = self.occurrences.get(query_occurrence)
        if query is None:
            return tuple()
        if isinstance(target_copy_id, tuple):
            copy_key = target_copy_id
        else:
            matches = [
                key for key in self.occurrences_by_copy
                if key[0] == query.get("family_id") and key[2] == target_copy_id
            ]
            if len(matches) != 1:
                return tuple()
            copy_key = matches[0]
        return tuple(
            {"target_occurrence_id": target, **self.evidence(query_occurrence, target)}
            for target in sorted(self.occurrences_by_copy.get(copy_key, ()))
        )

    def _competing_occurrences(
        self, query_key, target_key, query_occurrence, target_occurrence,
        query_positions0, target_positions0,
    ):
        """Find alternative occurrence placements supported by compatible paths."""
        competitors = set()
        query_copy = query_key[:3]
        target_copy = target_key[:3]
        target_transcripts = [
            key for key, transcript in self.transcripts.items()
            if key[:3] == target_copy and not transcript.unavailable_reason
        ]
        query_transcripts = [
            key for key, transcript in self.transcripts.items()
            if key[:3] == query_copy and not transcript.unavailable_reason
        ]
        for alternative_target_key in target_transcripts:
            pairs, _columns, _aligned_query, _aligned_target, _inverted = self._pair(
                query_key, alternative_target_key,
            )
            for (other_query, other_target), record in pairs.items():
                if other_query != query_occurrence or other_target == target_occurrence:
                    continue
                other_query_positions0 = {
                    query0 for query0, _target0 in record["positions0"]
                }
                if query_positions0 & other_query_positions0:
                    competitors.add(f"target:{other_target}")
        for alternative_query_key in query_transcripts:
            pairs, _columns, _aligned_query, _aligned_target, _inverted = self._pair(
                alternative_query_key, target_key,
            )
            for (other_query, other_target), record in pairs.items():
                if other_target != target_occurrence or other_query == query_occurrence:
                    continue
                other_target_positions0 = {
                    target0 for _query0, target0 in record["positions0"]
                }
                if target_positions0 & other_target_positions0:
                    competitors.add(f"query:{other_query}")
        return tuple(sorted(competitors))

    def evidence(self, query_occurrence, target_occurrence):
        query_keys = self.by_occurrence.get(query_occurrence, [])
        target_keys = self.by_occurrence.get(target_occurrence, [])
        result = {
            "protein_status": "unavailable",
            "protein_mapping_status": "uncovered",
            "protein_unavailable_reason": "no_CDS_transcript_path",
            "protein_membership_eligible": False,
            "protein_position_eligible": False,
            "protein_hard_observation_eligible": False,
            "protein_candidate_evidence_available": False,
        }
        unavailable, candidates = set(), []
        for query_key in query_keys:
            query = self.transcripts[query_key]
            for target_key in target_keys:
                target = self.transcripts[target_key]
                for side, transcript in (("query", query), ("target", target)):
                    if transcript.unavailable_reason:
                        unavailable.add(f"{side}:{transcript.key[-1]}:{transcript.unavailable_reason}")
                if query.unavailable_reason or target.unavailable_reason:
                    continue
                if not query.coding_lengths.get(query_occurrence) or not target.coding_lengths.get(target_occurrence):
                    unavailable.add("no_CDS_in_requested_occurrence")
                    continue
                result["protein_status"] = "no_aligned_CDS"
                pairs, known_columns, aligned_query, aligned_target, _inverted = self._pair(
                    query_key, target_key,
                )
                record = pairs.get((query_occurrence, target_occurrence))
                if not record:
                    continue
                positions = set(record["positions0"])
                identity = record["aa_matches"] / record["aa_pairs"]
                query_coverage = len(positions) / query.coding_lengths[query_occurrence]
                target_coverage = len(positions) / target.coding_lengths[target_occurrence]
                anchors = _anchor_metrics(
                    record, known_columns, aligned_query, aligned_target,
                    query_occurrence, target_occurrence,
                )
                ordered = sorted(positions)
                position_monotonic = all(
                    right[0] > left[0] and right[1] > left[1]
                    for left, right in zip(ordered, ordered[1:])
                )
                query_positions0 = {query0 for query0, _target0 in positions}
                target_positions0 = {target0 for _query0, target0 in positions}
                candidates.append(CodingProjectionCandidate(
                    query_transcript_key=query_key,
                    target_transcript_key=target_key,
                    coordinate_blocks=_coordinate_blocks0(positions),
                    aa_identity=identity,
                    query_cds_coverage=query_coverage,
                    target_cds_coverage=target_coverage,
                    known_aa_pairs=record["aa_pairs"],
                    blosum62_score=record["blosum62_score"],
                    gap_fraction=anchors["gap_fraction"],
                    left_anchor_pairs=anchors["left_pairs"],
                    right_anchor_pairs=anchors["right_pairs"],
                    left_anchor_score=anchors["left_score"],
                    right_anchor_score=anchors["right_score"],
                    left_anchor_supported=anchors["left_supported"],
                    right_anchor_supported=anchors["right_supported"],
                    terminal_side=anchors["terminal_side"],
                    msa_column_interval=anchors["column_interval"],
                    anchor_resolved=anchors["resolved"],
                    position_monotonic=position_monotonic,
                    competing_occurrences=self._competing_occurrences(
                        query_key, target_key, query_occurrence, target_occurrence,
                        query_positions0, target_positions0,
                    ),
                ))
        result["protein_unavailable_reason"] = ";".join(sorted(unavailable)) or ("NA" if candidates or result["protein_status"] != "unavailable" else "no_CDS_transcript_path")
        if not candidates:
            return result
        candidates = tuple(sorted(
            candidates,
            key=lambda candidate: (
                candidate.query_transcript_key, candidate.target_transcript_key,
                candidate.coordinate_blocks,
            ),
        ))
        candidate_set = CodingProjectionCandidateSet(candidates)
        best = max(candidates, key=lambda candidate: (
            candidate.position_eligible, candidate.anchor_resolved,
            candidate.known_aa_pairs, candidate.blosum62_score,
            candidate.aa_identity,
            max(candidate.query_cds_coverage, candidate.target_cds_coverage),
        ))
        ambiguous = (
            not candidate_set.coordinate_consensus
            or bool(candidate_set.competing_occurrences)
            or any(not candidate.position_monotonic for candidate in candidates)
        )
        mapping_status = (
            "ambiguous_mapping" if ambiguous
            else "resolved_local" if candidate_set.position_eligible
            else "supported_unanchored"
        )
        query_sources = self.transcripts[best.query_transcript_key].source_ids
        target_sources = self.transcripts[best.target_transcript_key].source_ids
        interval = best.msa_column_interval
        # ``protein_status`` remains available for the existing caller. New
        # state construction must use mapping_status and hard eligibility.
        result.update(
            protein_status="ambiguous_transcript_projection" if ambiguous else "supported",
            protein_mapping_status=mapping_status,
            protein_membership_eligible=candidate_set.membership_eligible,
            protein_position_eligible=candidate_set.position_eligible,
            protein_hard_observation_eligible=candidate_set.position_eligible,
            protein_candidate_evidence_available=candidate_set.candidate_evidence_available,
            protein_aa_identity=best.aa_identity,
            protein_query_cds_coverage=best.query_cds_coverage,
            protein_target_cds_coverage=best.target_cds_coverage,
            protein_known_aa_pairs=best.known_aa_pairs,
            protein_blosum62_score=best.blosum62_score,
            protein_gap_fraction=best.gap_fraction,
            protein_left_anchor_pairs=best.left_anchor_pairs,
            protein_right_anchor_pairs=best.right_anchor_pairs,
            protein_left_anchor_score=best.left_anchor_score,
            protein_right_anchor_score=best.right_anchor_score,
            protein_left_anchor_supported=best.left_anchor_supported,
            protein_right_anchor_supported=best.right_anchor_supported,
            protein_terminal_side=best.terminal_side,
            protein_msa_column_start0=interval.start0,
            protein_msa_column_end0=interval.end0,
            protein_msa_column_interval=f"{interval.start0}:{interval.end0}",
            protein_msa_mode=self.msa_mode,
            protein_candidate_mapping_count=len(candidates),
            protein_candidate_coordinate_consensus=candidate_set.coordinate_consensus,
            protein_competing_occurrences=";".join(candidate_set.competing_occurrences) or "NA",
            protein_candidate_details=json.dumps(
                [_candidate_json(candidate) for candidate in candidates],
                sort_keys=True, separators=(",", ":"),
            ),
            protein_candidate_set=candidate_set,
            protein_query_source_features=";".join(query_sources) or "NA",
            protein_target_source_features=";".join(target_sources) or "NA",
            protein_metrics_scope="best_transcript_pair",
            protein_best_query_transcript=best.query_transcript_id,
            protein_best_target_transcript=best.target_transcript_id,
            protein_supporting_transcripts=";".join(
                f"{candidate.query_transcript_id}>{candidate.target_transcript_id}"
                for candidate in candidates
            ),
        )
        if ambiguous:
            return result
        blocks = candidates[0].coordinate_blocks
        result["protein_projected_coordinate_blocks"] = blocks
        result["protein_projected_blocks"] = blocks
        return result

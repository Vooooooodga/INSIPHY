"""coding / transcripts: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from collections import defaultdict
from intraphy.preparation.features import FeatureHierarchy
from intraphy.coding.types import CodingBaseProjection
from intraphy.coding.types import CodingTranscript
from intraphy.coding.types import _EMPTY
from intraphy.coding.types import _STANDARD_CODE
from intraphy.coordinates import Interval0
from urllib.parse import unquote


def _copy_key(row):
    return row.get("family_id", ""), row.get("species", ""), row.get("gene_copy_id", "")


def _transcript_features(rows, transcript_id, *, feature_index=None):
    index = feature_index if feature_index is not None else FeatureHierarchy(rows)
    return index.transcript_features(transcript_id)


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


def build_coding_transcript(key, paths, occurrences, sequences, raw_features=(), *, feature_index=None):
    """Use path CDS intervals; raw feature attrs, when provided, are parsed dicts."""
    lengths = defaultdict(int)

    def unavailable(reason):
        return CodingTranscript(key, "", (), dict(lengths), reason)

    if raw_features and key[-1] not in (feature_index.by_id if feature_index is not None else {row.get("id") for row in raw_features}):
        features = _synthetic_transcript_features(raw_features, key[-1], paths, occurrences)
        if not features:
            return unavailable("unresolved_raw_transcript_ownership")
    else:
        features = _transcript_features(raw_features, key[-1], feature_index=feature_index)
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

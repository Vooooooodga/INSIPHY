"""Project external protein alignments onto annotated, transcript-specific CDS.

The fixed 0.70 amino-acid identity and 0.60 directional CDS coverage gates are
numerical homology heuristics, not probabilities. One directional coverage can
suffice for a complementary split; repeated reference coverage is resolved by
the existing occurrence graph. No protein evidence assigns an exonic role.

Directional coverage uses all annotated CDS bases of the requested occurrence
in that transcript as its denominator, including partial/unknown codons and
annotated stop codons. Those bases contribute no mapped support when excluded
from the protein alignment. Identity and coverage describe the best transcript
pair; accepted coordinates combine compatible, independently supported pairs.
"""

from collections import OrderedDict, defaultdict
from dataclasses import dataclass
from urllib.parse import unquote

from Bio.Data import CodonTable

from . import alignment


MIN_AA_IDENTITY = 0.70
MIN_CDS_COVERAGE = 0.60
_STANDARD_CODE = CodonTable.unambiguous_dna_by_id[1]
_EMPTY = {"", ".", "NA", None}


@dataclass
class CodingTranscript:
    key: tuple
    protein: str
    codons: tuple
    coding_lengths: dict
    unavailable_reason: str = ""


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
    for feature in features:
        if feature.get("type", "").lower() == "cds":
            location = (feature.get("seqid"), int(feature["start"]), int(feature["end"]), feature.get("strand"))
            raw_phases[location].add(str(feature.get("phase", ".")))

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
            segments.append((occurrence, start, end, int(phase)))
            lengths[occurrence["occurrence_id"]] += end - start + 1
    if not segments:
        return unavailable("no_annotated_CDS")
    locations = {(row.get("contig"), row.get("strand")) for row, _start, _end, _phase in segments}
    if len(locations) != 1:
        return unavailable("inconsistent_CDS_locus")
    strand = segments[0][0]["strand"]
    segments.sort(key=lambda segment: segment[1], reverse=strand == "-")
    ordered_intervals = sorted((start, end) for _row, start, end, _phase in segments)
    if any(left[1] >= right[0] for left, right in zip(ordered_intervals, ordered_intervals[1:])):
        return unavailable("overlapping_CDS_in_transcript")

    bases, positions = [], []
    initial_phase = segments[0][3]
    for occurrence, start, end, phase in segments:
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
        bases.extend(sequence[relative_start - 1:relative_end])
        positions.extend((occurrence_id, position) for position in range(relative_start, relative_end + 1))

    protein, codons = [], []
    # Initial partial codons are skipped once. Internal exon phases retain all bases.
    for offset in range(initial_phase, len(bases) - 2, 3):
        codon = "".join(bases[offset:offset + 3])
        amino_acid = "*" if codon in _STANDARD_CODE.stop_codons else _STANDARD_CODE.forward_table.get(codon, "X")
        protein.append(amino_acid)
        codons.append(tuple(positions[offset:offset + 3]) if amino_acid not in {"X", "*"} else None)
    if protein and protein[-1] == "*":
        protein.pop()
        codons.pop()
    protein = "".join(protein).replace("*", "X")
    if not protein or all(codon is None for codon in codons):
        return unavailable("no_complete_standard_codons")
    return CodingTranscript(key, protein, tuple(codons), dict(lengths))


def _aligned_occurrence_pairs(query, target, aligned_query, aligned_target):
    """Index each aligned codon once, including codons spanning exon junctions."""
    pairs = {}
    query_index = target_index = 0
    for query_aa, target_aa in zip(aligned_query, aligned_target):
        query_codon = query.codons[query_index] if query_aa != "-" else None
        target_codon = target.codons[target_index] if target_aa != "-" else None
        query_index += query_aa != "-"
        target_index += target_aa != "-"
        if query_codon is None or target_codon is None:
            continue
        touched = set()
        for (query_occ, query_pos), (target_occ, target_pos) in zip(query_codon, target_codon):
            key = (query_occ, target_occ)
            record = pairs.setdefault(key, {"positions": [], "aa_pairs": 0, "aa_matches": 0})
            record["positions"].append((query_pos, target_pos))
            touched.add(key)
        for key in touched:
            pairs[key]["aa_pairs"] += 1
            pairs[key]["aa_matches"] += query_aa == target_aa
    return pairs


def _coordinate_blocks(positions):
    blocks = []
    for query, target in sorted(positions):
        if blocks and query == blocks[-1][1] + 1 and target == blocks[-1][3] + 1:
            blocks[-1] = (blocks[-1][0], query, blocks[-1][2], target)
        else:
            blocks.append((query, query, target, target))
    return tuple(blocks)


class CodingProjectionIndex:
    """Run-local transcript alignments, bounded by both pair count and mapped bases."""

    MAX_CACHED_PAIRS = 32
    MAX_CACHED_BASE_PAIRS = 1_000_000

    def __init__(self, occurrences, sequences, transcript_paths, raw_features=(), threads=1):
        self.threads = threads
        self.cache = OrderedDict()
        self.cached_bases = 0
        self.transcripts = {}
        self.by_occurrence = defaultdict(list)
        by_copy = defaultdict(list)
        for row in raw_features:
            if row.get("ownership") == "target_gene_descendant":
                by_copy[_copy_key(row)].append(row)
        by_transcript = defaultdict(list)
        for row in transcript_paths:
            by_transcript[(*_copy_key(row), row["transcript_id"])].append(row)
        by_id = {row["occurrence_id"]: row for row in occurrences}
        for key, paths in sorted(by_transcript.items()):
            transcript = build_coding_transcript(key, paths, by_id, sequences, by_copy[key[:-1]])
            self.transcripts[key] = transcript
            for occurrence_id in {row["occurrence_id"] for row in paths}:
                self.by_occurrence[occurrence_id].append(key)

    def _pair(self, query_key, target_key):
        key = tuple(sorted((query_key, target_key)))
        if key in self.cache:
            self.cache.move_to_end(key)
            pairs = self.cache[key][0]
        else:
            query, target = (self.transcripts[item] for item in key)
            aligned_query, aligned_target = alignment.protein_pair_alignment(query.protein, target.protein, threads=self.threads)
            pairs = _aligned_occurrence_pairs(query, target, aligned_query, aligned_target)
            size = sum(len(record["positions"]) for record in pairs.values())
            if size <= self.MAX_CACHED_BASE_PAIRS:
                while self.cache and (len(self.cache) >= self.MAX_CACHED_PAIRS or self.cached_bases + size > self.MAX_CACHED_BASE_PAIRS):
                    _old_key, (_old_pairs, old_size) = self.cache.popitem(last=False)
                    self.cached_bases -= old_size
                self.cache[key] = (pairs, size)
                self.cached_bases += size
        return pairs, query_key != key[0]

    def evidence(self, query_occurrence, target_occurrence):
        query_keys = self.by_occurrence.get(query_occurrence, [])
        target_keys = self.by_occurrence.get(target_occurrence, [])
        result = {"protein_status": "unavailable", "protein_unavailable_reason": "no_CDS_transcript_path"}
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
                pairs, inverted = self._pair(query_key, target_key)
                pair_key = (target_occurrence, query_occurrence) if inverted else (query_occurrence, target_occurrence)
                record = pairs.get(pair_key)
                if not record:
                    continue
                positions = {(target, query) for query, target in record["positions"]} if inverted else set(record["positions"])
                identity = record["aa_matches"] / record["aa_pairs"]
                query_coverage = len(positions) / query.coding_lengths[query_occurrence]
                target_coverage = len(positions) / target.coding_lengths[target_occurrence]
                candidates.append((identity, query_coverage, target_coverage, positions, query_key[-1], target_key[-1]))
        result["protein_unavailable_reason"] = ";".join(sorted(unavailable)) or ("NA" if candidates or result["protein_status"] != "unavailable" else "no_CDS_transcript_path")
        if not candidates:
            return result
        supported = [item for item in candidates if item[0] >= MIN_AA_IDENTITY and max(item[1], item[2]) >= MIN_CDS_COVERAGE]
        best = max(supported or candidates, key=lambda item: (item[0], max(item[1], item[2]), len(item[3])))
        result.update(
            protein_status="supported" if supported else "low_protein_similarity_or_coverage",
            protein_aa_identity=best[0], protein_query_cds_coverage=best[1],
            protein_target_cds_coverage=best[2],
            protein_metrics_scope="best_transcript_pair",
            protein_best_query_transcript=best[4], protein_best_target_transcript=best[5],
            protein_supporting_transcripts=";".join(f"{item[4]}>{item[5]}" for item in supported) or "NA",
        )
        if not supported:
            return result
        positions = set().union(*(item[3] for item in supported))
        ordered = sorted(positions)
        if any(right[0] <= left[0] or right[1] <= left[1] for left, right in zip(ordered, ordered[1:])):
            result["protein_status"] = "ambiguous_transcript_projection"
            return result
        result["protein_projected_blocks"] = _coordinate_blocks(positions)
        return result

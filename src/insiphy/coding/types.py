"""coding / types: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from Bio.Align import substitution_matrices
from Bio.Data import CodonTable
from dataclasses import dataclass
from insiphy.coordinates import Interval0


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

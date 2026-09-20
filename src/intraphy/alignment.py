"""Compatibility exports. Implementations live in the documented submodules."""


# Public entry points; implementations have a single owner.
from intraphy.aligners.columns import (
    revcomp,
    ungapped_identity,
    _compress_ops,
    _parse_cigar,
    _nt,
    _known_base,
    _empty_stats,
    _stats_from_alignment_columns,
    _columns_from_cigar,
)
from intraphy.aligners.external import (
    _external_minimap2_stats,
    _external_miniprot_stats,
    _mafft_pair_alignment,
    protein_pair_alignment,
    protein_multiple_alignment,
    _external_mafft_stats,
    _external_lastz_stats,
)
from intraphy.aligners.formats import (
    _write_temp_fasta,
    _write_temp_fasta_records,
    _parse_fasta_records,
    _parse_paf_tags,
    _paired_bases_from_cigar,
)
from intraphy.aligners.legacy import (
    _candidate_from_legacy_stats,
    _legacy_coordinate_tuple,
    _legacy_stats_from_candidate,
    _candidate_as_legacy_hit,
)
from intraphy.aligners.pairwise import (
    _trim_terminal_overhangs,
    overlap_alignment_stats,
    _pairwise_aligner_stats,
    global_alignment_stats,
    local_alignment_stats,
    best_ungapped_hit,
)
from intraphy.aligners.projection import (
    _parse_gff_attributes,
    protein_locus_exons,
)
from intraphy.aligners.runner import (
    _backend_version,
    available_alignment_backends,
)
from intraphy.aligners.scoring import (
    splice_motif_score,
    phase_compatibility,
)
from intraphy.aligners.short import (
    _normalized_nucleotide_sequence,
    _pairwise_alignment_columns,
    _gap_blocks_from_coordinates,
    _normalized_search_interval,
    anchored_short_alignment,
)
from intraphy.aligners.types import (
    DNA_COMPLEMENT,
    KNOWN_NT,
    MAX_INTERNAL_DP_CELLS,
    MAX_MAFFT_PAIR_CELLS,
    MAX_OPTIMAL_ALIGNMENTS,
    NT_BLASTN_V1,
    NT_BLASTN_V1_MATCH,
    NT_BLASTN_V1_MISMATCH,
    NT_BLASTN_V1_GAP_OPEN,
    NT_BLASTN_V1_GAP_EXTEND,
    AlignmentStats,
    AlignmentGap,
    AlignmentCandidate,
    AlignmentCandidateSet,
    AlignmentBackendError,
    _gap_as_dict,
)
import re
import shutil
import subprocess
import tempfile

import tempfile

import unittest

from contextlib import redirect_stdout, nullcontext

from io import StringIO

from pathlib import Path

from intraphy.candidate_chain import (
    DEFAULT_CHAIN_CONFIGURATION,
    ChainCandidate,
    ChainPathMembership,
    classify_reference_coverage,
    ordered_candidate_chain,
)

from intraphy.coordinates import (
    ClosedInterval1,
    CoordinateBlock,
    Interval0,
    format_legacy_blocks,
    genome_interval_to_local,
    local_interval_to_genome,
    parse_legacy_blocks,
)

from intraphy.io import (
    LEGACY_STRUCTURAL_SITE_SCHEMA_VERSION,
    normalize_structural_site_row,
    read_structural_site_matrix,
    validate_structural_site_tip_rows,
    write_structural_site_matrix,
)

from unittest.mock import patch

from intraphy.cli import main

from intraphy.io import read_tsv

from intraphy.orthofinder import _resolve_members_to_locus, import_orthofinder


class InterfaceTestsSupport:
    def _orthofinder_root(self, root, header, row):
        og = root / "orthofinder" / "Orthogroups"
        og.mkdir(parents=True)
        (og / "Orthogroups.tsv").write_text(header + "\n" + row + "\n")
        return root / "orthofinder"

    def _orthofinder_txt_root(self, root, line):
        og = root / "orthofinder" / "Orthogroups"
        og.mkdir(parents=True)
        (og / "Orthogroups.txt").write_text(line + "\n")
        return root / "orthofinder"

    def _manifest(self, root, rows):
        path = root / "genomes.tsv"
        path.write_text(
            "species\tgenome_fasta\tannotation_file\n"
            + "\n".join(f"{species}\t{root / (species + '.fa')}\t{gff}" for species, gff in rows)
            + "\n"
        )
        return path

import shutil

import tempfile

import unittest

from pathlib import Path

from unittest.mock import patch

from intraphy.alignment import AlignmentBackendError, AlignmentStats, global_alignment_stats, phase_compatibility

from intraphy.annotation import _overlapping_annotation_role

from intraphy.cli import run_all

from intraphy.correspondence import simple_identity

from intraphy.io import read_tsv

from intraphy.benchmark import benchmark_events

from intraphy.calibration import calibrate_simulations

from intraphy.case import build_case, inspect_annotation, scan_hidden_segments

from intraphy.orthofinder import import_orthofinder

from intraphy.preprocess import derive_tables, extract_gene, graph_components

from intraphy.simulate import simulate_dataset

from intraphy.structural_phylogeny import build_structural_site_matrix, fit_model, infer_single_copy_phylogeny

from intraphy.tree import SpeciesTree

from intraphy.visualize import visualize_results

ROOT = Path(__file__).resolve().parents[1]


class RealCasePreparationTestsSupport:
    def write_case_files(self, tmp, species):
        genome = tmp / f"{species}.fa"
        annot = tmp / f"{species}.gff3"
        genome.write_text(">chr1\n" + "ACGT" * 120 + "\n")
        annot.write_text(
            "\n".join(
                [
                    f"chr1\tIntraPhy\tgene\t10\t160\t.\t+\t.\tID={species}_geneA;Name=GeneA;Alias=jgw,jingwei",
                    f"chr1\tIntraPhy\tmRNA\t10\t160\t.\t+\t.\tID={species}_txA;Parent={species}_geneA",
                    f"chr1\tIntraPhy\tCDS\t20\t50\t.\t+\t0\tID={species}_cds1;Parent={species}_txA",
                    f"chr1\tIntraPhy\tCDS\t100\t140\t.\t+\t1\tID={species}_cds2;Parent={species}_txA",
                ]
            )
            + "\n"
        )
        return genome, annot

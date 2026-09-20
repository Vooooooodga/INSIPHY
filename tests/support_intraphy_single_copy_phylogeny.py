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


class SingleCopyPhylogenyTestsSupport:
    def write_structural_case(self, root):
        input_dir = root / "input"
        result_dir = root / "result"
        input_dir.mkdir()
        result_dir.mkdir()
        (input_dir / "species_tree.tsv").write_text(
            "node_id\tparent_id\tlabel\tbranch_length\n"
            "root\t\troot\t0\n"
            "ab\troot\tab\t0.5\n"
            "a\tab\tA\t0.2\n"
            "b\tab\tB\t0.2\n"
            "cd\troot\tcd\t0.5\n"
            "c\tcd\tC\t0.2\n"
            "d\tcd\tD\t0.2\n"
        )
        occurrence_rows = [
            "occurrence_id\tfamily_id\tspecies\tgene_copy_id\ttranscript_id\tsegment_id\tcontig\tstart\tend\tstrand\trole\tpresence_status\tboundary_state\tevidence"
        ]
        element_rows = [
            "element_id\tfamily_id\thomology_id\toccurrence_id\tspecies\tgene_copy_id\telement_class\tdisplay_role\tsource_label\tsupport_type\tconfidence\tmembership_score\tmembership_call"
        ]
        path_rows = [
            "path_id\tfamily_id\tspecies\tgene_copy_id\ttranscript_id\tpath_rank\toccurrence_id\trole\tcontig\tstart\tend\tstrand\tphase\tpath_status"
        ]
        for species in "ABCD":
            roles = [
                "CDS",
                "CDS" if species in "AB" else "intron",
                "CDS" if species in "AC" else "intron",
                "CDS" if species in "AD" else "intron",
            ]
            for index, role in enumerate(roles, start=1):
                occurrence = f"{species}_e{index}"
                occurrence_rows.append(
                    f"{occurrence}\tfam\t{species}\t{species}_gene\t{species}_tx\te{index}\tchr1\t{index * 100}\t{index * 100 + 50}\t+\t{role}\tpresent\tconserved\ttest_fixture"
                )
                element_class = "exon_like" if role == "CDS" else "candidate_source"
                element_rows.append(
                    f"EG_{index}\tfam\tH_{index}\t{occurrence}\t{species}\t{species}_gene\t{element_class}\t{role}\tfam\ttest_fixture\thigh\t0.95\tcore_member"
                )
                path_rows.append(
                    f"{species}_p{index}\tfam\t{species}\t{species}_gene\t{species}_tx\t{index}\t{occurrence}\t{role}\tchr1\t{index * 100}\t{index * 100 + 50}\t+\t0\tcanonical_transcript_path"
                )
        (input_dir / "segment_occurrences.tsv").write_text("\n".join(occurrence_rows) + "\n")
        (input_dir / "transcript_paths.tsv").write_text("\n".join(path_rows) + "\n")
        (result_dir / "element_correspondence.tsv").write_text("\n".join(element_rows) + "\n")
        return input_dir, result_dir

import re

import tempfile

import unittest

from pathlib import Path

from xml.etree import ElementTree

from intraphy.io import read_tsv, write_tsv

from intraphy.visualize import draw_integrated_phylo_synteny, draw_phylogeny, draw_synteny, visualize_results


class VisualizationSemanticTestsSupport:
    def write_matches(self, directory, rows):
        write_tsv(
            directory / "segment_matches.tsv",
            rows,
            ["match_id", "query_occurrence_id", "subject_occurrence_id", "match_status",
             "alignment_strand", "projected_reference_blocks", "alignment_backend",
             "correspondence_basis", "protein_projected_blocks", "matched_blocks"],
        )

    def write_branch_events(self, directory, rows):
        write_tsv(
            directory / "branch_structural_events.tsv",
            rows,
            ["family_id", "layer", "site_id", "branch_scope", "event_type",
             "placement_status", "call_scope"],
        )

    def write_case(self, root):
        input_dir = root / "input"
        result_dir = root / "result"
        output_dir = root / "fig"
        input_dir.mkdir()
        result_dir.mkdir()
        output_dir.mkdir()
        (input_dir / "species_tree.tsv").write_text(
            "node_id\tparent_id\tlabel\tbranch_length\n"
            "root\t\troot\t0\n"
            "a\troot\tA\t0.1\n"
            "b\troot\tB\t0.1\n"
        )
        (input_dir / "segment_occurrences.tsv").write_text(
            "occurrence_id\tfamily_id\tspecies\tgene_copy_id\ttranscript_id\tsegment_id\tcontig\tstart\tend\tstrand\trole\tpresence_status\tboundary_state\tevidence\n"
            "A_e1\tfam\tA\tA_gene\tA_tx1\te1\tchr1\t100\t150\t+\tCDS\tpresent\tconserved\tannotated\n"
            "A_e2\tfam\tA\tA_gene\tA_tx2\te2\tchr1\t220\t260\t+\tCDS\tpresent\tconserved\tannotated\n"
            "A_e3\tfam\tA\tA_gene\tA_tx3\te3\tchr1\t320\t360\t+\tCDS\tpresent\tconserved\tannotated\n"
            "A_e4\tfam\tA\tA_gene\tA_tx3\te4\tchr1\t400\t440\t+\tCDS\tpresent\tconserved\tannotated\n"
            "B_e1\tfam\tB\tB_gene\tB_tx1\te1\tchr1\t100\t150\t+\tCDS\tpresent\tconserved\tannotated\n"
            "B_cand\tfam\tB\tB_gene\tB_tx1\tcand\tchr1\t220\t260\t+\tintron\tpresent\tunknown\tsequence_candidate\n"
            "B_pred\tfam\tB\tB_gene\tB_tx2\tpred\tchr1\t320\t360\t+\tCDS\tpresent\tpredicted\tsupports_hidden_segment\n"
            "B_unknown\tfam\tB\tB_gene\tB_tx2\te4\tchr1\t400\t440\t+\tunknown\tpresent\tconserved\tsequence_supported\n"
        )
        (input_dir / "transcript_paths.tsv").write_text(
            "path_id\tfamily_id\tspecies\tgene_copy_id\ttranscript_id\tpath_rank\toccurrence_id\trole\tcontig\tstart\tend\tstrand\tphase\tpath_status\n"
            "A_tx1_p1\tfam\tA\tA_gene\tA_tx1\t1\tA_e1\tCDS\tchr1\t100\t150\t+\t0\tobserved_transcript_path\n"
            "A_tx2_p0\tfam\tA\tA_gene\tA_tx2\t1\tA_e1\tCDS\tchr1\t100\t150\t+\t0\tobserved_transcript_path\n"
            "A_tx2_p1\tfam\tA\tA_gene\tA_tx2\t1\tA_e2\tCDS\tchr1\t220\t260\t+\t0\tobserved_transcript_path\n"
            "A_tx3_p1\tfam\tA\tA_gene\tA_tx3\t1\tA_e3\tCDS\tchr1\t320\t360\t+\t0\tobserved_transcript_path\n"
            "A_tx3_p2\tfam\tA\tA_gene\tA_tx3\t2\tA_e4\tCDS\tchr1\t400\t440\t+\t0\tobserved_transcript_path\n"
            "B_tx1_p1\tfam\tB\tB_gene\tB_tx1\t1\tB_e1\tCDS\tchr1\t100\t150\t+\t0\tobserved_transcript_path\n"
            "B_tx1_p2\tfam\tB\tB_gene\tB_tx1\t2\tB_cand\tintron\tchr1\t220\t260\t+\t.\tobserved_transcript_path\n"
            "B_tx2_p1\tfam\tB\tB_gene\tB_tx2\t1\tB_pred\tCDS\tchr1\t320\t360\t+\t0\tpredicted_transcript_path\n"
            "B_tx2_p2\tfam\tB\tB_gene\tB_tx2\t2\tB_unknown\tunknown\tchr1\t400\t440\t+\t.\tobserved_transcript_path\n"
        )
        (result_dir / "element_correspondence.tsv").write_text(
            "element_id\tfamily_id\thomology_id\toccurrence_id\tspecies\tgene_copy_id\telement_class\tdisplay_role\tsource_label\tsupport_type\tconfidence\tmembership_score\tmembership_call\tinferred_role\tpredicted_role\n"
            "EG_1\tfam\tH_1\tA_e1\tA\tA_gene\texon_like\tCDS\tfam\tannotation\t1\t1\tcore_member\tCDS\tCDS\n"
            "EG_1\tfam\tH_1\tB_e1\tB\tB_gene\texon_like\tCDS\tfam\tannotation\t1\t1\tcore_member\tCDS\tCDS\n"
            "EG_2\tfam\tH_2\tA_e2\tA\tA_gene\texon_like\tCDS\tfam\tannotation\t1\t1\tcore_member\tCDS\tCDS\n"
            "EG_2\tfam\tH_2\tB_cand\tB\tB_gene\tcandidate_source\tintron\tfam\tsequence_candidate\t0.8\t0.8\tambiguous_member\tnon_exonic_source\tNA\n"
            "EG_3\tfam\tH_3\tA_e3\tA\tA_gene\texon_like\tCDS\tfam\tannotation\t1\t1\tcore_member\tCDS\tCDS\n"
            "EG_3\tfam\tH_3\tB_pred\tB\tB_gene\texon_like\tCDS\tfam\tpredicted_exon_candidate\t0.8\t0.8\tambiguous_member\tpredicted_CDS\tCDS\n"
            "EG_4\tfam\tH_4\tA_e4\tA\tA_gene\texon_like\tCDS\tfam\tannotation\t1\t1\tcore_member\tCDS\tCDS\n"
            "EG_4\tfam\tH_4\tB_unknown\tB\tB_gene\texon_like\tunknown\tfam\tsequence_supported\t1\t1\tcore_member\tNA\tNA\n"
        )
        self.write_matches(input_dir, [
            {
                "match_id": f"m{i}", "query_occurrence_id": f"A_e{i}",
                "subject_occurrence_id": target, "match_status": "mapped",
                "alignment_strand": "+", "projected_reference_blocks": f"1-{length}:1-{length}",
                "alignment_backend": "mafft_overlap", "matched_blocks": f"1-{length}:1-{length}",
            }
            for i, target, length in ((1, "B_e1", 51), (2, "B_cand", 41), (3, "B_pred", 41), (4, "B_unknown", 41))
        ])
        return input_dir, result_dir, output_dir

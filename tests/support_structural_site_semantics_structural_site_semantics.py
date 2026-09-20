import json

import tempfile

import unittest

from pathlib import Path

from intraphy.io import normalize_structural_site_row, read_tsv

from intraphy.structural_sites import (
    _complete_tree_tip_observations,
    _continuous_reference_block_covers,
    _element_site_rows,
    _junction_site_rows,
    _mapped_reference_coordinate,
    _merge_catalogue_observations,
    _parse_projected_reference_blocks,
)


class StructuralSiteSemanticsTestsSupport:
    def _base_element(self):
        occurrences = [{
            "occurrence_id": "A_ex", "family_id": "fam", "species": "A",
            "gene_copy_id": "gA", "role": "CDS", "presence_status": "present",
            "contig": "chrA", "start": "1", "end": "30", "strand": "+",
            "source_feature_id": "exA",
        }]
        elements = [{
            "element_id": "EG1", "homology_id": "H1", "occurrence_id": "A_ex",
            "family_id": "fam", "species": "A", "gene_copy_id": "gA",
            "element_class": "exon_like", "membership_call": "core_member",
            "membership_edge_eligible": "1", "position_edge_eligible": "1",
            "correspondence_status": "resolved",
            "genomic_matched_blocks": "chrA:1-30:+", "parent_feature_ids": "exA",
        }]
        return occurrences, elements

    @staticmethod
    def _by_layer(rows, species):
        return {row["layer"]: row for row in rows if row["species"] == species}

    def _prediction_overlap_case(self, prediction_start, prediction_end, prediction_strand="+"):
        occurrences, elements = self._base_element()
        occurrences.append({
            "occurrence_id": "B_in", "family_id": "fam", "species": "B",
            "gene_copy_id": "gB", "role": "intron", "presence_status": "present",
            "contig": "chrB", "start": "50", "end": "250", "strand": "+",
        })
        elements.append({
            "element_id": "EG1", "homology_id": "H1", "occurrence_id": "B_in",
            "element_class": "candidate_source", "membership_call": "core_member",
            "membership_edge_eligible": "1", "position_edge_eligible": "1",
            "correspondence_status": "resolved", "genomic_matched_blocks": "chrB:100-120:+",
        })
        completion = [{
            "evidence_id": "predB", "homology_id": "H1", "family_id": "fam",
            "species": "B", "gene_copy_id": "gB", "homologous_dna_presence": "present",
            "homologous_dna_evidence": "protein_coding_projection",
            "candidate_resolution_status": "resolved", "predicted_role": "CDS",
            "predicted_role_blocks": json.dumps([{
                "block_id": "p", "target_contig": "chrB", "target_start": prediction_start,
                "target_end": prediction_end, "target_strand": prediction_strand,
            }]),
        }]
        rows = _element_site_rows(occurrences, elements, completion, {"fam"}, {"fam": {"A", "B"}})
        return self._by_layer(rows, "B")["exon_role"]

    def _junction_fixture(self, ambiguous_middle=False, missing_phase=False):
        root = Path(self.tempdir.name)
        input_dir, output_dir = root / "input", root / "output"
        input_dir.mkdir(); output_dir.mkdir()
        phase = "." if missing_phase else "0"
        paths = [
            ("R", "gR", "ref_tx", 1, "ref", "CDS", "0"),
            ("S", "gS", "split", 1, "left", "CDS", phase),
            ("S", "gS", "split", 2, "i1", "intron", "."),
            ("S", "gS", "split", 3, "middle", "CDS", phase),
            ("S", "gS", "split", 4, "i2", "intron", "."),
            ("S", "gS", "split", 5, "right", "CDS", phase),
            ("C", "gC", "continuous", 1, "continuous", "CDS", "0"),
            ("P", "gP", "partial", 1, "partial", "CDS", "0"),
            ("D", "gD", "repeat", 1, "repeat", "CDS", "0"),
        ]
        header = "path_id\tfamily_id\tspecies\tgene_copy_id\ttranscript_id\tpath_rank\toccurrence_id\trole\tphase\tpath_status\n"
        lines = [header] + [f"p{index}\tfam\t{sp}\t{gene}\t{tx}\t{rank}\t{occ}\t{role}\t{ph}\tannotated_transcript_path\n"
                            for index, (sp, gene, tx, rank, occ, role, ph) in enumerate(paths)]
        (input_dir / "transcript_paths.tsv").write_text("".join(lines))
        match_rows = [
            ("m1", "left", "ref", "1-30:1-30", "1", "resolved"),
            ("m2", "middle", "ref", "1-30:31-60", "0" if ambiguous_middle else "1",
             "ambiguous" if ambiguous_middle else "resolved"),
            ("m3", "right", "ref", "1-30:61-90", "1", "resolved"),
            ("m4", "continuous", "ref", "1-90:1-90;1-5:1-5", "1", "resolved"),
            ("m5", "partial", "ref", "1-30:1-30", "1", "resolved"),
            ("m6", "repeat", "ref", "1-90:1-90", "1", "resolved"),
        ]
        (output_dir / "segment_matches.tsv").write_text(
            "match_id\tquery_occurrence_id\tsubject_occurrence_id\tmatch_status\tmatched_blocks\talignment_strand\tposition_edge_eligible\tcandidate_resolution\n" +
            "".join(f"{mid}\t{query}\t{subject}\tmapped\t{blocks}\t+\t{eligible}\t{resolution}\n"
                    for mid, query, subject, blocks, eligible, resolution in match_rows))
        occurrences = []
        coordinates = {"ref": (1, 90), "left": (1, 30), "i1": (31, 40),
                       "middle": (41, 70), "i2": (71, 80), "right": (81, 110),
                       "continuous": (1, 90), "partial": (1, 30), "repeat": (1, 90)}
        species = {"ref": ("R", "gR"), "left": ("S", "gS"), "i1": ("S", "gS"),
                   "middle": ("S", "gS"), "i2": ("S", "gS"), "right": ("S", "gS"),
                   "continuous": ("C", "gC"), "partial": ("P", "gP"), "repeat": ("D", "gD")}
        for occurrence_id, (start, end) in coordinates.items():
            sp, gene = species[occurrence_id]
            occurrences.append({
                "occurrence_id": occurrence_id, "family_id": "fam", "species": sp,
                "gene_copy_id": gene, "role": "intron" if occurrence_id.startswith("i") else "CDS",
                "presence_status": "present", "contig": f"chr{sp}", "start": str(start),
                "end": str(end), "strand": "+", "phase": phase if sp == "S" else "0",
            })
        elements = []
        for occurrence_id in ("ref", "left", "middle", "right", "continuous", "partial", "repeat"):
            relation = "repeated_overlap" if occurrence_id == "repeat" else (
                "complementary" if occurrence_id in {"left", "middle", "right"} else "single")
            elements.append({
                "element_id": "EG1", "occurrence_id": occurrence_id,
                "element_class": "exon_like", "membership_call": "core_member",
                "position_edge_eligible": "0" if occurrence_id == "middle" and ambiguous_middle else "1",
                "correspondence_status": "ambiguous" if occurrence_id == "middle" and ambiguous_middle else "resolved",
                "reference_occurrence_id": "ref", "reference_coverage_relation": relation,
            })
        return input_dir, output_dir, occurrences, elements

    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tempdir.cleanup()

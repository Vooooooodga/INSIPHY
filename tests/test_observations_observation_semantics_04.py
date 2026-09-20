import json

import tempfile

import unittest

from pathlib import Path

from unittest.mock import patch

from intraphy.correspondence import _membership_match_details, infer_correspondence, match_total, tree_distances

from intraphy.coordinates import parse_legacy_blocks

from intraphy.elements import collect_element_profiles

from intraphy.io import read_tsv, write_tsv

from intraphy.cli import main as cli_main

from intraphy.preprocess import MATCH_FIELDS, _apply_ordered_candidate_chains, _species_tree_distances, assess_short_candidate_thresholds, cheap_match_evidence, cluster_segments, copy_order_context, extract_gene, graph_components, introns_from_path, match_evidence, read_annotation_for_gene, transcript_cds_length

from intraphy.alignment import AlignmentStats, phase_compatibility

from intraphy.structural_sites import (
    _completion_presence,
    _element_site_rows,
    _genomically_contiguous,
    _junction_site_rows,
)

from support_observations_observation_semantics import ObservationSemanticsTestsSupport

class ObservationSemanticsTests(ObservationSemanticsTestsSupport, unittest.TestCase):
    def test_membership_coverage_distinguishes_complementary_repeat_and_partial(self):
        occurrences = [
            {"occurrence_id": "ref", "family_id": "fam", "species": "A", "gene_copy_id": "gA", "transcript_id": "txA", "source_feature_id": "ref_feature", "start": "1", "end": "100"},
            {"occurrence_id": "b1", "family_id": "fam", "species": "B", "gene_copy_id": "gB", "transcript_id": "txB", "source_feature_id": "b1_feature", "start": "1", "end": "30", "path_roles": "exonic", "coding_roles": "CDS", "position_roles": "first", "annotation_source": "fixture", "original_attributes": "ID=b1", "partial_start": "0", "partial_end": "0", "path_role_records": '[{"transcript_id":"txB","path_role":"exonic"}]'},
            {"occurrence_id": "b2", "family_id": "fam", "species": "B", "gene_copy_id": "gB", "transcript_id": "txB", "source_feature_id": "b2_feature", "start": "31", "end": "60"},
            {"occurrence_id": "b3", "family_id": "fam", "species": "B", "gene_copy_id": "gB", "transcript_id": "txB", "source_feature_id": "b3_feature", "start": "61", "end": "100"},
            {"occurrence_id": "c1", "family_id": "fam", "species": "C", "gene_copy_id": "gC", "transcript_id": "txC", "source_feature_id": "c1_feature", "start": "1", "end": "40"},
            {"occurrence_id": "c2", "family_id": "fam", "species": "C", "gene_copy_id": "gC", "transcript_id": "txC", "source_feature_id": "c2_feature", "start": "101", "end": "140"},
            {"occurrence_id": "d1", "family_id": "fam", "species": "D", "gene_copy_id": "gD", "transcript_id": "txD", "source_feature_id": "d1_feature", "start": "1", "end": "40"},
        ]
        # A full-length independent match links the complete reference axis;
        # disjoint local memberships without this evidence are tested separately.
        occurrences.append({"occurrence_id": "e_ref", "family_id": "fam", "species": "E",
                            "gene_copy_id": "gE", "transcript_id": "txE",
                            "source_feature_id": "e_ref_feature", "start": "1", "end": "100"})
        element_rows = [
            {"element_id": "EG1", "occurrence_id": row["occurrence_id"], "species": row["species"], "gene_copy_id": row["gene_copy_id"]}
            for row in occurrences
        ]

        def match(match_id, query, blocks):
            parsed_blocks = parse_legacy_blocks(blocks)
            own = next(row for row in occurrences if row["occurrence_id"] == query)
            offset = int(own["start"]) - 1
            query_blocks = ";".join(
                f"chr{own['species']}:{offset + block.query.start0 + 1}-{offset + block.query.end0}:+"
                for block in parsed_blocks)
            reference_blocks = ";".join(
                f"chrA:{block.target.start0 + 1}-{block.target.end0}:+" for block in parsed_blocks)
            return {
                "match_id": match_id,
                "query_occurrence_id": query,
                "subject_occurrence_id": "ref",
                "match_status": "mapped",
                "candidate_resolution": "resolved",
                "candidate_ids": f"{match_id}.candidate_001",
                "retained_candidate_ids": f"{match_id}.candidate_001",
                "matched_blocks": blocks,
                "query_genomic_matched_blocks": query_blocks,
                "subject_genomic_matched_blocks": reference_blocks,
                "query_length": max(block.query.end0 for block in parsed_blocks),
                "target_length": 100,
            }

        scored = [
            match("m_e", "e_ref", "1-100:1-100"),
            match("m_b1", "b1", "1-30:1-30"),
            match("m_b2", "b2", "1-30:31-60"),
            match("m_b3", "b3", "1-40:61-100"),
            match("m_c1", "c1", "1-40:1-40"),
            match("m_c2", "c2", "1-40:1-40"),
            match("m_d1", "d1", "1-40:1-40"),
        ]
        _membership_match_details(element_rows, occurrences, scored)
        by_occurrence = {row["occurrence_id"]: row for row in element_rows}

        self.assertEqual(by_occurrence["b1"]["reference_coverage_relation"], "complementary_complete")
        self.assertEqual(by_occurrence["b2"]["reference_uncovered_bases"], 0)
        self.assertEqual(by_occurrence["b3"]["reference_covered_bases"], 100)
        self.assertEqual(by_occurrence["b1"]["reference_length"], 100)
        self.assertEqual(by_occurrence["b1"]["reference_length_source"], "aligned_reference_coordinates")
        self.assertEqual(by_occurrence["c1"]["reference_coverage_relation"], "repeated_overlap")
        self.assertNotEqual(by_occurrence["c1"]["repeat_instance_id"], "NA")
        self.assertEqual(by_occurrence["d1"]["reference_coverage_relation"], "single_partial")
        self.assertEqual(by_occurrence["d1"]["reference_uncovered_bases"], 60)
        self.assertEqual(by_occurrence["b1"]["parent_feature_ids"], "b1_feature")
        self.assertEqual(by_occurrence["b1"]["transcript_ids"], "txB")
        self.assertEqual(by_occurrence["b1"]["path_roles"], "exonic")
        self.assertEqual(by_occurrence["b1"]["coding_roles"], "CDS")
        self.assertEqual(by_occurrence["b1"]["position_roles"], "first")
        self.assertEqual(json.loads(by_occurrence["b1"]["path_role_records"])[0]["transcript_id"], "txB")

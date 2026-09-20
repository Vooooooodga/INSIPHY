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

from support_structural_site_semantics_structural_site_semantics import StructuralSiteSemanticsTestsSupport

class StructuralSiteSemanticsTests(StructuralSiteSemanticsTestsSupport, unittest.TestCase):
    def test_missing_coding_phase_is_recorded_without_erasing_exact_junction(self):
        fixture = self._junction_fixture(missing_phase=True)
        rows = _junction_site_rows(*fixture, {"fam"}, {"fam": {"R", "S", "C", "P", "D"}})
        split_rows = [row for row in rows if row["species"] == "S"]
        self.assertEqual({row["state"] for row in split_rows}, {"present"})
        self.assertTrue(all("junction_phase_unknown" in row["evidence"] for row in split_rows))

    def test_canonical_view_excludes_noncanonical_junction_paths(self):
        fixture = self._junction_fixture()
        path = fixture[0] / "transcript_paths.tsv"
        lines = path.read_text().splitlines()
        header = lines[0].split("\t")
        status_index = header.index("path_status")
        transcript_index = header.index("transcript_id")
        rewritten = [lines[0]]
        for line in lines[1:]:
            fields = line.split("\t")
            if fields[transcript_index] == "split":
                fields[status_index] = "canonical_transcript_path"
            rewritten.append("\t".join(fields))
        path.write_text("\n".join(rewritten) + "\n")

        rows = _junction_site_rows(
            *fixture,
            {"fam"},
            {"fam": {"R", "S", "C", "P", "D"}},
            annotation_view="canonical",
        )
        self.assertTrue(rows)
        self.assertEqual(
            {row["state"] for row in rows if row["species"] == "S"}, {"present"}
        )
        self.assertEqual(
            {row["state"] for row in rows if row["species"] == "C"}, {"unknown"}
        )
        self.assertEqual(
            {row["transcript_scope"] for row in rows},
            {"explicit_canonical_transcript"},
        )

    def test_site_only_complete_universe_entry_does_not_create_all_zero_state(self):
        rows = _merge_catalogue_observations(
            [],
            [{"family_id": "fam", "layer": "exon_presence", "site_id": "S_only"},
             {"family_id": "fam", "layer": "exon_presence", "site_id": "S_explicit",
              "species": "A", "state": "absent", "discovery_rule": "independent_catalogue",
              "observation_mask": "explicit_all_zero"},
             {"family_id": "fam", "layer": "exon_presence", "site_id": "S_generated",
              "species": "B", "state": "absent", "observation_mask": "generated_all_zero"}],
            {"fam"}, {"fam": {"A", "B"}}, annotation_view="repertoire",
        )
        by_site = {row["site_id"]: row for row in rows}
        self.assertEqual(set(by_site), {"S_explicit", "S_generated"})
        self.assertEqual(by_site["S_explicit"]["discovery_rule"], "independent_catalogue")
        self.assertEqual(by_site["S_generated"]["discovery_rule"], "independent_catalogue")
        self.assertEqual(by_site["S_explicit"]["observation_mask"], "explicit_all_zero")
        self.assertEqual(by_site["S_generated"]["observation_mask"], "generated_all_zero")

    def test_catalogue_role_view_must_match_analysis_and_is_preserved(self):
        base = {
            "family_id": "fam", "layer": "exon_role", "site_id": "S1",
            "species": "A", "state": "not_exonic", "state_0": "not_exonic",
            "state_1": "exonic", "annotation_view": "canonical",
            "transcript_scope": "explicit_canonical_transcript",
        }
        with self.assertRaisesRegex(ValueError, "conflicts with this analysis"):
            _merge_catalogue_observations(
                [], [base], {"fam"}, {"fam": {"A"}}, annotation_view="repertoire"
            )

        repertoire = {
            **base,
            "annotation_view": "repertoire",
            "transcript_scope": "annotated_transcript_repertoire",
        }
        with self.assertRaisesRegex(ValueError, "conflicts with this analysis"):
            _merge_catalogue_observations(
                [], [repertoire], {"fam"}, {"fam": {"A"}}, annotation_view="canonical"
            )

        rows = _merge_catalogue_observations(
            [], [base], {"fam"}, {"fam": {"A"}}, annotation_view="canonical"
        )
        self.assertEqual(rows[0]["annotation_view"], "canonical")

        conflicting_scope = {
            **base,
            "annotation_view": "repertoire",
        }
        with self.assertRaisesRegex(ValueError, "conflicts with transcript_scope"):
            _merge_catalogue_observations(
                [], [conflicting_scope], {"fam"}, {"fam": {"A"}},
                annotation_view="repertoire",
            )

    def test_catalogue_role_without_explicit_or_inferable_view_is_rejected(self):
        row = {
            "family_id": "fam", "layer": "splice_junction", "site_id": "J1",
            "state": "present", "state_0": "absent",
            "state_1": "present", "transcript_scope": "catalogue_defined",
        }
        with self.assertRaisesRegex(ValueError, "require an explicit or inferable"):
            _merge_catalogue_observations(
                [], [row], {"fam"}, {"fam": {"A"}}, annotation_view="canonical"
            )

    def test_catalogue_unknown_mask_is_not_replaced_with_observed(self):
        row = {
            "family_id": "fam", "layer": "exon_role", "site_id": "S1",
            "species": "A", "state": "exonic", "state_0": "not_exonic",
            "state_1": "exonic", "annotation_view": "repertoire",
            "transcript_scope": "annotated_transcript_repertoire",
            "observation_mask": "ambiguous",
        }
        merged = _merge_catalogue_observations(
            [], [row], {"fam"}, {"fam": {"A"}}, annotation_view="repertoire"
        )
        with self.assertRaisesRegex(ValueError, "invalid observation_mask"):
            normalize_structural_site_row(merged[0])

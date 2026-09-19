import json
import unittest
from unittest.mock import patch

from insiphy.alignment import AlignmentBackendError
from insiphy.coding_correspondence import CodingProjectionIndex, build_coding_transcript
from insiphy.correspondence import match_total
from insiphy.preprocess import cheap_match_evidence, cluster_segments
from insiphy.structural_sites import _mapped_reference_coordinate, _within_exon_boundary


def coding_rows(name, species, start, sequence, *, strand="+", phase="0", cds=None, transcript="tx"):
    end = start + len(sequence) - 1
    cds_left, cds_right = cds or (1, len(sequence))
    if strand == "-":
        coding_start, coding_end = end - cds_right + 1, end - cds_left + 1
    else:
        coding_start, coding_end = start + cds_left - 1, start + cds_right - 1
    occurrence = {
        "occurrence_id": name, "family_id": "fam", "species": species,
        "gene_copy_id": f"g{species}", "transcript_id": transcript,
        "role": "exon", "presence_status": "present", "contig": "chr",
        "start": str(start), "end": str(end), "strand": strand, "phase": phase,
        "transcript_order": "1", "splice_motif_score": "0.5",
    }
    path = {
        **occurrence, "path_id": f"{species}:{transcript}:{name}",
        "cds_phase": phase, "cds_length": str(coding_end - coding_start + 1),
        "cds_intervals": f"{coding_start}-{coding_end}",
        "path_status": "annotated_transcript_path", "coding_status": "coding",
    }
    return occurrence, path


def unchanged_protein_alignment(records, mode="linsi", threads=1):
    sequences = dict(records)
    if len({len(sequence) for sequence in sequences.values()}) > 1:
        raise AssertionError("This fixture requires equal-length aligned proteins")
    return sequences


def block_coordinates0(blocks):
    return tuple(
        (
            block.query.start0, block.query.end0,
            block.target.start0, block.target.end0,
        )
        for block in blocks
    )


def weak_dna_evidence(left, right, sequences, context, **kwargs):
    evidence = cheap_match_evidence(left, right, context, alignment_backend="mafft")
    evidence.update(
        alignment_score=0.4, coverage_score=1.0, sequence_score=0.58,
        total_score=evidence["total_score"] + 0.34 * 0.4 + 0.14,
        alignment_strand="+", alignment_cigar="3M", alignment_mode="overlap",
        alignment_meaning="deterministic DNA alignment fixture",
        alignment_requested_backend="mafft", query_alignment_start=1,
        query_alignment_end=3, target_alignment_start=1, target_alignment_end=3,
        projected_reference_blocks="1-3:1-3",
    )
    return evidence


class CodingCorrespondenceTests(unittest.TestCase):
    def split_fixture(self, strand="+"):
        left_start, right_start = (101, 201) if strand == "+" else (301, 201)
        left, left_path = coding_rows("left", "A", left_start, "ATGG", strand=strand)
        right, right_path = coding_rows("right", "A", right_start, "AAGCTTAATTTT", strand=strand, phase="2", cds=(1, 8))
        reference, reference_path = coding_rows("ref", "B", 501, "CCATGGAAGCTTAAGG", cds=(3, 14))
        occurrences = [left, right, reference]
        paths = [left_path, right_path, reference_path]
        sequences = {"left": "ATGG", "right": "AAGCTTAATTTT", "ref": "CCATGGAAGCTTAAGG"}
        return occurrences, paths, sequences

    def test_split_codon_utr_and_negative_strand_exact_projection(self):
        for strand in ("+", "-"):
            with self.subTest(strand=strand), patch(
                "insiphy.coding_correspondence.alignment.protein_multiple_alignment",
                side_effect=unchanged_protein_alignment,
            ) as aligner:
                occurrences, paths, sequences = self.split_fixture(strand)
                index = CodingProjectionIndex(occurrences, sequences, paths)
                left = index.evidence("left", "ref")
                right = index.evidence("right", "ref")
                reverse = index.evidence("ref", "right")
                self.assertEqual(index.transcripts[("fam", "A", "gA", "tx")].protein, "MEA")
                self.assertEqual(
                    block_coordinates0(left["protein_projected_blocks"]),
                    ((0, 4, 2, 6),),
                )
                self.assertEqual(
                    block_coordinates0(right["protein_projected_blocks"]),
                    ((0, 5, 6, 11),),
                )
                self.assertEqual(
                    block_coordinates0(reverse["protein_projected_blocks"]),
                    ((6, 11, 0, 5),),
                )
                block = right["protein_projected_coordinate_blocks"][0]
                self.assertEqual((block.query.start0, block.query.end0), (0, 5))
                self.assertEqual((block.target.start0, block.target.end0), (6, 11))
                self.assertEqual(right["protein_status"], "supported")
                self.assertEqual(right["protein_query_cds_coverage"], 5 / 8)
                self.assertEqual(right["protein_target_cds_coverage"], 5 / 12)
                self.assertEqual(reverse["protein_target_cds_coverage"], 5 / 8)
                right_bases = [
                    base
                    for residue in index.project("right")
                    for base in residue.cds_bases
                    if base.occurrence_id == "right" and base.occurrence_interval.start0 == 0
                ]
                self.assertTrue(right_bases)
                expected_genome_start0 = 200 if strand == "+" else 211
                self.assertEqual(right_bases[0].genome_interval.start0, expected_genome_start0)
                self.assertEqual(aligner.call_count, 1)
                self.assertEqual(aligner.call_args.args[0], (("coding_000001", "MEA"),))
                self.assertEqual(aligner.call_args.kwargs, {"mode": "linsi", "threads": 1})

    def test_initial_partial_codon_skipped_once_internal_phase_retained(self):
        first, first_path = coding_rows("a", "A", 1, "AATGG", phase="1")
        second, second_path = coding_rows("b", "A", 30, "AAGCT", phase="2")
        transcript = build_coding_transcript(
            ("fam", "A", "gA", "tx"), [first_path, second_path],
            {"a": first, "b": second}, {"a": "AATGG", "b": "AAGCT"},
        )
        self.assertEqual(transcript.protein, "MEA")
        self.assertEqual(transcript.codons[0], (("a", 2), ("a", 3), ("a", 4)))
        self.assertEqual(transcript.codons[1], (("a", 5), ("b", 1), ("b", 2)))

    def test_internal_stops_unknown_codons_and_protein_gaps_are_unmapped(self):
        query, query_path = coding_rows("q", "A", 1, "ATGNNNTAAGCTGCTTAA")
        target, target_path = coding_rows("t", "B", 1, "ATGGAAGCTGCTTAA")
        index = CodingProjectionIndex([query, target], {"q": "ATGNNNTAAGCTGCTTAA", "t": "ATGGAAGCTGCTTAA"}, [query_path, target_path])
        def gapped_family_alignment(records, mode="linsi", threads=1):
            return {
                record_id: {"MXXAA": "MXXAA", "MEAA": "ME-AA"}[sequence]
                for record_id, sequence in records
            }

        with patch(
            "insiphy.coding_correspondence.alignment.protein_multiple_alignment",
            side_effect=gapped_family_alignment,
        ):
            evidence = index.evidence("q", "t")
        self.assertEqual(index.transcripts[("fam", "A", "gA", "tx")].protein, "MXXAA")
        self.assertEqual(evidence["protein_status"], "supported")
        self.assertEqual(
            block_coordinates0(evidence["protein_projected_blocks"]),
            ((0, 3, 0, 3), (9, 15, 6, 12)),
        )
        self.assertEqual(evidence["protein_query_cds_coverage"], 9 / 18)
        self.assertEqual(evidence["protein_target_cds_coverage"], 9 / 15)

    def test_paths_have_their_own_cds_extents_and_unknown_phase_is_unavailable(self):
        occurrence, first_path = coding_rows("q", "A", 1, "ATGGCTGCTGCT", cds=(1, 9), transcript="tx1")
        second_path = {**first_path, "transcript_id": "tx2", "cds_intervals": "4-12", "cds_phase": "1"}
        third_path = {**first_path, "transcript_id": "tx_unknown", "cds_phase": "."}
        occurrence.update(cds_intervals="1-12", cds_length="12", cds_phase="0")
        index = CodingProjectionIndex([occurrence], {"q": "ATGGCTGCTGCT"}, [first_path, second_path, third_path])
        tx1 = index.transcripts[("fam", "A", "gA", "tx1")]
        tx2 = index.transcripts[("fam", "A", "gA", "tx2")]
        unknown = index.transcripts[("fam", "A", "gA", "tx_unknown")]
        self.assertEqual(tx1.protein, "MAA")
        self.assertEqual(tx1.coding_lengths["q"], 9)
        self.assertEqual(tx2.protein, "LL")
        self.assertEqual(tx2.codons[0], (("q", 5), ("q", 6), ("q", 7)))
        self.assertEqual(unknown.unavailable_reason, "unknown_CDS_phase")

    def test_raw_gene_metadata_and_per_interval_cds_phases(self):
        occurrence, path = coding_rows("q", "A", 1, "ATGGCCAAGCT", cds=(1, 4))
        path.update(cds_intervals="1-4;7-11", cds_length="9")
        raw = [
            {"id": "gA", "parent": "", "type": "gene", "attrs": {"transl_table": "1"}},
            {"id": "tx", "parent": "gA", "type": "mRNA", "attrs": {}},
            {"id": "c1", "parent": "tx", "type": "CDS", "seqid": "chr", "start": 1, "end": 4, "strand": "+", "phase": "0", "attrs": {}},
            {"id": "c2", "parent": "tx", "type": "CDS", "seqid": "chr", "start": 7, "end": 11, "strand": "+", "phase": "2", "attrs": {}},
        ]
        args = (("fam", "A", "gA", "tx"), [path], {"q": occurrence}, {"q": "ATGGCCAAGCT"})
        transcript = build_coding_transcript(*args, raw)
        self.assertEqual(transcript.protein, "MEA")
        self.assertEqual(transcript.codons[1], (("q", 4), ("q", 7), ("q", 8)))
        self.assertEqual(transcript.codon_sources[0][0].source_feature_ids, ("c1",))
        self.assertEqual(transcript.codon_sources[1][1].source_feature_ids, ("c2",))
        self.assertEqual(transcript.codon_sources[1][1].occurrence_interval.start0, 6)
        self.assertEqual(transcript.codon_sources[1][1].genome_interval.start0, 6)
        for attrs, expected in (
            ({"transl_table": "2"}, "unsupported_translation_table"),
            ({"transl_except": "(pos:4..6,aa:Sec)"}, "translation_exception"),
            ({"exception": "ribosomal slippage"}, "translation_exception"),
        ):
            with self.subTest(attrs=attrs):
                records = [{**raw[0], "attrs": attrs}, *raw[1:]]
                self.assertEqual(build_coding_transcript(*args, records).unavailable_reason, expected)
        conflicting = [*raw[:-1], {**raw[-1], "phase": "0"}]
        self.assertEqual(build_coding_transcript(*args, conflicting).unavailable_reason, "inconsistent_CDS_phase")

    def test_disagreeing_supported_isoform_coordinates_remain_ambiguous(self):
        query, first_path = coding_rows("q", "A", 1, "GCTGCTGCTGCT", cds=(1, 9), transcript="a")
        second_path = {**first_path, "transcript_id": "b", "cds_intervals": "4-12"}
        target, target_path = coding_rows("t", "B", 1, "GCTGCTGCT")
        index = CodingProjectionIndex([query, target], {"q": "GCTGCTGCTGCT", "t": "GCTGCTGCT"}, [first_path, second_path, target_path])
        with patch("insiphy.coding_correspondence.alignment.protein_multiple_alignment", side_effect=unchanged_protein_alignment):
            evidence = index.evidence("q", "t")
        self.assertEqual(evidence["protein_status"], "ambiguous_transcript_projection")
        self.assertEqual(evidence["protein_mapping_status"], "ambiguous_mapping")
        self.assertTrue(evidence["protein_candidate_evidence_available"])
        self.assertFalse(evidence["protein_membership_eligible"])
        self.assertFalse(evidence["protein_position_eligible"])
        self.assertFalse(evidence["protein_hard_observation_eligible"])
        self.assertNotIn("protein_projected_blocks", evidence)
        candidates = json.loads(evidence["protein_candidate_details"])
        self.assertEqual(
            [(row["query_transcript_id"], row["target_transcript_id"]) for row in candidates],
            [("a", "tx"), ("b", "tx")],
        )
        self.assertEqual(
            [[(block["query_start0"], block["query_end0"]) for block in row["projected_blocks"]]
             for row in candidates],
            [[(0, 9)], [(3, 12)]],
        )

    def test_anonymous_cds_are_selected_only_by_their_explicit_parent(self):
        occurrence, path = coding_rows("q", "A", 1, "ATGGCTGCT")
        raw = [
            {"id": "gA", "parent": "NA", "type": "gene", "attrs": {"transl_table": "1"}},
            {"id": "tx", "parent": "gA", "type": "mRNA", "attrs": {}},
            {"id": "other", "parent": "gA", "type": "mRNA", "attrs": {}},
            {"id": "NA", "parent": "tx", "type": "CDS", "seqid": "chr", "start": 1, "end": 9, "strand": "+", "phase": "0", "attrs": {}},
            {"id": "NA", "parent": "other", "type": "CDS", "seqid": "chr", "start": 1, "end": 9, "strand": "+", "phase": "1", "attrs": {"transl_table": "2"}},
            {"id": "NA", "parent": "NA", "type": "CDS", "seqid": "chr", "start": 1, "end": 9, "strand": "+", "phase": "2", "attrs": {"transl_except": "unsupported"}},
        ]
        transcript = build_coding_transcript(
            ("fam", "A", "gA", "tx"), [path], {"q": occurrence}, {"q": "ATGGCTGCT"}, raw,
        )
        self.assertEqual(transcript.unavailable_reason, "")
        self.assertEqual(transcript.protein, "MAA")
        other = build_coding_transcript(
            ("fam", "A", "gA", "other"), [{**path, "transcript_id": "other"}],
            {"q": occurrence}, {"q": "ATGGCTGCT"}, raw,
        )
        self.assertEqual(other.unavailable_reason, "unsupported_translation_table")

    def test_local_mapping_reports_identity_without_a_fixed_identity_gate(self):
        query, query_path = coding_rows("q", "A", 1, "ATGGCTGCT")
        target, target_path = coding_rows("t", "B", 1, "ATGGAAGCT")
        index = CodingProjectionIndex([query, target], {"q": "ATGGCTGCT", "t": "ATGGAAGCT"}, [query_path, target_path])
        with patch("insiphy.coding_correspondence.alignment.protein_multiple_alignment", side_effect=unchanged_protein_alignment):
            evidence = index.evidence("q", "t")
        self.assertEqual(evidence["protein_aa_identity"], 2 / 3)
        self.assertEqual(evidence["protein_query_cds_coverage"], 1.0)
        self.assertEqual(evidence["protein_status"], "supported")
        self.assertEqual(evidence["protein_mapping_status"], "supported_unanchored")
        self.assertFalse(evidence["protein_hard_observation_eligible"])

    def test_all_annotated_cds_bases_remain_in_coverage_denominator(self):
        for sequence, phase, expected in (("ATGTAA", "0", 3 / 6), ("AATGTAA", "1", 3 / 7), ("ATGNNN", "0", 3 / 6)):
            with self.subTest(sequence=sequence), patch(
                "insiphy.coding_correspondence.alignment.protein_multiple_alignment",
                side_effect=unchanged_protein_alignment,
            ):
                query, query_path = coding_rows("q", "A", 1, sequence, phase=phase)
                target, target_path = coding_rows("t", "B", 1, sequence, phase=phase)
                index = CodingProjectionIndex([query, target], {"q": sequence, "t": sequence}, [query_path, target_path])
                evidence = index.evidence("q", "t")
                self.assertEqual(evidence["protein_aa_identity"], 1.0)
                self.assertEqual(evidence["protein_query_cds_coverage"], expected)
                self.assertEqual(evidence["protein_target_cds_coverage"], expected)
                self.assertEqual(evidence["protein_status"], "supported")
                self.assertEqual(evidence["protein_mapping_status"], "supported_unanchored")
                self.assertIn("protein_projected_blocks", evidence)

    def test_synthetic_path_preserves_gene_and_exact_cds_translation_metadata(self):
        occurrence, path = coding_rows("q", "A", 1, "ATGGCTGCT", transcript="gA.synthetic_tx")
        owned = {"ownership": "target_gene_descendant"}
        gene = {**owned, "id": "gA", "parent": "NA", "type": "gene", "attrs": {"transl_table": "1"}}
        exon = {**owned, "id": "ex1", "parent": "gA", "type": "exon", "seqid": "chr", "start": 1, "end": 9, "strand": "+", "attrs": {}}
        cds = {**owned, "id": "NA", "parent": "gA", "type": "CDS", "seqid": "chr", "start": 1, "end": 9, "strand": "+", "phase": "0", "attrs": {"transl_table": "1"}}
        args = (("fam", "A", "gA", "gA.synthetic_tx"), [path], {"q": occurrence}, {"q": "ATGGCTGCT"})
        for parent in ("gA", "ex1"):
            with self.subTest(cds_parent=parent):
                transcript = build_coding_transcript(*args, [gene, exon, {**cds, "parent": parent}])
                self.assertEqual(transcript.unavailable_reason, "")
                self.assertEqual(transcript.protein, "MAA")
        cases = (
            ([{**gene, "attrs": {"transl_table": "2"}}, exon, cds], "unsupported_translation_table"),
            ([gene, exon, {**cds, "attrs": {"transl_except": "(pos:4..6,aa:Sec)"}}], "translation_exception"),
            ([gene, exon, {**cds, "phase": "."}], "unknown_CDS_phase"),
            ([gene, exon, {**cds, "end": 6}], "unresolved_raw_transcript_ownership"),
            ([gene, exon, {**cds, "parent": "other_tx"}], "unresolved_raw_transcript_ownership"),
            ([gene, exon, {**cds, "ownership": "overlapping_context"}], "unresolved_raw_transcript_ownership"),
            ([gene, exon, cds, {**owned, "id": "other_tx", "parent": "gA", "type": "mRNA", "attrs": {}}], "unresolved_raw_transcript_ownership"),
        )
        for raw, expected in cases:
            with self.subTest(expected=expected, raw=raw):
                self.assertEqual(build_coding_transcript(*args, raw).unavailable_reason, expected)

    def test_mutually_exclusive_isoform_coordinates_are_not_combined(self):
        query, first_path = coding_rows("q", "A", 1, "ATGGCTGCTGCT", cds=(1, 9), transcript="a")
        second_path = {**first_path, "transcript_id": "b", "cds_intervals": "4-12"}
        target, target_path = coding_rows("t", "B", 1, "ATGGCTGCTGCT", transcript="ref_tx")
        index = CodingProjectionIndex([query, target], {"q": "ATGGCTGCTGCT", "t": "ATGGCTGCTGCT"}, [first_path, second_path, target_path])

        def alignment_by_transcript(records, mode="linsi", threads=1):
            return {
                record_id: {"MAA": "MAA-", "AAA": "-AAA", "MAAA": "MAAA"}[sequence]
                for record_id, sequence in records
            }

        with patch("insiphy.coding_correspondence.alignment.protein_multiple_alignment", side_effect=alignment_by_transcript):
            evidence = index.evidence("q", "t")
        self.assertEqual(evidence["protein_mapping_status"], "ambiguous_mapping")
        self.assertTrue(evidence["protein_candidate_evidence_available"])
        self.assertFalse(evidence["protein_membership_eligible"])
        self.assertFalse(evidence["protein_position_eligible"])
        self.assertFalse(evidence["protein_hard_observation_eligible"])
        self.assertNotIn("protein_projected_blocks", evidence)
        self.assertEqual(evidence["protein_metrics_scope"], "best_transcript_pair")
        self.assertEqual(evidence["protein_best_query_transcript"], "a")
        self.assertEqual(evidence["protein_best_target_transcript"], "ref_tx")
        self.assertEqual(evidence["protein_query_cds_coverage"], 1.0)
        self.assertEqual(evidence["protein_target_cds_coverage"], 0.75)
        self.assertEqual(evidence["protein_supporting_transcripts"], "a>ref_tx;b>ref_tx")
        candidate_blocks = {
            row["query_transcript_id"]: row["projected_blocks"]
            for row in json.loads(evidence["protein_candidate_details"])
        }
        self.assertEqual(
            candidate_blocks,
            {
                "a": [{"query_start0": 0, "query_end0": 9, "target_start0": 0, "target_end0": 9}],
                "b": [{"query_start0": 3, "query_end0": 12, "target_start0": 3, "target_end0": 12}],
            },
        )

    def test_family_msa_deduplicates_identical_proteins_and_preserves_aliases(self):
        occurrences, paths, sequences = [], [], {}
        for species, offset in (("A", 0), ("B", 1000)):
            for label, amino_acids, start in (("left", 9, 1), ("middle", 2, 101), ("right", 9, 201)):
                name = f"{species}_{label}"
                sequence = "GCT" * amino_acids
                occurrence, path = coding_rows(name, species, start + offset, sequence)
                occurrences.append(occurrence)
                paths.append(path)
                sequences[name] = sequence

        with patch(
            "insiphy.coding_correspondence.alignment.protein_multiple_alignment",
            side_effect=unchanged_protein_alignment,
        ) as aligner:
            index = CodingProjectionIndex(occurrences, sequences, paths, threads=6)
            evidence = index.evidence("A_middle", "B_middle")
            projection = index.family_projection("fam")
            residues = index.project("A_middle")

        self.assertEqual(aligner.call_count, 1)
        self.assertEqual(aligner.call_args.kwargs, {"mode": "linsi", "threads": 6})
        self.assertEqual(len(projection.aligned_records), 1)
        aliases = next(iter(projection.aliases_by_record.values()))
        self.assertEqual(aliases, (("fam", "A", "gA", "tx"), ("fam", "B", "gB", "tx")))
        self.assertEqual(evidence["protein_mapping_status"], "resolved_local")
        self.assertTrue(evidence["protein_candidate_evidence_available"])
        self.assertTrue(evidence["protein_membership_eligible"])
        self.assertTrue(evidence["protein_position_eligible"])
        self.assertTrue(evidence["protein_hard_observation_eligible"])
        self.assertEqual(evidence["protein_left_anchor_pairs"], 9)
        self.assertEqual(evidence["protein_right_anchor_pairs"], 9)
        self.assertEqual(evidence["protein_msa_column_interval"], "9:11")
        self.assertEqual(evidence["protein_msa_column_start0"], 9)
        self.assertEqual(evidence["protein_msa_column_end0"], 11)
        self.assertEqual([item.residue_index0 for item in residues], [9, 10])
        self.assertEqual([item.msa_column0 for item in residues], [9, 10])
        first_base = residues[0].cds_bases[0]
        self.assertEqual((first_base.occurrence_interval.start0, first_base.occurrence_interval.end0), (0, 1))
        self.assertEqual((first_base.genome_interval.start0, first_base.genome_interval.end0), (100, 101))
        self.assertEqual(first_base.source_path_id, "A:tx:A_middle")

    def test_identical_protein_aliases_keep_paths_and_accept_identical_blocks(self):
        occurrences, paths, sequences = [], [], {}
        for species, offset in (("A", 0), ("B", 1000)):
            species_paths = []
            for label, amino_acids, start in (("left", 9, 1), ("middle", 2, 101), ("right", 9, 201)):
                name = f"{species}_{label}"
                sequence = "GCT" * amino_acids
                occurrence, path = coding_rows(name, species, start + offset, sequence)
                occurrences.append(occurrence)
                species_paths.append(path)
                sequences[name] = sequence
            paths.extend(species_paths)
            paths.extend({
                **path,
                "transcript_id": "tx_alias",
                "path_id": path["path_id"].replace(":tx:", ":tx_alias:"),
            } for path in species_paths)

        with patch(
            "insiphy.coding_correspondence.alignment.protein_multiple_alignment",
            side_effect=unchanged_protein_alignment,
        ):
            index = CodingProjectionIndex(occurrences, sequences, paths)
            evidence = index.evidence("A_middle", "B_middle")
            projection = index.family_projection("fam")

        aliases = next(iter(projection.aliases_by_record.values()))
        self.assertEqual(len(aliases), 4)
        self.assertEqual(evidence["protein_candidate_mapping_count"], 4)
        self.assertTrue(evidence["protein_candidate_coordinate_consensus"])
        self.assertEqual(evidence["protein_mapping_status"], "resolved_local")
        self.assertTrue(evidence["protein_position_eligible"])
        self.assertEqual(
            block_coordinates0(evidence["protein_projected_blocks"]),
            ((0, 6, 0, 6),),
        )

    def test_focal_occurrence_cannot_supply_its_own_anchor(self):
        query, query_path = coding_rows("q", "A", 1, "GCT" * 11)
        target_focal, focal_path = coding_rows("t_focal", "B", 1, "GCT" * 2)
        target_other, other_path = coding_rows("t_other", "B", 101, "GCT" * 9)
        index = CodingProjectionIndex(
            [query, target_focal, target_other],
            {"q": "GCT" * 11, "t_focal": "GCT" * 2, "t_other": "GCT" * 9},
            [query_path, focal_path, other_path],
        )
        with patch(
            "insiphy.coding_correspondence.alignment.protein_multiple_alignment",
            side_effect=unchanged_protein_alignment,
        ):
            evidence = index.evidence("q", "t_focal")

        self.assertEqual(evidence["protein_right_anchor_pairs"], 0)
        self.assertEqual(evidence["protein_terminal_side"], "left")
        self.assertEqual(evidence["protein_mapping_status"], "supported_unanchored")
        self.assertTrue(evidence["protein_candidate_evidence_available"])
        self.assertFalse(evidence["protein_membership_eligible"])
        self.assertFalse(evidence["protein_position_eligible"])
        self.assertFalse(evidence["protein_hard_observation_eligible"])

    def test_uncovered_family_columns_are_reported_without_coordinates(self):
        query, query_path = coding_rows("q", "A", 1, "GCTGCTGCT")
        target, target_path = coding_rows("t", "B", 1, "GAAGAAGAA")

        def disjoint_alignment(records, mode="linsi", threads=1):
            return {
                record_id: {"AAA": "AAA---", "EEE": "---EEE"}[sequence]
                for record_id, sequence in records
            }

        index = CodingProjectionIndex(
            [query, target], {"q": "GCTGCTGCT", "t": "GAAGAAGAA"},
            [query_path, target_path],
        )
        with patch(
            "insiphy.coding_correspondence.alignment.protein_multiple_alignment",
            side_effect=disjoint_alignment,
        ):
            evidence = index.evidence("q", "t")
        self.assertEqual(evidence["protein_mapping_status"], "uncovered")
        self.assertEqual(evidence["protein_status"], "no_aligned_CDS")
        self.assertNotIn("protein_projected_blocks", evidence)

    def test_protein_blocks_enter_graph_and_junction_preserving_raw_dna(self):
        for strand in ("+", "-"):
            with self.subTest(strand=strand), patch(
                "insiphy.preprocess.match_evidence", side_effect=weak_dna_evidence,
            ), patch(
                "insiphy.coding_correspondence.alignment.protein_multiple_alignment",
                side_effect=unchanged_protein_alignment,
            ):
                occurrences, paths, sequences = self.split_fixture(strand)
                homology, matches = cluster_segments(
                    occurrences, sequences, transcript_paths=paths,
                )
                by_occurrence = {row["occurrence_id"]: row["homology_id"] for row in homology}
                self.assertEqual(len(set(by_occurrence.values())), 1)
                for match in matches:
                    self.assertEqual(match["alignment_score"], "0.4")
                    self.assertEqual(match["projected_reference_blocks"], "1-3:1-3")
                    self.assertEqual(match["alignment_cigar"], "3M")
                    self.assertEqual(match["dna_match_status"], "low_similarity")
                    self.assertEqual(match["match_status"], "mapped")
                    self.assertEqual(match["correspondence_basis"], "annotated_CDS_protein")
                    self.assertEqual(match["protein_metrics_scope"], "best_transcript_pair")
                    self.assertEqual(match["protein_best_query_transcript"], "tx")
                    self.assertEqual(match["protein_best_target_transcript"], "tx")
                    self.assertGreater(match_total(match), float(match["total_score"]))
                boundary = _within_exon_boundary(
                    "fam", "EG_1", "left", "right", {row["occurrence_id"]: row for row in occurrences},
                    {"EG_1": "ref"}, matches,
                )
                self.assertIsNotNone(boundary)
                self.assertEqual((boundary["donor_projection"], boundary["acceptor_projection"]), (6, 7))
                self.assertIsNone(_mapped_reference_coordinate(matches, "right", "ref", 12))

    def test_repeated_reference_coverage_does_not_merge_distinct_exons(self):
        first, first_path = coding_rows("first", "A", 1, "ATGGCTGCT", transcript="one")
        second, second_path = coding_rows("second", "A", 101, "ATGGCTGCT", transcript="two")
        reference, reference_path = coding_rows("ref", "B", 1, "ATGGCTGCT")
        with patch("insiphy.preprocess.match_evidence", side_effect=weak_dna_evidence), patch(
            "insiphy.coding_correspondence.alignment.protein_multiple_alignment", side_effect=unchanged_protein_alignment,
        ):
            index = CodingProjectionIndex(
                [first, second, reference],
                {"first": "ATGGCTGCT", "second": "ATGGCTGCT", "ref": "ATGGCTGCT"},
                [first_path, second_path, reference_path],
            )
            evidence = index.evidence("first", "ref")
            homology, _matches = cluster_segments(
                [first, second, reference], {"first": "ATGGCTGCT", "second": "ATGGCTGCT", "ref": "ATGGCTGCT"},
                transcript_paths=[first_path, second_path, reference_path],
            )
        self.assertEqual(evidence["protein_mapping_status"], "ambiguous_mapping")
        self.assertFalse(evidence["protein_position_eligible"])
        self.assertEqual(evidence["protein_competing_occurrences"], "query:second")
        by_occurrence = {row["occurrence_id"]: row["homology_id"] for row in homology}
        self.assertNotEqual(by_occurrence["first"], by_occurrence["second"])

    def test_unavailable_metadata_preserves_dna_and_internal_backend_never_calls_protein(self):
        occurrences, paths, sequences = self.split_fixture()
        paths[0]["cds_phase"] = "."
        for backend, expected in (("mafft", "unavailable"), ("internal", "disabled_backend")):
            with self.subTest(backend=backend), patch(
                "insiphy.preprocess.match_evidence", side_effect=weak_dna_evidence,
            ), patch("insiphy.coding_correspondence.alignment.protein_multiple_alignment") as aligner:
                _homology, matches = cluster_segments(occurrences, sequences, transcript_paths=paths, aligner=backend)
                aligner.assert_not_called()
                self.assertTrue(all(row["protein_status"] == expected for row in matches))
                self.assertTrue(all(row["alignment_score"] == "0.4" and row["match_status"] == "low_similarity" for row in matches))
                if backend == "mafft":
                    self.assertTrue(all("unknown_CDS_phase" in row["protein_unavailable_reason"] for row in matches))

    def test_external_aligner_error_propagates(self):
        occurrences, paths, sequences = self.split_fixture()
        with patch("insiphy.preprocess.match_evidence", side_effect=weak_dna_evidence), patch(
            "insiphy.coding_correspondence.alignment.protein_multiple_alignment", side_effect=AlignmentBackendError("explicit MAFFT failure"),
        ):
            with self.assertRaisesRegex(AlignmentBackendError, "explicit MAFFT failure"):
                cluster_segments(occurrences, sequences, transcript_paths=paths)

    def test_cache_is_run_local_and_bounded_by_mapped_bases(self):
        occurrences, paths, sequences = self.split_fixture()
        index = CodingProjectionIndex(occurrences, sequences, paths)
        index.MAX_CACHED_BASE_PAIRS = 1
        with patch("insiphy.coding_correspondence.alignment.protein_multiple_alignment", side_effect=unchanged_protein_alignment):
            self.assertEqual(index.evidence("left", "ref")["protein_status"], "supported")
        self.assertFalse(index.cache)
        self.assertEqual(index.cached_bases, 0)
        self.assertFalse(CodingProjectionIndex(occurrences, sequences, paths).cache)


if __name__ == "__main__":
    unittest.main()

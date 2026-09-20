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


class ObservationSemanticsTestsSupport:
    def _intronic_cds_prediction_case(self):
        occurrences = [
            {
                "occurrence_id": "Bign_exon", "family_id": "OG0006454",
                "species": "Bombus_ignitus", "gene_copy_id": "gene:ENSBIGG00000010441",
                "role": "exon", "presence_status": "present", "contig": "LG5",
                "start": "1947147", "end": "1947446", "strand": "-",
            },
            {
                "occurrence_id": "Bter_intron", "family_id": "OG0006454",
                "species": "Bombus_terrestris", "gene_copy_id": "gene-LOC100645322",
                "role": "intron", "presence_status": "present", "contig": "NC_063273.1",
                "start": "2472297", "end": "2507265", "strand": "-",
            },
        ]
        elements = [
            {"element_id": "EG_0029", "homology_id": "HC_0029", "occurrence_id": "Bign_exon", "element_class": "exon_like", "membership_call": "core_member", "genomic_matched_blocks": "LG5:1947147-1947446:-"},
            {"element_id": "EG_0029", "homology_id": "HC_0029", "occurrence_id": "Bter_intron", "element_class": "candidate_source", "membership_call": "core_member", "genomic_matched_blocks": "NC_063273.1:2472297-2507265:-"},
        ]
        candidate = {
            "family_id": "OG0006454", "homology_id": "HC_0029",
            "species": "Bombus_terrestris", "gene_copy_id": "gene-LOC100645322",
            "interval": "NC_063273.1:2476398-2476441:-",
            "completion_call": "predicted_exon_candidate",
            "evidence_status": "supports_hidden_segment",
            "annotation_status": "protein_projection_supports_missing_cds",
            "inferred_event": "protein_cds_projection", "inferred_role": "predicted_CDS",
            "predicted_role": "CDS", "support_score": "0.8667",
            "homologous_dna_presence": "present",
            "homologous_dna_evidence": "resolved_protein_coding_projection",
            "candidate_resolution_status": "resolved",
            "primary_mapping_status": "protein_projection", "correspondence_status": "resolved",
            "interval_scope": "projected_cds_container", "confidence_flag": "high",
        }
        return occurrences, elements, candidate

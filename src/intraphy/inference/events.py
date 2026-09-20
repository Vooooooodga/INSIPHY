"""inference / events: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations




def _transition_event_type(layer, source_index, target_index):
    if source_index == target_index:
        return "no_change"
    if layer == "splice_junction":
        return "intron_gain" if (source_index, target_index) == (0, 1) else "intron_loss"
    if layer == "exon_presence":
        return "sequence_gain" if (source_index, target_index) == (0, 1) else "sequence_loss"
    if layer == "exon_role":
        return "exon_role_gain" if (source_index, target_index) == (0, 1) else "exon_role_loss"
    return "state_gain" if (source_index, target_index) == (0, 1) else "state_loss"


def _transition_structural_relation(layer, site_kind, source_index, target_index):
    if layer != "splice_junction" or site_kind not in {
        "within_element_junction",
        "within_exon_boundary",
    }:
        return "NA"
    if (source_index, target_index) == (0, 1):
        return "exon_split"
    if (source_index, target_index) == (1, 0):
        return "exon_fusion"
    return "NA"

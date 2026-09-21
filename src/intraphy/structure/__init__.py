"""Exon configurations and elementary structural edits (V19).

Alignment coordinates describe a correspondence hypothesis, not independent
characters. A configuration contains complete, ordered exon spans and the
material required to interpret insertions and irreversible deletions.
"""
from .types import ExonConfiguration, ExonSpan, Material, ObservationEvidence

MODEL_VERSION = "exon_configuration_v2"
SCHEMA_VERSION = "intraphy.exon-configurations/2"
__all__ = ["ExonConfiguration", "ExonSpan", "Material", "ObservationEvidence"]

"""aligners / scoring: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations




def splice_motif_score(oriented_intron_seq: str) -> tuple[float, str, str]:
    seq = (oriented_intron_seq or "").upper()
    if len(seq) < 4:
        return 0.0, "NA", "NA"
    donor = seq[:2]
    acceptor = seq[-2:]
    pair = f"{donor}-{acceptor}"
    if pair == "GT-AG":
        return 1.0, donor, acceptor
    if pair in {"GC-AG", "AT-AC"}:
        return 0.8, donor, acceptor
    if donor in {"GT", "GC", "AT"} or acceptor in {"AG", "AC"}:
        return 0.45, donor, acceptor
    return 0.1, donor, acceptor


def phase_compatibility(left_phase: str, right_phase: str, left_cds_length=None) -> str:
    if left_phase in {"", ".", "NA", None} or right_phase in {"", ".", "NA", None} or left_cds_length is None:
        return "unknown"
    try:
        left_phase = int(left_phase)
        right_phase = int(right_phase)
        coding_bases = int(left_cds_length) - left_phase
    except (TypeError, ValueError):
        return "unknown"
    if coding_bases < 0:
        return "incompatible"
    expected_right_phase = (3 - (coding_bases % 3)) % 3
    return "compatible" if right_phase == expected_right_phase else "incompatible"

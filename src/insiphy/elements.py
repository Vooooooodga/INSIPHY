"""Shared exon-like element semantics for INSIPHY."""

from collections import defaultdict

from .io import norm_state, to_float


EXON_LIKE_ROLES = {"CDS", "exon", "UTR", "noncoding_exon"}
NONCODING_ROLES = {"intron", "intergenic"}
HIDDEN_COMPLETION_CALLS = {
    "hidden_segment_candidate",
    "shifted_splice_site_candidate",
    "joined_exon_candidate",
    "hidden_segment_with_frame_disruption",
    "predicted_exon_candidate",
}
HIDDEN_EVIDENCE_STATUS = {"supports_hidden_segment", "conflicts_annotation"}
ROLE_ALIASES = {
    "predicted_CDS": "CDS",
    "predicted_cds": "CDS",
}


def element_id_from_homology(homology_id):
    """Return the user-facing exon-like element id for an internal component id."""
    if homology_id.startswith("HC_"):
        return f"EG_{homology_id.split('_', 1)[1]}"
    legacy_prefix = "HS" + "G_"
    if homology_id.startswith(legacy_prefix):
        return f"EG_{homology_id.split('_', 1)[1]}"
    if homology_id.startswith("H_"):
        return f"EG_{homology_id.split('_', 1)[1]}"
    return f"EG_{homology_id}"


def role_bucket(role):
    role = ROLE_ALIASES.get(role, role)
    if role == "CDS":
        return "CDS"
    if role in {"exon", "UTR", "noncoding_exon"}:
        return "exon_or_UTR"
    if role in NONCODING_ROLES:
        return "non_exonic_source"
    return "unknown"


def evidence_role(row):
    predicted = row.get("predicted_role", "")
    inferred = row.get("inferred_role", "")
    return ROLE_ALIASES.get(predicted, predicted) or ROLE_ALIASES.get(inferred, inferred)


def evidence_promotes_element(row):
    status = row.get("evidence_status", "")
    completion = row.get("completion_call", "")
    candidate_role = evidence_role(row)
    return (
        candidate_role in EXON_LIKE_ROLES
        or (completion in HIDDEN_COMPLETION_CALLS and candidate_role in EXON_LIKE_ROLES)
        or (status in HIDDEN_EVIDENCE_STATUS and candidate_role in EXON_LIKE_ROLES)
    )


def collect_element_profiles(homology, occurrences, evidence_rows=None):
    """Classify internal homology groups into public EGs and context groups."""
    occ_by_id = {row["occurrence_id"]: row for row in occurrences}
    profiles = defaultdict(lambda: {"roles": set(), "families": set(), "candidate": False})
    for row in homology:
        hid = row["homology_id"]
        occ = occ_by_id.get(row["occurrence_id"], {})
        role = occ.get("role", "unknown")
        profiles[hid]["roles"].add(role)
        if occ.get("family_id"):
            profiles[hid]["families"].add(occ["family_id"])
    for row in evidence_rows or []:
        hid = row.get("homology_id")
        if not hid:
            continue
        profiles[hid]["roles"].add(evidence_role(row) or row.get("inferred_role", "unknown"))
        if row.get("family_id"):
            profiles[hid]["families"].add(row["family_id"])
        if evidence_promotes_element(row):
            profiles[hid]["candidate"] = True

    element_by_homology = {}
    for hid, profile in profiles.items():
        if profile["roles"] & EXON_LIKE_ROLES or profile["candidate"]:
            element_by_homology[hid] = element_id_from_homology(hid)
    return profiles, element_by_homology


def element_class_for_occurrence(occ, homology_id, profiles, element_by_homology):
    if homology_id not in element_by_homology:
        return "context"
    role = ROLE_ALIASES.get(occ.get("role", "unknown"), occ.get("role", "unknown"))
    if norm_state(occ.get("presence_status")) == "absent":
        return "absent"
    if role in EXON_LIKE_ROLES:
        return "exon_like"
    return "candidate_source"


def element_role_from_occurrences(rows):
    roles = {ROLE_ALIASES.get(row.get("role", "unknown"), row.get("role", "unknown")) for row in rows if norm_state(row.get("presence_status")) == "present"}
    if not roles:
        return "absent"
    if "CDS" in roles:
        return "CDS"
    if roles & {"UTR", "noncoding_exon", "exon"}:
        return "exon_or_UTR"
    if roles & NONCODING_ROLES:
        return "non_exonic_source"
    return "unknown"


def confidence_score(confidence):
    text = str(confidence or "").strip().lower()
    if text in {"high", "strong", "curated"}:
        return 0.95
    if text in {"medium", "moderate"}:
        return 0.75
    if text in {"low", "weak"}:
        return 0.45
    return to_float(confidence, 0.5)


def membership_call(confidence, degree=0):
    base = confidence_score(confidence)
    adjusted = min(1.0, 0.75 * base + 0.25 * min(1.0, degree / 2))
    call = "core_member" if adjusted >= 0.7 else "ambiguous_member"
    return adjusted, call

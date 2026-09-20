"""mapping / ordered paths: explicit implementation ownership."""
from __future__ import annotations

from intraphy.mapping.chain_graph import _path_membership_for_context
from intraphy.mapping.chain_graph import _solve_context
from intraphy.mapping.chain_types import CandidateChainResult
from intraphy.mapping.chain_types import ChainCandidate
from intraphy.mapping.chain_types import ChainContextSummary
from intraphy.mapping.chain_types import DEFAULT_CHAIN_CONFIGURATION
from typing import Iterable
from typing import Mapping
from typing import Optional
from typing import Tuple


def ordered_candidate_chain(
    candidates: Iterable[ChainCandidate],
    score_delta: float,
    start_ids: Optional[Iterable[str]] = None,
    end_ids: Optional[Iterable[str]] = None,
    context_anchor_ids: Optional[
        Mapping[tuple, Tuple[Iterable[str], Iterable[str]]]
    ] = None,
    configuration_name: str = DEFAULT_CHAIN_CONFIGURATION.name,
):
    """Return candidates participating in near-optimal monotone chains.

    Coordinates must already be normalized to the transcriptional direction.
    Scores are comparable only within one named scoring scheme. Supplying start
    or end IDs requests an anchored path; an unreachable anchor is an error in
    the caller's candidate construction rather than evidence of sequence loss.
    """

    candidates = tuple(
        sorted(
            candidates,
            key=lambda item: (
                item.query.start0,
                item.query.end0,
                item.target.start0,
                item.target.end0,
                item.candidate_id,
            ),
        )
    )
    if not candidates:
        return CandidateChainResult(
            0.0,
            frozenset(),
            frozenset(),
            float(score_delta),
            True,
            configuration_name=configuration_name,
        )
    if score_delta < 0:
        raise ValueError("score_delta must be non-negative")
    schemes = {item.score_scheme for item in candidates}
    if len(schemes) != 1:
        raise ValueError("candidate-chain scores must use one scoring scheme")
    identifiers = [item.candidate_id for item in candidates]
    if len(set(identifiers)) != len(identifiers):
        raise ValueError("candidate IDs must be unique within a chain")
    start_ids = frozenset(start_ids or ())
    end_ids = frozenset(end_ids or ())
    unknown_anchors = (start_ids | end_ids) - set(identifiers)
    if unknown_anchors:
        raise ValueError("unknown chain anchors: " + ", ".join(sorted(unknown_anchors)))

    has_path_memberships = any(candidate.path_memberships for candidate in candidates)
    contexts = sorted(
        {
            membership.context
            for candidate in candidates
            for membership in candidate.path_memberships
        }
    ) if has_path_memberships else [None]
    anchors_by_context = {}
    for context in contexts:
        if context_anchor_ids is not None and context not in context_anchor_ids:
            continue
        context_start, context_end = start_ids, end_ids
        has_explicit_context_anchors = (
            context_anchor_ids is not None and context in context_anchor_ids
        )
        if has_explicit_context_anchors:
            raw_start, raw_end = context_anchor_ids[context]
            context_start = frozenset(raw_start or ())
            context_end = frozenset(raw_end or ())
        context_ids = {
            candidate.candidate_id
            for candidate in candidates
            if context is None or _path_membership_for_context(candidate, context) is not None
        }
        unknown = (context_start | context_end) - context_ids
        if unknown:
            if not has_explicit_context_anchors:
                continue
            raise ValueError(
                "unknown chain anchors for transcript context: "
                + ", ".join(sorted(unknown))
            )
        anchors_by_context[context] = (context_start, context_end)
    contexts = [context for context in contexts if context in anchors_by_context]

    exact_solutions = []
    for context in contexts:
        context_start, context_end = anchors_by_context[context]
        solution = _solve_context(
            candidates, context, context_start, context_end, 0.0,
        )
        if solution is not None:
            exact_solutions.append((context, solution))
    if not exact_solutions:
        # An anchor-bounded search can legitimately have no connecting path:
        # the intervening sequence may be absent, rearranged, or unresolved.
        # Preserve every reported candidate as descriptive membership evidence
        # and let the caller keep coordinates/absence in the unresolved state.
        # Raising here aborts the complete observation table and loses that
        # distinction.
        all_ids = frozenset(item.candidate_id for item in candidates)
        used_start_ids = frozenset(
            candidate_id
            for start, _end in anchors_by_context.values()
            for candidate_id in start
        )
        used_end_ids = frozenset(
            candidate_id
            for _start, end in anchors_by_context.values()
            for candidate_id in end
        )
        return CandidateChainResult(
            best_score=0.0,
            retained_ids=all_ids,
            best_path_member_ids=frozenset(),
            score_delta=float(score_delta),
            local_mode=not used_start_ids and not used_end_ids,
            retained_edges=frozenset(),
            start_anchor_ids=used_start_ids,
            end_anchor_ids=used_end_ids,
            ambiguity_status="no_connecting_chain",
            configuration_name=configuration_name,
            ambiguous_ids=all_ids,
        )
    best_score = max(solution["best_score"] for _context, solution in exact_solutions)
    solutions = []
    for context, exact_solution in exact_solutions:
        context_best = exact_solution["best_score"]
        remaining_delta = float(score_delta) - (best_score - context_best)
        if remaining_delta < 0:
            continue
        context_start, context_end = anchors_by_context[context]
        solutions.append(
            (
                context,
                _solve_context(
                    candidates, context, context_start, context_end, remaining_delta,
                    prepared=exact_solution["_prepared"],
                ),
            )
        )
    tolerance = 1e-12 * max(1.0, abs(best_score))
    retained = set()
    best_path_members = set()
    retained_edges = set()
    context_summaries = []
    ambiguous_ids = set()
    for context, solution in solutions:
        context_score = solution["best_score"]
        context_retained = solution["retained"]
        context_best = solution["best"]
        context_edges = solution["retained_edges"]
        if context_score + float(score_delta) + tolerance < best_score:
            continue
        retained.update(context_retained)
        retained_edges.update(context_edges)
        if abs(context_score - best_score) <= tolerance:
            best_path_members.update(context_best)
        context_start, context_end = anchors_by_context[context]
        ambiguous_ids.update(solution["ambiguous_ids"])
        context_summaries.append(ChainContextSummary(
            context=context,
            best_score=context_score,
            retained_ids=frozenset(context_retained),
            best_path_member_ids=frozenset(context_best),
            retained_edges=frozenset(context_edges),
            start_anchor_ids=frozenset(context_start),
            end_anchor_ids=frozenset(context_end),
            best_path_count_capped=solution["best_path_count_capped"],
            near_optimal_path_count_capped=solution["near_optimal_path_count_capped"],
            ambiguous_ids=frozenset(solution["ambiguous_ids"]),
        ))
    # A candidate must occur in every retained context, not merely be
    # obligatory in the context where it happens to occur.
    if context_summaries:
        mandatory = set.intersection(*(set(item.retained_ids - item.ambiguous_ids)
                                       for item in context_summaries))
        ambiguous_ids = retained - mandatory
    ambiguity = (
        "multiple_near_optimal_chains"
        if ambiguous_ids
        else "unique_within_reported_candidates"
    )
    used_start_ids = frozenset(
        candidate_id
        for summary in context_summaries
        for candidate_id in summary.start_anchor_ids
    )
    used_end_ids = frozenset(
        candidate_id
        for summary in context_summaries
        for candidate_id in summary.end_anchor_ids
    )
    return CandidateChainResult(
        best_score=best_score,
        retained_ids=frozenset(retained),
        best_path_member_ids=frozenset(best_path_members),
        score_delta=float(score_delta),
        local_mode=not used_start_ids and not used_end_ids,
        retained_edges=frozenset(retained_edges),
        start_anchor_ids=used_start_ids,
        end_anchor_ids=used_end_ids,
        ambiguity_status=ambiguity,
        configuration_name=configuration_name,
        best_path_count_capped=max(
            (summary.best_path_count_capped for summary in context_summaries),
            default=0,
        ),
        near_optimal_path_count_capped=max(
            (summary.near_optimal_path_count_capped for summary in context_summaries),
            default=0,
        ),
        ambiguous_ids=frozenset(ambiguous_ids),
        context_summaries=tuple(context_summaries),
    )


def genomic_candidate_chain(candidates, score_delta, start_ids=None, end_ids=None):
    """Physical DNA channel, independent of annotation/transcript membership.

    Requires one locus pair and one score scheme, with coordinates already
    oriented by the caller. It does not construct an RNA path or validate exons.
    Non-collinear candidates remain descriptive, never an absence call.
    """
    from dataclasses import replace
    physical = tuple(replace(candidate, path_memberships=()) for candidate in candidates)
    return ordered_candidate_chain(physical, score_delta, start_ids, end_ids,
                                   configuration_name="genomic_DNA_channel_v017")

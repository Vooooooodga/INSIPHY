"""Ordered candidate chains for intragenic sequence correspondence."""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from typing import FrozenSet, Iterable, Mapping, Optional, Tuple

from .coordinates import Interval0


@dataclass(frozen=True)
class ChainConfiguration:
    name: str = "insiphy_nt_chain_delta_v1"
    nucleotide_match_reward: float = 2.0
    relative_delta: float = 0.05

    def score_delta(self, best_score: float, score_scheme: str) -> float:
        if str(score_scheme).startswith("nt_"):
            return max(
                2.0 * self.nucleotide_match_reward,
                self.relative_delta * abs(float(best_score)),
            )
        return 0.0

    @property
    def delta_rule(self):
        return "max(2*match_reward,relative_delta*abs(best_score))"


DEFAULT_CHAIN_CONFIGURATION = ChainConfiguration()


@dataclass(frozen=True)
class ChainPathMembership:
    """One compatible query/target transcript-path pairing for a candidate."""

    query_path_id: str
    target_path_id: str
    query_order: int
    target_order: int
    query_contig: str
    target_contig: str
    query_strand: str
    target_strand: str

    def __post_init__(self):
        if not self.query_path_id or not self.target_path_id:
            raise ValueError("chain path membership requires both path IDs")
        if self.query_strand not in {"+", "-"} or self.target_strand not in {"+", "-"}:
            raise ValueError("chain path membership requires explicit gene strands")

    @property
    def context(self):
        return (
            self.query_path_id,
            self.target_path_id,
            self.query_contig,
            self.target_contig,
            self.query_strand,
            self.target_strand,
        )


@dataclass(frozen=True)
class ChainCandidate:
    candidate_id: str
    query: Interval0
    target: Interval0
    score: float
    score_scheme: str
    relative_strand: str = "+"
    path_memberships: Tuple[ChainPathMembership, ...] = ()

    def __post_init__(self):
        if not self.candidate_id:
            raise ValueError("candidate_id is required")
        if not isfinite(float(self.score)):
            raise ValueError(f"candidate score must be finite: {self.candidate_id}")
        if not self.score_scheme:
            raise ValueError(f"score_scheme is required: {self.candidate_id}")
        if self.relative_strand not in {"+", "-"}:
            raise ValueError(f"invalid relative strand: {self.relative_strand}")


@dataclass(frozen=True)
class ChainContextSummary:
    context: Optional[tuple]
    best_score: float
    retained_ids: FrozenSet[str]
    best_path_member_ids: FrozenSet[str]
    retained_edges: FrozenSet[Tuple[str, str]]
    start_anchor_ids: FrozenSet[str]
    end_anchor_ids: FrozenSet[str]
    best_path_count_capped: int
    near_optimal_path_count_capped: int
    ambiguous_ids: FrozenSet[str] = frozenset()


@dataclass(frozen=True)
class CandidateChainResult:
    best_score: float
    retained_ids: FrozenSet[str]
    best_path_member_ids: FrozenSet[str]
    score_delta: float
    local_mode: bool
    retained_edges: FrozenSet[Tuple[str, str]] = frozenset()
    start_anchor_ids: FrozenSet[str] = frozenset()
    end_anchor_ids: FrozenSet[str] = frozenset()
    ambiguity_status: str = "unique_within_reported_candidates"
    configuration_name: str = DEFAULT_CHAIN_CONFIGURATION.name
    best_path_count_capped: int = 0
    near_optimal_path_count_capped: int = 0
    ambiguous_ids: FrozenSet[str] = frozenset()
    context_summaries: Tuple[ChainContextSummary, ...] = ()


@dataclass(frozen=True)
class CoverageClassification:
    relation: str
    covered_bases: int
    overlap_bases: int
    uncovered_bases: Optional[int]


def _ordered_before(left, right, left_path=None, right_path=None):
    path_compatible = True
    if left_path is not None or right_path is not None:
        path_compatible = (
            left_path is not None
            and right_path is not None
            and left_path.context == right_path.context
            and left_path.query_order < right_path.query_order
            and left_path.target_order < right_path.target_order
        )
    return (
        path_compatible
        and left.relative_strand == right.relative_strand == "+"
        and left.query.end0 <= right.query.start0
        and left.target.end0 <= right.target.start0
    )


def _path_membership_for_context(candidate, context):
    if context is None:
        return None
    return next(
        (membership for membership in candidate.path_memberships if membership.context == context),
        None,
    )


def _top_complete_paths(candidates, candidate_indices, predecessors, terminal_indices, start_ids):
    """Return the two highest-scoring distinct complete paths in one DAG."""

    def preference(item):
        score, predecessor = item
        predecessor_key = predecessor if predecessor is not None else (-1, -1)
        return -score, predecessor_key

    def retain_two(items):
        retained = []
        seen = set()
        for item in items:
            if item[1] in seen:
                continue
            seen.add(item[1])
            retained.append(item)
            retained.sort(key=preference)
            if len(retained) > 2:
                retained.pop()
        return retained

    paths_by_index = {index: [] for index in candidate_indices}
    path_predecessor = {}
    for index in candidate_indices:
        candidate = candidates[index]
        options = []
        if not start_ids or candidate.candidate_id in start_ids:
            options.append((float(candidate.score), None))
        for predecessor in predecessors[index]:
            for score, path_id in paths_by_index[predecessor]:
                options.append((score + float(candidate.score), path_id))
        selected = retain_two(options)
        for rank, (score, predecessor_id) in enumerate(selected):
            path_id = (index, rank)
            path_predecessor[path_id] = predecessor_id
            paths_by_index[index].append((score, path_id))

    complete = []
    for index in terminal_indices:
        complete.extend(paths_by_index[index])
    selected_complete = retain_two(complete)

    def materialize(path_id):
        path = []
        while path_id is not None:
            candidate_index, _rank = path_id
            path.append(candidates[candidate_index].candidate_id)
            path_id = path_predecessor[path_id]
        return tuple(reversed(path))

    return tuple(
        (score, materialize(path_id))
        for score, path_id in selected_complete
    )


def _solve_context(candidates, context, start_ids, end_ids, score_delta, *, prepared=None):
    negative_infinity = float("-inf")
    if prepared is None:
        candidate_indices = [
            index
            for index, candidate in enumerate(candidates)
            if context is None or _path_membership_for_context(candidate, context) is not None
        ]
        if start_ids and not any(candidates[index].candidate_id in start_ids for index in candidate_indices):
            return None
        if end_ids and not any(candidates[index].candidate_id in end_ids for index in candidate_indices):
            return None

        predecessors = {index: [] for index in candidate_indices}
        successors = {index: [] for index in candidate_indices}
        for offset, left_index in enumerate(candidate_indices):
            left = candidates[left_index]
            left_path = _path_membership_for_context(left, context)
            for right_index in candidate_indices[offset + 1 :]:
                right = candidates[right_index]
                right_path = _path_membership_for_context(right, context)
                if _ordered_before(left, right, left_path, right_path):
                    predecessors[right_index].append(left_index)
                    successors[left_index].append(right_index)

        negative_infinity = float("-inf")
        forward = {index: negative_infinity for index in candidate_indices}
        for index in candidate_indices:
            candidate = candidates[index]
            starts_here = not start_ids or candidate.candidate_id in start_ids
            previous = max((forward[item] for item in predecessors[index]), default=negative_infinity)
            if starts_here:
                previous = max(0.0, previous)
            if previous != negative_infinity:
                forward[index] = float(candidate.score) + previous

        terminal_indices = [
            index for index in candidate_indices
            if not end_ids or candidates[index].candidate_id in end_ids
        ]
        reachable_terminal = [index for index in terminal_indices if forward[index] != negative_infinity]
        if not reachable_terminal:
            return None
        best_score = max(forward[index] for index in reachable_terminal)

        backward = {index: negative_infinity for index in candidate_indices}
        for index in reversed(candidate_indices):
            candidate = candidates[index]
            ends_here = not end_ids or candidate.candidate_id in end_ids
            following = max((backward[item] for item in successors[index]), default=negative_infinity)
            if ends_here:
                following = max(0.0, following)
            if following != negative_infinity:
                backward[index] = float(candidate.score) + following

    else:
        candidate_indices, predecessors, successors, forward, backward, reachable_terminal, best_score = prepared

    tolerance = 1e-12 * max(1.0, abs(best_score))
    threshold = best_score - float(score_delta) - tolerance
    through = {
        index: forward[index] + backward[index] - float(candidates[index].score)
        if forward[index] != negative_infinity and backward[index] != negative_infinity
        else negative_infinity
        for index in candidate_indices
    }
    retained = {
        candidates[index].candidate_id for index, value in through.items() if value >= threshold
    }
    best = {
        candidates[index].candidate_id
        for index, value in through.items()
        if abs(value - best_score) <= tolerance
    }
    retained_edges = {
        (candidates[left].candidate_id, candidates[right].candidate_id)
        for left in candidate_indices
        for right in successors[left]
        if forward[left] != negative_infinity
        and backward[right] != negative_infinity
        and forward[left] + backward[right] >= threshold
    }
    complete_paths = _top_complete_paths(
        candidates,
        candidate_indices,
        predecessors,
        reachable_terminal,
        start_ids,
    )
    best_path_count = sum(
        1 for score, _path in complete_paths if abs(score - best_score) <= tolerance
    )
    near_path_count = sum(1 for score, _path in complete_paths if score >= threshold)
    # A node is mandatory iff removing it eliminates *all* near-optimal paths.
    # Work on the existing DAG: do not introduce new edges across a removed node.
    shared_near_ids = set()
    for excluded in candidate_indices:
        candidate_id = candidates[excluded].candidate_id
        if candidate_id not in retained:
            continue
        without = {}
        for index in candidate_indices:
            if index == excluded:
                without[index] = negative_infinity
                continue
            previous = max((without.get(item, negative_infinity) for item in predecessors[index]),
                           default=negative_infinity)
            if not start_ids or candidates[index].candidate_id in start_ids:
                previous = max(0.0, previous)
            without[index] = previous + float(candidates[index].score)
        alternative = max((without.get(index, negative_infinity) for index in reachable_terminal),
                          default=negative_infinity)
        if alternative < threshold:
            shared_near_ids.add(candidate_id)
    return {
        "_prepared": (candidate_indices, predecessors, successors, forward, backward, reachable_terminal, best_score),
        "best_score": best_score,
        "retained": retained,
        "best": best,
        "retained_edges": retained_edges,
        "complete_paths": complete_paths,
        "best_path_count_capped": min(2, best_path_count),
        "near_optimal_path_count_capped": min(2, near_path_count),
        "ambiguous_ids": retained - shared_near_ids,
    }


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


def classify_reference_coverage(
    intervals: Iterable[Interval0],
    parent_interval: Optional[Interval0] = None,
):
    """Describe whether several projections are complementary or repetitive."""

    intervals = tuple(sorted(intervals))
    if not intervals:
        uncovered = parent_interval.length if parent_interval is not None else None
        return CoverageClassification("uncovered", 0, 0, uncovered)
    if parent_interval is not None and any(
        interval.start0 < parent_interval.start0 or interval.end0 > parent_interval.end0
        for interval in intervals
    ):
        raise ValueError("reference projection lies outside the parent interval")

    events = []
    for interval in intervals:
        events.append((interval.start0, 1))
        events.append((interval.end0, -1))
    covered_bases = 0
    overlap_bases = 0
    depth = 0
    previous = events[0][0]
    for coordinate, delta in sorted(events, key=lambda item: (item[0], item[1])):
        span = coordinate - previous
        if depth > 0:
            covered_bases += span
        if depth > 1:
            overlap_bases += span
        depth += delta
        previous = coordinate
    uncovered = None if parent_interval is None else parent_interval.length - covered_bases

    if len(intervals) == 1:
        if parent_interval is None:
            relation = "single_projection"
        else:
            relation = "single_complete" if uncovered == 0 else "single_partial"
    elif overlap_bases:
        relation = "repeated_overlap"
    else:
        if parent_interval is None:
            relation = "complementary_disjoint"
        else:
            relation = "complementary_complete" if uncovered == 0 else "complementary_partial"
    return CoverageClassification(relation, covered_bases, overlap_bases, uncovered)

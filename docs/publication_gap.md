# Publication Gap

## Current status

Status date: 2026-09-18.

Latest v0.14.0 status: run `20260918_121100_insiphy` completed all ten re-inference/visualization tasks and five thread-comparison tasks from unchanged cases/evidence. **136 selected formal regression tests passed** under Slurm job `61625` in 17.063 s ([report](/data/projects/intragenic_structure/results/20260918_121100_insiphy/regression_tests.txt)), including four new role-conflict tests. All six core tables agree in every case pair. The scheduled execution and targeted checks are complete; scientific limitations remain as recorded below.

Historical baseline `20260918_111929_insiphy` completed after recovery at 12:10:30 on 2026-09-18 without repeating alignments. Its 132-test result predates the dsx role-conflict correction; its three role contrasts and conditional ER posteriors are not the final result. The final run changes Bter EG_0029 to unknown role while retaining sequence presence, leaves two observation-level role contrasts, and reports the role ER fit as nonidentifiable with no node/branch posteriors.

Final focal observations retain the independently audited Hdac3 `D470_A471` contrast (1 present, 3 absent, 0 unknown) at the exact 62 bp intron, with Dana split and melanogaster-subgroup fusion equally parsimonious. rec8 `D1277_A1317` encloses the exact 39 bp interval (1 present, 1 absent, 2 unknown), with Cryo fusion and Octo split equally possible. rec8 is **partial recovery**: direct intron DNA correspondence and role contrast remain missing. spo5's focal 55 bp difference is **not recovered**. Linked final matrices, branch tables and statistical outputs are in the [completed assessment](real_data_benchmark.md#final-biological-assessment).

The earlier iteration passed 97 selected tests under job `61583` and completed ten real-data runs with exit status zero. Saved comparison reports showed identical observations in six core tables for each thread-count pair. Hdac3, spo5 and rec8 failed to recover the expected positive differences; statistical review identified invalid fitted posteriors without observed state contrast. Those historical findings remain in the [first-iteration benchmark](real_data_benchmark.md#v014-first-iteration-provisional-results).

Current test statuses refer to 136 selected formal regression tests; the earlier 132-test result is preserved as the intermediate baseline. Four new role-conflict tests passed. Legacy experimental and toy test classes were intentionally excluded. Coverage includes observation coding, exact ID mapping, transcript repertoire, predicted roles, tree validity, probability calculations, coding projections, phase handling, observed-species summaries and SVG semantics. Passing these tests establishes the checked behavior; empirical accuracy, small-sample statistical performance and publication readiness still require separate evidence.

The current documentation therefore supports a development snapshot and review plan. It does not claim publication readiness or completed biological benchmarking.

## Minimum formal claim after v0.14 acceptance

After tests and real runs pass, a first manuscript can reasonably claim:

- a formal representation of gene-internal synteny as homologous sequence units, exonic roles and splice junctions;
- sequence-assisted recovery of missing or incomplete annotation inside known orthologous genes;
- branch placement of qualitative structural changes under maximum parsimony;
- optional fixed-tree likelihood comparisons for gain/loss rate structure when data are estimable;
- visualization of confirmed exon correspondence with tree context.

This first claim should stay within single-copy genes.

## v0.14 unified repair ledger

Implementation and empirical scope are tracked separately in the status column.
Passing regression tests establishes the stated code behavior. Remaining
real-case limitations include incomplete sequence/role recovery, unknown
observations, ambiguous event direction and unavailable rate-model comparisons.
The current scientific outcomes do not mark every implemented repair as a
fully accepted biological capability.

| Theme | Source | v0.14 measure | Status |
|---|---|---|---|
| Annotation search bounds too narrow | previous 9 issues | Preserve original bounds, linked features and expanded search window | Implemented; recovery beyond truncated gene bounds remains quantitatively unassessed |
| DNA presence confused with exon role | previous 9 issues; external report | Separate sequence presence, confirmed role, predicted role and non-exonic source | Regression passed; final dsx EG_0029 retains DNA presence while its conflicting role is unknown; rec8 role recovery remains partial |
| Unsupported absence calls | previous 9 issues; external report | Require flank/context evidence and keep unresolved cases unknown | Targeted deletion/flank tests passed; empirical sequence-deletion error rate remains unmeasured |
| Split/fusion limited to simple cases | previous 9 issues | Allow ordered `1:n` and `n:1` complementary projections | Regression passed; exact Hdac3 contrast and partial rec8 boundary recovered, both direction-ambiguous; spo5 unrecovered |
| Splice-boundary heuristic coordinates | previous 9 issues; external report | Use projected donor/acceptor coordinates and exact boundary evidence | Regression passed; exact Hdac3 and rec8 focal boundaries independently audited; broader coverage remains unmeasured |
| Strand/orientation loss | previous 9 issues | Preserve relative orientation and reject incompatible splice projection | Negative-strand tests passed; rec8 focal negative-strand mapping audited; broader orientation cases remain unmeasured |
| Alignment metric inconsistency | previous 9 issues | Use common alignment-statistics adapter and explicit backend reporting | Tests passed; independent Hdac3 review confirmed preserved DNA metrics and separately labelled amino-acid identity, CDS coverage and selected coding blocks; broader accuracy remains unmeasured |
| Ascertainment mismatch | previous 9 issues; external report | Align site inclusion rules with conditional likelihood mode | Three-mode and enumerated-pattern tests passed; statistical performance review pending |
| Understated uncertainty | previous 9 issues; external report | Mark invalid fits and unestimated sensitivity directly; do not invent ranges | Regression passed; final dsx role ER is nonidentifiable and emits no posteriors; uncertainty coverage remains unestablished |
| Known non-exonic/absent records dropped | old review 12; external report | Keep non-exonic source and explicit absence as separate evidence classes | Retention tests passed; explicit non-exonic observations remain in two final dsx role contrasts; sequence-absence accuracy remains unmeasured |
| Sequence hit upgraded to hard exon | old review 12; external report | `predicted_exon_candidate` remains predicted evidence | Tests passed; final EG_0029 conflict becomes unknown role, with DNA presence retained; no hard predicted exon added |
| Longest-CDS representative bug | old review 12 | Correct CDS-length interpretation in representative selection | CDS-metadata union regression passed; real transcript review pending |
| Splice boundary merge too broad | old review 12 | Remove length-ratio fallback for junction homology | Implemented; Hdac3 and rec8 focal coordinate mappings audited; broader boundary-shift recovery remains unmeasured |
| Disconnected or malformed tree accepted | old review 12 | Validate rooted, connected, acyclic tree with unique labels | Invalid-topology regression tests passed |
| Invalid model used for posteriors | old review 12; external report | Posterior tables and figures require valid selected fit and finite probabilities | Regression passed; published repaired RpL32 outputs suppress no-contrast fitted posteriors; broad statistical validation remains incomplete |
| LRT eligibility confused with posterior validity | external report; parent visual review | Gate probability display on selected-fit validity, not ER/ARD P-value availability | Regression passed; historical baseline ER conditional posteriors were distinct from unavailable LRTs; final role ER is invalid and emits no posteriors |
| Site universe and ascertainment inconsistent | old review 12 | Require explicit complete-universe catalogue for complete-universe mode | Catalogue and mode-selection tests passed; real input review pending |
| Site dependence overstated as event count | old review 12 | Document site-level histories and avoid independent molecular-event counts | Documentation updated |
| Internal repeats confused with split | old review 12 | Distinguish complementary projections from repeated full-overlap hits | Ordered-fragment, repeat and alternative-exon regressions passed; real positive cases pending |
| Multi-transcript conflict treated as unknown only | old review 12; user review | Preserve all transcript paths as repertoire evidence | Repertoire, CDS-phase and separate-lane tests passed; real transcript review pending |
| Extant summaries named as ancestral | external report | Rename to `observed_element_tree_coverage.tsv` and `observed_intragenic_paths.tsv` | Output-name regression passed; real cases produce the renamed extant summaries; no joint ancestral structure claim |
| Same-species fragments inflate conservation calls | correspondence review, 2026-09-18 | Set `segment_conservation.conservation_call` to `observed_single_species`, `observed_in_multiple_species` or `no_species_presence_observed` using distinct species with observed presence; exclude absence records and count same-species fragments once; retain mean DNA identity as descriptive | Observed-species summary regression passed; real acceptance pending |
| Candidate/predicted visualized as confirmed homology | external report; parent visual review | Confirmed ribbons exclude predicted, unknown and candidate-source boxes | SVG semantic tests passed; rendered real-figure inspection pending |
| Same-species transcript lanes linked as homology | parent visual review | Use facing lanes for shared occurrences in cross-species ribbons and retain distinct split members | Shared-isoform ribbon tests passed; rendered real-figure inspection pending |
| Missing branch lengths interrupt drawing | parent dependency review | Use all-unit figure depths when any non-root length is missing; preserve all supplied lengths including zero when complete | Missing/root-missing/zero-length drawing tests passed; CTMC requirements remain strict |
| Raw feature diversity underdocumented | user review | Preserve original features and ownership in `raw_gene_features.tsv` | Retention tests passed and first real runs produced the table; additional statistical layers remain future work |
| Whole-exon DNA threshold excludes a CDS/UTR segment | real-data review, 2026-09-18 | Map CDS-derived proteins in full-transcript MAFFT alignments to occurrence-local coding intervals, retaining original DNA evidence | Regression passed; audited Hdac3 focal junction recovered; rec8 remains partial and spo5 unrecovered |
| Minus-strand or missing phase affects coding correspondence | real-data repair review, 2026-09-18 | Use transcript-oriented coding coordinates and mark missing phase unavailable | Negative-strand and unknown-phase tests passed; rec8 focal mapping audited; broader incomplete-CDS coverage remains unmeasured |
| No observed contrast still produces fitted posteriors | statistical review after first real run | Classify no-contrast layers as unavailable for fitted branch histories | Tests passed; RpL32 and dsx no-contrast layers report no fitted posteriors; final dsx role layer is also withheld for nonidentifiability |
| Whole-box ribbons overstate partial coding correspondence | real-figure review, 2026-09-18 | Clip ribbons to direct accepted blocks, dispatch by correspondence basis, retain full annotation boxes | SVG tests passed; written Hdac3 and rec8 example endpoints match selected protein blocks; rendered review pending |
| Predicted CDS overlaps a hard non-exonic observation | dsx biological-source audit, EG_0029 in Bter | Resolve supported overlapping prediction/non-exonic evidence as role uncertainty while preserving confirmed exonic usage and unresolved sequence absence; unrelated predictions do not alter roles | Fix implemented, four new regressions passed, final Bter EG_0029 role is unknown while presence remains present; all five core-table comparisons agree across threads |
| Legacy calibration evidence attributed to the formal model | external report | CLI labels the retained simulator/calibration commands as experimental multi-copy; formal single-copy evidence comes from its own probability tests and real-case runs | Legacy code retained; no simulation study performed and no claim that those legacy results validate single-copy P values |

Targeted test definitions are in [observations](../tests/test_observations.py), [interfaces](../tests/test_interfaces.py), [sequence evidence](../tests/test_evidence.py), [coding correspondence](../tests/test_coding_correspondence.py), [statistics](../tests/test_statistics.py), and [visualization](../tests/test_visualization.py). The current total of 136 counts selected formal regression tests; real biological cases are counted separately.

## Real-data work still required

The final five-case re-inference and all core-table thread comparisons are complete. Historical full-analysis timings remain separate: the [baseline record](real_data_benchmark.md#completed-baseline-comparisons-and-runtime) reports dsx at 1,306.460 s and 666.822 s, with 2.14 GiB peak RSS. These are alignment-inclusive baseline measurements, not final inference-only task times.

Scientific limits remain after the completed fixes: Hdac3 direction is ambiguous; rec8 lacks direct intron DNA/role recovery; spo5 is unrecovered; RpL32 and dsx have incomplete observations. The two remaining dsx role contrasts each have five unknown species and support alternative parsimony placements, with no valid fitted branch posterior in the final run. These limits are distinct from the implemented, tested and observed EG_0029 conflict correction.

1. **Curated positive single-copy cases**
   Select real genes with literature-supported intron loss, exonization, splice-boundary shift, split/fusion or role transition. Reconfirm modern assembly coordinates and single-copy status.

2. **Conserved negative controls**
   Analyze multiple stable single-copy genes across the same taxa. Report unsupported changes, unknown fraction and correspondence ambiguity.

3. **Divergence gradient**
   Include close and moderately diverged species sets to measure where sequence-assisted internal correspondence starts to fail.

4. **Comparator runs**
   Compare annotation-only, sequence-only and full structure-aware correspondence. These are method-ablation comparisons, not pathway analyses.

5. **Manual biological audit**
   Inspect alignments, splice boundaries, and generated SVG figures for every showcase case.

6. **Runtime profile**
   Report elapsed time and memory as species number, exon count and candidate pair count grow.

## Statistical limitations to report

- Parsimony conditions on equal costs and the supplied root.
- ER/ARD P values are asymptotic model-comparison P values.
- A valid branch-history posterior can exist without a valid ER/ARD LRT.
- A small gene may contain too few informative sites for rate asymmetry.
- Neighboring exons and junctions can be coupled by one biological change.
- Tree topology, branch lengths, orthology and correspondence are conditioned upon.
- Current posterior intervals do not propagate full joint parameter uncertainty.

## Multi-copy work deferred

Duplicated genes need gene-tree reconciliation, copy-aware orthology/paralogy, and internal-structure histories that can differ among copies. Existing multi-copy routines remain experimental. They should not be used for formal v0.14 performance claims.

## First-paper readiness gate

The current 136 selected formal regression tests passed under job `61625`. Before manuscript submission, the project still needs:

- broader acceptance evidence beyond the completed five-case assessment;
- resolution or measured characterization of the remaining rec8, spo5 and missing-observation limitations;
- additional curated positives and conserved negative controls;
- visual review of synteny and integrated tree figures;
- quantitative comparison against appropriate alternatives and coverage/error assessment;
- a clear statement that mechanism interpretation requires external biological evidence.

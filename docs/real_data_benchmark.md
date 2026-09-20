# Real-Data Demonstration

## v0.15.0 formal assessment: pending

The v0.15.0 five-case assessment has not been run. No v0.15 recovery count,
site count, likelihood result, thread comparison, or figure assessment is
available yet. The planned roles are:

| Case | v0.15 assessment role | Upstream qualification |
|---|---|---|
| RpL32 | Qualitative conservation control | Supplied single-copy grouping must be recorded with the run |
| spo5 / `spbc29a10_02` | Positive observable benchmark | Curated grouping only; no independent genome-wide single-copy qualification yet |
| rec8 / `spog_00055` | Positive observable benchmark | Curated grouping only; no independent genome-wide single-copy qualification yet |
| Hdac3 | Positive observable benchmark | Curated grouping only; no independent genome-wide single-copy qualification yet |
| dsx | Descriptive case for annotation coverage and unresolved structure | No positive-event or mechanism claim |

The formal run must first create and freeze one schema-v3 repertoire matrix per
case. Parsimony, ER/ARD, and foreground analyses must consume that same file.
Canonical analyses, when run, use a separate view-specific sensitivity matrix.
After completion, this section will record the run identifier, exact input and
execution provenance, upstream single-copy qualification, matrix source and
view, total and CTMC-included site counts, focal-site recovery, parsimony
alternatives, model-test availability, thread comparison, and rendered-figure
review. Until those fields are filled, the v0.15 status remains `pending`.

All numerical results below belong to pre-v0.15 software and are historical
records. They do not describe v0.15 behavior.

## Historical v0.14.0 completed assessment

Run `20260918_121100_insiphy` started at 12:15 on 2026-09-18. Its **136 selected
formal regression tests passed** under Slurm job `61625` in **17.063 s**
([test report](/data/projects/intragenic_structure/results/20260918_121100_insiphy/regression_tests.txt)).
This included four new role-conflict tests. All ten re-inference/visualization
tasks and five thread-comparison tasks completed. The run covers parsimony,
ER/ARD, foreground analysis and SVG output for five biological cases at 1 and
16 threads, reusing existing prepared cases and sequence evidence without
repeating alignments. This completed the scheduled v0.14 software assessment;
the scientific limitations below remain.

### v0.14 outputs and thread agreement

Each comparison report records `same_observations=true` for all six core
tables. Counts exclude headers; branch rows represent endpoint assignments,
including unchanged or ambiguous histories, and are not independent events.

| Case and comparison report | Occurrence (= homology) rows | Correspondence rows | Completion rows | Species-state rows | Branch rows |
|---|---:|---:|---:|---:|---:|
| [RpL32](/data/projects/intragenic_structure/results/20260918_121100_insiphy/rpl32_control/thread_comparison.tsv) | 41 | 24 | 10 | 100 | 160 |
| [dsx](/data/projects/intragenic_structure/results/20260918_121100_insiphy/dsx_bees/thread_comparison.tsv) | 163 | 97 | 109 | 546 | 961 |
| [Hdac3](/data/projects/intragenic_structure/results/20260918_121100_insiphy/hdac3/thread_comparison.tsv) | 22 | 13 | 0 | 36 | 56 |
| [spo5](/data/projects/intragenic_structure/results/20260918_121100_insiphy/spbc29a10_02/thread_comparison.tsv) | 8 | 6 | 14 | 40 | 60 |
| [rec8](/data/projects/intragenic_structure/results/20260918_121100_insiphy/spog_00055/thread_comparison.tsv) | 34 | 19 | 22 | 96 | 150 |

The six tables are `case/segment_occurrences.tsv`, `case/segment_homology.tsv`,
`results_parsimony/element_correspondence.tsv`,
`results_parsimony/annotation_completion_candidates.tsv`,
`results_parsimony/structural_site_matrix.tsv`, and
`results_parsimony/branch_structural_events.tsv`. These comparisons establish
agreement for the named tables only.

### v0.14 biological assessment

- **RpL32:** 20 sites (5 presence, 5 role, 10 junction), with no observed
  contrast and zero minimum changes. Unknown observations remain; this
  supports conservation of the observed structure with incomplete coverage.
  [Site summary](/data/projects/intragenic_structure/results/20260918_121100_insiphy/rpl32_control/threads_1/analysis/results_parsimony/structural_site_summary.tsv).
- **Hdac3:** `D470_A471` recovers the audited 62 bp focal intron contrast,
  1 present / 3 absent / 0 unknown. Dana split and melanogaster-subgroup fusion
  remain equally parsimonious, each `possible`.
  [Site summary](/data/projects/intragenic_structure/results/20260918_121100_insiphy/hdac3/threads_1/analysis/results_parsimony/structural_site_summary.tsv),
  [branch table](/data/projects/intragenic_structure/results/20260918_121100_insiphy/hdac3/threads_1/analysis/results_parsimony/branch_structural_events.tsv).
- **rec8: partial recovery.** `D1277_A1317` encloses the exact 39 bp focal
  interval, with 1 present / 1 absent / 2 unknown. Cryo fusion and Octo split
  remain equally possible. Direct intron DNA correspondence and the associated
  exon/intron role contrast remain unrecovered.
  [Site summary](/data/projects/intragenic_structure/results/20260918_121100_insiphy/spog_00055/threads_16/analysis/results_parsimony/structural_site_summary.tsv),
  [branch table](/data/projects/intragenic_structure/results/20260918_121100_insiphy/spog_00055/threads_16/analysis/results_parsimony/branch_structural_events.tsv).
- **spo5: unrecovered.** The focal 55 bp difference still has no corresponding
  junction in the output; 5 presence and 5 role sites remain without an
  observed contrast. The reused protein-assisted evidence fails its existing
  similarity/coverage criteria.
  [Site summary](/data/projects/intragenic_structure/results/20260918_121100_insiphy/spbc29a10_02/threads_1/analysis/results_parsimony/structural_site_summary.tsv).
- **dsx:** the conflicting Bter role at `EG_0029` is now `unknown`, while
  sequence presence is retained. Two observation-level role contrasts remain:
  `EG_0013` (Acer/Amel) and `EG_0055` (Bpas/Bter), each with five unknown species.
  They do not uniquely determine a historical event or its mechanism.
  [Site summary](/data/projects/intragenic_structure/results/20260918_121100_insiphy/dsx_bees/threads_16/analysis/results_parsimony/structural_site_summary.tsv),
  [observation matrix](/data/projects/intragenic_structure/results/20260918_121100_insiphy/dsx_bees/threads_16/analysis/results_parsimony/structural_site_matrix.tsv).

Final dsx has 24 presence sites (74 known / 94 unknown observations), 24 role
sites (63 known: 61 exonic and 2 non-exonic / 105 unknown), and 30 junction
sites (65 known / 145 unknown). EG_0029 retains Bign exonic usage and has zero
minimum changes. EG_0013 and EG_0055 each require one minimum change overall,
with respectively two and three possible placements and no unique branch.

### v0.14 dsx likelihood status

The final [ER/ARD tests](/data/projects/intragenic_structure/results/20260918_121100_insiphy/dsx_bees/threads_16/analysis/results_er-ard/model_tests.tsv)
and [foreground tests](/data/projects/intragenic_structure/results/20260918_121100_insiphy/dsx_bees/threads_16/analysis/results_foreground/model_tests.tsv)
report all six P/Q values as `NA`, with `parameters_not_estimable`. Eligible
role-fit sites decrease from the baseline's 18 to 17, and informative role
patterns from 3 to 2. ER/ARD and foreground role LRT statistics are 3.06583 and
0.224956, respectively; they are diagnostic values without valid P values. Presence
and junction layers have no observed contrast. The role-layer ER fit is now
`nonidentifiable`; role ARD and foreground fits are `boundary_limited`.
[ER/ARD fit diagnostics](/data/projects/intragenic_structure/results/20260918_121100_insiphy/dsx_bees/threads_16/analysis/results_er-ard/model_fits.tsv).
The final ER/ARD and foreground node/branch posterior tables contain headers
only. The intermediate baseline's successful role ER fit and its 468/432
posterior rows are historical and must not be used to interpret this final run.
Default parsimony retains the qualitative alternatives. No significant rate
comparison or validated branch-probability claim follows from this assessment.

### v0.14 figure assessment scope

SVG semantic regression tests passed. Read-only checks of the final integrated
SVGs confirmed `protein_projected_blocks` endpoints `1-597:471-1067` for
[Hdac3](/data/projects/intragenic_structure/results/20260918_121100_insiphy/hdac3/threads_1/analysis/figures/integrated_phylo_synteny.svg)
and `1-77:1317-1393` for
[rec8](/data/projects/intragenic_structure/results/20260918_121100_insiphy/spog_00055/threads_16/analysis/figures/integrated_phylo_synteny.svg),
matching the retained coding-projection evidence. Full annotated exon boxes
remain intact. No rendered SVG inspection was performed; visual readability
and the entire figure set still need review.

## Historical v0.14 intermediate baseline

<details>
<summary>132-test baseline and observations before the final role-conflict repair</summary>

The following records describe `20260918_111929_insiphy`. Pending statements
and the three dsx role contrasts in this historical snapshot were superseded
by the completed assessment above.

Status date: **2026-09-18**. The repaired iteration passed **132 selected formal
regression tests** under Slurm job `61609`: 17.717 s in the
[test report](/data/projects/intragenic_structure/results/20260918_111929_insiphy/regression_tests.txt).
The task trace reports 19.7 s elapsed time and 129.5 MB peak RSS for this
regression task. These resource figures describe the test task only.

Run `20260918_111929_insiphy` completed after recovery at **12:10:30 on
2026-09-18**, retaining all five cases at 1 and 16 threads as an intermediate
alignment/evidence baseline. The coordinator's earlier SIGTERM (exit 143) was
recovered without repeating alignments. All five saved reports show
`same_observations=true` for the six core tables. Execution completion and these
agreements are separate from biological acceptance.

A dsx audit identified conflicting predicted-CDS and hard non-exonic evidence.
Run `20260918_121100_insiphy` reuses the exact prepared cases and evidence for
the scoped role-conflict repair; its regression and execution status are
recorded above. The baseline's 132 passing tests precede that fix. Final method
outcomes remain on hold until revised inference is complete and reviewed.
Independent focal-site review confirms the baseline recovery grades and branch
ambiguity below. The first-iteration findings remain separate historical records.

### Completed baseline comparisons and runtime

The six compared tables are `case/segment_occurrences.tsv`,
`case/segment_homology.tsv`, `results_parsimony/element_correspondence.tsv`,
`results_parsimony/annotation_completion_candidates.tsv`,
`results_parsimony/structural_site_matrix.tsv`, and
`results_parsimony/branch_structural_events.tsv`. All six agree in each saved
report: [RpL32](/data/projects/intragenic_structure/results/20260918_111929_insiphy/rpl32_control/thread_comparison.tsv),
[dsx](/data/projects/intragenic_structure/results/20260918_111929_insiphy/dsx_bees/thread_comparison.tsv),
[Hdac3](/data/projects/intragenic_structure/results/20260918_111929_insiphy/hdac3/thread_comparison.tsv),
[spo5](/data/projects/intragenic_structure/results/20260918_111929_insiphy/spbc29a10_02/thread_comparison.tsv),
and [rec8](/data/projects/intragenic_structure/results/20260918_111929_insiphy/spog_00055/thread_comparison.tsv).
This comparison does not cover every likelihood output or rendered figure.

Baseline dsx elapsed times were **1,306.460 s (21 min 46.460 s)** with 1 thread
and **666.822 s (11 min 6.822 s)** with 16 threads; peak RSS was **2.14 GiB** in
both tasks. These are the completed baseline's measurements, separate from
the older first iteration and forthcoming inference-only run.

### RpL32 from the repaired run

Read-only counts from both published thread settings under
`20260918_111929_insiphy/rpl32_control` give the following data rows, excluding
headers. The saved baseline comparison report also confirms matching contents
for these six tables.

| Table | 1 thread | 16 threads |
|---|---:|---:|
| `segment_occurrences.tsv` | 41 | 41 |
| `segment_homology.tsv` | 41 | 41 |
| `element_correspondence.tsv` | 24 | 24 |
| `annotation_completion_candidates.tsv` | 10 | 10 |
| `structural_site_matrix.tsv` | 100 | 100 |
| `branch_structural_events.tsv` | 160 | 160 |

The [1-thread site summary](/data/projects/intragenic_structure/results/20260918_111929_insiphy/rpl32_control/threads_1/analysis/results_parsimony/structural_site_summary.tsv)
contains 5 sequence-presence sites, 5 role sites and 10 junction sites across
five species. All 20 sites have `no_observed_contrast` and `min_changes=0`.
Unknown observations remain in each layer; this supports conserved observed
structure with incomplete coverage. The 160 branch-table rows describe
endpoint assignments. The [completion summary](/data/projects/intragenic_structure/results/20260918_111929_insiphy/rpl32_control/threads_1/analysis/results_parsimony/annotation_completion_summary.tsv)
contains 9 `ambiguous_evidence` records and 1 `predicted_exon_candidate`.

All three [ER/ARD tests](/data/projects/intragenic_structure/results/20260918_111929_insiphy/rpl32_control/threads_1/analysis/results_er-ard/model_tests.tsv)
and three [foreground tests](/data/projects/intragenic_structure/results/20260918_111929_insiphy/rpl32_control/threads_1/analysis/results_foreground/model_tests.tsv)
have `P=NA / parameters_not_estimable`. The ER/ARD
[structural-change table](/data/projects/intragenic_structure/results/20260918_111929_insiphy/rpl32_control/threads_1/analysis/results_er-ard/structural_changes.tsv)
records `posterior_not_reported`, `fit_status=not_estimable` and
`inference_status=no_observed_contrast`. This is real-case evidence for the
specific no-contrast output repair.

### Bee dsx: 16-thread statistical audit

The completed statistical audit reports the following baseline matrix counts
for seven species. Biological-source review subsequently found a conflicting
evidence defect at `EG_0029`: Bter's 44 bp predicted CDS candidate overlaps a
hard `not_exonic` observation from an annotated intron. Explicit annotated
introns support the recorded non-exonic sources, but this overlap still requires
matrix-level conflict handling. Counts and statistical outputs below describe
the baseline and await revised inference; they do not establish three confirmed
role changes or a final post-repair count.

| Layer | Sites | Known observations | Unknown observations | Sites with observed contrast |
|---|---:|---:|---:|---:|
| Sequence presence | 24 | 74 / 168 | 94 | 0 |
| Exonic role | 24 | 64 / 168 | 104 | 3 |
| Splice junction | 30 | 65 / 210 | 145 | 0 |

Role observations comprise 61 exonic and 3 non-exonic states. The three role
contrasts are `EG_0013` (Acer/Amel), `EG_0029` (Bign/Bter), and `EG_0055`
(Bpas/Bter). Each has five unknown species and respectively 2, 2 and 3 possible
parsimony branch placements in the baseline. The EG_0029 evidence conflict is
a correctness defect repaired after this baseline; the new targeted tests
passed, while revised real-data outputs remain pending.

All six model-comparison P/Q values are `NA`. Sequence-presence and junction
layers have `no_observed_contrast` and emit no fitted posteriors. The role-layer
ER fit is numerically successful with an interior estimate and **does emit
conditional-MLE posteriors**: 468 node-posterior rows and 432 branch-posterior
rows. Its profile optimization failed and sensitivity limits are `NA`; these
conditional probabilities do not include an estimated parameter-uncertainty
range. ARD and foreground fits are boundary solutions, and no foreground
posteriors are emitted. Unavailable LRTs and availability of a selected ER
conditional posterior are reported separately.

### Provisional recovered observations

Counts refer to species with a splice junction at the indicated homologous
position. The `D..._A...` coordinates are local to the reference occurrence;
they are not assembly coordinates. Independent review matched the recovered
boundaries to the modern focal intervals in the [external site audit](real_positive_cases.md#hdac3-fig-s4).
Both recovered contrasts retain two equally parsimonious directional histories.

| Case and inspected run | Reference-local boundary | Present | Absent | Unknown | Audited outcome |
|---|---|---:|---:|---:|---|
| Hdac3, 1 thread | `D470_A471` | 1 | 3 | 0 | Exact 62 bp focal intron contrast recovered in all four species; direction ambiguous |
| rec8, 16 threads | `D1277_A1317` | 1 | 1 | 2 | Partial recovery: correct boundary around the 39 bp focal interval; direct DNA correspondence and role contrast missing |
| spo5, 1 thread | No corresponding junction recovered | NA | NA | NA | Focal 55 bp difference not recovered |

The [Hdac3 observation matrix](/data/projects/intragenic_structure/results/20260918_111929_insiphy/hdac3/threads_1/analysis/results_parsimony/structural_site_matrix.tsv)
records the intron in D. ananassae and continuous sequence in the three other
species. The focal 62 bp intron is `NC_057927.1:3929894-3929955:+`.
Its [branch table](/data/projects/intragenic_structure/results/20260918_111929_insiphy/hdac3/threads_1/analysis/results_parsimony/branch_structural_events.tsv)
reports two alternative one-change histories: `root->Drosophila_ananassae`
split, or `root->melanogaster_subgroup` fusion. Both are `possible`;
neither direction is required across all minimum-change histories.

The [rec8 matrix](/data/projects/intragenic_structure/results/20260918_111929_insiphy/spog_00055/threads_16/analysis/results_parsimony/structural_site_matrix.tsv)
records S. octosporus as present, S. cryophilus as absent, and S. pombe and
S. japonicus as `not_covered`. The boundary encloses reference-exon positions
1278-1316, exactly the 39 bp coding interval
`NW_013185626.1:118741-118779:-` confirmed in the [external audit](real_positive_cases.md#rec8-fig-3b).
Direct homology between this DNA and the corresponding introns, together with
their exon/intron role contrast, remains unrecovered. This case is therefore
**partial recovery**. Its [branch table](/data/projects/intragenic_structure/results/20260918_111929_insiphy/spog_00055/threads_16/analysis/results_parsimony/branch_structural_events.tsv)
retains two alternative one-change histories: `cry_octo->Schizosaccharomyces_cryophilus`
fusion or `cry_octo->Schizosaccharomyces_octosporus` split, each `possible`.

The [spo5 summary](/data/projects/intragenic_structure/results/20260918_111929_insiphy/spbc29a10_02/threads_1/analysis/results_parsimony/structural_site_summary.tsv)
still has five sequence-presence sites, five role sites and no junction site.
Its [match table](/data/projects/intragenic_structure/results/20260918_111929_insiphy/spbc29a10_02/threads_1/analysis/case/segment_matches.tsv)
labels all attempted protein-assisted pairs `low_protein_similarity_or_coverage`.
The S. octosporus/S. cryophilus exon pair remains accepted on DNA evidence.
The focal 55 bp S. pombe intron, `NC_003423.3:2536843-2536897:+`, and its
cross-species role difference remain unrecovered under the current acceptance
criteria. [External focal-site evidence](real_positive_cases.md#spo5-fig-3a)
is retained separately. No thresholds were changed for this result assessment.

### Likelihood status

Reviewed ER/ARD comparisons in [Hdac3](/data/projects/intragenic_structure/results/20260918_111929_insiphy/hdac3/threads_1/analysis/results_er-ard/model_tests.tsv),
[rec8](/data/projects/intragenic_structure/results/20260918_111929_insiphy/spog_00055/threads_16/analysis/results_er-ard/model_tests.tsv)
and [spo5](/data/projects/intragenic_structure/results/20260918_111929_insiphy/spbc29a10_02/threads_1/analysis/results_er-ard/model_tests.tsv)
all report `p_value=NA` and `test_status=parameters_not_estimable`. Hdac3 and
rec8 each have one informative junction pattern. Their qualitative branch
alternatives remain available under parsimony; no likelihood significance
or uniquely supported direction is established. The 132 passing regression
tests do not establish statistical performance across real biological cases.

### SVG provenance inspection

Read-only inspection of the written SVG and match tables found matching
protein-projection endpoints in two concrete examples:

- Hdac3 `match_00083`: the
  [integrated SVG](/data/projects/intragenic_structure/results/20260918_111929_insiphy/hdac3/threads_1/analysis/figures/integrated_phylo_synteny.svg)
  uses `1-597:471-1067` from `protein_projected_blocks`. The full source and
  target annotation boxes represent 859 and 1,148 bases respectively. The
  separate DNA alignment reaches source position 859 and target position 1,148;
  those longer spans are not used for this protein-supported ribbon.
- rec8 `match_00298`: the
  [integrated SVG](/data/projects/intragenic_structure/results/20260918_111929_insiphy/spog_00055/threads_16/analysis/figures/integrated_phylo_synteny.svg)
  uses `1-77:1317-1393`, matching the accepted coding projection. Its separate
  DNA blocks extend to source position 254 and target position 1,525.

Both paths record `correspondence_basis=annotated_CDS_protein` and
`projection_field=protein_projected_blocks`. This inspection confirms these
stored examples only. No renderer was run; layout readability and the complete
figure set still require visual review. Overall positive-case acceptance and
statistical validation remain incomplete.

</details>

## v0.14 first iteration: provisional results

Assessment date: **2026-09-18**. The first complete v0.14 run finished with
exit status zero for all ten case tasks: five biological cases, each analyzed
with 1 and 16 threads. The 97 selected formal regression tests passed under
Slurm job `61583`; legacy experimental and toy test classes were excluded.
**Biological acceptance remains pending.** Hdac3, spo5 and rec8 all failed to
recover the expected positive structural differences. A further statistical
defect was identified in fitted models with no observed state contrast.

The completed first-iteration records are:

- [Selected regression test output](/data/projects/intragenic_structure/results/20260916_214900_insiphy/regression_tests.txt).
- [Execution trace and resource measurements](/data/projects/intragenic_structure/logs/20260916_214900_insiphy/nextflow/trace.txt).
- [Run summary](/data/projects/intragenic_structure/logs/20260916_214900_insiphy/run_summary.md), with the completed invocation beginning at 10:59:55 and ending at 11:18:33 on 2026-09-18.
- Result root: `/data/projects/intragenic_structure/results/20260916_214900_insiphy`.

The date in this reused run-directory name predates its execution. Each case
has `threads_1/analysis` and `threads_16/analysis` outputs. These records describe
the first iteration; acceptance of the next repair requires new results.

### Thread-count comparison

All five saved `thread_comparison.tsv` reports record `same_observations=true`
for each of six tables: `case/segment_occurrences.tsv`,
`case/segment_homology.tsv`, `results_parsimony/element_correspondence.tsv`,
`results_parsimony/annotation_completion_candidates.tsv`,
`results_parsimony/structural_site_matrix.tsv`, and
`results_parsimony/branch_structural_events.tsv`.

| Case and comparison report | Occurrence rows (= homology rows) | Correspondence rows | Completion candidate rows | Site observation rows | Branch endpoint rows |
|---|---:|---:|---:|---:|---:|
| [RpL32](/data/projects/intragenic_structure/results/20260916_214900_insiphy/rpl32_control/thread_comparison.tsv) | 41 | 24 | 9 | 95 | 152 |
| [Bee dsx](/data/projects/intragenic_structure/results/20260916_214900_insiphy/dsx_bees/thread_comparison.tsv) | 163 | 97 | 116 | 574 | 1021 |
| [Hdac3](/data/projects/intragenic_structure/results/20260916_214900_insiphy/hdac3/thread_comparison.tsv) | 22 | 13 | 3 | 40 | 60 |
| [spo5 (`spbc29a10_02`)](/data/projects/intragenic_structure/results/20260916_214900_insiphy/spbc29a10_02/thread_comparison.tsv) | 8 | 6 | 14 | 40 | 60 |
| [rec8 (`spog_00055`)](/data/projects/intragenic_structure/results/20260916_214900_insiphy/spog_00055/thread_comparison.tsv) | 34 | 19 | 33 | 120 | 180 |

Counts above are data rows excluding headers. A site-observation row records
one site in one species; branch endpoint rows include unchanged and unresolved
assignments. Biological event recovery is assessed separately. The comparison
reports establish agreement for these six tables; likelihood outputs and
rendered figures require their own review.

### Runtime

These are case-task elapsed times and peak resident memory reported in the
execution trace. Each task completed with exit status zero.

| Case | 1 thread | 16 threads | Peak RSS, 1 thread | Peak RSS, 16 threads |
|---|---:|---:|---:|---:|
| RpL32 | 2 min 33 s | 2 min | 108.6 MB | 110.1 MB |
| Bee dsx | 18 min 13 s | 8 min 1 s | 2.1 GB | 2.1 GB |
| Hdac3 | 1 min 40 s | 1 min 30 s | 108.6 MB | 109.5 MB |
| spo5 | 19.6 s | 16.8 s | 115.4 MB | 108.4 MB |
| rec8 | 53.1 s | 31.3 s | 108 MB | 108.7 MB |

### Biological and statistical findings

- **Hdac3:** the expected positive difference was not recovered. The current
  [site summary](/data/projects/intragenic_structure/results/20260916_214900_insiphy/hdac3/threads_1/analysis/results_parsimony/structural_site_summary.tsv)
  contains four presence sites, four exonic-role sites and two junction sites,
  each observed as state 1 in all four species. The whole-exon DNA threshold
  rejects the right-hand segment containing CDS and UTR, leaving the expected
  additional junction out of the structural observations.
- **spo5:** the expected positive difference was not recovered. Its
  [site summary](/data/projects/intragenic_structure/results/20260916_214900_insiphy/spbc29a10_02/threads_1/analysis/results_parsimony/structural_site_summary.tsv)
  contains five presence and five role sites, with no junction site. One
  sequence group is observed in two species; the other groups each have one
  observed species. All recorded observations are state 1 with remaining cells
  unknown, leaving no observed structural contrast.
- **rec8:** the expected positive difference was not recovered. The
  [ER/ARD test table](/data/projects/intragenic_structure/results/20260916_214900_insiphy/spog_00055/threads_1/analysis/results_er-ard/model_tests.tsv)
  reports zero informative patterns in all three fitted layers. The selected
  annotation's partial 5-prime end and strand/phase handling remain relevant
  to the correspondence review.
- **RpL32 and bee dsx:** retain their roles as a conserved-control case and a
  transcript-repertoire case. Their first-iteration results require review
  after the correspondence and statistical repairs. All recorded ER/ARD tests
  across the five cases have `p_value=NA`; this does not resolve the posterior
  defect described below.

The statistical review identified a remaining validity error: layers whose
observed states are all 1 can still be labelled `fit_status=success` and emit
fitted posteriors despite lacking an observed 0/1 contrast. RpL32
[model fits](/data/projects/intragenic_structure/results/20260916_214900_insiphy/rpl32_control/threads_1/analysis/results_er-ard/model_fits.tsv)
and [structural changes](/data/projects/intragenic_structure/results/20260916_214900_insiphy/rpl32_control/threads_1/analysis/results_er-ard/structural_changes.tsv)
provide a concrete example. Those probabilities remain unaccepted. The 97
passing tests did not cover this failure adequately; it is an open repair in
the [publication ledger](publication_gap.md#v014-unified-repair-ledger).

### Authorized next repair

The next iteration will map CDS-derived proteins in full-transcript MAFFT
alignments back to genomic coding intervals, repair minus-strand phase
handling, and classify unestimable models without observed state contrast
before selecting fits or emitting posteriors. These changes require targeted
regression tests and new real-data runs. Recovery of the published candidate
site, its modern sequence correspondence and its directional branch placement
must each be evaluated; the restricted species trees may leave direction
ambiguous. [Case definitions and biological references](real_positive_cases.md)
document the input accessions and outstanding site-mapping questions.

This first iteration establishes completion and six-table agreement between
thread counts. Positive-case recovery and valid fitted posterior reporting
remain required for v0.14 acceptance.

## Historical v0.13 demonstration

<details>
<summary>Archived 2026-09-16 results and commands</summary>

All statements and counts below describe the archived v0.13 run. Current
v0.14 first-iteration findings are recorded above.

## RpL32 conserved control

The first formal single-copy demonstration uses RpL32 orthologs from:

- Drosophila melanogaster, GCF_000001215.4;
- Drosophila simulans, GCF_016746395.2;
- Drosophila erecta, GCF_003286155.1;
- Drosophila yakuba, GCF_016746365.2;
- Drosophila teissieri, GCF_016746235.2.

Genome FASTA and GFF3 files are stored outside the repository under
`/data/db/genome`. The repository contains accession-level manifests and the
species tree.

## Complete user path

```bash
insiphy build-case \
  --manifest examples/real_cases/rpl32_control/manifest.tsv \
  --species-tree examples/real_cases/rpl32_control/species_tree.tsv \
  --output-dir work/rpl32_control \
  --aligner mafft \
  --context-aligner minimap2 \
  --flank 1000 \
  --max-extension 10000 \
  --threads 4

insiphy run \
  --input-dir work/rpl32_control \
  --output-dir results/rpl32_control \
  --analysis-scope single-copy \
  --model parsimony \
  --evidence-aligner miniprot \
  --threads 4

insiphy visualize \
  --input-dir work/rpl32_control \
  --result-dir results/rpl32_control \
  --output-dir figures/rpl32_control
```

These commands show the package's user-facing stages. The server schedules them
with Slurm and an external Nextflow workflow. The package has no Nextflow
dependency. Version 0.13.0 uses the existing dependency container with an explicit
source overlay and four CPUs per case. Default parsimony uses the rooted tree;
optional ER/ARD analysis uses unit branch lengths for RpL32 and supplied lengths
for dsx, with estimated root frequency. No simulations or checksum runs were used.

Completed application run `20260916_202436_insiphy` is stored under
`/data/projects/intragenic_structure/results/20260916_202436_insiphy`.
It reuses the prepared cases archived in run `20260916_164410_insiphy`; genome
extraction and exon comparisons were completed in `20260916_162953_insiphy`.
The current run recomputes annotation evidence, both phylogenetic analyses, and
three SVG figures after repairing deletion-spanning genomic alignment.

## RpL32: conserved control

Three homologous exon groups and two splice boundaries have observations in
all five species. Sequence presence, exon role and both boundaries are conserved
in this matrix. Equal-cost parsimony requires zero changes at these sites.
No annotation-completion query remains for this group.

All three optional ER/ARD comparisons report `parameters_not_estimable` and
`p_value=NA`. With no observed state contrast, these data do not identify
gain/loss asymmetry. The embedded-null optimization start gives the same fitted
likelihood under ER and ARD. Earlier boundary patterns and P values from v0.12
are historical outputs based on different correspondence coding.

This control does not measure sensitivity to known exon gains or losses.

## Bee dsx: alternative exon organization

Seven bee orthologs from OrthoFinder group OG0006454 retain all annotated
transcripts: Apis cerana, Apis mellifera, Bombus ignitus, Bombus pascuorum,
Bombus terrestris, Frieseomelitta varia and Tetragonisca angustula. Genome and
annotation accessions are in the case manifest. No expression or sex-specific
transcript data enter this analysis.

In run `20260916_202436_insiphy`, the raw matrix contains 42 exon-associated
sequence groups, their 42 role observations, and 57 projected junction sites.
These group counts include alternative annotated intervals; they are not the
number of exons in a single transcript. Unknown observations are frequent:

| Layer | Sites | Observed cells | Unknown cells | Sites observed in all seven species |
|---|---:|---:|---:|---:|
| Sequence presence | 42 | 95 | 199 | 1 |
| Exonic role | 42 | 90 | 204 | 0 |
| Splice junction | 57 | 84 | 315 | 0 |

No site has a supported state-0/state-1 contrast in this run. Parsimony requires
zero changes under the observed constraints, with fully missing sites labeled
uninformative. This leaves dsx evolutionary history unresolved; the large number
of unknowns precludes a conclusion of whole-gene structural conservation.

Of 152 supplementary comparisons, 142 remain ambiguous, five suggest a boundary
conflict, four support sequence within already annotated exons, and one predicts
a missing coding interval in B. terrestris (NC_063273.1:2476398-2476441, minus
strand). This 44 bp interval is a protein-projection candidate, not a confirmed
transcribed exon or a dated exon-gain event. Its support score summarizes
alignment identity/coverage and has no posterior-probability interpretation.

Only sites with at least two observed species enter optional likelihood fitting:
24 presence sites, 24 role sites and 18 junction sites. Each layer has zero
observed state contrasts; all three model tests report parameters not estimable
and P=NA. The previously reported P=0.4524 used the older junction coding and
does not describe these repaired observations.

These gene-level likelihoods are conditional on the supplied tree, inferred
homologous elements and annotated isoforms. Adjacent exon and junction sites
are correlated; tests pooling them assume conditional independence. Missing
annotations and unknown states limit inference from the apparent absence of
an exon. Default figures show correspondence and parsimony placement categories.
Optional fitted probabilities condition on estimated parameters; their joint
parameter-sensitivity ranges remain unestimated.

## Deletion evidence and computational limits

The repaired deletion search uses minimap2 on the genomic interval spanning
homologous flanking exons. Both flanks must occur in one supported alignment,
with a query-only gap covering the expected exon and contiguous target bases
beside that gap. Source flanks follow one annotated transcript path. The search
can run without a protein reference, including for noncoding exons.

In the current 152 dsx comparisons, 72 lack a source flanking pair and 72 lack
an observed target homologous flank. Eight comparisons reach genomic alignment:
six have insufficient alignment identity, one lacks joint projection of both
flanks, and one has no gap spanning the expected exon. None supports a deletion.
The earlier internal alignment-size errors are absent from these results.
The structural matrix and annotation-completion counts remain unchanged.

Positive annotation completion in this run uses miniprot. For 122 comparisons,
the chosen reference has no usable protein context; this run does not perform
a supplementary DNA-presence search for those queries. Their unresolved status
therefore reflects both available evidence and this analysis choice. Noncoding
sequence recovery remains an important real-data evaluation requirement.

Each case received four CPUs and 12 GB RAM. Nextflow reports task execution
times of 9.9 s for RpL32 and 64 s for dsx, with peak resident memory of 107.1 MB
and 140.8 MB respectively. These measurements cover annotation completion,
inference and drawing from prepared cases; they exclude the earlier extraction
and pairwise alignments. Both Slurm jobs completed with exit status zero.

Figures in each case's `figures/` directory are `intragenic_synteny.svg`,
`phylogenetic_event_map.svg`, and `integrated_phylo_synteny.svg`. The default
color encoding includes homologous exon connections; the event panel displays
only the changes permitted by the reported reconstruction.

## Next real-data set

The next benchmark should contain single-copy genes with independently
documented structural changes. Candidate cases must satisfy:

- one ortholog per species;
- assembly and annotation versions available;
- enough species to distinguish alternative branch placements;
- sequence-level evidence for the altered exon or junction;
- a published history that can be reviewed independently of IntraPhy.

Jingwei and Sdic remain useful for future multi-copy development. Their
duplication histories place them outside the current formal single-copy
benchmark.

</details>

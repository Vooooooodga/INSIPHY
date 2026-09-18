# Publication Readiness Review

Date: 2026-09-18

Final v0.14.0 assessment: **136 selected formal regression tests passed** under Slurm `61625` in 17.063 s ([report](/data/projects/intragenic_structure/results/20260918_121100_insiphy/regression_tests.txt)), including four new role-conflict tests. All ten re-inference/visualization tasks and five thread-comparison tasks completed. All six core tables agree within every pair. Existing cases/evidence were reused without repeated alignments. Scheduled execution and targeted checks are complete; scientific limitations remain below.

Historical intermediate baseline `20260918_111929_insiphy` completed at 12:10:30 on 2026-09-18 after recovery. Its 132-test result, three dsx role contrasts and conditional ER posteriors predate the final conflict repair. They are retained as historical records, separate from the completed final assessment.

Final focal audit: Hdac3 recovers the exact 62 bp intron contrast (1/3/0 species present/absent/unknown), with Dana split and melanogaster-subgroup fusion equally parsimonious. rec8 recovers the exact 39 bp interval boundary (1/1/2), with Cryo fusion and Octo split equally possible; direct intron DNA correspondence and role contrast remain missing, so recovery is partial. spo5's 55 bp focal difference is unrecovered. Final dsx has two observed role contrasts; ER is nonidentifiable and all node/branch posterior tables are header-only. Written SVG endpoints in the inspected Hdac3/rec8 matches agree with selected coding blocks; no rendered inspection was performed. See the [final observations](real_data_benchmark.md#final-biological-assessment).

## Decision

INSIPHY v0.14.0 completes the scheduled single-copy correctness assessment with 136 selected tests and five paired real-case runs. The dsx conflict repair has targeted-test and final-output evidence. Earlier 97- and 132-test iterations remain historical. Scientific recovery is uneven, observations are incomplete, and likelihood parameters remain insufficiently determined in these cases. The project is not yet publication-ready.

## Test evidence received

The [final report](/data/projects/intragenic_structure/results/20260918_121100_insiphy/regression_tests.txt) records 136 selected formal regression tests. Legacy experimental and toy test classes were intentionally excluded. Targeted checks cover:

- observation states, predicted roles, transcript repertoire, exact ID mapping and raw-feature retention;
- tree validity, zero-length transitions, pruning and posterior calculations against enumeration, expected transition counts, ascertainment and unsuccessful fits;
- SVG semantics for valid selected-fit posteriors, unknown/predicted roles, cross-species ribbons, separate transcript lanes and missing-length layout.
- full-transcript coding projections, missing/negative-strand phase, observed-species conservation summaries, no-contrast fitted-posterior exclusion, and partial-CDS/1:n ribbon endpoints.

The [repair ledger](publication_gap.md#v014-unified-repair-ledger) records the specific tested behaviors and remaining empirical work. SVG semantic tests establish the checked drawing behavior; readability of the real figures still requires inspection.

## Historical first-iteration findings

In the first iteration, Hdac3, spo5 and rec8 did not recover their focal positive differences. The Hdac3 whole-exon DNA threshold rejected the right-hand segment containing CDS and UTR, leaving the expected junction out of the structural observations. These historical findings motivated full-transcript CDS/protein mapping and phase repairs. The repaired Hdac3 result now recovers the focal contrast; rec8 remains partial and spo5 unrecovered.

The first-iteration statistical review found all-1/no-observed-contrast layers labelled as successful fits and used to produce posteriors. The repair passed targeted regression, and the current published RpL32 tables now record `posterior_not_reported` in all three layers. The earlier 97-test run did not cover this failure adequately. The figure check uses the selected fit's status, with validity determined in the statistical calculation.

[First-iteration records, table counts and runtime measurements](real_data_benchmark.md#v014-first-iteration-provisional-results) document that historical run. They remain separate from intermediate-baseline full-analysis measurements and the final inference-only assessment.

## Current statistical assessment

The repaired RpL32 run contains 20 sites (5 presence, 5 role, 10 junction), all with zero minimum changes and no observed contrast; coverage remains incomplete. Its six likelihood comparisons report unavailable P/Q values and no fitted posteriors. The [current RpL32 record](real_data_benchmark.md#rpl32-from-the-repaired-run) uses the newly published files.

In final dsx outputs, Bter EG_0029 role is unknown and Bign exonic usage is retained; that site's minimum changes decrease from 1 to 0. Two contrasts remain, EG_0013 and EG_0055, each with five unknown species and respectively two and three possible branch placements. Presence has 74 known/94 unknown observations, role 63 known (61 exonic, 2 non-exonic)/105 unknown, and junction 65 known/145 unknown. Eligible role-fit sites decrease from 18 to 17. Role ER is now nonidentifiable; ARD and foreground fits are boundary-limited. All six P/Q values are NA and all four node/branch posterior tables are header-only. Role LRT statistics 3.06583 and 0.224956 are diagnostics only. The historical baseline's valid conditional ER posteriors do not describe this final matrix.

## What is scientifically coherent now

- The method starts from an upstream single-copy ortholog set.
- The gene-internal observations are explicit: sequence presence, exonic role and splice junction.
- Predicted roles remain prediction evidence.
- Missing annotation and unresolved correspondence remain unknown.
- Default analysis is qualitative parsimony branch placement.
- Optional likelihood analyses are model comparisons on a fixed tree.
- P values are asymptotic model-comparison quantities, not event-truth probabilities.
- Mechanism interpretation is left to downstream biological analysis.

## What changed in v0.14

- all-transcript repertoire is treated as the main structure representation;
- predicted CDS evidence is separated from confirmed CDS homology;
- figures require selected-fit success and finite endpoint probabilities; the no-contrast repair passed regression and is observed in the repaired RpL32 output;
- confirmed homology ribbons exclude predicted, unknown and candidate-source members; block-clipping tests for partial CDS and complementary split projections passed;
- extant path and coverage output names no longer imply ancestral graph reconstruction;
- code policy now centers on interface validation, adapters around mature tools, and explicit failure states.

## Evidence still missing

- broader biological assessment beyond these five completed case pairs;
- direct DNA correspondence and role recovery for the partially recovered rec8 focal interval, and recovery of spo5;
- broader curated positive-case coverage beyond the audited focal examples;
- conserved controls beyond one showcase gene;
- visual inspection of generated SVGs;
- runtime and memory scaling beyond the recorded five-case workload;
- written comparison against annotation-only and sequence-only ablations.

## Acceptance criteria for the first manuscript

1. Every showcase case starts from raw genome and annotation files.
2. The selected genes satisfy the single-copy input condition.
3. Structural matrices, correspondence tables and figures are manually reviewed.
4. At least one positive case demonstrates a known internal structural change.
5. At least one conserved-control panel estimates the unsupported-change rate.
6. Optional likelihood results are reported only when `fit_status`, identifiability and endpoint probabilities are valid.
7. The manuscript states that multi-copy analysis is deferred.

## Current risk assessment

The no-contrast posterior correction now has regression and current RpL32 output evidence. Remaining risks include partial rec8 recovery, unrecovered spo5, uncertain dsx role evidence, and incomplete estimation of parameter uncertainty. The selected tests and focal audits support their specific checked behaviors. Broader sequence correspondence, state construction and branch inference performance remain unmeasured. Known single-copy genes enter as homologous loci; internal sequence units and splice junctions become structural sites; branch histories are reconstructed under defined parsimony or likelihood assumptions.

Broader empirical comparisons and resolution or characterization of the remaining scientific limitations are required before general performance or biological discovery claims.

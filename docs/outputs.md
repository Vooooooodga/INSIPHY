> **V18 legacy documentation / retained reference.** The default V19 model and
> commands are described in [the README](../README.md) and
> [exon_structure_model.md](exon_structure_model.md). Do not interpret the old
> independent-layer analyses as exon configuration inference.

# Output interpretation

## Evidence and characters

`element_correspondence.tsv` records local element memberships, actual matched
intervals and eligibility. `splice_boundary_correspondence.tsv` records exact
projected positions. Unresolved candidates are evidence records, not hard
phylogenetic observations.

`character_catalogue.tsv` defines one row per family/layer/site. `count_unit_id`
identifies the character whose state transitions are counted. `coordinates_available`
is an evidence-availability flag, not an accuracy score. `character_coordinates.tsv`
contains actual genomic intervals (1-based closed) or explicitly named reference
splice coordinates. `character_dependencies.tsv` records linked character IDs
and the basis; `shared_mutation_inferred` is false.

`structural_site_matrix.tsv` is the complete effective matrix.
`analysis_structural_site_matrix.tsv` is the requested subset. Input values,
effective masks and reasons appear in `structural_observation_eligibility.tsv`.
Per-character coverage includes called species, MRCA and representation of root
subtrees. The native annotation background is unchanged by selection.

## Parsimony

`node_structural_states.tsv` and `branch_structural_events.tsv` retain all-optimum
sets. Event IDs identify a character, branch and directed endpoint change.
`possible` placements can be mutually exclusive; never count table rows as a
historical event total. A split/fusion descriptor does not add an event.

`structural_site_summary.tsv` gives each character's minimum.
`gene_change_summary.tsv` sums within each layer and explicitly leaves mutation
event counts unestimated. `minimum_change_history.tsv` gives one compatible
minimum-change reconstruction with probability NA. `ancestral_state_consistency.tsv`
flags conflicts across separately inferred sequence-presence and exon-role layers.
No complete ancestral transcript is reconstructed.

## Likelihood

`model_fits.tsv` contains fit status and available parameters, profiles and AIC.
`model_tests.tsv` distinguishes valid tests, unidentifiable parameters, linked
characters, failed optimization and other unavailable outcomes. `p_value` and
`q_value` concern rate-model comparisons only.

`node_state_posteriors.tsv`, `branch_transition_posteriors.tsv` and
`structural_changes.tsv` describe conditional probabilities and expected counts
only when an eligible fit exists. Known dependent-character blocks produce no
formal fit or posterior. Unavailable fields are NA, not zero and not evidence
for conservation.

## Visualization

`index.html` indexes per-target structure, correspondence and history SVGs.
Only the target's actual matched blocks receive correspondence ribbons; native
annotation remains the full background. Unknown or inapplicable observations are
not displayed as absence. A history view describes the selected target's
conditional reconstruction, not a complete ancestral gene model.

## Execution provenance

Primary commands write `execution.json`, `environment.json`, `intraphy.log` and
where relevant `run_result.json`. A failed run cannot be represented as a previous
successful run. Raw-input synthetic expected files are evaluation references and
are not consumed by the analysis.

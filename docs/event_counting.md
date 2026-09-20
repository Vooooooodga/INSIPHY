# What is counted?

## Four distinct objects

A homologous sequence interval is a proposed correspondence. A structural
character is a property of that correspondence. A character-state transition is
a change under a phylogenetic model. A molecular mutation is a DNA process.
These objects are not interchangeable.

A single contiguous deletion may affect two exons and several junctions. Coding
those observations separately does not establish several independent mutations.
Conversely, one splice-position character may gain and lose introns repeatedly
across a tree. Representing a whole gene as one multistate character would still
require a transition model; it would not restrict the gene to one event.

## Elementary changes in v18

The formal output includes sequence-presence gain/loss, exon-role gain/loss and
intron/junction gain/loss. These are elementary changes of the stated binary
characters. `exon_split` or `exon_fusion` can describe a junction change without
creating another event. Multiple cutpoints remain multiple characters. The
program does not infer compound events or simultaneous mutations.

Each directed branch change has a unique family/character/branch/direction ID.
Different views of the same change share that identity. A `possible` placement
is an alternative compatible with at least one optimal history; possible
placements cannot be added to obtain a historical event count.

The gene summary is stratified by layer. `minimum_character_changes` adds the
site-specific minima under the stated model. `mutation_event_count` is always
`not_estimated`. The deterministic history is one compatible witness; it has no
assigned probability. All-optimum sets remain the primary parsimony result.

## Why a character catalogue is retained

A drawing may split one exon into several visual rectangles. This does not
increase the character count. Actual local correspondence refinement can create
distinct biological hypotheses, so a changed correspondence algorithm may change
the catalogue. Such runs must not be presented as identical character universes.

An independent interval-based indel catalogue could code a supported contiguous
deletion as one character. V18 does not automatically infer minimal mutation
intervals from a collection of missing exon observations. Indel coding requires
explicit boundary and applicability rules (see GapCoder in the references).

Dependence metadata remains necessary without compound-event interpretation.
A linked group records shared observation structure, not proof of a shared
mutation. Never delete linkage merely to enable an independent-character LRT.

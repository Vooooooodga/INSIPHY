# IntraPhy maintenance rules

Use the `intraphy` package and command only. Active code, messages and documentation
are English. Preserve established biological terms and define technical fields.

Scientific observations belong in `observations/`. Inference consumes a frozen
matrix and must not change homology or invent negative states. DNA presence,
annotation-conditional exon identity and splice positions have distinct semantics.
Unknown, absent, inapplicable and predicted-only evidence remain distinct.

Count structural character changes, not asserted mutation events. Do not add
possible placements or double-count split/fusion descriptions. Do not infer
compound events. Keep dependence metadata and block all independent-character
likelihood outputs for known linked included sites.

Use actual matched intervals; do not infer complete intronic DNA homology from
flanking exon correspondence. References do not imply ancestors. Repertoire
means supplied annotated paths, not all biological transcript usage.

Keep production modules at or below 500 physical lines. Prefer cohesive explicit
owner modules and ordinary imports. No runtime source extraction, forwarding old
namespace, content-hash pipeline or unnecessary framework. Patch tests at the
actual implementation owner. Never relax evidence requirements to repair a test.

Before a release, run source-layout checks, the complete test suite with real
MAFFT/minimap2, clean wheel smoke tests and raw-input integration. Archive actual
versions and failures. Separate unit/integration success from unperformed
biological benchmarking and statistical calibration. Preserve historical logs.

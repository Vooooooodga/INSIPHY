# V19 implementation and limits

| Plan item | Implemented code | Verification / limitation |
|---|---|---|
| FASTA/GFF/tree, no manual manifest | `commands/exons.py`, existing `inputs/` | `analyze`, direct loci and whole-genome ortholog selectors |
| Complete native exons, not alignment-fragment events | `structure/native.py`, `types.py` | Source IDs and local comparison spans both retained |
| Protein + genomic evidence | existing coding projection + `structure/alignment.py`, `corroboration.py` | Real MAFFT/minimap2; alignment conflicts unresolved |
| Predicted-only exon/boundary alternatives | `structure/build.py`, `observations.py` | Evidence and annotation views separate; no tree-driven repair |
| True coexistence vs uncertainty | Typed observations and conditional scenarios | No usage estimates, no CTMC for an unmodeled repertoire |
| Local configurations and material provenance | `structure/material.py`, `space.py` | Complete enumeration only within declared finite boundaries |
| Elementary split/fusion/shifts/interval indels | `structure/edits.py` | Shared registry; consequences not extra events |
| Reversible property vs irreversible source loss | `structure/origins.py`, `edits.py` | Explicit single introduction-opportunity model, not universal Dollo |
| All-optimal parsimony and compatible witness | `configuration_dp.py`, `configuration_history.py` | Exhaustive small-tree comparisons and edit-count counterexamples |
| CTMC and ancestor/branch quantities | `configuration_ctmc.py`, `configuration_model.py` | Matrix, enumeration, no-jump and marked-count numerical tests |
| Shared rates and foreground contrast | `exon_rates.py`, `exon_resampling.py` | Fixed relative rates + fitted scalar; validity/discovery gates |
| Gene-level bootstrap | `exon_rates.py` | Whole genes sampled, failed draws retained |
| CDS phase and UTR handling | `structure/consequences.py` | Phase 0/1/2 exact-fusion raw tests; exceptions explicit |
| CESAR2 optional adapter | `aligners/cesar_adapter.py`, `realign-exons` | Input/command interface tests only; not automatically used as truth |
| Copy, inversion and complex sources | `structure/alignment.py`, unresolved scope | Detection/abstention, NOT copy/rearrangement history |
| Raw validation infrastructure | `verification/exon_cases.py`, `tools/validate_v19.py` | 20 synthetic scenarios; not empirical biological accuracy |
| Read-only actual result diagrams | `reporting/exon_results.py` | Tree + native exon structures + IDs from saved event rows |

The design's real-case benchmarks (FDPS, MAMSTR, SLC7A6, KANK1, EF-1alpha),
comparison-method benchmark runs and broad discovery-aware statistical calibration
are not fabricated as completed work. V19 supplies the implementation and test
infrastructure; these remain empirical acceptance work for a methods paper.

State/source/candidate caps are observable failure states, not hidden truncation.
Annotation evidence and costs are not calibrated likelihoods. The core engine does
not claim every candidate alignment or molecular mechanism has been resolved.

A gene-only selected locus can be retained as unknown. For staged `build-case`, use
`--allow-unannotated-loci` explicitly; `analyze` retains such loci by default. No
exon is synthesized from the span. Entirely missing/unlocated gene loci still require
upstream locus identification rather than an invented negative observation.

Optional CESAR2 example (tool and taxon-appropriate profiles supplied by the user):

```bash
intraphy realign-exons --input-dir result/prepared_inputs --family-id family \
    --reference-species Species_A --query-species Species_B \
    --reference-transcript transcript_ID --cesar /path/to/cesar \
    --profile-dir /path/to/chosen/profiles --codon-matrix /path/to/codon_matrix.txt \
    --output-dir cesar_prediction
```

The complete standard-code reference requirements apply to this adapter, not to
whether native UTR/frameshift-containing structures are retained by IntraPhy.
Predictions remain separate; there is no automatic GFF replacement or evidence
promotion. The optional adapter does not ship CESAR source, profiles or benchmarks.

# Quick start: ordinary genomic files

Install IntraPhy, MAFFT and minimap2 as described in [installation](installation.md).
An ordinary run uses **genomic FASTA + GFF/GTF + a supplied rooted species tree**.
It does not require a user-written manifest. Read [input details](inputs.md) for
whole-genome selection, combined locus FASTA, AGAT and coordinate conventions.

```bash
intraphy example --output-dir example
intraphy check --fasta example --gff example --species-tree example/species_tree.nwk
intraphy build-case --fasta example --gff example \
  --species-tree example/species_tree.nwk --output-dir work/example --threads 2
intraphy run --input-dir work/example --output-dir results/example --threads 2
intraphy visualize --input-dir work/example --result-dir results/example \
  --output-dir figures/example
```

Open `figures/example/index.html`. The linked methods guide uses labelled
synthetic examples; the target panels below it use the actual saved run results.
The standalone guide needs no input data or alignment tools:

```bash
intraphy explain --output-dir method-guide
```

For whole-genome annotations add `--orthologs OG0001.fa` to `check`/`build-case`.
FASTA member IDs select exact annotated loci; their sequences are not a substitute
for genomic DNA. File stems match tree tips. All supplied target transcript paths
are retained by default. No separate upstream/downstream files are required.
An unavailable CTMC fit is not evidence that no event occurred. Read the
[counting rules](event_counting.md) before totaling branch rows.

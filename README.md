# IntraPhy

**Phylogenetic analysis of gene structure · version 0.18.0**

Which homologous regions of a gene are retained? Which remain exonic, and
where have intron positions changed? IntraPhy compares **single-copy orthologous
genes** using genomic sequences, gene annotations and a rooted species phylogeny.
It reconstructs local structural character changes, not physical mutation counts.

![From homologous regions to ancestral gene-structure changes](docs/figures/method_overview.svg)

*One gene family, three steps: establish local sequence correspondence; define
comparable DNA, exonic-status and intron-position characters; reconstruct their
evolution on the species tree. The example includes missing annotation, DNA loss
and intron loss, not only exon splitting. Drawings are synthetic; branch labels
illustrate minimum-change placements. [Detailed algorithm](docs/figures/homology_inference.svg)
· [Probability model](docs/figures/phylogenetic_model.svg)
· [Mathematical explanation and limits](docs/model_bridge.md).*

## Install

Python 3.10 or later, MAFFT and minimap2 are required for the standard workflow.
On Ubuntu 24.04:

```bash
sudo apt-get update
sudo apt-get install -y mafft minimap2
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -c requirements-ci.txt .
intraphy --version
intraphy inspect-aligners
```

The executable, Python package and module entry point are all `intraphy`.
There is no previous-name runtime alias. AGAT and miniprot are optional;
[installation details](docs/installation.md) describe their separate roles.

## Your inputs: FASTA + GFF + rooted tree

**No user-written manifest is required.** The standard input layout is:

```text
genomes/                 annotations/              orthologs/
  Species_A.fa             Species_A.gff3            OG0001.fa
  Species_B.fa             Species_B.gff3            OG0002.fa
  Species_C.fa             Species_C.gff3
species_tree.nwk
```

FASTA/GFF file stems match the species-tree tip names. GFF sequence IDs and
coordinates refer to the matching genomic FASTA. Each optional ortholog FASTA
contains the selected upstream gene-family members, identified by exact gene,
transcript or protein IDs present in the annotation. AGAT-style `gene=` metadata
is accepted. Multiple isoforms of one locus are allowed; paralogous loci are not.

```bash
intraphy check \
  --fasta genomes/ --gff annotations/ --orthologs orthologs/ \
  --species-tree species_tree.nwk

intraphy build-case \
  --fasta genomes/ --gff annotations/ --orthologs orthologs/ \
  --species-tree species_tree.nwk --output-dir work/genes --threads 8

intraphy run --input-dir work/genes --output-dir results/genes --threads 8

intraphy visualize \
  --input-dir work/genes --result-dir results/genes --output-dir figures/genes
```

If every species GFF already contains exactly one selected gene, omit
`--orthologs`. A **combined genomic locus FASTA** is also accepted with separate
species GFFs, provided different species use distinct FASTA record IDs and the
GFF coordinates refer to those records:

```bash
intraphy build-case \
  --fasta orthologous_genomic_loci.fa --gff locus_annotations/ \
  --species-tree species_tree.nwk --output-dir work/one_gene
```

Spliced CDS, transcripts and proteins **cannot provide introns or genomic
flanks**. Such a FASTA can select members through `--orthologs`, but genomic DNA
is still needed for structure analysis. IntraPhy does not infer gene orthology
from FASTA similarity. [Input formats, ID matching and coordinate requirements](docs/inputs.md).

## Flanking regions are extracted automatically

`build-case` reads the target gene and available flanks from genomic FASTA.
There is no separate flank FASTA to prepare. `--flank 1000` and
`--max-extension 10000` specify search limits, not homology or accuracy thresholds.
Sequence beyond the supplied contig/crop is unavailable and is never invented.

For reusable small inputs, an optional export writes matching locus FASTA/GFF
pairs with **rebased coordinates**, retaining strand and CDS phase:

```bash
intraphy extract-loci \
  --fasta genomes/ --gff annotations/ --orthologs orthologs/ \
  --species-tree species_tree.nwk --flank 1000 --output-dir loci
```

Each family directory can then be supplied directly to `build-case`. A source
coordinate map records the actual available flanks. AGAT normalization is also
available as an explicit, separately logged operation:

```bash
intraphy normalize-annotation --gff annotations/ --output-dir normalized_annotations
```

This command requires `agat_convert_sp_gxf2gxf.pl`. Annotation normalization may
change IDs, boundaries or features; it is not independent biological evidence.
[AGAT usage and limitations](docs/inputs.md#optional-agat-preparation).

## Run the synthetic raw-input example

```bash
intraphy example --output-dir example
intraphy check --fasta example/ --gff example/ --species-tree example/species_tree.nwk
intraphy build-case \
  --fasta example/ --gff example/ --species-tree example/species_tree.nwk \
  --output-dir work/example --threads 2
intraphy run --input-dir work/example --output-dir results/example --threads 2
intraphy visualize \
  --input-dir work/example --result-dir results/example --output-dir figures/example
```

Open `figures/example/index.html`. No inferred homology or observation matrix is
supplied to this example. `--scenario conserved` and `--scenario
annotation_dropout` provide additional controls. The example is an implementation
test; it does not establish biological accuracy or statistical calibration.

## See what each stage means

[How sequence alignment establishes structural correspondence](docs/figures/homology_inference.svg) ·
[How gene structures become a phylogenetic model](docs/figures/phylogenetic_model.svg) ·
[Figure captions and literature](docs/figure.md)

```bash
# No genomic data or aligners needed for the teaching gallery.
intraphy explain --output-dir method-guide
```

Every default result gallery links to a separate synthetic methods guide.
Data-derived target figures are never replaced by teaching examples.

The [software architecture](docs/architecture.md) is documented separately.
The full character matrix retains unknown and inapplicable observations. Sequence
correspondence is evaluated before ancestral reconstruction; the reference species
is not assumed to be ancestral.


## Principal results and interpretation

| Output | Meaning |
|---|---|
| `input_targets.tsv` | Automatically resolved genes and source files; an output, not a required input table |
| `character_catalogue.tsv` | Character identity, state definition, counting unit and dependence metadata |
| `character_coordinates.tsv` | Actual aligned intervals or projected splice positions |
| `structural_site_matrix.tsv` | Complete structural character matrix, including unknown and inapplicable observations |
| `structural_character_eligibility.tsv` | Known 0/1 observations and their phylogenetic coverage |
| `branch_structural_events.tsv` | Elementary changes in all or some minimum-change histories |
| `gene_change_summary.tsv` | Minimum character changes per gene family and layer; no mutation count |
| `minimum_change_history.tsv` | One globally compatible optimum, not a probability sample |
| `ancestral_state_consistency.tsv` | Applicability conflicts across separately reconstructed layers |

A single mutation can affect several characters; a single character can change
repeatedly. Alternative possible placements are not added as separate changes.
There is no formal compound-event summary. [Counting rules](docs/event_counting.md).

The default retains all characters. Optional coverage filtering does not trim
DNA or connect previously separated exons. Seven present plus three explicitly
absent observations are 100% callable, not 70%.

## Optional likelihood analysis

ER/ARD and foreground CTMC analyses use the same structural character matrix:

```bash
intraphy infer-phylogeny \
  --input-dir work/genes \
  --structural-site-matrix results/genes/structural_site_matrix.tsv \
  --output-dir results/genes_erard --model er-ard
```

Known linked characters jointly included in a family-layer block the independent
fit, AIC, intervals, test and posterior. Other invalid fits remain unavailable
with explicit reasons. P values compare rate models; they do not establish a
named historical event. Finite-sample calibration remains unassessed.
[Statistical assumptions](docs/statistical_model.md).

## Limits and reproducibility

An unaligned interval is not automatically absent. Corresponding flanking exons
do not establish homology of every intronic nucleotide. Supplied transcript
repertoire does not measure tissue-specific usage. Formal inference excludes
complete ancestral transcripts, exon-shuffling mechanisms, complex rearrangements,
gene-duplication histories, selection and phenotypic causes.
[Scope and the 28 cases](docs/scope_policy.md).

Primary commands write `intraphy.log`, `environment.json` and `execution.json`.
Nonempty output directories are refused. `--force` preserves an identifiable old
IntraPhy output in a timestamped backup. Run `intraphy COMMAND --help` for options.

[Method](docs/method.md) · [Architecture](docs/architecture.md) ·
[CLI](docs/cli.md) · [Outputs](docs/outputs.md) ·
[Migration](docs/MIGRATION_0.18.md) · [Validation](docs/validation.md) ·
[References](docs/references.md)

```bash
python tools/check_source_layout.py
python -m unittest discover -s tests -v
python tools/render_method_figures.py --png
python -m build
```

Historical v16/v17 evidence remains labelled as historical. Current verification
distinguishes unit tests, real alignment-tool integration, synthetic raw inputs,
and unperformed biological/statistical calibration.

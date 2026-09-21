"""Raw genomic/GFF examples with structural changes OR observation damage.

The truth is a separate evaluation file. Preparation and inference never read it.
These fixtures validate contracts, not genome-wide operating characteristics.
"""
from __future__ import annotations
import json
from pathlib import Path
import random
from ..aligners.columns import revcomp

SCENARIOS = ("conserved", "split_insertion", "intronization", "fusion_phase0", "fusion_phase1", "fusion_phase2",
             "exon_deletion", "multi_exon_deletion", "annotation_dropout", "assembly_gap",
             "donor_shift", "acceptor_shift", "neutral_upstream_indel", "duplication", "inversion",
             "negative_strand", "coexisting", "utr", "microexon", "noncanonical")
TAXA = ("Species_A", "Species_B", "Species_C", "Species_D")


def write_exon_example(output_dir, scenario="split_insertion", seed=19):
    if scenario not in SCENARIOS:
        raise ValueError("Unknown V19 raw scenario")
    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    if any(root.iterdir()):
        raise ValueError("Example output must be empty")
    rng = random.Random(seed)
    dna = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    codons = [a+b+c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a+b+c not in {"TGA", "TAA", "TAG"}]
    phase = int(scenario[-1]) if scenario.startswith("fusion_phase") else 0
    lengths = [180+phase, 9 if scenario == "microexon" else 180, 180-phase, 180]
    coding = "ATG"+"".join(rng.choice(codons) for _ in range((sum(lengths)-6)//3))+"TAA"
    exons, offset = [], 0
    for length in lengths:
        exons.append(coding[offset:offset+length])
        offset += length
    introns = [("GC" if scenario == "noncanonical" else "GT")+dna(76)+"AG" for _ in range(3)]
    prefix, suffix = dna(90), dna(90)
    ancestral, positions = prefix, []
    for i, exon in enumerate(exons):
        positions.append((len(ancestral), len(ancestral)+len(exon)))
        ancestral += exon
        if i < len(introns):
            ancestral += introns[i]
    ancestral += suffix
    truth = {"scenario": scenario, "seed": seed, "biological_accuracy_benchmark": False,
             "changed_tip": "Species_D", "molecular_mutation_count_claimed": False}
    selectors = []
    for species in TAXA:
        seq, spans = ancestral, list(positions)
        if species == "Species_D":
            def delete(a, b):
                nonlocal seq, spans
                seq = seq[:a]+seq[b:]
                adjusted = []
                for x, y in spans:
                    if y <= a:
                        adjusted.append((x, y))
                    elif x >= b:
                        adjusted.append((x-(b-a), y-(b-a)))
                    else:
                        left, right = min(x, a), max(a, y-(b-a))
                        if left < right:
                            adjusted.append((left, right))
                spans = adjusted
            if scenario in {"split_insertion", "negative_strand"}:
                a, b = spans[1]
                cut, insert = a+90, "GT"+dna(76)+"AG"
                seq = seq[:cut]+insert+seq[cut:]
                spans = [spans[0], (a, cut), (cut+80, b+80), *((x+80, y+80) for x, y in spans[2:])]
            elif scenario == "intronization":
                a, b = spans[1]
                spans = [spans[0], (a, a+60), (a+90, b), *spans[2:]]
            elif scenario.startswith("fusion_phase"):
                start, end = spans[0][1], spans[1][0]
                delete(start, end)
                spans = [(spans[0][0], spans[1][1]), *spans[2:]]
            elif scenario in {"exon_deletion", "multi_exon_deletion"}:
                a, b = spans[1][0]-20, spans[2 if scenario == "multi_exon_deletion" else 1][1]+20
                delete(a, b)
            elif scenario == "annotation_dropout":
                spans.pop(1)
            elif scenario == "assembly_gap":
                a, b = spans[1]
                seq = seq[:a]+"N"*(b-a)+seq[b:]
            elif scenario == "donor_shift":
                spans[1] = (spans[1][0], spans[1][1]-3)
            elif scenario == "acceptor_shift":
                spans[1] = (spans[1][0]+3, spans[1][1])
            elif scenario == "neutral_upstream_indel":
                seq = seq[:40]+dna(17)+seq[40:]
                spans = [(a+17, b+17) for a, b in spans]
            elif scenario == "duplication":
                a, b = spans[1]
                insert = "GT"+dna(76)+"AG"+seq[a:b]
                seq = seq[:b]+insert+seq[b:]
                spans = [spans[0], spans[1], (b+80, b+len(insert)), *((x+len(insert), y+len(insert)) for x, y in spans[2:])]
            elif scenario == "inversion":
                a, b = spans[1]
                seq = seq[:a]+revcomp(seq[a:b])+seq[b:]
        contig, gene, tx = species+"_chr", species+"_gene", species+"_tx"
        negative = scenario == "negative_strand"
        strand = "-" if negative else "+"
        genomic = revcomp(seq) if negative else seq
        def coords(a, b):
            return (len(seq)-b+1, len(seq)-a) if negative else (a+1, b)
        left, right = coords(spans[0][0], spans[-1][1])
        lines = [f"{contig}\tfixture\tgene\t{left}\t{right}\t.\t{strand}\t.\tID={gene}"]
        paths = [(tx, spans)]
        if scenario == "coexisting" and species == "Species_D":
            changed = list(spans)
            changed[1] = (changed[1][0], changed[1][1]-3)
            paths.append((tx+"_alternative", changed))
        for transcript, path in paths:
            lines.append(f"{contig}\tfixture\tmRNA\t{left}\t{right}\t.\t{strand}\t.\tID={transcript};Parent={gene}")
            coding_length = 0
            for i, (a, b) in enumerate(path):
                x, y = coords(a, b)
                lines.append(f"{contig}\tfixture\texon\t{x}\t{y}\t.\t{strand}\t.\tID={transcript}_e{i+1};Parent={transcript}")
                if scenario != "utr":
                    cds_phase = (3-coding_length % 3) % 3
                    lines.append(f"{contig}\tfixture\tCDS\t{x}\t{y}\t.\t{strand}\t{cds_phase}\tID={transcript}_c{i+1};Parent={transcript}")
                    coding_length += b-a
        (root/f"{species}.fa").write_text(f">{contig}\n{genomic}\n")
        (root/f"{species}.gff3").write_text("##gff-version 3\n"+"\n".join(lines)+"\n")
        selectors.append(f">{tx} gene={gene}\n{''.join(seq[a:b] for a,b in spans)}\n")
    (root/"orthologs").mkdir()
    (root/"orthologs/example.fa").write_text("".join(selectors))
    (root/"species_tree.nwk").write_text("((Species_A:0.2,Species_B:0.2):0.2,(Species_C:0.2,Species_D:0.2):0.2)root;\n")
    (root/"truth.json").write_text(json.dumps(truth, indent=2)+"\n")
    return root

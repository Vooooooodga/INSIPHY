"""Small synthetic FASTA/GFF3 inputs for integration tests, not calibration data."""
from pathlib import Path
import json
import random

from ..storage.tabular import write_tsv

SPECIES = ("Taxon_A", "Taxon_B", "Taxon_C", "Taxon_D")


def build_native_example(output_dir, scenario="splice_difference", seed=18):
    """Create raw inputs only; no inferred homology or character table is supplied.

    The truth file is for independent evaluation. The analysis never reads it.
    Annotation dropout changes the supplied annotation; DNA remains unchanged.
    """
    if scenario not in {"conserved", "splice_difference", "annotation_dropout"}:
        raise ValueError("Unsupported native example scenario")
    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    if any(root.iterdir()):
        raise ValueError("Example output directory must be empty")
    rng = random.Random(seed)
    codons = [a+b+c for a in "ACGT" for b in "ACGT" for c in "ACGT"
              if a+b+c not in {"TAA", "TAG", "TGA"}]
    exons = ["".join(rng.choice(codons) for _ in range(60)) for _ in range(3)]
    exons[0] = "ATG" + exons[0][3:]
    exons[-1] = exons[-1][:-3] + "TAA"
    def random_dna(length):
        return "".join(rng.choice("ACGT") for _ in range(length))
    manifests = []
    member_records = []
    for species in SPECIES:
        segments = list(exons)
        if species == "Taxon_D" and scenario == "splice_difference":
            segments = [exons[0], exons[1][:90], exons[1][90:], exons[2]]
        sequence, positions = random_dna(100), []
        for segment in segments:
            start = len(sequence)+1
            sequence += segment
            positions.append((start, len(sequence)))
            sequence += "GT" + random_dna(76) + "AG"
        sequence += random_dna(100)
        contig, gene, transcript = species + "_chr", species + "_gene", species + "_transcript"
        start, end = positions[0][0], positions[-1][1]
        records = [f"{contig}\tsynthetic\tgene\t{start}\t{end}\t.\t+\t.\tID={gene}",
                   f"{contig}\tsynthetic\tmRNA\t{start}\t{end}\t.\t+\t.\tID={transcript};Parent={gene}"]
        for index, (left, right) in enumerate(positions, 1):
            if scenario == "annotation_dropout" and species == "Taxon_D" and index == 2:
                continue
            records.extend([
                f"{contig}\tsynthetic\texon\t{left}\t{right}\t.\t+\t.\tID={transcript}_exon{index};Parent={transcript}",
                f"{contig}\tsynthetic\tCDS\t{left}\t{right}\t.\t+\t0\tID={transcript}_cds{index};Parent={transcript}"])
        (root / f"{species}.fa").write_text(f">{contig}\n{sequence}\n")
        (root / f"{species}.gff3").write_text("##gff-version 3\n" + "\n".join(records) + "\n")
        member_records.append(f">{transcript} gene={gene}\n{''.join(segments)}\n")
        manifests.append({"species": species, "family_id": "example_gene", "gene_id": gene,
                          "genome_fasta": f"{species}.fa", "annotation_file": f"{species}.gff3"})
    (root / "orthologs").mkdir()
    (root / "orthologs/example_gene.fa").write_text("".join(member_records))
    write_tsv(root / 'manifest.tsv', manifests,
              ['species', 'family_id', 'gene_id', 'genome_fasta', 'annotation_file'])
    (root / 'species_tree.nwk').write_text('((Taxon_A:0.1,Taxon_B:0.1):0.1,(Taxon_C:0.1,Taxon_D:0.1):0.1)root;\n')
    (root / 'expected.json').write_text(json.dumps({"scenario": scenario, "seed": seed,
        "purpose": "raw_input_integration_test", "taxa": list(SPECIES),
        "known_process": "one_added_intron_in_Taxon_D" if scenario == "splice_difference" else "no_DNA_structural_change",
        "annotation_modified": scenario == "annotation_dropout",
        "biological_accuracy_calibration": False}, indent=2) + '\n')
    return root / 'manifest.tsv'

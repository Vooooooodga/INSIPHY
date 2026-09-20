"""Export paired genomic FASTA/GFF loci with an explicit coordinate translation."""
from __future__ import annotations

from pathlib import Path
from shutil import copyfile

from ..preparation.annotation_index import read_annotation_for_gene
from ..preparation.features import FeatureHierarchy
from ..storage.fasta import fasta_record_length, read_fasta_interval
from ..storage.tabular import write_tsv
from .selection import InputSelection


def _gff_line(feature, contig: str, offset: int) -> str:
    attrs = dict(feature.get("attrs", {}))
    if feature.get("id"):
        attrs["ID"] = feature["id"]
    parents = FeatureHierarchy.parents(feature)
    if parents:
        attrs["Parent"] = ",".join(parents)
    # Values from the GFF parser already preserve GFF escaping. Do not encode
    # them a second time, or gene/transcript identifiers would silently change.
    attributes = ";".join(f"{key}={value}" for key, value in attrs.items()) or "."
    return "\t".join(map(str, (contig, feature.get("source", "provided"), feature["type"],
        int(feature["start"])-offset, int(feature["end"])-offset,
        feature.get("score", "."), feature["strand"], feature.get("phase", "."), attributes)))


def export_loci(selection: InputSelection, species_tree: str, output_dir: str,
                flank: int = 1000) -> list[dict]:
    """Export forward-genomic sequence; preserve negative-strand annotations.

    Each family has one genomic FASTA and one GFF3 per species. The exported GFF
    contains only the selected gene and its descendants. Neighbouring genes are
    not deleted from the source; they remain available in whole-genome runs.
    """
    if flank < 0:
        raise ValueError("--flank must be nonnegative")
    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    mapping = []
    selection.write(root)
    for row in selection.rows:
        directory = root / row["family_id"]
        directory.mkdir(exist_ok=True)
        features, gene, _, bounds = read_annotation_for_gene(row["annotation_file"], row["gene_id"])
        start, end = int(bounds["linked_start"]), int(bounds["linked_end"])
        length = fasta_record_length(row["genome_fasta"], gene["seqid"])
        left, right = max(1, start-flank), min(length, end+flank)
        sequence = read_fasta_interval(row["genome_fasta"], gene["seqid"], left, right)
        contig = row["species"] + "_locus"
        (directory / f"{row['species']}.fa").write_text(
            f">{contig}\n" + "\n".join(sequence[i:i+80] for i in range(0, len(sequence), 80)) + "\n")
        target_features = FeatureHierarchy(features).transcript_features(gene["id"])
        if not any(r.get("type") == "gene" for r in target_features):
            target_features.insert(0, gene)
        seen, lines = set(), []
        for feature in target_features:
            line = _gff_line(feature, contig, left-1)
            if line not in seen:
                seen.add(line)
                lines.append(line)
        (directory / f"{row['species']}.gff3").write_text(
            f"##gff-version 3\n##sequence-region {contig} 1 {len(sequence)}\n" + "\n".join(lines) + "\n")
        source_tree = Path(species_tree)
        copyfile(source_tree, directory / ("species_tree.tsv" if source_tree.suffix == ".tsv" else "species_tree.nwk"))
        upstream = start-left if gene["strand"] == "+" else right-end
        downstream = right-end if gene["strand"] == "+" else start-left
        mapping.append({"family_id": row["family_id"], "species": row["species"], "gene_id": row["gene_id"],
            "source_fasta": row["genome_fasta"], "source_gff": row["annotation_file"],
            "source_contig": gene["seqid"], "source_start": left, "source_end": right,
            "export_contig": contig, "export_start": 1, "export_end": len(sequence),
            "original_strand": gene["strand"], "sequence_orientation": "forward_genomic",
            "coordinate_rule": "source_position=export_position+source_start-1",
            "requested_flank_bp": flank, "upstream_available_bp": upstream,
            "downstream_available_bp": downstream,
            "flank_status": "truncated_at_input_boundary" if min(upstream, downstream) < flank else "complete"})
    write_tsv(root / "locus_coordinate_map.tsv", mapping, list(mapping[0]))
    return mapping

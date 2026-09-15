"""Small simulation fixtures for INSIPHY benchmarking."""

import random
from pathlib import Path

from .io import read_tsv, write_tsv


def mutate(seq, rng, rate=0.05):
    alphabet = "ACGT"
    out = []
    for base in seq:
        if rng.random() < rate:
            out.append(rng.choice([b for b in alphabet if b != base]))
        else:
            out.append(base)
    return "".join(out)


def simulate_negative_dataset(output_dir, seed=7, hidden_dropout=False):
    rng = random.Random(seed)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    species_tree = [
        {"node_id": "root", "parent_id": "", "label": "root", "branch_length": "1.0"},
        {"node_id": "sp1", "parent_id": "root", "label": "Sp1", "branch_length": "1.0"},
        {"node_id": "sp2", "parent_id": "root", "label": "Sp2", "branch_length": "1.0"},
        {"node_id": "clade34", "parent_id": "root", "label": "clade34", "branch_length": "1.0"},
        {"node_id": "sp3", "parent_id": "clade34", "label": "Sp3", "branch_length": "1.0"},
        {"node_id": "sp4", "parent_id": "clade34", "label": "Sp4", "branch_length": "1.0"},
    ]
    base = {
        "A": "ATGGCCGATGCCGATGCCGATGCCGATGCC",
        "B": "GTACCGTACCGTACCGTACCGTACCGTACC",
    }
    rows = []
    homology = []
    fasta_records = {}
    adj = []
    copy_context = []
    evidence = []
    matches = []
    family = "sim_negative"
    for species in ["Sp1", "Sp2", "Sp3", "Sp4"]:
        copy_context.append({"family_id": family, "species": species, "gene_copy_id": "copy1", "copy_class": "single_copy", "copy_subclass": "single_copy", "copy_span": f"{species}_chr1:100-232"})
        prev = None
        for idx, seg in enumerate(["A", "B"], start=1):
            occ_id = f"{species}_copy1_{seg}"
            rows.append(
                {
                    "occurrence_id": occ_id,
                    "family_id": family,
                    "species": species,
                    "gene_copy_id": "copy1",
                    "transcript_id": f"{species}_tx1",
                    "role": "CDS",
                    "role_set": "CDS",
                    "presence_status": "present",
                    "contig": f"{species}_chr1",
                    "start": 100 * idx,
                    "end": 100 * idx + len(base[seg]) - 1,
                    "strand": "+",
                    "phase": "0",
                    "source_feature_id": f"{seg}_copy1",
                    "boundary_class": "annotated_segment",
                    "splice_motif_score": "0.8",
                    "splice_donor": "GT",
                    "splice_acceptor": "AG",
                    "frame_status": "coding_frame_annotated",
                }
            )
            homology.append({"homology_id": f"H_neg_{seg}", "occurrence_id": occ_id, "support_type": "simulation_truth", "confidence": "1.0", "source_label": "ancestral_source"})
            fasta_records[occ_id] = mutate(base[seg], rng, 0.02)
            dropout = hidden_dropout and species == "Sp3" and seg == "B"
            evidence.append(
                {
                    "evidence_id": f"ev_{occ_id}",
                    "family_id": family,
                    "species": species,
                    "gene_copy_id": "copy1",
                    "homology_id": f"H_neg_{seg}",
                    "annotation_status": "sequence_completed" if dropout else "annotated",
                    "evidence_status": "supports_hidden_segment" if dropout else "supports_annotation",
                    "inferred_role": "CDS",
                    "contig": f"{species}_chr1",
                    "start": 100 * idx,
                    "end": 100 * idx + len(base[seg]) - 1,
                    "strand": "+",
                    "sequence_score": "0.95",
                    "left_synteny_score": "0.9",
                    "right_synteny_score": "0.9",
                    "splice_motif_score": "0.8",
                    "phase_compatibility": "compatible",
                    "inferred_event": "annotation_dropout_control" if dropout else "annotated_segment",
                    "frame_status": "coding_frame_preserved",
                }
            )
            if prev:
                adj.append({"adjacency_id": f"{species}_copy1_{prev}_{seg}", "family_id": family, "species": species, "gene_copy_id": "copy1", "left_occurrence_id": f"{species}_copy1_{prev}", "right_occurrence_id": occ_id, "adjacency_status": "present"})
            prev = seg
    for i, left in enumerate(homology):
        for right in homology[i + 1 :]:
            if left["homology_id"] == right["homology_id"]:
                matches.append(
                    {
                        "match_id": f"neg_match_{len(matches) + 1:04d}",
                        "query_occurrence_id": left["occurrence_id"],
                        "subject_occurrence_id": right["occurrence_id"],
                        "alignment_score": "0.95",
                        "coverage_score": "1.0",
                        "left_context_score": "0.9",
                        "right_context_score": "0.9",
                        "boundary_score": "1.0",
                        "phase_score": "1.0",
                        "order_score": "1.0",
                        "strand_score": "1.0",
                        "splice_score": "1.0",
                        "size_ratio": "1.0",
                        "total_score": "0.95",
                        "match_status": "mapped",
                    }
                )
    write_tsv(output_dir / "species_tree.tsv", species_tree, ["node_id", "parent_id", "label", "branch_length"])
    write_tsv(output_dir / "segment_occurrences.tsv", rows, ["occurrence_id", "family_id", "species", "gene_copy_id", "transcript_id", "role", "role_set", "presence_status", "contig", "start", "end", "strand", "phase", "source_feature_id", "boundary_class", "splice_motif_score", "splice_donor", "splice_acceptor", "frame_status"])
    write_tsv(output_dir / "segment_homology.tsv", homology, ["homology_id", "occurrence_id", "support_type", "confidence", "source_label"])
    write_tsv(output_dir / "physical_adjacencies.tsv", adj, ["adjacency_id", "family_id", "species", "gene_copy_id", "left_occurrence_id", "right_occurrence_id", "adjacency_status"])
    write_tsv(output_dir / "copy_context.tsv", copy_context, ["family_id", "species", "gene_copy_id", "copy_class", "copy_subclass", "copy_span"])
    write_tsv(output_dir / "copy_relationships.tsv", [], ["family_id", "species", "query_copy_id", "subject_copy_id", "relationship_class", "synteny_score", "distance_bp", "evidence"])
    write_tsv(output_dir / "sequence_synteny_evidence.tsv", evidence, ["evidence_id", "family_id", "species", "gene_copy_id", "homology_id", "annotation_status", "evidence_status", "inferred_role", "contig", "start", "end", "strand", "sequence_score", "left_synteny_score", "right_synteny_score", "splice_motif_score", "phase_compatibility", "inferred_event", "frame_status"])
    write_tsv(output_dir / "segment_matches.tsv", matches, ["match_id", "query_occurrence_id", "subject_occurrence_id", "alignment_score", "coverage_score", "left_context_score", "right_context_score", "boundary_score", "phase_score", "order_score", "strand_score", "splice_score", "size_ratio", "total_score", "match_status"])
    write_tsv(output_dir / "truth_events.tsv", [], ["family_id", "event_class", "branch_scope", "object_id", "notes"])
    write_tsv(output_dir / "simulation_truth_segments.tsv", [{"family_id": family, "homology_id": f"H_neg_{seg}", "truth_state": "conserved_present"} for seg in ["A", "B"]], ["family_id", "homology_id", "truth_state"])
    with (output_dir / "segment_sequences.fasta").open("w") as handle:
        for name, seq in sorted(fasta_records.items()):
            handle.write(f">{name}\n{seq}\n")
    return []


SCENARIO_TRUTH_FILTERS = {
    "exonization": {"exonization_candidate"},
    "source_join": {"segment_fusion_or_new_adjacency"},
    "tandem_duplication": {"copy_duplication_or_expansion"},
    "segment_split_fusion": {"segment_fusion_or_new_adjacency", "segment_split_or_adjacency_loss"},
}


def filter_truth_events(output_dir, allowed_classes):
    output_dir = Path(output_dir)
    truth = read_tsv(output_dir / "truth_events.tsv", optional=True)
    kept = [row for row in truth if row.get("event_class") in allowed_classes]
    write_tsv(output_dir / "truth_events.tsv", kept, ["family_id", "event_class", "branch_scope", "object_id", "notes"])
    return kept


def simulate_dataset(output_dir, seed=7, scenario="compound"):
    if scenario in {"negative_control", "annotation_dropout"}:
        return simulate_negative_dataset(output_dir, seed, hidden_dropout=(scenario == "annotation_dropout"))
    rng = random.Random(seed)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    species_tree = [
        {"node_id": "root", "parent_id": "", "label": "root"},
        {"node_id": "sp1", "parent_id": "root", "label": "Sp1"},
        {"node_id": "sp2", "parent_id": "root", "label": "Sp2"},
        {"node_id": "clade34", "parent_id": "root", "label": "clade34"},
        {"node_id": "sp3", "parent_id": "clade34", "label": "Sp3"},
        {"node_id": "sp4", "parent_id": "clade34", "label": "Sp4"},
    ]
    base = {
        "A": "ATGGCCGATGCCGATGCCGATGCCGATGCC",
        "B": "GTACCGTACCGTACCGTACCGTACCGTACC",
        "C": "TTTAAACCCGGGTTTAAACCCGGGTTTAAA",
    }
    rows = []
    homology = []
    fasta_records = {}
    adj = []
    copy_context = []
    copy_relationships = []
    evidence = []
    matches = []
    family = "sim_gene"
    copy_map = {"Sp1": ["copy1"], "Sp2": ["copy1"], "Sp3": ["copy1", "copy2"], "Sp4": ["copy1", "copy2"]}
    role_by_species = {"Sp1": {"A": "CDS", "B": "CDS"}, "Sp2": {"A": "CDS", "B": "CDS"}, "Sp3": {"A": "CDS", "B": "CDS", "C": "CDS"}, "Sp4": {"A": "CDS", "B": "CDS", "C": "CDS"}}
    for species, copies in copy_map.items():
        copy_class = "single_copy" if len(copies) == 1 else "tandem_multi_copy"
        for copy in copies:
            copy_context.append({"family_id": family, "species": species, "gene_copy_id": copy, "copy_class": copy_class, "copy_subclass": copy_class, "copy_span": f"{species}_chr1:100-330"})
        if len(copies) > 1:
            copy_relationships.append({"family_id": family, "species": species, "query_copy_id": copies[0], "subject_copy_id": copies[1], "relationship_class": "tandem_duplication_candidate", "synteny_score": "0.9", "distance_bp": "1000", "evidence": "simulation_truth"})
        for copy in copies:
            segments = ["A", "B"]
            if species in {"Sp3", "Sp4"} and copy == "copy1":
                segments = ["A", "C", "B"]
            prev = None
            for idx, seg in enumerate(segments, start=1):
                occ_id = f"{species}_{copy}_{seg}"
                role = role_by_species[species].get(seg, "intron")
                rows.append(
                    {
                        "occurrence_id": occ_id,
                        "family_id": family,
                        "species": species,
                        "gene_copy_id": copy,
                        "transcript_id": f"{species}_{copy}_tx1",
                        "role": role,
                        "role_set": role,
                        "presence_status": "present",
                        "contig": f"{species}_chr1",
                        "start": 100 * idx,
                        "end": 100 * idx + len(base[seg]) - 1,
                        "strand": "+",
                        "phase": "0",
                        "source_feature_id": f"{seg}_{copy}",
                        "boundary_class": "annotated_segment",
                        "splice_motif_score": "0.8",
                        "splice_donor": "GT",
                        "splice_acceptor": "AG",
                        "frame_status": "coding_frame_annotated",
                    }
                )
                source = "ancestral_source" if seg in {"A", "B"} else "intronic_source"
                homology.append({"homology_id": f"H_sim_{seg}", "occurrence_id": occ_id, "support_type": "simulation_truth", "confidence": "1.0", "source_label": source})
                fasta_records[occ_id] = mutate(base[seg], rng)
                evidence.append(
                    {
                        "evidence_id": f"ev_{occ_id}",
                        "family_id": family,
                        "species": species,
                        "gene_copy_id": copy,
                        "homology_id": f"H_sim_{seg}",
                        "annotation_status": "annotated" if seg != "C" else "sequence_completed",
                        "evidence_status": "supports_annotation" if seg != "C" else "supports_hidden_segment",
                        "inferred_role": role,
                        "contig": f"{species}_chr1",
                        "start": 100 * idx,
                        "end": 100 * idx + len(base[seg]) - 1,
                        "strand": "+",
                        "sequence_score": "0.95",
                        "left_synteny_score": "0.85",
                        "right_synteny_score": "0.85",
                        "splice_motif_score": "0.7",
                        "phase_compatibility": "compatible",
                    }
                )
                if prev:
                    adj.append({"adjacency_id": f"{species}_{copy}_{prev}_{seg}", "family_id": family, "species": species, "gene_copy_id": copy, "left_occurrence_id": f"{species}_{copy}_{prev}", "right_occurrence_id": occ_id, "adjacency_status": "present"})
                prev = seg
            if species in {"Sp1", "Sp2"} and copy == "copy1":
                occ_id = f"{species}_{copy}_C_absent"
                rows.append(
                    {
                        "occurrence_id": occ_id,
                        "family_id": family,
                        "species": species,
                        "gene_copy_id": copy,
                        "transcript_id": f"{species}_{copy}_tx1",
                        "role": "intron",
                        "role_set": "intron",
                        "presence_status": "absent",
                        "contig": f"{species}_chr1",
                        "start": 300,
                        "end": 300 + len(base["C"]) - 1,
                        "strand": "+",
                        "phase": "0",
                        "source_feature_id": "C_absent",
                        "boundary_class": "simulated_absence",
                        "splice_motif_score": "0.8",
                        "splice_donor": "GT",
                        "splice_acceptor": "AG",
                        "frame_status": "not_coding_or_unknown",
                    }
                )
                homology.append({"homology_id": "H_sim_C", "occurrence_id": occ_id, "support_type": "simulation_absence", "confidence": "1.0", "source_label": "intronic_source"})
                adj.append({"adjacency_id": f"{species}_{copy}_A_C_absent", "family_id": family, "species": species, "gene_copy_id": copy, "left_occurrence_id": f"{species}_{copy}_A", "right_occurrence_id": occ_id, "adjacency_status": "absent"})
                adj.append({"adjacency_id": f"{species}_{copy}_C_B_absent", "family_id": family, "species": species, "gene_copy_id": copy, "left_occurrence_id": occ_id, "right_occurrence_id": f"{species}_{copy}_B", "adjacency_status": "absent"})
    for i, left in enumerate(homology):
        for right in homology[i + 1 :]:
            if left["homology_id"] == right["homology_id"]:
                matches.append(
                    {
                        "match_id": f"sim_match_{len(matches) + 1:04d}",
                        "query_occurrence_id": left["occurrence_id"],
                        "subject_occurrence_id": right["occurrence_id"],
                        "alignment_score": "0.9",
                        "coverage_score": "1.0",
                        "left_context_score": "0.8",
                        "right_context_score": "0.8",
                        "boundary_score": "1.0",
                        "phase_score": "1.0",
                        "order_score": "0.9",
                        "strand_score": "1.0",
                        "splice_score": "1.0",
                        "size_ratio": "1.0",
                        "total_score": "0.91",
                        "match_status": "mapped",
                    }
                )
    truth = [
        {"family_id": family, "event_class": "exonization_candidate", "branch_scope": "root->clade34", "object_id": "H_sim_C", "notes": "simulated hidden intronic-source segment becomes CDS"},
        {"family_id": family, "event_class": "segment_fusion_or_new_adjacency", "branch_scope": "root->clade34", "object_id": "H_sim_A__H_sim_C", "notes": "simulated new internal adjacency"},
        {"family_id": family, "event_class": "copy_duplication_or_expansion", "branch_scope": "root->clade34", "object_id": family, "notes": "simulated extra copy in the derived clade"},
    ]
    write_tsv(output_dir / "species_tree.tsv", species_tree, ["node_id", "parent_id", "label"])
    write_tsv(output_dir / "segment_occurrences.tsv", rows, ["occurrence_id", "family_id", "species", "gene_copy_id", "transcript_id", "role", "role_set", "presence_status", "contig", "start", "end", "strand", "phase", "source_feature_id", "boundary_class", "splice_motif_score", "splice_donor", "splice_acceptor", "frame_status"])
    write_tsv(output_dir / "segment_homology.tsv", homology, ["homology_id", "occurrence_id", "support_type", "confidence", "source_label"])
    write_tsv(output_dir / "physical_adjacencies.tsv", adj, ["adjacency_id", "family_id", "species", "gene_copy_id", "left_occurrence_id", "right_occurrence_id", "adjacency_status"])
    write_tsv(output_dir / "copy_context.tsv", copy_context, ["family_id", "species", "gene_copy_id", "copy_class", "copy_subclass", "copy_span"])
    write_tsv(output_dir / "copy_relationships.tsv", copy_relationships, ["family_id", "species", "query_copy_id", "subject_copy_id", "relationship_class", "synteny_score", "distance_bp", "evidence"])
    write_tsv(output_dir / "sequence_synteny_evidence.tsv", evidence, ["evidence_id", "family_id", "species", "gene_copy_id", "homology_id", "annotation_status", "evidence_status", "inferred_role", "contig", "start", "end", "strand", "sequence_score", "left_synteny_score", "right_synteny_score", "splice_motif_score", "phase_compatibility"])
    write_tsv(output_dir / "segment_matches.tsv", matches, ["match_id", "query_occurrence_id", "subject_occurrence_id", "alignment_score", "coverage_score", "left_context_score", "right_context_score", "boundary_score", "phase_score", "order_score", "strand_score", "splice_score", "size_ratio", "total_score", "match_status"])
    write_tsv(output_dir / "truth_events.tsv", truth, ["family_id", "event_class", "branch_scope", "object_id", "notes"])
    write_tsv(output_dir / "simulation_truth_segments.tsv", [{"family_id": family, "homology_id": f"H_sim_{seg}", "truth_state": "derived_present" if seg == "C" else "ancestral_present"} for seg in ["A", "B", "C"]], ["family_id", "homology_id", "truth_state"])
    with (output_dir / "segment_sequences.fasta").open("w") as handle:
        for name, seq in sorted(fasta_records.items()):
            handle.write(f">{name}\n{seq}\n")
    if scenario in SCENARIO_TRUTH_FILTERS:
        return filter_truth_events(output_dir, SCENARIO_TRUTH_FILTERS[scenario])
    return truth

"""Real-case preparation utilities for INSIPHY."""

from collections import Counter, defaultdict
from pathlib import Path

from .alignment import local_alignment_stats, revcomp as reverse_complement, splice_motif_score
from .io import parse_fasta, read_tsv, write_tsv
from .preprocess import extract_gene, derive_tables, read_annotation, translate_cds


def split_aliases(value):
    tokens = []
    for part in str(value or "").replace(",", ";").split(";"):
        part = part.strip()
        if part:
            tokens.append(part)
    return tokens


def infer_manifest_source_label(row):
    explicit = row.get("source_label", "")
    if explicit and explicit not in {"NA", "unknown", "unknown_source"}:
        return explicit
    text = " ".join([row.get("role_hint", ""), row.get("gene_symbol", ""), row.get("gene_copy_id", "")]).lower()
    if "derived" in text or "sdic" in text or "jgw" in text or "jingwei" in text:
        return "unknown_source"
    if "anxb10" in text or "annexin" in text:
        return "AnxB10"
    if "short_wing" in text or text.endswith(" sw") or "_sw" in text or " sw_" in text:
        return "sw"
    if "ymp" in text or "yande" in text or "yellow_emperor" in text or "yellow emperor" in text:
        return "ymp"
    if "adh" in text:
        return "Adh"
    return "unknown_source"


def infer_manifest_copy_role(row):
    explicit = row.get("copy_role", "")
    if explicit and explicit != "NA":
        return explicit
    text = " ".join([row.get("role_hint", ""), row.get("gene_symbol", ""), row.get("gene_copy_id", "")]).lower()
    if "derived" in text or "sdic" in text or "jgw" in text or "jingwei" in text:
        return "derived"
    if "source" in text:
        return "source"
    if "background" in text or "ortholog" in text:
        return "background"
    return "candidate"


def feature_tokens(feature):
    attrs = feature.get("attrs", {})
    values = [feature.get("id", ""), feature.get("name", ""), feature.get("parent", "")]
    for key in ["ID", "Name", "Alias", "gene_id", "gene_name", "transcript_id", "Parent", "Dbxref", "description"]:
        values.extend(split_aliases(attrs.get(key, "")))
    tokens = set()
    for value in values:
        for token in split_aliases(value):
            tokens.add(token)
            if ":" in token:
                tokens.add(token.split(":")[-1])
    return {token for token in tokens if token}


def collect_queries(queries=None, alias_file=None):
    out = []
    for query in queries or []:
        out.extend(split_aliases(query))
    if alias_file:
        for row in read_tsv(alias_file, optional=True):
            values = []
            for field in ["query", "gene_id", "gene_symbol", "target_id", "aliases"]:
                values.extend(split_aliases(row.get(field, "")))
            out.extend(values)
    seen = set()
    unique = []
    for query in out:
        key = query.lower()
        if key not in seen:
            seen.add(key)
            unique.append(query)
    return unique


def child_counts(features, gene):
    gene_ids = feature_tokens(gene)
    counts = Counter()
    for feature in features:
        if feature is gene:
            continue
        parent_tokens = set(split_aliases(feature.get("parent", "")))
        overlaps = feature["seqid"] == gene["seqid"] and feature["start"] <= gene["end"] and feature["end"] >= gene["start"]
        if parent_tokens & gene_ids or (overlaps and feature["type"].lower() in {"mrna", "transcript", "exon", "cds"}):
            counts[feature["type"].lower()] += 1
    return counts


def inspect_annotation(annotation, output_dir, queries=None, alias_file=None, species="NA", case_id="case"):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    features = read_annotation(annotation)
    query_list = collect_queries(queries, alias_file)
    rows = []
    for query in query_list:
        q = query.lower()
        for feature in features:
            tokens = feature_tokens(feature)
            token_lut = {token.lower(): token for token in tokens}
            text = " ".join(sorted(tokens)).lower()
            if q not in token_lut and q not in text:
                continue
            counts = child_counts(features, feature)
            rows.append(
                {
                    "case_id": case_id,
                    "species": species,
                    "query": query,
                    "matched_value": token_lut.get(q, query),
                    "feature_id": feature.get("id", "NA") or "NA",
                    "feature_name": feature.get("name", "NA") or "NA",
                    "feature_type": feature.get("type", "NA"),
                    "contig": feature.get("seqid", "NA"),
                    "start": feature.get("start", "NA"),
                    "end": feature.get("end", "NA"),
                    "strand": feature.get("strand", "NA"),
                    "parent": feature.get("parent", "NA") or "NA",
                    "transcript_count": counts.get("mrna", 0) + counts.get("transcript", 0),
                    "exon_count": counts.get("exon", 0),
                    "cds_count": counts.get("cds", 0),
                    "child_feature_count": sum(counts.values()),
                    "match_rank": "exact_token" if q in token_lut else "text_contains",
                }
            )
    fields = [
        "case_id",
        "species",
        "query",
        "matched_value",
        "feature_id",
        "feature_name",
        "feature_type",
        "contig",
        "start",
        "end",
        "strand",
        "parent",
        "transcript_count",
        "exon_count",
        "cds_count",
        "child_feature_count",
        "match_rank",
    ]
    write_tsv(output_dir / "gene_candidate_report.tsv", rows, fields)
    return rows


def write_provenance(manifest_rows, output_dir):
    fields = ["case_id", "species", "assembly", "annotation", "source_url", "release", "genome_fasta", "annotation_file", "role_hint", "source_label", "copy_role", "notes"]
    rows = []
    for row in manifest_rows:
        rows.append({field: row.get(field, "NA") for field in fields})
    write_tsv(Path(output_dir) / "case_provenance.tsv", rows, fields)


def copy_optional_tree(tree_path, output_dir, output_name):
    if not tree_path:
        return
    path = Path(tree_path)
    if path.exists():
        (Path(output_dir) / output_name).write_text(path.read_text())


def build_case(
    manifest,
    output_dir,
    identity_threshold=0.7,
    species_tree=None,
    transcript_policy="canonical",
    canonical_rule="longest_cds",
    aligner="auto",
    threads=1,
    min_size_ratio=0.25,
    copy_tree=None,
    gene_tree=None,
):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = read_tsv(manifest, ["case_id", "species", "family_id", "gene_id", "gene_copy_id", "genome_fasta", "annotation_file"])
    write_provenance(rows, output_dir)
    report = []
    appended = False
    for row in rows:
        required_paths = [row.get("genome_fasta", ""), row.get("annotation_file", "")]
        missing_path = [path for path in required_paths if not path or path == "TBD" or not Path(path).exists()]
        if missing_path or row.get("gene_id", "") in {"", "TBD"}:
            message = ";".join(missing_path) or "missing_gene_id"
            raise SystemExit(f"manifest row has incomplete required input for {row.get('gene_copy_id', row.get('gene_id', 'unknown'))}: {message}")
        before = len(read_tsv(output_dir / "segment_occurrences.tsv", optional=True))
        extract_gene(
            row["genome_fasta"],
            row["annotation_file"],
            row["gene_id"],
            row["family_id"],
            row["species"],
            row["gene_copy_id"],
            output_dir,
            append=appended,
            transcript_policy=transcript_policy,
            canonical_rule=canonical_rule,
            source_label=infer_manifest_source_label(row),
            copy_role=infer_manifest_copy_role(row),
        )
        appended = True
        after = len(read_tsv(output_dir / "segment_occurrences.tsv", optional=True))
        report.append(
            {
                "case_id": row.get("case_id", "NA"),
                "species": row.get("species", "NA"),
                "gene_id": row.get("gene_id", "NA"),
                "gene_copy_id": row.get("gene_copy_id", "NA"),
                "status": "extracted",
                "segment_count": after - before,
                "message": "ok",
            }
        )
    copy_optional_tree(species_tree, output_dir, "species_tree.tsv")
    copy_optional_tree(copy_tree, output_dir, "copy_tree.tsv")
    copy_optional_tree(gene_tree, output_dir, "gene_tree.tsv")
    if appended:
        derive_tables(output_dir, output_dir, identity_threshold, aligner=aligner, threads=threads, min_size_ratio=min_size_ratio)
    write_tsv(output_dir / "case_build_report.tsv", report, ["case_id", "species", "gene_id", "gene_copy_id", "status", "segment_count", "message"])
    return report


def revcomp(seq):
    return reverse_complement(seq)


def best_ungapped_hit(query, target):
    query = query.upper()
    target = target.upper()
    if not query or not target:
        return {"start": 0, "end": 0, "identity": 0.0, "coverage": 0.0}
    if len(target) < len(query):
        span = len(target)
        matches = sum(1 for a, b in zip(query[:span], target) if a == b and a != "N" and b != "N")
        return {"start": 1, "end": span, "identity": matches / max(1, span), "coverage": span / len(query)}
    best = {"start": 1, "end": len(query), "identity": -1.0, "coverage": 1.0}
    qlen = len(query)
    for offset in range(0, len(target) - qlen + 1):
        window = target[offset : offset + qlen]
        matches = sum(1 for a, b in zip(query, window) if a == b and a != "N" and b != "N")
        identity = matches / qlen
        if identity > best["identity"]:
            best = {"start": offset + 1, "end": offset + qlen, "identity": identity, "coverage": 1.0}
    return best


def scan_hidden_segments(source_fasta, target_fasta, output_dir, family_id="NA", species="NA", gene_copy_id="NA", min_identity=0.75, min_coverage=0.5, aligner="internal", threads=1):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    sources = parse_fasta(source_fasta)
    targets = parse_fasta(target_fasta)
    rows = []
    for query_id, query_seq in sorted(sources.items()):
        if aligner == "miniprot" and set(query_seq.upper()) <= set("ACGTUN-"):
            query_seq = translate_cds(query_seq)
        best = None
        for target_id, target_seq in sorted(targets.items()):
            orientations = [("+", query_seq)] if aligner == "miniprot" else [("+", query_seq), ("-", revcomp(query_seq))]
            for strand, qseq in orientations:
                hit = local_alignment_stats(qseq, target_seq, backend=aligner, threads=threads)
                target_fragment = target_seq[max(0, hit.target_start - 1) : hit.target_end]
                motif, donor, acceptor = splice_motif_score(target_fragment)
                frame_status = "coding_frame_preserved" if (hit.query_end - hit.query_start + 1) % 3 == 0 else "frameshift_or_stop_risk"
                row = {
                    "family_id": family_id,
                    "species": species,
                    "gene_copy_id": gene_copy_id,
                    "query_id": query_id,
                    "target_id": target_id,
                    "start": hit.target_start,
                    "end": hit.target_end,
                    "strand": strand,
                    "identity": hit.identity,
                    "coverage": hit.query_coverage or hit.coverage,
                    "alignment_score": hit.score,
                    "alignment_cigar": hit.cigar,
                    "alignment_backend": hit.backend,
                    "splice_motif_score": motif,
                    "splice_donor": donor,
                    "splice_acceptor": acceptor,
                    "frame_status": frame_status,
                }
                if best is None or (row["identity"], row["coverage"]) > (best["identity"], best["coverage"]):
                    best = row
        if best is None:
            continue
        best["support_call"] = "hidden_segment_candidate" if best["identity"] >= min_identity and best["coverage"] >= min_coverage else "low_support"
        best["identity"] = f"{best['identity']:.6g}"
        best["coverage"] = f"{best['coverage']:.6g}"
        best["alignment_score"] = f"{best['alignment_score']:.6g}"
        best["splice_motif_score"] = f"{best['splice_motif_score']:.6g}"
        rows.append(best)
    write_tsv(
        output_dir / "hidden_segment_scan.tsv",
        rows,
        ["family_id", "species", "gene_copy_id", "query_id", "target_id", "start", "end", "strand", "identity", "coverage", "alignment_score", "alignment_cigar", "alignment_backend", "splice_motif_score", "splice_donor", "splice_acceptor", "frame_status", "support_call"],
    )
    return rows

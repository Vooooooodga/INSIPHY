"""preparation / annotation_index: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from functools import lru_cache
from bisect import bisect_left, bisect_right
from collections import defaultdict
from .features import FeatureHierarchy
from intraphy.storage.tabular import open_text
from pathlib import Path


def parse_attributes(raw):
    attrs = {}
    for part in raw.strip().strip(";").split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
        elif " " in part:
            key, value = part.split(" ", 1)
            value = value.strip().strip('"')
        else:
            continue
        attrs[key.strip()] = value.strip().strip('"')
    return attrs


def split_ids(value):
    out = []
    for part in str(value or "").replace(";", ",").split(","):
        part = part.strip().strip('"')
        if part:
            out.append(part)
    return out


def iter_annotation(path):
    with open_text(path) as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            attrs = parse_attributes(parts[8])
            ftype = parts[2].lower()
            transcript_types = {"mrna", "transcript", "lnc_rna", "ncrna", "rrna", "trna"}
            if ftype == "gene":
                feature_id = attrs.get("ID") or attrs.get("gene_id") or attrs.get("Name") or ""
                parent = attrs.get("Parent") or ""
            elif ftype in transcript_types:
                feature_id = attrs.get("ID") or attrs.get("transcript_id") or attrs.get("Name") or ""
                parent = attrs.get("Parent") or attrs.get("gene_id") or ""
            else:
                feature_id = attrs.get("ID") or attrs.get("exon_id") or attrs.get("protein_id") or attrs.get("transcript_id") or attrs.get("gene_id") or ""
                parent = attrs.get("Parent") or attrs.get("transcript_id") or attrs.get("gene_id") or ""
            yield {
                    "seqid": parts[0],
                    "source": parts[1],
                    "type": parts[2],
                    "start": int(parts[3]),
                    "end": int(parts[4]),
                    "score": parts[5],
                    "strand": parts[6],
                    "phase": parts[7],
                    "attrs": attrs,
                    "id": feature_id,
                    "parent": parent,
                    "parents": split_ids(parent),
                    "name": attrs.get("Name") or attrs.get("gene_name") or "",
                }


def read_annotation(path):
    return list(iter_annotation(path))


def _annotation_row_key(row):
    return (
        row.get("seqid", ""),
        row.get("type", ""),
        int(row.get("start", 0)),
        int(row.get("end", 0)),
        row.get("strand", "."),
        row.get("id", ""),
        row.get("parent", ""),
    )


def _copy_annotation_row(row):
    copied = dict(row)
    copied["attrs"] = dict(row.get("attrs", {}))
    copied["parents"] = list(row.get("parents", []))
    return copied


class AnnotationIndex(FeatureHierarchy):
    """One parsed GFF/GTF resource with ID, Parent, token and genomic indexes.

    Spatial queries return source-order records, retain antisense overlaps, and
    use the parser's 1-based closed coordinates. No presence/role is inferred.
    """

    def __init__(self, rows):
        super().__init__(rows)
        self.by_token = defaultdict(list)
        grouped = defaultdict(list)
        for ordinal, row in enumerate(self.rows):
            for token in feature_tokens(row):
                self.by_token[token].append(row)
            grouped[row["seqid"]].append((row["start"], row["end"], ordinal))
        self._spatial = {}
        for seqid, values in grouped.items():
            values.sort()
            starts, max_ends, end_so_far = [], [], -1
            for start, end, _ in values:
                starts.append(start)
                end_so_far = max(end_so_far, end)
                max_ends.append(end_so_far)
            self._spatial[seqid] = (values, starts, max_ends)

    def overlap(self, seqid, start, end):
        if start > end:
            return []
        values, starts, max_ends = self._spatial.get(seqid, ((), (), ()))
        upper = bisect_right(starts, end)
        lower = bisect_left(max_ends, start, 0, upper)
        ordinals = sorted(ordinal for _, stop, ordinal in values[lower:upper]
                          if stop >= start)
        return [self.rows[ordinal] for ordinal in ordinals]

    def for_gene(self, gene_id):
        matched = self.by_token.get(gene_id, ())
        genes = [row for row in matched if row["type"].lower() == "gene"]
        exact = [row for row in genes if gene_id in {
            row.get("id", ""), row.get("attrs", {}).get("ID", ""),
            row.get("attrs", {}).get("gene_id", ""),
        }]
        if len(exact) > 1:
            raise SystemExit(f"gene_id maps to multiple exact gene loci: {gene_id}")
        if not exact and len(genes) > 1:
            loci = ",".join(f"{r['seqid']}:{r['start']}-{r['end']}:{r.get('id')}" for r in genes)
            raise SystemExit(f"gene_id ambiguously matches multiple gene loci: {gene_id} ({loci})")
        if exact or genes:
            gene = (exact or genes)[0]
            matching_children = []
        else:
            matching_children = [row for row in matched if row["type"].lower() != "gene"]
            if not matching_children:
                raise SystemExit(f"gene_id not found in annotation: {gene_id}")
            gene, _ = locate_gene(matching_children, gene_id)
        gene_ids = {v for v in (gene.get("id"), gene.get("attrs", {}).get("ID"),
                               gene.get("attrs", {}).get("gene_id")) if v}
        linked = {_annotation_row_key(gene): gene}
        link_ids = set(gene_ids)
        for row in matching_children:
            if row["seqid"] == gene["seqid"]:
                linked[_annotation_row_key(row)] = row
                if row.get("id"):
                    link_ids.add(row["id"])
        pending = sorted(link_ids)
        visited = set()
        while pending:
            identifier = pending.pop()
            if identifier in visited:
                continue
            visited.add(identifier)
            for row in self.children.get(identifier, ()):
                if row["seqid"] != gene["seqid"] or row["type"].lower() == "gene":
                    continue
                linked[_annotation_row_key(row)] = row
                if row.get("id") and row["id"] not in link_ids:
                    link_ids.add(row["id"])
                    pending.append(row["id"])
        linked_start = min(row["start"] for row in linked.values())
        linked_end = max(row["end"] for row in linked.values())
        region_rows = self.overlap(gene["seqid"], linked_start, linked_end)
        deduped = {_annotation_row_key(row): row for row in region_rows}
        for row in linked.values():
            deduped[_annotation_row_key(row)] = row
        bounds = dict(annotation_start=gene["start"], annotation_end=gene["end"],
                      linked_start=linked_start, linked_end=linked_end)
        return tuple(deduped.values()), gene, frozenset(gene_ids | link_ids), bounds


@lru_cache(maxsize=2)
def _load_annotation_index_cached(path, modified_ns, size):
    return AnnotationIndex(iter_annotation(path))


def load_annotation_index(path):
    """Bounded process-local cache. Resource changes invalidate by stat, not hashes.

    Workflows should treat resource files as immutable for the duration of a run.
    Call clear_annotation_cache after replacing a resource while preserving stat.
    """
    resource = Path(path).resolve()
    stat = resource.stat()
    return _load_annotation_index_cached(str(resource), stat.st_mtime_ns, stat.st_size)


def clear_annotation_cache():
    _load_annotation_index_cached.cache_clear()


def _read_annotation_for_gene_cached(path, gene_id):
    """Compatibility lookup; the parsed resource rather than each gene is cached."""
    return load_annotation_index(path).for_gene(gene_id)


def read_annotation_for_gene(path, gene_id):
    rows, gene, gene_ids, bounds = _read_annotation_for_gene_cached(str(Path(path)), gene_id)
    return (
        [_copy_annotation_row(row) for row in rows],
        _copy_annotation_row(gene),
        set(gene_ids),
        dict(bounds),
    )


def feature_tokens(feature):
    attrs = feature.get("attrs", {})
    values = [feature.get("id", ""), feature.get("name", ""), feature.get("parent", "")]
    for key in ["ID", "Name", "Alias", "gene_id", "gene_name", "transcript_id", "Parent", "Dbxref"]:
        values.extend(split_ids(attrs.get(key, "")))
    tokens = set()
    for value in values:
        for token in split_ids(value):
            tokens.add(token)
            if ":" in token:
                tokens.add(token.split(":")[-1])
    return {token for token in tokens if token}


def feature_matches_gene(feature, gene_ids):
    return bool(feature_tokens(feature) & set(gene_ids))


def locate_gene(features, gene_id):
    gene_ids = {gene_id}
    genes = [row for row in features if row["type"].lower() == "gene" and feature_matches_gene(row, gene_ids)]
    if genes:
        exact_genes = [
            row for row in genes
            if gene_id in {row.get("id", ""), row.get("attrs", {}).get("ID", ""), row.get("attrs", {}).get("gene_id", "")}
        ]
        if len(exact_genes) > 1:
            raise SystemExit(f"gene_id maps to multiple exact gene loci: {gene_id}")
        if len(exact_genes) == 1:
            gene = exact_genes[0]
            return gene, gene_ids | feature_tokens(gene)
        if len(genes) > 1:
            raise SystemExit(f"gene_id ambiguously matches multiple gene loci: {gene_id}")
        gene = genes[0]
        return gene, gene_ids | feature_tokens(gene)
    children = [row for row in features if feature_matches_gene(row, gene_ids)]
    if not children:
        raise SystemExit(f"gene_id not found in annotation: {gene_id}")
    seqids = {row["seqid"] for row in children}
    if len(seqids) != 1:
        raise SystemExit(f"gene_id maps to multiple contigs: {gene_id}")
    strand = children[0]["strand"]
    gene = {
        "seqid": children[0]["seqid"],
        "type": "gene",
        "start": min(row["start"] for row in children),
        "end": max(row["end"] for row in children),
        "strand": strand,
        "phase": ".",
        "attrs": {"ID": gene_id, "Name": gene_id},
        "id": gene_id,
        "parent": "",
        "parents": [],
        "name": gene_id,
    }
    return gene, gene_ids


def overlaps_gene(feature, gene):
    return feature["seqid"] == gene["seqid"] and feature["start"] <= gene["end"] and feature["end"] >= gene["start"]

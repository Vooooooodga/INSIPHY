"""Optional CESAR2 gene-mode evidence export; never silently revises observations.

Input syntax follows hillerlab/CESAR2.0 README (direct gene mode). Profiles and
codon matrix are explicit; human splice priors are never assumed for other taxa.
No third-party source is copied or redistributed by this adapter.
"""
from __future__ import annotations
from pathlib import Path
import shutil
from ..structure.native import load_native
from ..structure.serialization import write_json
from .recorded import run_recorded


def gene_mode_input(reference, query, transcript, profiles: Path) -> str:
    if reference.translation_exceptions:
        raise ValueError("CESAR adapter does not reinterpret translation exceptions")
    if reference.partial:
        raise ValueError("CESAR full-gene profile requires a complete reference CDS; native structure analysis remains available")
    if transcript not in reference.paths:
        raise ValueError("Reference transcript is absent from the prepared annotation")
    by_id = {e.id: e for e in reference.exons}
    parts = reference.coding_by_transcript.get(transcript, ())
    if not parts or parts[0][2] != 0:
        raise ValueError("CESAR gene mode requires an explicitly phased complete reference CDS")
    chunks = []
    for identifier in reference.paths[transcript]:
        exon = by_id[identifier]
        exonic = [p for p in parts if exon.start <= p[0] < p[1] <= exon.end]
        if not exonic:
            continue  # UTR-only exons remain in native records; not coding input.
        sequence = ""
        for a, b, phase in exonic:
            left, right = (a-reference.search_start+1, b-reference.search_start+1) if reference.strand == "+" else (reference.search_end-b, reference.search_end-a)
            sequence += reference.sequence[left:right]
        chunks.append(sequence)
    if not chunks or not chunks[0].startswith("ATG") or chunks[-1][-3:] not in {"TAA", "TAG", "TGA"} or sum(map(len, chunks)) % 3:
        raise ValueError("CESAR adapter needs a complete standard-code reference; this is not an evolutionary filter")
    names = ("firstCodon_profile.txt", "lastCodon_profile.txt", "acc_profile.txt", "do_profile.txt")
    for name in names:
        if not (profiles/name).is_file():
            raise FileNotFoundError(f"Explicit CESAR profile missing: {profiles/name}")
    records, completed = [], 0
    for i, sequence in enumerate(chunks):
        leading = (3-completed % 3) % 3 if i else 0
        trailing = (len(sequence)-leading) % 3
        body_end = len(sequence)-trailing
        masked = sequence[:leading].lower()+sequence[leading:body_end]+sequence[body_end:].lower()
        acceptor = profiles/(names[0] if i == 0 else names[2])
        donor = profiles/(names[1] if i == len(chunks)-1 else names[3])
        records.append(f">reference_exon_{i+1}\t{acceptor.resolve()}\t{donor.resolve()}\n{masked}\n")
        completed += len(sequence)
    return "".join(records)+"####\n>query_genomic_locus\n"+query.sequence+"\n"


def realign(args):
    loci = [l for l in load_native(args.input_dir) if l.family == args.family_id]
    def select(species):
        matches = [l for l in loci if l.species == species]
        if len(matches) != 1:
            raise ValueError(f"Expected one selected locus for {species}")
        return matches[0]
    reference, query = select(args.reference_species), select(args.query_species)
    if reference.species == query.species:
        raise ValueError("Reference and query must be distinct species")
    transcript = args.reference_transcript
    if not transcript and len(reference.paths) == 1:
        transcript = next(iter(reference.paths))
    if not transcript:
        raise ValueError("Multiple reference annotations: specify --reference-transcript")
    out = Path(args.output_dir)
    content = gene_mode_input(reference, query, transcript, Path(args.profile_dir))
    source = out/"cesar_gene_mode.fa"
    source.write_text(content)
    executable = shutil.which(args.cesar)
    if not executable:
        raise FileNotFoundError("CESAR executable unavailable; no internal replacement was used")
    matrix = Path(args.codon_matrix).resolve()
    if not matrix.is_file():
        raise FileNotFoundError("The explicit CESAR codon matrix is missing")
    result = run_recorded([executable, str(source.resolve()), "--matrix", str(matrix)], out, "cesar", args.timeout)
    if not result.read_text().strip():
        raise ValueError("CESAR returned empty evidence")
    write_json(out/"cesar_evidence.json", {"status": "prediction_only", "reference_species": reference.species,
        "query_species": query.species, "reference_transcript": transcript, "profile_dir": str(Path(args.profile_dir).resolve()),
        "codon_matrix": str(matrix), "input": source.name, "raw_prediction": result.name,
        "updates_original_annotation": False, "enters_primary_history_automatically": False,
        "scope": "explicit_optional_realigner_evidence; not independent transcript validation"})

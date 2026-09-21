"""Small explicit structural examples evaluated by the production model."""
from __future__ import annotations
from pathlib import Path
import json
from ..structure.types import (Catalogue, ExonSpan as E, ExonConfiguration as C,
    Material, ObservationEvidence as O, ConfigurationAlternative)
from ..structure.edits import EDIT_KINDS
from ..structure.serialization import write_catalogues, write_json
from ..storage.tabular import write_tsv
from ..inference.configuration_run import infer_configurations


def build_examples(output_dir):
    out = Path(output_dir); raw = out/"example_inputs"; raw.mkdir(parents=True, exist_ok=True)
    full, split, shifted = C((E(0, 180),)), C((E(0, 60), E(90, 180))), C((E(0, 168),))
    catalogues = []
    for number, (name, base, changed, taxon) in enumerate((
        ("split", full, split, "D"), ("fusion", split, full, "B"), ("boundary", full, shifted, "C"))):
        spans = tuple(sorted(set((*base.exons, *changed.exons))))
        obs = tuple(O(s, (changed if s == taxon else base,)) for s in "ABCD")
        junctions = ((60, 90),) if name != "boundary" else ()
        catalogues.append(Catalogue("Example_gene", name, 180, spans, junctions,
            observations=obs, boundary_candidates=spans, alignment_offset=number*250))
    obs = tuple(O(s, (C((), (0,)) if s == "D" else C(full.exons, (1,)),),
                  material_presence=(0 if s == "D" else 1,)) for s in "ABCD")
    catalogues.append(Catalogue("Example_gene", "deletion", 180, full.exons, (),
        (Material("source_exon", 0, 180),), obs, boundary_candidates=full.exons, alignment_offset=750))
    alternative = ConfigurationAlternative(full, split.key, "A", "A_transcript", "synthetic_identical_sequence", 1.)
    obs = tuple(O(s, (split if s == "D" else full,), "partial" if s == "D" else "observed",
                  alternatives=(alternative,) if s == "D" else ()) for s in "ABCD")
    catalogues.append(Catalogue("Annotation_control", "uncertainty", 180,
        (E(0, 180), E(0, 60), E(90, 180)), ((60, 90),), observations=obs,
        boundary_candidates=(E(0, 180), E(0, 60), E(90, 180))))
    write_catalogues(raw/"catalogue.jsonl", catalogues)
    rows = [{"node_id": "root", "parent_id": "", "label": "root", "branch_length": "NA"},
        {"node_id": "AB", "parent_id": "root", "label": "AB", "branch_length": .5},
        {"node_id": "CD", "parent_id": "root", "label": "CD", "branch_length": .5}]
    rows += [{"node_id": s, "parent_id": "AB" if s in "AB" else "CD", "label": s, "branch_length": .5} for s in "ABCD"]
    write_tsv(raw/"species_tree.tsv", rows, ["node_id", "parent_id", "label", "branch_length"])
    rates = {k: .7 if k in {"split", "fusion"} else .04 for k in EDIT_KINDS}
    write_json(raw/"rates.json", {"schema": "intraphy.exon-rates/1", "rates": rates,
        "provenance": "Illustrative fixed parameters for software-generated teaching figures; not fitted to biology"})
    summary = infer_configurations(raw, out/"example_results", model="exon-ctmc",
        configurations=raw/"catalogue.jsonl", rates=raw/"rates.json", expected_edits=True,
        origin_root_sensitivity=(.25, 4.))
    data = json.loads((out/"example_results/exon_history.json").read_text())
    data["synthetic"] = True
    data["fixed_teaching_rates"] = rates
    data["summary"] = summary
    data["biological_accuracy_validation"] = False
    write_json(out/"illustrative_model.json", data)
    return data

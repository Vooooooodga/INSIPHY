"""Read-only exon structures with distinct parsimony and CTMC result figures."""
from __future__ import annotations
from html import escape
import json
from pathlib import Path
import re
import shutil
import numpy as np
from ..run_result import result_model
from ..structure import MODEL_VERSION
from ..topology import SpeciesTree
from ..storage.tabular import read_tsv
from ..storage.fasta import parse_fasta
from .exon_drawing import text, line, rect, exon_track, paragraphs, document, TEAL, GREY


def _alternatives(observation):
    unique = {}
    for item in observation.get("alternatives", ()):
        config = item["configuration"]
        key = tuple((e["start"], e["end"]) for e in config["exons"])
        unique[key] = config
    native = {tuple((e["start"], e["end"]) for e in c["exons"])
              for c in observation.get("configurations", ())}
    return [unique[k] for k in sorted(unique) if k not in native]


def _positions(tree, observed):
    ys, depth, cursor = {}, {tree.root: 0}, 175
    # Preorder gives an actual planar embedding, even if the input TSV is shuffled.
    for leaf in (n for n in tree.preorder() if n in tree.leaves):
        o = observed.get(tree.label[leaf], {})
        ys[leaf] = cursor
        lanes = max(1, len(o.get("configurations", ())))+min(3, len(_alternatives(o)))
        cursor += max(100, lanes*29+48)
    for n in tree.postorder():
        if tree.children.get(n): ys[n] = sum(ys[c] for c in tree.children[n])/len(tree.children[n])
    for n in tree.preorder():
        for c in tree.children.get(n, ()): depth[c] = depth[n]+1
    xs = {n: 155+130*depth[n] for n in depth}
    return xs, ys, cursor


def _posterior(unit, tree):
    nodes = unit.get("ctmc", {}).get("nodes", {})
    for node, values in nodes.items():
        a = np.asarray(values, float)
        if (node not in tree.parent or a.shape != (len(unit["states"]),) or
                not np.isfinite(a).all() or (a < -1e-10).any() or abs(a.sum()-1) > 1e-6):
            raise ValueError("Invalid saved ancestral probabilities; renderer will not normalize them")
    for branch in unit.get("ctmc", {}).get("branches", ()):
        if (branch.get("parent"), branch.get("child")) not in tree.edges():
            raise ValueError("Saved branch is absent from the result tree")
        for key in ("probability_different_endpoints", "probability_at_least_one_edit"):
            value = branch.get(key)
            if value is not None and (not np.isfinite(value) or not -1e-10 <= value <= 1+1e-10):
                raise ValueError("Invalid saved branch probability")
        for key, value in branch.items():
            if key.startswith("expected_") and (not np.isfinite(value) or value < -1e-10):
                raise ValueError("Invalid saved expected edit count")
    return nodes


def _ancestor_glyph(unit, node, x, y, mode, observation_view, nodes):
    states = unit.get("states", ())
    if not states:
        return []
    if mode == "ctmc":
        if node not in nodes:
            return [text(x, y-20, "probability unavailable", 10)]
        values = np.asarray(nodes[node], float)
        index = int(np.argmax(values)); probability = float(values[index])
        caption = f'{states[index]["state_id"]}  P={probability:.3f}'
        attrs = f'data-node="{escape(node, quote=True)}" data-state-index="{index}" data-probability="{probability:.17g}"'
        body = [f'<g {attrs}>', rect(x-128, y-64, 116, 57, "white")]
        body += exon_track(x-122, y-34, 84, unit["catalogue"]["length"], states[index]["exons"], height=11)
        body += [text(x-122, y-49, caption, 10), rect(x-122, y-18, 84, 4, GREY),
                 rect(x-122, y-18, 84*probability, 4, TEAL), '</g>']
        return body
    possibilities = set()
    for h in unit.get("views", {}).get(observation_view, {}).get("histories", ()):
        possibilities.update(h.get("nodes", {}).get(node, ()))
    if not possibilities:
        return [text(x-122, y-18, "unresolved", 10)]
    index = min(possibilities)
    body = [rect(x-128, y-61, 116, 51, "white")]
    body += exon_track(x-122, y-29, 84, unit["catalogue"]["length"], states[index]["exons"],
                      dash="4,2" if len(possibilities) > 1 else "", height=11)
    body.append(text(x-122, y-45, f'{states[index]["state_id"]}'+
        (f" / {len(possibilities)} optimal states" if len(possibilities) > 1 else ""), 10))
    return body


def _ribbons(c, observed, tree, ys, alignment, left, width):
    elements = []
    labels = [tree.label[n] for n in tree.preorder() if n in tree.leaves]
    for sa, sb in zip(labels, labels[1:]):
        oa, ob = observed.get(sa, {}), observed.get(sb, {})
        ca, cb = oa.get("configurations", ()), ob.get("configurations", ())
        if (sa not in alignment or sb not in alignment or len(ca) != 1 or len(cb) != 1 or
            any(o.get("kind") not in {"observed", "partial"} for o in (oa, ob))):
            continue
        offset = c["alignment_offset"]
        for ea in ca[0]["exons"]:
            for eb in cb[0]["exons"]:
                a, b = max(ea["start"], eb["start"]), min(ea["end"], eb["end"])
                valid = "".join("x" if alignment[sa][offset+i] in "ACGT" and alignment[sb][offset+i] in "ACGT"
                                else "." for i in range(a, b))
                for run in re.finditer("x+", valid):
                    l, r = left+width*(a+run.start())/c["length"], left+width*(a+run.end())/c["length"]
                    ya, yb = ys[tree.leaf_by_label[sa]]+9, ys[tree.leaf_by_label[sb]]-9
                    elements.append(f'<path d="M{l},{ya} L{r},{ya} L{r},{yb} L{l},{yb} Z" '
                                    'fill="#AAC7CF" opacity="0.25" data-role="paired-base-ribbon"/>')
    return elements


def draw_unit(unit, tree, events, alignment=None, path=None, *, mode=None, observation_view="evidence"):
    """All probabilities and edit IDs are read from saved results, never refitted."""
    mode = mode or ("ctmc" if unit.get("ctmc", {}).get("nodes") else "parsimony")
    if mode not in {"ctmc", "parsimony"}: raise ValueError("Unknown figure inference mode")
    c = unit["catalogue"]
    observed = {o["species"]: o for o in c["observations"]}
    xs, ys, bottom = _positions(tree, observed)
    nodes = _posterior(unit, tree) if mode == "ctmc" else {}
    title = f'{c["family"]} / {c["unit"]} — '+("conditional CTMC" if mode == "ctmc" else "maximum parsimony")
    elements = [text(30, 32, title, 22),
        text(30, 58, f"Observation conditions: {observation_view}. Each box is an exon, not an independent alignment fragment.", 14),
        text(30, 81, "Gene structures use homologous coordinates; ribbons show paired bases only. Branch lengths are drawn schematically.", 13),
        text(30, 105, "Most probable local ancestor shown; bar = its probability, grey = all other states. No joint ancestor is implied." if mode == "ctmc" else
             "Internal glyphs show one compatible state; dashed glyphs mark multiple optima. Required is conditional, not biological confirmation.", 12)]
    branches = {(b["parent"], b["child"]): b for b in unit.get("ctmc", {}).get("branches", ())}
    for parent, child in tree.edges():
        elements += [line(xs[parent], ys[parent], xs[parent], ys[child]), line(xs[parent], ys[child], xs[child], ys[child])]
        if mode == "ctmc":
            b = branches.get((parent, child), {})
            if "probability_at_least_one_edit" in b:
                p = b["probability_at_least_one_edit"]
                elements.append(text((xs[parent]+xs[child])/2, ys[child]+18, f"p={p:.3f}", 10, "middle",
                    f'data-branch="{escape(child, quote=True)}" data-any-edit-probability="{p:.17g}"'))
        else:
            labels = [e["event_id"].split("/")[-1]+("?" if e["support"] != "required" else "")
                      for e in events if e["parent"] == parent and e["child"] == child]
            if labels: elements.append(text((xs[parent]+xs[child])/2, ys[child]+19, ",".join(labels), 10, "middle"))
    for node in tree.preorder():
        if tree.children.get(node):
            elements += [f'<circle cx="{xs[node]}" cy="{ys[node]}" r="3" fill="#20313D"/>',
                         text(xs[node]+6, ys[node]+12, node, 10)]
            elements += _ancestor_glyph(unit, node, xs[node], ys[node], mode, observation_view, nodes)
    label_x = max(xs.values())+30
    left, track_width = max(xs.values())+240., 920.
    if alignment: elements += _ribbons(c, observed, tree, ys, alignment, left, track_width)
    for leaf in (n for n in tree.preorder() if n in tree.leaves):
        species = tree.label[leaf]; o = observed.get(species, {})
        elements += [text(label_x, ys[leaf]+5, species, 14), text(label_x, ys[leaf]+24, o.get("kind", "unknown"), 11)]
        configs = o.get("configurations", ()) or [{"exons": []}]
        for lane, config in enumerate(configs):
            yy = ys[leaf]+lane*29
            elements += exon_track(left, yy, track_width, c["length"], config["exons"],
                                   GREY if o.get("kind") in {"unknown", "excluded"} else TEAL)
            if not config["exons"]:
                absent = any(p == 0 for p in o.get("material_presence", ()))
                elements.append(text(left+8, yy-14, "DNA tract absent: see genomic evidence" if absent else
                                     "No exon annotated here; absence is not established", 12))
            for material, present in zip(c.get("material", ()), o.get("material_presence", ())):
                if present == 0:
                    x = left+track_width*material["start"]/c["length"]
                    w = track_width*(material["end"]-material["start"])/c["length"]
                    elements += [rect(x, yy-5, w, 10, "white", "2,2"), text(x+w/2, yy+4, "//", 10, "middle")]
        alternatives = _alternatives(o)
        for i, config in enumerate(alternatives[:3]):
            yy = ys[leaf]+(len(configs)+i)*29
            elements += exon_track(left, yy, track_width, c["length"], config["exons"], "none", "4,3", 12)
            elements.append(text(left-115, yy+4, "candidate", 10))
        if len(alternatives) > 3:
            elements.append(text(left, ys[leaf]+(len(configs)+3)*29, f"+{len(alternatives)-3} candidates in the saved catalogue", 10))
    yy = bottom+8
    if mode == "ctmc":
        elements += [text(30, yy, "Branch p = probability of at least one edit; not the probability of different endpoints.", 14)]
        yy += 27
        if not nodes:
            elements.append(text(30, yy, "CTMC probabilities unavailable: inspect model_diagnostics.json; no parsimony substitute is drawn.", 13))
            yy += 25
        for b in branches.values():
            value = (f'{b["parent"]} → {b["child"]}: P(different endpoints)={b["probability_different_endpoints"]:.4f}; '
                     f'P(at least one edit)={b.get("probability_at_least_one_edit", float("nan")):.4f}')
            if "expected_edits" in b: value += f'; expected edits={b["expected_edits"]:.4f}'
            elements.append(text(30, yy, value, 12)); yy += 23
        elements.append(text(30, yy+5, "Fixed rates, root/origin prior, tree, candidate catalogue and correspondence are conditions, not jointly validated truth.", 12)); yy += 30
    else:
        elements.append(text(30, yy, "Elementary edits — IDs match structural_history.tsv; possible placements must not be added.", 14)); yy += 28
        for e in events:
            event_id = e["event_id"]
            value = f'{event_id} | {e["operation"]} | {e["support"]} | {e["parent"]} → {e["child"]} | {e.get("consequences", "")}'
            rows, used = paragraphs(30, yy, value, 170, 12, 18)
            elements += [f'<a xlink:href="structural_history.tsv" data-event-id="{escape(event_id, quote=True)}">', *rows, '</a>']; yy += used+7
        if not events:
            elements.append(text(30, yy, "No edit listed. Unknown or unanalysed input is not evidence of conservation.", 13)); yy += 25
    svg = document(int(max(1600, left+track_width+30)), int(yy+35), elements, title, "exon_configuration_v2:"+mode)
    if path is not None: Path(path).write_text(svg)
    return svg


def render_exon_results(result_dir, output_dir):
    result, out = Path(result_dir), Path(output_dir)
    if result_model(result) not in {"exon-parsimony", "exon-ctmc"}: raise ValueError("Expected completed exon-configuration result")
    data = json.loads((result/"exon_history.json").read_text())
    if data.get("model") != MODEL_VERSION:
        raise ValueError("Saved history is not the current configuration model; do not relabel old results")
    tree = SpeciesTree(data["tree"])
    events = read_tsv(result/"structural_history.tsv")
    preparation = {r["family_id"]: r for r in read_tsv(result/"exon_preparation_summary.tsv", optional=True)}
    out.mkdir(parents=True, exist_ok=True)
    for name in ("structural_history.tsv", "exon_ancestral_states.tsv", "exon_branch_posteriors.tsv", "model_diagnostics.json"):
        if (result/name).is_file() and (result/name).resolve() != (out/name).resolve(): shutil.copyfile(result/name, out/name)
    figures, sections = [], []
    for index, unit in enumerate(data["units"], 1):
        relevant = [r for r in events if r["family_id"] == unit["family_id"] and r["unit_id"] == unit["unit_id"]]
        record = preparation.get(unit["family_id"], {})
        source = result/record.get("evidence_directory", "")/"genomic_alignment.fa"
        alignment = parse_fasta(source) if source.is_file() else None
        for mode in (("parsimony", "ctmc") if data["engine"] == "exon-ctmc" else ("parsimony",)):
            filename = f'exon_structure_{index:05d}'+("_ctmc" if mode == "ctmc" else "")+".svg"
            draw_unit(unit, tree, relevant, alignment, out/filename, mode=mode, observation_view=data["observation_view"])
            figures.append({"family": unit["family_id"], "unit": unit["unit_id"], "svg": filename, "mode": mode,
                "computed_from": "exon_history.json; structural_history.tsv", "inference_performed_by_renderer": False})
            sections.append(f'<h2>{escape(unit["family_id"])} / {escape(unit["unit_id"])} — {mode}</h2><object data="{filename}" type="image/svg+xml" style="width:100%"></object>')
    (out/"figure_manifest.json").write_text(json.dumps(figures, indent=2)+"\n")
    (out/"index.html").write_text('<!doctype html><html lang="en"><meta charset="utf-8"><title>IntraPhy exon structures</title>'
        '<style>body{font:16px Arial;max-width:1600px;margin:30px auto;padding:20px}object{margin:10px 0 35px}</style>'
        '<h1>Exon structural evolution</h1><p>Solid boxes: annotation. Dashed boxes: sequence-supported alternatives. '
        'Multiple solid lanes: coexisting annotated structures, not usage measurements. // denotes absent genomic material in common coordinates.</p>'
        '<p>Parsimony and CTMC are separate views. Full probabilities: <a href="exon_ancestral_states.tsv">ancestral states</a>; '
        '<a href="exon_branch_posteriors.tsv">branches</a>. Read <a href="model_diagnostics.json">diagnostics</a> before interpreting missing events.</p>'
        +''.join(sections)+'</html>\n')
    return figures

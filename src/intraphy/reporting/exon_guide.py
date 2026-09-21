"""Complete current-model guide with reproducible, explicitly synthetic examples."""
from __future__ import annotations
from html import escape
from pathlib import Path
import json
from .exon_guide_examples import build_examples
from .exon_guide_plates import PLATES
from .exon_results import render_exon_results
from ..structure.serialization import write_json


def render_guide(output_dir, *, png=False):
    out=Path(output_dir);out.mkdir(parents=True,exist_ok=True)
    data=build_examples(out)
    results=render_exon_results(out/"example_results",out/"example_figures")
    records=[]
    for name,(draw,caption) in PLATES.items():
        path=out/(name+".svg");path.write_text(draw(data))
        records.append({"svg":path.name,"caption":caption,
            "purpose":"synthetic_teaching_schematic","data_result":False,
            "model":data["model"],"computed_from":"illustrative_model.json",
            "algorithm_version":"0.19.1","inference_performed_by_renderer":False})
    if png:
        try:
            import cairosvg
        except ImportError as exc:
            raise RuntimeError("PNG rendering requires optional CairoSVG") from exc
        for record in records:
            path=out/record["svg"]
            cairosvg.svg2png(url=str(path),write_to=str(path.with_suffix(".png")))
            record["png"]=path.with_suffix(".png").name
    write_json(out/"figure_manifest.json",records)
    sections=[]
    for record in records:
        title=record["caption"]
        sections.append(f'<h2>{escape(title)}</h2><object data="{record["svg"]}" type="image/svg+xml" style="width:100%"></object>')
    (out/"index.html").write_text('<!doctype html><html lang="en"><meta charset="utf-8">'
        '<title>IntraPhy 0.19.1 model guide</title><style>body{font:17px Arial;max-width:1600px;'
        'margin:30px auto;padding:20px;line-height:1.5}object{margin-bottom:35px}h2{margin-top:45px}</style>'
        '<h1>V19 exon configuration model · 0.19.1</h1>'
        '<p><strong>Purpose:</strong> exon structural evolution from genomic FASTA, GFF and a rooted species tree. '
        'These examples require no transcriptome and make no statement about exon usage.</p>'
        '<p>Original exons are compared as local configurations. Whole-locus sequence supports correspondence and '
        'annotation alternatives; the model then applies elementary structural edits to these configurations.</p>'
        '<p><strong>Synthetic, not biological validation:</strong> probabilities and minimum edit counts are computed '
        'by this version. Rates, input structures and outputs are saved in '
        '<a href="illustrative_model.json">illustrative_model.json</a>. '
        '<a href="example_figures/index.html">Inspect the separate parsimony and CTMC result figures</a>.</p>'
        +''.join(sections)+'<h2>Limits</h2><p>The model is conditional on an explicit finite candidate catalogue, '
        'sequence correspondence, annotation and tree. It does not reconstruct copy genealogies, inversions, '
        'exon shuffling, all extinct exons, complete transcript repertoires or molecular repair mechanisms. '
        'A candidate is not confirmed by a simpler evolutionary history.</p></html>\n')
    return records

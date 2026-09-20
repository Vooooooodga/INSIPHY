"""Three-tier biological guide, separate from data-derived result graphics."""
from __future__ import annotations
from html import escape
import json
from pathlib import Path
from . import overview, homology, model
from .example import probability_example

FIGURES=(
    ('method_overview','Overview: homologous regions to ancestral changes',overview),
    ('homology_inference','Algorithm: alignments, genomic context and incomplete annotation',homology),
    ('phylogenetic_model','Model: character coding and ancestral-state probabilities',model),
)
CAPTIONS={
 'method_overview': 'The same four-species gene family is followed throughout. Teal marks a homologous DNA segment; purple marks exonic sequence around a separate intron position. Dashed exons are annotation candidates, not confirmed transcription. Branch labels illustrate minimum-change placements of two different characters. Optional likelihood inference estimates ancestral states and rates when its assumptions and diagnostics permit. No complete ancestral transcript or physical mutation mechanism is inferred.',
 'homology_inference': 'Protein alignments are projected through coding coordinates; genomic alignments search DNA beyond existing exon annotations. Ordered-chain dynamic programming distinguishes complementary interval coverage from competing copies. Flank-bounded short alignments are evaluated with affine gaps. A DNA match without an informative transcript annotation is coded present for DNA but unknown for exonic status. A covering annotated transcript can instead support a non-exonic state, conditional on that annotation. Supported deletion, unresolved alignment and assembly gaps are not interchangeable. Illustrations explain rules rather than calibrated detection performance.',
 'phylogenetic_model': 'Rows are species; columns refer to the highlighted DNA interval, its exonic status, and the purple intron position. The tree is an input phylogeny, not an inferred gene tree. The probability example conditions on fixed rates and branch lengths and uses the intron column only. Pie fractions are computed by summing ancestral-state assignments and cross-checked against the production pruning algorithm; they are not arbitrary percentages or fitted biological results. Unknown and inapplicable observations have distinct meanings, though each is marginalized by the separate binary model. Rate fitting, AIC and posteriors are unavailable for known dependent included characters; small-sample calibration remains unassessed.',
}

def render_guide(output_dir, *, png=False):
    rasterizer=None
    if png:
        try:import cairosvg as rasterizer
        except ImportError as exc:raise RuntimeError('PNG export requires optional CairoSVG; SVG needs only the standard library.') from exc
    out=Path(output_dir);out.mkdir(parents=True,exist_ok=True)
    records=[];sections=[]
    for name,title,module in FIGURES:
        path=module.draw(out/f'{name}.svg')
        if rasterizer:rasterizer.svg2png(url=str(path),write_to=str(path.with_suffix('.png')),scale=1.5)
        records.append({'figure':name,'svg':path.name,'title':title,'purpose':'synthetic_teaching_schematic',
                        'data_result':False,'source_module':module.__name__,'caption':CAPTIONS[name]})
        sections.append(f'<section id="{name}"><h2>{escape(title)}</h2><a href="{path.name}">Editable SVG</a>'
            f'<object type="image/svg+xml" data="{path.name}" style="width:100%"></object><p>{escape(CAPTIONS[name])}</p></section>')
    (out/'figure_manifest.json').write_text(json.dumps(records,indent=2)+'\n',encoding='utf-8')
    (out/'illustrative_model.json').write_text(json.dumps(probability_example(),indent=2)+'\n',encoding='utf-8')
    (out/'index.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><title>IntraPhy method</title>'
        '<style>body{font:17px system-ui;margin:30px auto;padding:0 25px;max-width:1600px;line-height:1.6;color:#202D38}'
        'section{margin:50px 0}a{color:inherit}p{max-width:1100px}</style><h1>How IntraPhy reconstructs gene-structure changes</h1>'
        '<p>Read the overview first, then the sequence-matching algorithm and probability model. All examples are synthetic.'
        ' Real-data graphics are generated separately.</p>'+''.join(sections)+'</html>\n',encoding='utf-8')
    return records

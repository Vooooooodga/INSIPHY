"""Original, schematic figures for the methods description."""
from pathlib import Path
from . import overview, correspondence, states, events, inference

FIGURES = (
    ("method_scheme", overview.draw),
    ("exon_correspondence", correspondence.draw),
    ("structural_states", states.draw),
    ("evolutionary_events", events.draw),
    ("phylogenetic_inference", inference.draw),
)


def render_publication(output_dir, *, png=False, pdf=False):
    """Write native SVG plates and optional CairoSVG exports; return written paths."""
    rasterizer = None
    if png or pdf:
        try:
            import cairosvg as rasterizer
        except ImportError as exc:
            formats = "PNG and PDF" if png and pdf else "PNG" if png else "PDF"
            raise RuntimeError(f"{formats} export requires optional CairoSVG") from exc
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    written = []
    for name, draw in FIGURES:
        svg = out / f"{name}.svg"
        svg.write_text(draw(), encoding="utf-8")
        written.append(svg)
        if png:
            target = svg.with_suffix(".png")
            rasterizer.svg2png(url=str(svg), write_to=str(target), output_width=2700)
            written.append(target)
        if pdf:
            target = svg.with_suffix(".pdf")
            rasterizer.svg2pdf(url=str(svg), write_to=str(target))
            written.append(target)
    return written

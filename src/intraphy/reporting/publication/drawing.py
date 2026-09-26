"""Compact, accessible SVG layout helpers for the publication schematics."""
from __future__ import annotations
from html import escape
from ..exon_drawing import text as svg_text, exon_track

INK = "#20313D"
TEAL = "#356B78"
GOLD = "#B87923"
PURPLE = "#72558A"
RED = "#A34F4F"
MUTED = "#637580"
PALE = "#F2F5F6"
RULE = "#CBD3D7"


def text(x, y, value, size=16, anchor="start", color=INK, bold=False):
    weight = ' font-weight="600"' if bold else ""
    # At a 180 mm plate width, 28 SVG units correspond to approximately 8 pt.
    plate_size = max(28, size * 1.6)
    return svg_text(x, y, value, plate_size, anchor,
                    f'fill="{color}" font-family="Arial, sans-serif"{weight}')


def math_text(x, y, segments, size=16, color=INK, bold=False):
    """Draw mathematical labels with true SVG subscripts and restored baseline."""
    plate_size = max(28, size * 1.6)
    weight = ' font-weight="600"' if bold else ""
    parts = []
    previous_offset = 0.0
    for value, subscript in segments:
        current_offset = plate_size * .30 if subscript else 0.0
        font_size = plate_size * .70 if subscript else plate_size
        delta = current_offset - previous_offset
        parts.append(f'<tspan font-size="{font_size:.1f}" dy="{delta:.1f}">{escape(value)}</tspan>')
        previous_offset = current_offset
    return (f'<text x="{x}" y="{y}" font-size="{plate_size}" fill="{color}" '
            f'font-family="Arial, sans-serif"{weight}>{"".join(parts)}</text>')


def line(x1, y1, x2, y2, color=INK, width=2, dash="", arrow=False):
    marker = ' marker-end="url(#arrow)"' if arrow else ""
    return (f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" '
            f'stroke-width="{width}" stroke-dasharray="{dash}"{marker}/>')


def box(x, y, width, height, fill="white", stroke=RULE, radius=8):
    return (f'<rect x="{x}" y="{y}" width="{width}" height="{height}" rx="{radius}" '
            f'fill="{fill}" stroke="{stroke}" stroke-width="1.5"/>')


def tree(x, y, scale=1, labels=True):
    """A rooted four-tip species tree with no numerical branch lengths."""
    xa, xb = x+70*scale, x+145*scale
    out = [line(x, y+30, x, y+150),
           line(x, y+30, xa, y+30), line(x, y+150, xa, y+150),
           line(xa, y, xa, y+60), line(xa, y+120, xa, y+180),
           line(xa, y, xb, y), line(xa, y+60, xb, y+60),
           line(xa, y+120, xb, y+120), line(xa, y+180, xb, y+180)]
    if labels:
        out.extend(text(xb+8, y+6+i*60, s, 13) for i,s in enumerate("ABCD"))
    return out


def document(width, height, title, body, description):
    """Return a self-contained SVG with a title, accessible palette and arrowhead."""
    physical_height = height * 180 / width
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="180mm" height="{physical_height:.2f}mm" '
            f'viewBox="0 0 {width} {height}"><title>{escape(title)}</title>'
            f'<desc>{escape(description)}</desc>'
            '<defs><marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" '
            'markerWidth="7" markerHeight="7" orient="auto"><path d="M0 0L10 5L0 10Z" '
            f'fill="{INK}"/></marker></defs>'
            '<rect width="100%" height="100%" fill="white"/>'
            f'<g font-family="Arial, sans-serif" fill="{INK}">{"".join(body)}</g></svg>\n')


def track(x, y, width, length, exons, color=TEAL, height=20, dash=""):
    return exon_track(x, y, width, length, exons, fill=color, height=height, dash=dash)


def bullet(x, y, mark, label, color=INK):
    return [text(x, y, mark, 18, color=color, bold=True), text(x+24, y, label, 15, color=color)]

"""Small deterministic SVG primitives for exon structures, never inference."""
from __future__ import annotations
from html import escape
import textwrap

INK = "#20313D"
TEAL = "#427788"
AMBER = "#AC641C"
GREY = "#E1E6E8"


def text(x, y, value, size=14, anchor="start", attrs=""):
    return f'<text x="{x:.2f}" y="{y:.2f}" font-size="{size}" text-anchor="{anchor}" {attrs}>{escape(str(value))}</text>'


def line(x1, y1, x2, y2, stroke="#647583", dash=""):
    return f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" stroke="{stroke}" stroke-width="1.5" stroke-dasharray="{dash}"/>'


def rect(x, y, width, height, fill=TEAL, dash="", attrs=""):
    return f'<rect x="{x:.2f}" y="{y:.2f}" width="{max(width, .5):.2f}" height="{height:.2f}" fill="{fill}" stroke="#294353" stroke-dasharray="{dash}" {attrs}/>'


def exon_track(x, y, width, length, exons, fill=TEAL, dash="", height=16):
    result = [line(x, y, x+width, y, "#A5B0B8")]
    for exon in exons:
        a, b = (exon["start"], exon["end"]) if isinstance(exon, dict) else exon
        result.append(rect(x+width*a/length, y-height/2, width*(b-a)/length, height, fill, dash))
    return result


def paragraphs(x, y, message, chars=135, size=13, leading=19):
    lines = textwrap.wrap(str(message), width=chars) or [""]
    return [text(x, y+i*leading, value, size) for i, value in enumerate(lines)], len(lines)*leading


def document(width, height, body, title, model="exon_configuration_v2"):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
            f'width="{width}" height="{height}" viewBox="0 0 {width} {height}" data-model="{escape(model, quote=True)}">'
            f'<title>{escape(title)}</title><rect width="100%" height="100%" fill="white"/>'
            f'<g font-family="Arial, sans-serif" fill="{INK}">'+"\n".join(body)+"</g></svg>\n")

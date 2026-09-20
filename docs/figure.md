# Methods architecture figure

The README embeds `figures/architecture.svg`, a vector diagram with selectable
text and explicit evidence, character and inference stages. It is independent of
model result plots and does not imply biological validation.

Regenerate the exact SVG with the standard library:

```bash
python tools/render_architecture.py
```

For a high-resolution PNG preview:

```bash
python -m pip install cairosvg
python tools/render_architecture.py --png
```

The SVG is the publication source. It can be edited in Inkscape or another SVG
editor, and exported to the format and dimensions required by a journal. Use the
SVG for vector artwork; the PNG is a preview. Rendering requires locally available
fonts; no font files are distributed. Labels use established gene-structure and
phylogenetic terminology. Do not remove the distinction between structural
character changes and mutation events when adapting the figure.

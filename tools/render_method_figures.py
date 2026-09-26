#!/usr/bin/env python3
"""Render original IntraPhy explanatory SVGs and optional PNG previews."""
from pathlib import Path
import argparse
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=None)
    parser.add_argument('--png', action='store_true', help='Also render PNG previews; requires CairoSVG.')
    parser.add_argument('--pdf', action='store_true', help='Also render PDF plates; requires CairoSVG.')
    parser.add_argument('--publication', action='store_true',
                        help='Render the five publication schematics instead of the interactive guide.')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    if args.publication:
        from intraphy.reporting.publication import render_publication
        output_dir = args.output_dir or root / 'docs/figures/publication'
        for path in render_publication(output_dir, png=args.png, pdf=args.pdf):
            print(path)
        return
    if args.pdf:
        parser.error('--pdf is available with --publication')
    from intraphy.reporting.exon_guide import render_guide
    output_dir = args.output_dir or root / 'docs/figures'
    for row in render_guide(output_dir, png=args.png):
        print(output_dir / row['svg'])


if __name__ == '__main__':
    main()

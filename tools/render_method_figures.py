#!/usr/bin/env python3
"""Render original IntraPhy explanatory SVGs and optional PNG previews."""
from pathlib import Path
import argparse
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from intraphy.reporting.methods import render_guide


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path,
                        default=Path(__file__).resolve().parents[1] / 'docs/figures')
    parser.add_argument('--png', action='store_true', help='Also render PNG previews; requires CairoSVG.')
    args = parser.parse_args()
    for row in render_guide(args.output_dir, png=args.png):
        print(args.output_dir / row['svg'])


if __name__ == '__main__':
    main()

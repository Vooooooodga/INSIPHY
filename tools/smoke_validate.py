#!/usr/bin/env python3
"""Run an installed-package CLI example and retain reproducible diagnostic outputs."""
from pathlib import Path
import argparse
import json
import os
import subprocess
import sys
import xml.etree.ElementTree as ET


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--python',default=sys.executable,help='Interpreter of the installed distribution.')
    args=parser.parse_args(); root=args.output_dir.resolve()
    if root.exists() and any(root.iterdir()): parser.error('Smoke output directory must be empty')
    root.mkdir(parents=True,exist_ok=True)
    commands=[]
    def run(*values):
        command=[args.python,'-I','-m','intraphy',*map(str,values)]
        result=subprocess.run(command,cwd=root,env={k:v for k,v in os.environ.items() if k!='PYTHONPATH'},
                              text=True,capture_output=True,timeout=240)
        commands.append(dict(command=command,returncode=result.returncode,
                             stdout=result.stdout,stderr=result.stderr))
        (root/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
        if result.returncode: raise RuntimeError(result.stderr or result.stdout)
    run('--version')
    run('inspect-aligners')
    run('example','--output-dir',root/'native')
    run('check','--manifest',root/'native/manifest.tsv','--species-tree',root/'native/species_tree.nwk')
    run('build-case','--manifest',root/'native/manifest.tsv','--species-tree',root/'native/species_tree.nwk',
        '--output-dir',root/'prepared','--threads',2)
    run('run','--input-dir',root/'prepared','--output-dir',root/'parsimony','--threads',2)
    run('infer-phylogeny','--input-dir',root/'prepared','--output-dir',root/'erard','--model','er-ard',
        '--structural-site-matrix',root/'parsimony/structural_site_matrix.tsv')
    run('visualize','--input-dir',root/'prepared','--result-dir',root/'parsimony','--output-dir',root/'figures')
    svg=list((root/'figures').rglob('*.svg'))
    if not svg: raise AssertionError('No result figures generated')
    for path in svg: ET.parse(path)
    print(json.dumps(dict(commands=len(commands),svg_files=len(svg),status='passed',
                         validation_scope='installed_package_synthetic_raw_input'),indent=2))

if __name__=='__main__':
    main()

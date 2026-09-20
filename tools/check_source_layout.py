#!/usr/bin/env python3
"""Check module size, parseability, active English text and the package namespace."""
from pathlib import Path
import ast
import re
import sys


def main():
    root=Path(__file__).resolve().parents[1]
    failures=[]; counts=[]
    for folder in ('src','tests','tools'):
        for path in (root/folder).rglob('*.py'):
            text=path.read_text(encoding='utf-8'); count=len(text.splitlines())
            counts.append((count,str(path.relative_to(root))))
            if count>500: failures.append(f'{path.relative_to(root)}: {count} lines exceeds 500')
            try: ast.parse(text,filename=str(path),feature_version=(3,10))
            except SyntaxError as exc: failures.append(str(exc))
            if re.search('[\u3400-\u9fff]',text): failures.append(f'{path.relative_to(root)}: active non-English text')
    if (root/'src'/('insi'+'phy')).exists(): failures.append('Removed package namespace is still present')
    print(f'{len(counts)} Python modules; maximum {max(counts)[0]} lines ({max(counts)[1]})')
    for failure in failures: print(failure,file=sys.stderr)
    return int(bool(failures))

if __name__=='__main__':
    raise SystemExit(main())

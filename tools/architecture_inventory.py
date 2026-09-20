#!/usr/bin/env python3
"""Read-only Python architecture inventory. Python >=3.9, stdlib only.

Run on a LOCAL checkout; does not import or execute the audited package.
No file checksums, source dumps, network or package installation. Writes only --out; audited source is not modified.
AST branch counts are review aids, not a certified complexity measure.
"""
from __future__ import annotations

import argparse
import ast
import copy
import csv
import json
import re
import subprocess
from collections import defaultdict
from pathlib import Path

HASH_WORDS = re.compile(r'\b(?:hashlib|sha256|sha1|md5|hexdigest|file_digest|checksum)\b|\bhash\s*\(')
SERIAL_CALLS = {'read_tsv', 'write_tsv', 'json.loads', 'json.dumps', 'read_text', 'write_text'}


def call_name(node):
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        base = call_name(node.value)
        return (base + '.' if base else '') + node.attr
    return ''


def write_table(path, rows, fields):
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, '') for key in fields})


def module_name(path, source_root):
    parts = list(path.relative_to(source_root).with_suffix('').parts)
    if parts[-1] == '__init__':
        parts.pop()
    return '.'.join(parts)


def function_signature(node):
    # Exact AST clone detection except outer function name, docstring, locations.
    # Identifiers in the body are NOT normalized: near-duplicates are not detected.
    clone = copy.deepcopy(node)
    clone.name = '__FUNCTION__'
    if clone.body and isinstance(clone.body[0], ast.Expr):
        expr = clone.body[0].value
        if isinstance(expr, ast.Constant) and isinstance(expr.value, str):
            clone.body = clone.body[1:]
    return ast.dump(clone, include_attributes=False)


class InventoryVisitor(ast.NodeVisitor):
    def __init__(self, path, module, known_modules, package):
        self.path, self.module = path, module
        self.known_modules, self.package = known_modules, package
        self.scope = []
        self.functions, self.imports, self.clones = [], [], []

    def visit_ClassDef(self, node):
        self.scope.append(node.name)
        self.generic_visit(node)
        self.scope.pop()

    def visit_FunctionDef(self, node):
        end = getattr(node, 'end_lineno', node.lineno)
        args = len(node.args.posonlyargs) + len(node.args.args) + len(node.args.kwonlyargs)
        args += int(node.args.vararg is not None) + int(node.args.kwarg is not None)
        # This includes nested control flow; report it as a proxy, not cyclomatic complexity.
        branches = sum(isinstance(x, (ast.If, ast.For, ast.AsyncFor, ast.While,
                                      ast.ExceptHandler, ast.IfExp)) for x in ast.walk(node))
        branches += sum(len(x.values) - 1 for x in ast.walk(node) if isinstance(x, ast.BoolOp))
        name = '.'.join(self.scope + [node.name])
        row = {'file': self.path, 'function': name, 'start': node.lineno, 'end': end,
               'span_lines': end - node.lineno + 1, 'argument_count': args,
               'branch_points_proxy': branches}
        self.functions.append(row)
        if row['span_lines'] >= 8:
            self.clones.append((function_signature(node), row))
        self.scope.append(node.name)
        self.generic_visit(node)
        self.scope.pop()

    visit_AsyncFunctionDef = visit_FunctionDef

    def _record(self, target, line, imported_name=''):
        if target not in self.known_modules:
            return
        self.imports.append({'source': self.module, 'target': target, 'line': line,
                             'scope': '.'.join(self.scope) or '<module>',
                             'imported_name': imported_name,
                             'private_name': imported_name.startswith('_')})

    def visit_Import(self, node):
        for item in node.names:
            self._record(item.name, node.lineno)

    def visit_ImportFrom(self, node):
        if node.level:
            bits = self.package.split('.')
            trim = node.level - 1
            if trim > len(bits):
                return
            base_bits = bits[:len(bits) - trim] if trim else bits
            base = '.'.join(base_bits + ([node.module] if node.module else []))
        else:
            base = node.module or ''
        for item in node.names:
            candidate_module = '.'.join(part for part in (base, item.name) if part)
            if candidate_module in self.known_modules:
                self._record(candidate_module, node.lineno, item.name)
            else:
                self._record(base, node.lineno, item.name)


def strongly_connected(graph):
    index = 0
    stack, on_stack, indices, lows, found = [], set(), {}, {}, []

    def visit(node):
        nonlocal index
        indices[node] = lows[node] = index
        index += 1
        stack.append(node)
        on_stack.add(node)
        for other in sorted(graph.get(node, ())):
            if other not in indices:
                visit(other)
                lows[node] = min(lows[node], lows[other])
            elif other in on_stack:
                lows[node] = min(lows[node], indices[other])
        if lows[node] == indices[node]:
            component = []
            while True:
                other = stack.pop()
                on_stack.remove(other)
                component.append(other)
                if other == node:
                    break
            if len(component) > 1 or node in graph.get(node, ()):
                found.append(sorted(component))

    for node in sorted(graph):
        if node not in indices:
            visit(node)
    return sorted(found)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('repo', type=Path)
    parser.add_argument('--source', default='src/intraphy')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    repo = args.repo.resolve()
    package_root = (repo / args.source).resolve()
    if not package_root.is_dir():
        parser.error(f'Source directory not found: {package_root}')
    paths = sorted(package_root.rglob('*.py'))
    if not paths:
        parser.error('No Python source files found')
    args.out.mkdir(parents=True, exist_ok=True)
    source_root = package_root.parent
    modules = {module_name(path, source_root) for path in paths}
    files, functions, imports, hashes, failures = [], [], [], [], []
    clones = defaultdict(list)
    graph = {module: set() for module in modules}
    for path in paths:
        relative = str(path.relative_to(repo))
        data = path.read_bytes()
        try:
            text = data.decode('utf-8-sig')
            tree = ast.parse(text, filename=relative)
        except (UnicodeError, SyntaxError) as error:
            failures.append({'file': relative, 'error': str(error)})
            continue
        module = module_name(path, source_root)
        package = module if path.name == '__init__.py' else module.rpartition('.')[0]
        visitor = InventoryVisitor(relative, module, modules, package)
        visitor.visit(tree)
        functions.extend(visitor.functions)
        imports.extend(visitor.imports)
        for signature, row in visitor.clones:
            clones[signature].append(row)
        for edge in visitor.imports:
            graph[edge['source']].add(edge['target'])
        lines = text.splitlines()
        for number, line in enumerate(lines, 1):
            if HASH_WORDS.search(line):
                hashes.append({'file': relative, 'line': number, 'text': line.strip()[:240],
                               'classification': 'text_mention_requires_manual_review'})
        calls = [call_name(x.func) for x in ast.walk(tree) if isinstance(x, ast.Call)]
        files.append({'file': relative, 'bytes': len(data), 'physical_lines': len(lines),
                      'nonblank_lines': sum(bool(line.strip()) for line in lines),
                      'functions_including_nested': len(visitor.functions),
                      'largest_function_span': max((r['span_lines'] for r in visitor.functions), default=0),
                      'internal_import_targets': len(graph[module]),
                      'dot_get_calls': sum(name.endswith('.get') for name in calls),
                      'serialization_call_sites': sum(name in SERIAL_CALLS or name.rsplit('.', 1)[-1] in SERIAL_CALLS for name in calls),
                      'bare_or_broad_except': sum(isinstance(x, ast.ExceptHandler) and
                           (x.type is None or isinstance(x.type, ast.Name) and x.type.id in {'Exception', 'BaseException'})
                           for x in ast.walk(tree))})
    duplicated = []
    for group in clones.values():
        if len(group) > 1:
            duplicated.append({'members': [{k: r[k] for k in ('file', 'function', 'start', 'end')} for r in group]})
    try:
        result = subprocess.run(['git', '-C', str(repo), 'rev-parse', 'HEAD'],
                                capture_output=True, text=True, timeout=5, check=False)
        revision = result.stdout.strip() if result.returncode == 0 else None
    except (OSError, subprocess.TimeoutExpired):
        revision = None
    summary = {'revision': revision, 'source': str(package_root), 'file_count': len(files),
               'source_bytes': sum(x['bytes'] for x in files),
               'physical_lines': sum(x['physical_lines'] for x in files),
               'static_import_cycles': strongly_connected(graph),
               'exact_ast_clone_groups': duplicated, 'parse_failures': failures,
               'scope': 'Static AST analysis; local/conditional imports included; dynamic imports not resolved. '
                        'Cycles are potential dependencies, not proof of runtime failure. '
                        'Hash mentions include comments/identifiers, not confirmed redundant checks. '
                        'No audited code is executed and no source checksums are computed.'}
    (args.out / 'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    for name, rows, fields in (
        ('files.tsv', files, ['file', 'bytes', 'physical_lines', 'nonblank_lines', 'functions_including_nested',
                             'largest_function_span', 'internal_import_targets', 'dot_get_calls',
                             'serialization_call_sites', 'bare_or_broad_except']),
        ('functions.tsv', functions, ['file', 'function', 'start', 'end', 'span_lines', 'argument_count', 'branch_points_proxy']),
        ('imports.tsv', imports, ['source', 'target', 'line', 'scope', 'imported_name', 'private_name']),
        ('hash_mentions.tsv', hashes, ['file', 'line', 'classification', 'text']),
    ):
        write_table(args.out / name, rows, fields)
    print(json.dumps({key: summary[key] for key in ('revision', 'file_count', 'source_bytes', 'physical_lines')}, indent=2))
    print('Largest files:')
    for row in sorted(files, key=lambda r: r['bytes'], reverse=True)[:8]:
        print(f"  {row['file']}: {row['physical_lines']} lines, {row['bytes']} bytes")
    print(f"Potential import-cycle groups: {len(summary['static_import_cycles'])}; AST clone groups: {len(duplicated)}")
    print(f'Reports: {args.out.resolve()}')
    if failures:
        raise SystemExit(2)


if __name__ == '__main__':
    main()

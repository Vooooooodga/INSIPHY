#!/usr/bin/env python3
"""Run source-checkout CLI smoke checks on a synthetic fixture.

Requires the source tree including tests. Replaces only validation/cli_fixture_v017.
It does not run biological demos or install external programs.
"""
from pathlib import Path
import ast, json, subprocess, sys, os, shutil, xml.etree.ElementTree as ET
root=Path(__file__).resolve().parents[1]
work=root/'validation'/'cli_fixture_v017'
if work.exists():shutil.rmtree(work)
work.mkdir(parents=True)
sys.path[:0]=[str(root/'src'),str(root/'tests')]
from test_insiphy import SingleCopyPhylogenyTests
from insiphy.structural_sites import build_structural_site_matrix
fixture=work/'prepared'
fixture.mkdir(exist_ok=True)
input_dir,result_dir=SingleCopyPhylogenyTests().write_structural_case(fixture)
rows,_=build_structural_site_matrix(input_dir,result_dir)
# Published local smoke fixture: no real species and no fitting truth is invented.
matrix=result_dir/'structural_site_matrix.tsv'
from insiphy.io import write_structural_site_matrix
write_structural_site_matrix(matrix,rows)
commands=[]
env=dict(os.environ,PYTHONPATH=str(root/'src'))
def run(args):
 p=subprocess.run([sys.executable,'-m','insiphy.cli']+list(map(str,args)),env=env,text=True,capture_output=True,timeout=40)
 commands.append({'args':list(map(str,args)),'exit_code':p.returncode,'stdout':p.stdout,'stderr':p.stderr})
 if p.returncode:raise RuntimeError(p.stderr or p.stdout)
 return p
run(['--version']);run(['--help'])

foreground=work/'foreground.tsv'
foreground.write_text('parent_id\tchild_id\nab\ta\n')
checks=[]
for model in ('parsimony','er-ard','foreground'):
    for threads in (1,2):
        out=work/(model+'_'+str(threads))
        args=['infer-phylogeny','--input-dir',input_dir,'--output-dir',out,
              '--model',model,'--structural-site-matrix',matrix,'--threads',threads]
        if model=='foreground':args+=['--foreground-branches',foreground]
        run(args)
        manifest=json.loads((out/'run_result.json').read_text())
        assert manifest['model']==model and manifest['status']=='completed'
        assert manifest['artifacts'] and all((out/name).is_file() for name in manifest['artifacts'])
    files = ['node_structural_states.tsv','branch_structural_events.tsv','structural_site_summary.tsv','compound_structural_events.tsv'] if model=='parsimony' else ['model_fits.tsv','model_tests.tsv','node_state_posteriors.tsv','branch_transition_posteriors.tsv','structural_changes.tsv']
    for name in files:
        left=(work/(model+'_1')/name).read_text()
        right=(work/(model+'_2')/name).read_text()
        assert left==right,(model,name)
    checks.append({'model':model,'threads':[1,2],'exactly_equal_result_tables':files})
    run(['visualize','--input-dir',input_dir,'--result-dir',work/(model+'_1'),
         '--output-dir',work/(model+'_figures')])
    svg=list((work/(model+'_figures')).rglob('*.svg'))
    assert svg
    for path in svg:ET.parse(path)
    checks[-1]['valid_svg_files']=len(svg)
for path in (root/'src').rglob('*.py'):
    ast.parse(path.read_text(),filename=str(path),feature_version=(3,9))
commands=[{**row,'args':[a.replace(str(root),'PROJECT') for a in row['args']]} for row in commands]
(root/'validation'/'cli_smoke_results_v017.json').write_text(json.dumps({'commands':commands,'checks':checks,'syntax_python39':True,'data':'synthetic four-species fixture, not biological validation'},indent=2)+'\n')
print(json.dumps(checks,indent=2))

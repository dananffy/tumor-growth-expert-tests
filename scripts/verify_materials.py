"""Check the code/input inventory. This does not run scientific experiments."""
import argparse,ast,hashlib,json,sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
GENERATED={'runs','results','outputs','demo_runs','working','__pycache__','.git','.venv','venv','.ipynb_checkpoints'}

def members():
    for path in ROOT.rglob('*'):
        rel=path.relative_to(ROOT)
        if not path.is_file():continue
        generated=set(rel.parts)&GENERATED
        if generated and not (path.name=='.gitkeep' and generated=={'results'}):continue
        if rel.as_posix()=='MANIFEST.json' or path.suffix in {'.pyc','.asv'}:continue
        yield path,rel.as_posix()

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--refresh',action='store_true',help='Explicitly regenerate the file inventory after an adopted change.')
    args=parser.parse_args()
    if args.refresh:
        files={rel:{'bytes':path.stat().st_size,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()} for path,rel in sorted(members(),key=lambda x:x[1])}
        (ROOT/'MANIFEST.json').write_text(json.dumps({'scope':'Source and required inputs; generated results excluded.','files':files},ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n')
        print(f'Regenerated MANIFEST.json: {len(files)} files.');return
    expected=json.loads((ROOT/'MANIFEST.json').read_text(encoding='utf-8'))['files']
    actual=dict((rel,path) for path,rel in members())
    errors=[]
    for name in sorted(set(expected)-set(actual)):errors.append('Missing: '+name)
    for name in sorted(set(actual)-set(expected)):errors.append('Unlisted: '+name)
    for name in sorted(set(actual)&set(expected)):
        data=actual[name].read_bytes()
        if len(data)!=expected[name]['bytes'] or hashlib.sha256(data).hexdigest()!=expected[name]['sha256']:
            errors.append('Changed: '+name)
    for path,rel in members():
        if path.suffix=='.py':
            try:ast.parse(path.read_text(encoding='utf-8'),filename=rel)
            except SyntaxError as err:errors.append(str(err))
        if path.suffix=='.ipynb':
            for index,cell in enumerate(json.loads(path.read_text(encoding='utf-8'))['cells']):
                if cell['cell_type']=='code':
                    try:ast.parse(''.join(cell['source']),filename=f'{rel}:cell{index+1}')
                    except SyntaxError as err:errors.append(str(err))
    tests=json.loads((ROOT/'TESTS.json').read_text(encoding='utf-8'))
    for test in tests:
        for rel in test['files']:
            if not (ROOT/rel).is_file():errors.append(test['id']+' dependency missing: '+rel)
    code=ROOT/'A/baseline/code'
    cases=json.loads((code/'cases.json').read_text(encoding='utf-8'))
    for name,case in cases.items():
        if case['kind']=='python' and not (code/case['script']).is_file():errors.append('Case script missing: '+name)
        if case['kind']=='cpp' and not (code/case['source']/'CMakeLists.txt').is_file():errors.append('CMake source missing: '+name)
        if case.get('prepare_working_source') and not (code/case['prepare_working_source']).is_dir():errors.append('Preparation source missing: '+name)
    if errors:
        print('\n'.join(errors));raise SystemExit(1)
    print(f'OK: {len(expected)} files, {len(tests)} test entries, {len(cases)} A-group configurations.')
    print('SHA256, declared inputs and Python/Notebook syntax checked; scientific results are not certified here.')

if __name__=='__main__':main()

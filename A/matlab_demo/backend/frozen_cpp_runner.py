"""Installed in a frozen packet's code/ directory; Ubuntu standard library."""
import argparse,datetime,hashlib,json,os,platform,shutil,subprocess,time,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--case',required=True);ap.add_argument('--build-only',action='store_true');ap.add_argument('--skip-build',action='store_true');a=ap.parse_args()
    cases=json.loads((ROOT/'cases.json').read_text(encoding='utf-8'));case=cases[a.case]
    source=ROOT/case['source'];hashes={str(p.relative_to(source)):sha(p) for p in source.rglob('*') if p.is_file()}
    key=hashlib.sha256(json.dumps(hashes,sort_keys=True).encode()).hexdigest()[:12]
    cache=Path.home()/'.cache/tumor_expert_test'/f"github_{a.case}_{key}";cache.mkdir(parents=True,exist_ok=True)
    src=cache/'source';build=cache/'build';exe=build/case['executable']
    run_id=datetime.datetime.now().astimezone().strftime('%Y%m%d_%H%M%S_%f')+'_frozen_'+a.case
    group=case['group'];run=ROOT/'runs'/group/run_id
    (run/'result/points').mkdir(parents=True);(run/('result/functions' if group=='T01' else 'result/solutions')).mkdir()
    config={'grid':case['grid'],'dt':case['dt'],'threads':2}
    if group=='T01':config['final_time']=case['end']
    else:config.update({'end':case['end'],'mode':case['mode']})
    meta={'run_id':run_id,'case_id':a.case,'source_origin':'A baseline 2026-11-09','config':config,'source_sha256':hashes,'platform':platform.platform(),'parameters':case.get('parameters',{}),'versions':{tool:subprocess.check_output([tool,'--version'],text=True).splitlines()[0] for tool in ['g++','cmake']}}
    start=time.perf_counter()
    with (run/'build.log').open('w') as log:
        if not a.skip_build:
            shutil.copytree(source,src,dirs_exist_ok=True)
            for cmd in [['cmake','-S',str(src),'-B',str(build),'-DCMAKE_BUILD_TYPE=Release'],['cmake','--build',str(build),'--parallel','2']]:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
        elif not exe.is_file():raise RuntimeError('No matching prebuilt executable; run --build-only first')
    meta['build_seconds']=time.perf_counter()-start;meta['executable_sha256']=sha(exe)
    if a.build_only:
        (run/'build_only.json').write_text(json.dumps(meta,indent=2));print(json.dumps({'case':a.case,'build_only':True,'cache':str(cache)}));return
    command=[str(exe),str(case['grid']),str(case['dt']),str(case['end'])]
    if group=='T01':command.append('2')
    meta['command']=command;start=time.perf_counter()
    with (run/'run.log').open('w') as log:result=subprocess.run(command,cwd=run,stdout=log,stderr=subprocess.STDOUT,timeout=7200)
    meta['run_seconds']=time.perf_counter()-start;meta['exit_code']=result.returncode
    (run/'config_and_environment.json').write_text(json.dumps(meta,indent=2))
    with zipfile.ZipFile(run/'solver_source_snapshot.zip','w',zipfile.ZIP_DEFLATED) as z:
        for p in source.rglob('*'):
            if p.is_file():z.write(p,p.relative_to(source))
    if result.returncode:raise RuntimeError(f'Failed; inspect {run}/run.log')
    pointer='latest_run.txt' if group=='T01' else f"latest_{case['mode']}.txt"
    (ROOT/'runs'/group/pointer).write_text(run_id)
    print(json.dumps({'run_directory':str(run),'run_seconds':meta['run_seconds'],'exit_code':0}))
if __name__=='__main__':main()

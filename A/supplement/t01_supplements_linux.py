"""Isolated supplementary runs, Linux stdlib; original freeze is read-only."""
import datetime,hashlib,json,os,platform,shutil,subprocess,time,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parent;CODE=ROOT.parent
SOURCE=ROOT/'sources/t01'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    hashes={p.relative_to(SOURCE).as_posix():sha(p) for p in SOURCE.rglob('*') if p.is_file()}
    cache=Path.home()/'.cache/tumor_expert_test/review_20261110_t01';cache.mkdir(parents=True,exist_ok=True)
    shutil.copytree(SOURCE,cache/'source',dirs_exist_ok=True)
    with (ROOT/'t01_build.log').open('w') as log:
        for cmd in [['cmake','-S',str(cache/'source'),'-B',str(cache/'build'),'-DCMAKE_BUILD_TYPE=Release'],['cmake','--build',str(cache/'build'),'--parallel','2']]:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
    records=[]
    for grid,dt in [(256,.02),(256,.005),(256,.0025),(128,.0025),(64,.0025)]:
        rid=datetime.datetime.now().strftime('%Y%m%d_%H%M%S_%f')+f'_n{grid}_dt{dt:g}'
        run=ROOT/'runs/T01'/rid;(run/'result/points').mkdir(parents=True);(run/'result/functions').mkdir()
        cmd=[str(cache/'build/t01'),str(grid),str(dt),'.2','2'];start=time.perf_counter()
        meta={'run_id':rid,'config':{'grid':grid,'dt':dt,'final_time':.2,'threads':2},'parameters':{'lambda':.5,'c_B':10,'G0':.2,'R0':.5,'n_ctrl':32,'box':[-1.5,1.5],'M':4*grid},'source_sha256':hashes,'source_note':'Exact frozen T01 sources; no numerical changes.','build_seconds':None,'command':cmd,'platform':platform.platform(),'executable_sha256':sha(cache/'build/t01')}
        (run/'config_and_environment.json').write_text(json.dumps(meta,indent=2))
        with (run/'run.log').open('w') as log:r=subprocess.run(cmd,cwd=run,stdout=log,stderr=subprocess.STDOUT,timeout=3600)
        meta.update({'run_seconds':time.perf_counter()-start,'exit_code':r.returncode});(run/'config_and_environment.json').write_text(json.dumps(meta,indent=2))
        with zipfile.ZipFile(run/'solver_source_snapshot.zip','w',zipfile.ZIP_DEFLATED) as z:
            for p in SOURCE.rglob('*'):
                if p.is_file():z.write(p,p.relative_to(SOURCE))
        assert r.returncode==0,(rid,r.returncode)
        records.append({'run_id':rid,'grid':grid,'dt':dt,'run_seconds':meta['run_seconds']});(ROOT/'t01_supplement_runs.json').write_text(json.dumps(records,indent=2));print(json.dumps(records[-1]),flush=True)
if __name__=='__main__':main()

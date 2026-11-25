"""Spatial onset check using the actual stage-1 nutrient/pressure solvers."""
import json,datetime,hashlib,shutil,subprocess,time,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def main():
    rid=datetime.datetime.now().astimezone().strftime('%Y%m%d_%H%M%S_%f')+'_onset_probe';run=ROOT/'runs/T03'/rid;run.mkdir(parents=True)
    source=ROOT/'working/t03_emergence';cache=Path('/tmp/tumor_t03_onset_probe');src=cache/'source';src.mkdir(parents=True,exist_ok=True)
    shutil.copytree(source,src,dirs_exist_ok=True);shutil.copy2(ROOT/'t03/onset_probe.cpp',src/'src/test.cpp')
    with zipfile.ZipFile(run/'solver_source_snapshot.zip','w',zipfile.ZIP_DEFLATED) as z:
        for p in src.rglob('*'):
            if p.is_file():z.write(p,p.relative_to(src))
        z.write(__file__,'workflow/run_onset_probe.py')
    with (run/'build.log').open('w') as log:
        for cmd in [['cmake','-S',str(src),'-B',str(cache/'build'),'-DCMAKE_BUILD_TYPE=Release'],['cmake','--build',str(cache/'build'),'--parallel','2']]:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
    values=[];roots=[];start=time.perf_counter()
    for grid in [128,256]:
        def evaluate(R):
            folder=run/f'n{grid}_sample{len(values):02d}';folder.mkdir()
            with (folder/'run.log').open('w') as log:subprocess.run([str(cache/'build/t03'),str(grid),repr(R)],cwd=folder,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=180)
            row=json.loads((folder/'probe.json').read_text());values.append(row);print(json.dumps(row),flush=True);return row['center_pressure']
        lo,hi=2.24,2.26;flo,fhi=evaluate(lo),evaluate(hi)
        assert flo>0 and fhi<0
        for _ in range(12):
            candidate=lo-flo*(hi-lo)/(fhi-flo);fc=evaluate(candidate)
            if abs(fc)<1e-9:break
            if fc>0:lo,flo=candidate,fc
            else:hi,fhi=candidate,fc
        assert abs(fc)<1e-7
        roots.append({'grid':grid,'center_pressure_zero_radius_estimate':candidate,'pressure_at_estimate':fc})
    summary={'run_id':rid,'samples':values,'roots':roots,'run_seconds':time.perf_counter()-start,'method':'Stationary circular stage-1 BVP with actual upstream solvers; safeguarded secant interpolation in radius. Separates spatial onset error from time-step and contour-resolution delay.'}
    (run/'onset_summary.json').write_text(json.dumps(summary,indent=2));(ROOT/'runs/T03/latest_onset_probe.txt').write_text(rid)
    print(json.dumps(roots,indent=2))
if __name__=='__main__':main()

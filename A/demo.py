"""T01 demonstration using the same baseline source and configuration."""
import argparse,datetime,json,platform,subprocess,sys,time
from pathlib import Path

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'demo_runs'

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    choice=parser.add_mutually_exclusive_group(required=True)
    choice.add_argument('--prepare',action='store_true')
    choice.add_argument('--run',action='store_true')
    choice.add_argument('--show-latest',action='store_true')
    args=parser.parse_args()
    if args.show_latest:
        print((OUT/'latest.json').read_text(encoding='utf-8'));return
    stamp=datetime.datetime.now().astimezone().strftime('%Y%m%d_%H%M%S_%f')
    folder=OUT/stamp;folder.mkdir(parents=True)
    flags=['--build-only'] if args.prepare else ['--skip-build']
    command=[sys.executable,str(ROOT/'run.py'),'--case','T01_radial_64',*flags]
    started=time.perf_counter()
    with (folder/'execution.log').open('w',encoding='utf-8') as log:
        result=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT)
    record={'stage':'build' if args.prepare else 'new_run_and_postprocessing',
            'date':datetime.datetime.now().astimezone().isoformat(),
            'exit_code':result.returncode,'wall_seconds':time.perf_counter()-started,
            'platform':platform.platform(),'python':sys.version,'command':command,
            'record_type':'Local execution; expert witnessing is recorded separately.'}
    if args.run and result.returncode==0:
        base=ROOT/'baseline/code/runs/T01'
        run_id=(base/'latest_run.txt').read_text().strip()
        record['run_id']=run_id
        record['result_path']=str((base/run_id).relative_to(ROOT))
        record['summary']=json.loads((base/run_id/'summary.json').read_text(encoding='utf-8'))
    text=json.dumps(record,ensure_ascii=False,indent=2)
    (folder/'record.json').write_text(text,encoding='utf-8')
    (OUT/'latest.json').write_text(text,encoding='utf-8')
    print(text)
    raise SystemExit(result.returncode)

if __name__=='__main__':main()

"""Run one selected A-group case from this repository, from any directory."""
import argparse,json,os,subprocess,sys
from pathlib import Path

HERE=Path(__file__).resolve().parent
CODE=HERE/'baseline/code'

def main():
    cases=json.loads((CODE/'cases.json').read_text(encoding='utf-8'))
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--list',action='store_true')
    parser.add_argument('--case',choices=list(cases))
    parser.add_argument('--build-only',action='store_true')
    parser.add_argument('--skip-build',action='store_true')
    args=parser.parse_args()
    if args.list:
        for key,value in cases.items():print(key+' : '+value['description'])
        return
    if not args.case:parser.error('Select --case or --list')
    for group in ['t01','t02','t03','t04']:(CODE/group/'results').mkdir(parents=True,exist_ok=True)
    flags=['--case',args.case]
    if args.build_only:flags.append('--build-only')
    if args.skip_build:flags.append('--skip-build')
    if os.name=='nt' and args.case=='T03_onset_probe':
        path=(CODE/'frozen_runner.py').resolve().as_posix()
        path='/mnt/'+path[0].lower()+path[2:]
        command=['wsl.exe','-d','Ubuntu-24.04','--','python3',path,*flags]
    else:command=[sys.executable,str(CODE/'frozen_runner.py'),*flags]
    subprocess.run(command,check=True)

if __name__=='__main__':main()

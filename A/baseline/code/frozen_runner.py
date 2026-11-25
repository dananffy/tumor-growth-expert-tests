"""Cross-platform front end installed in a frozen packet's code/ directory."""
import argparse,json,os,subprocess,sys,shutil
from pathlib import Path
ROOT=Path(__file__).resolve().parent
def linux_path(p):
    p=p.resolve();s=p.as_posix()
    if p.drive:return '/mnt/'+p.drive[0].lower()+s[2:]
    return str(p)
def main():
    cases=json.loads((ROOT/'cases.json').read_text(encoding='utf-8'));ap=argparse.ArgumentParser()
    ap.add_argument('--case',choices=list(cases));ap.add_argument('--list',action='store_true');ap.add_argument('--build-only',action='store_true');ap.add_argument('--skip-build',action='store_true');a=ap.parse_args()
    if a.list:
        for key,c in cases.items():print(key,'—',c['description'])
        return
    if not a.case:ap.error('Choose --case or --list')
    c=cases[a.case]
    if c['kind']=='cpp':
        flags=['--case',a.case]+(['--build-only'] if a.build_only else [])+(['--skip-build'] if a.skip_build else [])
        if os.name=='nt':command=['wsl.exe','-d','Ubuntu-24.04','--','python3',linux_path(ROOT/'frozen_cpp_runner.py'),*flags]
        else:command=[sys.executable,str(ROOT/'frozen_cpp_runner.py'),*flags]
        subprocess.run(command,check=True)
        if a.build_only:return
        post=ROOT/('t01/postprocess.py' if c['group']=='T01' else 't03/postprocess.py')
        command=[sys.executable,str(post)]+(['--mode',c['mode']] if c['group']=='T03' else [])
        subprocess.run(command,check=True)
    else:
        if a.build_only or a.skip_build:ap.error('Build flags apply only to C++ cases')
        if c.get('prepare_working_source'):
            shutil.copytree(ROOT/c['prepare_working_source'],ROOT/'working/t03_emergence',dirs_exist_ok=True)
        subprocess.run([sys.executable,str(ROOT/c['script']),*c['arguments']],check=True)
if __name__=='__main__':main()

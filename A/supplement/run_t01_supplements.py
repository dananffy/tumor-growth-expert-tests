"""Windows orchestrator: isolated C++ computations and frozen postprocessing."""
from pathlib import Path
import json,os,shutil,subprocess,sys
ROOT=Path(__file__).resolve().parent;CODE=ROOT.parent
def main():
    env=os.environ.copy();env['WSL_UTF8']='1';env['PYTHONIOENCODING']='utf-8'
    p=ROOT/'t01_supplements_linux.py'
    if os.name=='nt':
        linux='/mnt/'+p.drive[0].lower()+p.as_posix()[2:]
        command=['wsl.exe','-d','Ubuntu-24.04','--','python3',linux]
    else:command=[sys.executable,str(p)]
    subprocess.run(command,env=env,check=True)
    assert (ROOT/'t01/postprocess.py').is_file()
    for row in json.loads((ROOT/'t01_supplement_runs.json').read_text()):
        subprocess.run([sys.executable,str(ROOT/'t01/postprocess.py'),'--run',row['run_id']],env=env,check=True)
if __name__=='__main__':main()

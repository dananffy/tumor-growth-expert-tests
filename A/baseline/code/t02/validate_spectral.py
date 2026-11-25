"""Independent circular-base linearized interface equations, mode ell=8.

From continuity/normal-flux matching, q=c1(R) per unit displacement is
-[(lambda-1)c0(R)+cB]/[sqrt(lambda) Iell'/Iell - Kell'/Kell].
The kinematic condition gives sigma=c0(R)-(ell+1)V0/R
+q/lambda*(sqrt(lambda) Iell'/Iell-ell/R), for G0=1.
These are evaluated independently of the nonlinear collocation solver.
"""
import sys,json,csv,hashlib,argparse
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
if sys.platform=='win32':sys.path.insert(0,str(ROOT/'.python_packages'))
import numpy as np
from scipy.special import iv,kv,ivp,kvp
from scipy.integrate import solve_ivp
def base(R,lam):
    s=np.sqrt(lam);c=100*iv(0,s*R)*kv(1,R)/(s*iv(1,s*R)*kv(0,R)+iv(0,s*R)*kv(1,R))
    v=c/s*iv(1,s*R)/iv(0,s*R)
    return c,v
def growth(R,lam,ell=8):
    c,v=base(R,lam);a=np.sqrt(lam)*ivp(ell,np.sqrt(lam)*R)/iv(ell,np.sqrt(lam)*R);b=kvp(ell,R)/kv(ell,R)
    q=-((lam-1)*c+100)/(a-b)
    return c-(ell+1)*v/R+q/lam*(a-ell/R)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run',required=True);a=ap.parse_args();run=ROOT/'runs/T02'/a.run
    summary=json.loads((run/'summary.json').read_text());validation=[]
    for s in summary['summaries']:
        lam=s['lambda'];d=np.load(run/f'lambda{lam:g}.npz');z=d['history'];times=z[:,0]
        ref=solve_ivp(lambda t,y:[base(y[0],lam)[1],growth(y[0],lam)*y[1]],(0,times[-1]),[4.,.01],t_eval=times,rtol=1e-11,atol=1e-13,method='DOP853')
        assert ref.success
        record={'lambda':lam,'initial_linear_growth_rate':growth(4,lam),'final_linear_growth_rate':growth(ref.y[0,-1],lam),'numerical_absolute_amplitude_ratio':s['absolute_amplitude_ratio'],'linear_reference_amplitude_ratio':float(ref.y[1,-1]/.01),'max_mean_radius_difference_vs_unperturbed_ODE':float(np.max(abs(z[:,1]-ref.y[0]))),'max_relative_amplitude_difference_vs_linear_reference':float(np.max(abs(z[:,2]-ref.y[1])/ref.y[1]))}
        validation.append(record)
        with (run/f'lambda{lam:g}_linear_comparison.csv').open('w',newline='') as f:
            w=csv.writer(f);w.writerow(['time','radius_numerical','radius_reference','amplitude_numerical','amplitude_linear_reference']);w.writerows(zip(times,z[:,1],ref.y[0],z[:,2],ref.y[1]))
    (run/'linear_validation.json').write_text(json.dumps(validation,indent=2))
    (run/'output_checksums.json').write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in run.iterdir() if p.is_file() and p.name!='output_checksums.json'},indent=2))
    out=ROOT/f't02/results/spectral_n{summary["summaries"][0]["nodes"]}_dr{summary["summaries"][0]["max_radial_increment"]:g}'
    import shutil
    shutil.copy2(run/'linear_validation.json',out/'linear_validation.json');print(json.dumps(validation,indent=2))
if __name__=='__main__':main()

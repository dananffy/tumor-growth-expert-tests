"""Planar traveling-front BVP + boundary-motion check (ZAMP Sec 2.2.1/2).

This solves the nutrient/pressure finite-difference systems and obtains front
speed from the numerical pressure derivative. Analytic profiles are references,
not substituted for the numerical elliptic solutions. The half-line traveling
profile is truncated to a 20-unit comoving interval on each applicable side.
"""
import sys,csv,json,hashlib,datetime
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
if sys.platform=='win32':sys.path.insert(0,str(ROOT/'.python_packages'))
import numpy as np
from scipy.linalg import solve_banded
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
def bvp(q,f,h,right):
    # Unknowns are nodes 0..N-1; node N is the prescribed right Dirichlet value.
    n=len(q);ab=np.zeros((3,n));ab[1]=2/h**2+q;ab[0,1:]=-1/h**2;ab[2,:-1]=-1/h**2
    ab[0,1]=-2/h**2 # reflected ghost at the left Neumann boundary
    rhs=f.copy();rhs[-1]+=right/h**2
    u=solve_banded((1,1),ab,rhs)
    return np.r_[u,right]
def main():
    rid=datetime.datetime.now().astimezone().strftime('%Y%m%d_%H%M%S_%f')+'_planar_wave';run=ROOT/'runs/T01'/rid;run.mkdir(parents=True)
    lam=.5;G=.2;cB=10.;L=20.;rows=[];profiles={}
    for regime in ['vitro','vivo']:
        A=cB if regime=='vitro' else cB/(1+np.sqrt(lam))
        exactspeed=G*A/np.sqrt(lam)
        for N in [100,200,400,800]:
            h=L/N
            if regime=='vitro':
                xi=np.linspace(-L,0,N+1);c=bvp(np.full(N,lam),np.zeros(N),h,cB)
                cref=cB*np.exp(np.sqrt(lam)*xi);cin=c
            else:
                xi=np.linspace(-L,L,2*N+1);q=np.where(xi[:-1]<0,lam,1.);source=np.where(xi[:-1]<0,0.,cB)
                q[N]=(lam+1)/2;source[N]=cB/2
                c=bvp(q,source,h,cB);cin=c[:N+1]
                cref=np.where(xi<=0,A*np.exp(np.sqrt(lam)*xi),cB-A*np.sqrt(lam)*np.exp(-xi))
            p=bvp(np.zeros(N),G*cin[:-1],h,0.);left=np.linspace(-L,0,N+1)
            pref=G*A/lam*(1-np.exp(np.sqrt(lam)*left))
            speed=(4*p[-2]-p[-3])/(2*h)
            rel=abs(speed-exactspeed)/exactspeed
            assert np.isfinite(c).all() and np.isfinite(p).all() and min(c)>0 and min(p)>=0
            times=np.linspace(0,.5,21);front=times*speed
            np.savez_compressed(run/f'{regime}_n{N}.npz',xi=xi,c=c,c_reference=cref,xi_pressure=left,p=p,p_reference=pref,time=times,front_numerical=front,front_reference=times*exactspeed)
            rows.append({'regime':regime,'intervals_inside':N,'dx':h,'speed_numerical':speed,'speed_reference':exactspeed,'relative_speed_error':rel,'nutrient_Linf_error':float(abs(c-cref).max()),'pressure_Linf_error':float(abs(p-pref).max()),'front_position_error_at_0_5':float(abs(front[-1]-.5*exactspeed))})
            if N==800:profiles[regime]=(xi,c,cref,left,p,pref,times,front,exactspeed)
    with (run/'wave_error_summary.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    fig,ax=plt.subplots(2,3,figsize=(13,7),layout='constrained')
    for i,(regime,data) in enumerate(profiles.items()):
        xi,c,cref,left,p,pref,times,front,exact=data
        ax[i,0].plot(xi,c,label='Numerical BVP');ax[i,0].plot(xi,cref,'--',label='Half-line reference');ax[i,0].set(title=f'{regime}: nutrient',xlabel='xi = x - front',xlim=(-8,6 if regime=='vivo' else 0))
        ax[i,1].plot(left,p,label='Numerical BVP');ax[i,1].plot(left,pref,'--',label='Half-line reference');ax[i,1].set(title='Pressure',xlabel='xi',xlim=(-8,0))
        ax[i,2].plot(times,front,label='Numerical front');ax[i,2].plot(times,times*exact,'--',label='Traveling-wave reference');ax[i,2].set(title='Front propagation',xlabel='Time',ylabel='Front position')
    for a in ax.flat:a.legend(fontsize=8);a.grid(alpha=.2)
    fig.suptitle('T01 planar traveling-front supplement | G0=0.2, cB=10, lambda=0.5')
    for ext in ['png','pdf']:fig.savefig(run/f'T01_traveling_wave.{ext}',dpi=170)
    plt.close(fig)
    s={'run_id':rid,'rows':rows,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'reference':'ZAMP Tumor boundary instability induced by nutrient consumption, Sec 2.2.1/2; half-line traveling plane.','scope':'One-dimensional planar free-boundary BVP with numerical boundary propagation; not a finite-density transient PDE or nonsymmetric-wave validation.','truncation_length':L,'parameters':{'lambda':lam,'G0':G,'c_B':cB}}
    (run/'wave_summary.json').write_text(json.dumps(s,indent=2));(run/'output_checksums.json').write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in run.iterdir() if p.is_file() and p.name!='output_checksums.json'},indent=2))
    import shutil
    out=ROOT/'t01/results/traveling_wave';out.mkdir(exist_ok=True)
    for name in ['wave_summary.json','wave_error_summary.csv','T01_traveling_wave.png','T01_traveling_wave.pdf']:shutil.copy2(run/name,out/name)
    print(json.dumps(rows,indent=2))
if __name__=='__main__':main()

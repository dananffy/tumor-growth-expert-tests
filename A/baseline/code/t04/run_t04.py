"""1D prediction-correction implementation, based on autophagy Appendix B.

Velocity correction follows u=-d_x p and p=gamma/(gamma-1)*n**(gamma-1).
The printed exponent in Eq. (52) is inconsistent with Eq. (6); this explicit
choice is recorded, rather than silently reproducing the inconsistent exponent.
"""
import sys,argparse,csv,datetime,hashlib,json,time,shutil
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
if sys.platform=='win32':sys.path.insert(0,str(ROOT/'.python_packages'))
import numpy as np
from scipy.linalg import solve_banded
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt

def pressure(n,gamma):return gamma/(gamma-1)*n**(gamma-1)
def velocity(n,gamma,h):return np.r_[0.,-np.diff(pressure(n,gamma))/h,0.]
def flux(n,u,h):
    d=np.diff(n)/h;s=np.zeros_like(n)
    same=d[:-1]*d[1:]>0
    s[1:-1]=np.where(same,np.sign(d[:-1])*np.minimum(np.abs(d[:-1]),np.abs(d[1:])),0)
    left=n[:-1]+h/2*s[:-1];right=n[1:]-h/2*s[1:]
    return np.r_[0.,np.maximum(u[1:-1],0)*left+np.minimum(u[1:-1],0)*right,0.]
def predict(n,u,S,gamma,h,dt):
    face=(n[:-1]+n[1:])/2;w=n**(gamma-2)
    main=1+dt*gamma/h**2*(w[:-1]+w[1:])*face
    lower=-dt*gamma/h**2*w[1:-1]*face[:-1]
    upper=-dt*gamma/h**2*w[1:-1]*face[1:]
    ab=np.zeros((3,len(face)));ab[1]=main;ab[0,1:]=upper;ab[2,:-1]=lower
    rhs=u[1:-1]-dt*gamma*np.diff(w*S)/h
    return np.r_[0.,solve_banded((1,1),ab,rhs),0.]
OCCUPANCY_CUTOFF=1e-12
def nutrient(n,n2,h,a):
    c=np.ones_like(n);occupied=np.flatnonzero(n>OCCUPANCY_CUTOFF)
    if len(occupied)<2:return c
    lo=max(1,occupied[0]);hi=min(len(n)-2,occupied[-1]);nn=n[lo:hi+1];N=len(nn)
    ab=np.zeros((3,N));ab[1]=2/h**2+nn;ab[0,1:]=-1/h**2;ab[2,:-1]=-1/h**2
    rhs=a*n2[lo:hi+1];rhs=rhs.copy();rhs[0]+=1/h**2;rhs[-1]+=1/h**2
    c[lo:hi+1]=solve_banded((1,1),ab,rhs);return c
def simulate(gamma,h,dtmax,end,autophagy):
    intervals=round(12/h)
    if intervals%2 or abs(intervals*h-12)>1e-10:raise ValueError('dx must divide [-6,6] into an even number of intervals')
    x=(np.arange(intervals+1)-intervals//2)*h;D=.3 if autophagy else 0.;a=.5 if autophagy else 0.
    mu=brentq(lambda q:D*q*(1-q)-q+(1-q),0,1) if autophagy else 1.
    # The exact support endpoints are zero, including floating-point values
    # within a roundoff-scaled distance of R. Tiny positive pressure at an
    # endpoint would otherwise become O(1) density when gamma is large.
    def pref(x,R):return np.maximum((1-(1-mu)*a)*(1-np.cosh(x)/np.cosh(R))+.5*(1-mu)*(a-D)*(R*R-x*x),0)*(np.abs(x)<R-32*np.finfo(float).eps*max(R,1))
    n=(np.maximum((gamma-1)/gamma*pref(x,1),0))**(1/(gamma-1))
    initial_symmetry=float(np.max(abs(n-n[::-1])))
    n1=mu*n;n2=(1-mu)*n;u=velocity(n,gamma,h);c=nutrient(n,n2,h,a)
    t=0.;steps=0;rejected=0;balance=0.;log=[];start=time.perf_counter()
    frames={0.:np.array([n1,n2,c,pressure(n,gamma)])}
    targets=[.25,.5,.75,end];targets=sorted(set(q for q in targets if q<=end))
    while t<end-1e-12:
        next_target=min(q for q in targets if q>t+1e-12)
        dt=min(dtmax,next_target-t)
        S=n1*c+n2*(c-D)
        for attempt in range(30):
            us=predict(n,u,S,gamma,h,dt)
            if np.max(np.abs(us))*dt/h>.45:dt*=.5;rejected+=1;continue
            z1=n1-dt*np.diff(flux(n1,us,h))/h;z2=n2-dt*np.diff(flux(n2,us,h))/h
            if autophagy:
                a11=1-dt*(c-1);a22=1-dt*(c-D-1);det=a11*a22-dt*dt
                q1=(a22*z1+dt*z2)/det;q2=(dt*z1+a11*z2)/det
            else:q1=z1/(1-dt*c);q2=np.zeros_like(q1)
            if not np.isfinite(q1+q2).all() or min(q1.min(),q2.min()) < -1e-12:
                dt*=.5;rejected+=1;continue
            break
        else:raise RuntimeError('Time step failed')
        q1=np.maximum(q1,0);q2=np.maximum(q2,0)
        defect=abs(np.sum(q1+q2-n)*h-dt*np.sum(c*q1+(c-D)*q2)*h)
        balance=max(balance,float(defect));n1,n2=q1,q2;n=n1+n2
        u=velocity(n,gamma,h);c=nutrient(n,n2,h,a);t+=dt;steps+=1
        assert n.max()<3 and min(c)>=-1e-10 and n[:3].max()<1e-12 and n[-3:].max()<1e-12
        radii=[max(abs(x[n>threshold]),default=0) for threshold in [1e-6,1e-4,.01]]
        log.append([t,dt,*radii,n.min(),n.max(),float(np.sum(n)*h),defect])
        if abs(t-next_target)<1e-10:frames[round(t,10)]=np.array([n1,n2,c,pressure(n,gamma)])
        if steps>200000:raise RuntimeError('Exceeded step limit')
    ref=solve_ivp(lambda t,R:(1-(1-mu)*a)*np.tanh(R)+(1-mu)*(a-D)*R,(0,end),[1.],rtol=1e-12,atol=1e-14,dense_output=True)
    R=float(ref.y[0,-1]);nref=(abs(x)<R).astype(float);p_ref=pref(x,R)
    massfrac=np.max(abs(n1[n>1e-6]/n[n>1e-6]-mu))
    return {'x':x,'n1':n1,'n2':n2,'n':n,'c':c,'p':pressure(n,gamma),'nref':nref,'pref':p_ref,'trajectory':np.array(log),'frames':frames}, {'gamma':gamma,'dx':h,'dt_max':dtmax,'end':end,'mu_star':mu,'R_reference':R,'density_L1_error':float(np.sum(abs(n-nref))*h),'pressure_L1_error':float(np.sum(abs(pressure(n,gamma)-p_ref))*h),'pressure_Linf_error':float(np.max(abs(pressure(n,gamma)-p_ref))),'radius_threshold_1e-4':float(log[-1][3]),'radius_abs_error_threshold_1e-4':abs(float(log[-1][3])-R),'steps':steps,'rejected':rejected,'max_mass_balance_defect':balance,'max_well_mixed_fraction_error':float(massfrac),'run_seconds':time.perf_counter()-start}
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--dx',type=float,default=.04);ap.add_argument('--dt',type=float,default=.002);ap.add_argument('--gammas',type=int,nargs='+',default=[2,8,16,32,64]);ap.add_argument('--mode',choices=['basic','autophagy'],default='autophagy');ap.add_argument('--end',type=float,default=1.);ap.add_argument('--occupancy-cutoff',type=float,default=1e-12);ap.add_argument('--label',default='');a=ap.parse_args()
    global OCCUPANCY_CUTOFF
    OCCUPANCY_CUTOFF=a.occupancy_cutoff
    rid=datetime.datetime.now().astimezone().strftime('%Y%m%d_%H%M%S_%f')+'_'+a.mode;run=ROOT/'runs/T04'/rid;run.mkdir(parents=True)
    (run/'source_snapshot').mkdir();shutil.copy2(__file__,run/'source_snapshot/run_t04.py')
    summaries=[];fig,ax=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    for g in a.gammas:
        data,s=simulate(g,a.dx,a.dt,a.end,a.mode=='autophagy')
        n0=data['frames'][0.][0]+data['frames'][0.][1]
        s['initial_density_mirror_defect']=float(abs(n0-n0[::-1]).max())
        s['final_density_mirror_defect']=float(abs(data['n']-data['n'][::-1]).max())
        s['final_pressure_mirror_defect']=float(abs(data['p']-data['p'][::-1]).max())
        s['occupancy_cutoff']=OCCUPANCY_CUTOFF
        assert s['initial_density_mirror_defect']<1e-12
        summaries.append(s)
        np.savez_compressed(run/f'gamma_{g}.npz',**{k:v for k,v in data.items() if k!='frames'})
        np.savez_compressed(run/f'frames_{g}.npz',**{f't_{t}':v for t,v in data['frames'].items()})
        with (run/f'trajectory_{g}.csv').open('w',newline='') as f:
            w=csv.writer(f);w.writerow(['time','dt','radius_1e-6','radius_1e-4','radius_1e-2','density_min','density_max','mass','mass_balance_defect']);w.writerows(data['trajectory'])
        ax[0].plot(data['x'],data['n'],label=f'gamma={g}');ax[1].plot(data['x'],data['p'],label=f'gamma={g}')
        print(json.dumps(s),flush=True)
    for a_,key,title in [(ax[0],'nref','Total density'),(ax[1],'pref','Pressure')]:
        a_.plot(data['x'],data[key],'k--',label='Free-boundary reference');a_.set(xlabel='x',title=title,xlim=(-3,3));a_.legend(fontsize=8);a_.grid(alpha=.2)
    fig.suptitle(f'T04 {a.mode} | t={a.end}, dx={a.dx}, max dt={a.dt}')
    for ext in ['png','pdf']:fig.savefig(run/f'T04_profiles.{ext}',dpi=180)
    plt.close(fig)
    with (run/'error_summary.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(summaries[0]));w.writeheader();w.writerows(summaries)
    result={'run_id':rid,'config':vars(a),'summaries':summaries,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'algorithm_version':'method-review-v2-symmetric-grid','velocity_correction':'Eq (6) p=gamma/(gamma-1)*n^(gamma-1); u=-grad(p). Printed exponent gamma in Eq (52) is not used.','nutrient':f'Quasi-static Eq (8) solved in the occupied interval, c_B=1 outside; occupancy cutoff {OCCUPANCY_CUTOFF:g}.','boundary':'Zero density flux at +/-6, checked distant from support.','scope':'New implementation and internal checks; method correspondence and expert approval pending.'}
    (run/'summary.json').write_text(json.dumps(result,indent=2));(run/'output_checksums.json').write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in run.iterdir() if p.is_file() and p.name!='output_checksums.json'},indent=2))
    suffix='' if a.occupancy_cutoff==1e-12 else f'_cutoff{a.occupancy_cutoff:g}'
    if a.label:suffix+='_'+a.label
    out=ROOT/f't04/results/{a.mode}_dx{a.dx:g}{suffix}';out.mkdir(parents=True,exist_ok=True)
    for name in ['summary.json','error_summary.csv','T04_profiles.png','T04_profiles.pdf']:shutil.copy2(run/name,out/name)
if __name__=='__main__':main()

"""T03 radial references use the upstream Bessel matching equations."""
import sys,json,csv,hashlib,re,argparse
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
if sys.platform=='win32':sys.path.insert(0,str(ROOT/'.python_packages'))
sys.path.insert(0,str(ROOT/'t01'))
from postprocess import spline_radius
import numpy as np
from scipy.special import iv,kv
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt

def coefficients(r,R):
    L=(iv(0,np.sqrt(.001)*r)*iv(1,r)-np.sqrt(.001)*iv(1,np.sqrt(.001)*r)*iv(0,r))/(iv(0,np.sqrt(.001)*r)*kv(1,r)+np.sqrt(.001)*iv(1,np.sqrt(.001)*r)*kv(0,r))
    b0=10/(iv(0,R)+L*kv(0,R));b1=b0*L
    A=r*(b0*iv(1,r)-b1*kv(1,r))-r*r*5/2
    B=b0*iv(0,r)+b1*kv(0,r)-r*r*5/4-A*np.log(r)
    return b0,b1,A,B
def condition(r,R):
    b0,b1,A,B=coefficients(r,R)
    return R*R*5/4+A*np.log(R)+B-(b0*iv(0,R)+b1*kv(0,R))
CRITICAL=brentq(lambda R:10*(1-1/iv(0,R))-5*R*R/4,2,3,xtol=1e-13)
def inner(R):
    if R<=CRITICAL:return 0.
    return brentq(condition,1e-6,.99*R,args=(R,),xtol=1e-12)
def rhs(t,y):
    R=float(y[0]);r=inner(R)
    if not r:return [10*iv(1,R)/iv(0,R)-5*R/2]
    b0,b1,A,B=coefficients(r,R)
    return [b0*iv(1,R)-b1*kv(1,R)-R*5/2-A/R]
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--mode',choices=['core','emergence'],default='core');ap.add_argument('--run');a=ap.parse_args()
    parent=ROOT/'runs/T03';run=parent/(a.run or (parent/f'latest_{a.mode}.txt').read_text().strip())
    meta=json.loads((run/'config_and_environment.json').read_text());cfg=meta['config'];mode=cfg['mode']
    if meta.get('assembly'):
        assert meta['assembly']['verified_overlap'] and meta['assembly']['continuation_exit_code']==0
    else:assert meta['exit_code']==0
    files=sorted((run/'result/points').glob('ctrl_points_1_*.txt'),key=lambda p:int(p.stem.split('_')[-1]))
    ids=[int(p.stem.split('_')[-1]) for p in files]
    offset=cfg.get('time_offset',0.)
    if mode=='core':times=offset+np.asarray(ids)*cfg['dt']
    else:
        logged=np.loadtxt(run/'result/time_log.txt');times=np.atleast_1d(logged)[ids]
    assert np.all(np.diff(times)>0) and abs(times[0]-offset)<1e-12
    Rinit=2.5 if mode=='core' else 2.2
    ref=solve_ivp(rhs,(0,offset+cfg['end']),[Rinit],t_eval=times,rtol=1e-10,atol=1e-12,method='DOP853',dense_output=True)
    assert ref.success
    Rout=ref.y[0];Rin=np.array([inner(R) for R in Rout])
    rn0=[];rn1=[]
    for i,p in zip(ids,files):
        outer=np.loadtxt(p);assert outer.shape==(64,2);assert np.isfinite(outer).all()
        rn1.append(spline_radius(outer))
        q=run/f'result/points/ctrl_points_0_{i}.txt'
        if not q.exists() or not q.stat().st_size:rn0.append(0.)
        else:rn0.append(spline_radius(np.loadtxt(q)))
    rn0=np.array(rn0);rn1=np.array(rn1)
    assert np.all(rn0<rn1)
    checked=0
    for p in (run/'result').rglob('*.txt'):
        if p.stat().st_size:
            assert np.isfinite(np.loadtxt(p)).all(),p
            checked+=1
    errors0=np.abs(rn0-Rin);errors1=np.abs(rn1-Rout)
    rows=list(zip(ids,times,rn0,Rin,errors0,rn1,Rout,errors1))
    with (run/'radii_comparison.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['step','time','inner_numerical','inner_reference','inner_abs_error','outer_numerical','outer_reference','outer_abs_error']);w.writerows(rows)
    fig,ax=plt.subplots(2,2,figsize=(10,8),layout='constrained')
    for a_,num,reference,title in [(ax[0,0],rn0,Rin,'Inner radius'),(ax[0,1],rn1,Rout,'Outer radius')]:
        a_.plot(times,reference,label='Bessel/ODE reference');a_.plot(times,num,'o--',ms=3,label='Numerical');a_.set(title=title,xlabel='Time',ylabel='Radius');a_.legend()
    ax[1,0].plot(times,errors0,label='Inner');ax[1,0].plot(times,errors1,label='Outer');ax[1,0].set(title='Absolute radius errors',xlabel='Time');ax[1,0].legend()
    theta=np.linspace(0,2*np.pi,300)
    for k in sorted(set([0,len(ids)//2,len(ids)-1])):
        ax[1,1].plot(Rout[k]*np.cos(theta),Rout[k]*np.sin(theta),lw=.8,label=f't={times[k]:.3f} outer ref')
        pts=np.loadtxt(files[k]);ax[1,1].plot(pts[:,0],pts[:,1],'.',ms=2)
        if rn0[k]>0:
            pts=np.loadtxt(run/f'result/points/ctrl_points_0_{ids[k]}.txt');ax[1,1].plot(np.r_[pts[:,0],pts[0,0]],np.r_[pts[:,1],pts[0,1]],lw=1)
    ax[1,1].set(title='Boundary snapshots',xlabel='x',ylabel='y',aspect='equal');ax[1,1].legend(fontsize=7)
    for a_ in ax.flat:a_.grid(alpha=.25)
    fig.suptitle(f'T03 {mode} | grid={cfg["grid"]}, dt={cfg["dt"]}')
    for ext in ['png','pdf']:fig.savefig(run/f'T03_overview.{ext}',dpi=170)
    plt.close(fig)
    # The emergence driver saves final boundaries but not final fields. Choose
    # the latest actual field snapshot, and label it with its matching time.
    field_ids=sorted(int(p.stem.split('_')[-1]) for p in (run/'result/solutions').glob('c_*.txt'))
    field_index=next(i for i in reversed(field_ids) if (run/f'result/solutions/p_{i}.txt').exists() and i in ids)
    field_time=float(times[ids.index(field_index)])
    Rfield=float(ref.sol(field_time)[0]);rfield=inner(Rfield)
    coords=np.linspace(-5,5,cfg['grid']+1);xx,yy=np.meshgrid(coords,coords,indexing='ij');rr=np.hypot(xx,yy)
    cnum=np.loadtxt(run/f'result/solutions/c_{field_index}.txt');pnum=np.loadtxt(run/f'result/solutions/p_{field_index}.txt')
    if rfield>0:
        b0,b1,A,B=coefficients(rfield,Rfield);rsafe=np.maximum(rr,rfield)
        cref=b0*iv(0,rsafe)+b1*kv(0,rsafe)
        ccore=(b0*iv(0,rfield)+b1*kv(0,rfield))*iv(0,np.sqrt(.001)*rr)/iv(0,np.sqrt(.001)*rfield)
        cref=np.where(rr<rfield,ccore,cref)
        pref=rsafe**2*5/4+A*np.log(rsafe)+B-(b0*iv(0,rsafe)+b1*kv(0,rsafe));pref=np.where(rr<rfield,0,pref)
    else:
        cref=10*iv(0,rr)/iv(0,Rfield);pref=10*(1-iv(0,rr)/iv(0,Rfield))-5*(Rfield**2-rr**2)/4
    actual_outer=np.loadtxt(run/f'result/points/ctrl_points_1_{field_index}.txt');outer_min=np.linalg.norm(actual_outer,axis=1).min()
    mask=rr<outer_min;common=rr<min(outer_min,Rfield)-10/cfg['grid']
    ff,aa=plt.subplots(2,2,figsize=(9,8),layout='constrained')
    for axis,values,title in zip(aa.flat,[np.where(mask,cnum,np.nan),np.where(mask,pnum,np.nan),np.where(common,abs(cnum-cref),np.nan),np.where(common,abs(pnum-pref),np.nan)],['Nutrient','Pressure','Nutrient error (common interior)','Pressure error (common interior)']):
        im=axis.pcolormesh(xx,yy,values,shading='auto',rasterized=True);axis.set(title=title,xlabel='x',ylabel='y',xlim=(-3.2,3.2),ylim=(-3.2,3.2),aspect='equal');ff.colorbar(im,ax=axis,shrink=.8)
    ff.suptitle(f'T03 {mode} fields | actual saved time={field_time:g}')
    for ext in ['png','pdf']:ff.savefig(run/f'T03_fields.{ext}',dpi=170)
    plt.close(ff)
    log=(run/'run.log').read_text();issues=[s for s in log.splitlines() if any(x in s.lower() for x in ['failed','warning','invalid','not inside'])]
    detected=re.findall(r'Necrotic core detected at time =\s*([\d.eE+-]+)',log)
    summary={'run_id':meta['run_id'],'mode':mode,'grid':cfg['grid'],'dt':cfg['dt'],'last_saved_time':float(times[-1]),'saved_states':len(times),'run_seconds':meta['run_seconds'],'inner_final':float(rn0[-1]),'inner_reference_final':float(Rin[-1]),'outer_final':float(rn1[-1]),'outer_reference_final':float(Rout[-1]),'inner_rmse':float(np.sqrt(np.mean(errors0**2))),'outer_rmse':float(np.sqrt(np.mean(errors1**2))),'inner_max_abs_error':float(errors0.max()),'outer_max_abs_error':float(errors1.max()),'critical_outer_radius':CRITICAL,'detected_time_logged':detected,'first_saved_nonzero_core_time':float(times[np.flatnonzero(rn0>0)[0]]) if np.any(rn0>0) else None,'checked_finite_files':checked,'solver_issue_lines':issues,'scope':'Numerical reproduction; no expert acceptance conclusion.'}
    (run/'summary.json').write_text(json.dumps(summary,indent=2))
    summary['last_saved_field_time']=field_time
    summary['execution_kind']='verified_checkpoint_assembly' if meta.get('assembly') else ('restart_segment' if offset else 'single_process')
    if Rinit<CRITICAL<float(Rout[-1]):
        summary['reference_core_nucleation_time']=float(brentq(lambda t:float(ref.sol(t)[0])-CRITICAL,0,cfg['end'],xtol=1e-12))
    summary['nutrient_max_abs_error_common_interior']=float(np.max(abs(cnum-cref)[common]))
    summary['pressure_max_abs_error_common_interior']=float(np.max(abs(pnum-pref)[common]))
    h=10/cfg['grid']
    lap_p=(pnum[2:,1:-1]+pnum[:-2,1:-1]+pnum[1:-1,2:]+pnum[1:-1,:-2]-4*pnum[1:-1,1:-1])/h**2
    lap_c=(cnum[2:,1:-1]+cnum[:-2,1:-1]+cnum[1:-1,2:]+cnum[1:-1,:-2]-4*cnum[1:-1,1:-1])/h**2
    checkmask=(rr<outer_min-2*h)[1:-1,1:-1]
    dual=-lap_p-(cnum[1:-1,1:-1]-5)
    pp=pnum[1:-1,1:-1]
    summary['obstacle_pressure_min_regular_interior']=float(pp[checkmask].min())
    summary['obstacle_dual_min_regular_interior']=float(dual[checkmask].min())
    summary['obstacle_complementarity_max_regular_interior']=float(abs(pp*dual)[checkmask].max())
    inner_file=run/f'result/points/ctrl_points_0_{field_index}.txt'
    inner_radius=spline_radius(np.loadtxt(inner_file)) if inner_file.exists() and inner_file.stat().st_size else 0.
    nutrient_mask=checkmask & (abs(rr[1:-1,1:-1]-inner_radius)>2*h)
    factor=np.where(rr[1:-1,1:-1]<inner_radius,.001,1.)
    nutrient_residual=-lap_c+factor*cnum[1:-1,1:-1]
    summary['nutrient_equation_max_abs_residual_away_from_interfaces']=float(abs(nutrient_residual)[nutrient_mask].max())
    (run/'summary.json').write_text(json.dumps(summary,indent=2))
    (run/'output_checksums.json').write_text(json.dumps({str(p.relative_to(run)):hashlib.sha256(p.read_bytes()).hexdigest() for p in run.rglob('*') if p.is_file() and p.name!='output_checksums.json'},indent=2))
    suffix='_restart' if offset else ''
    out=ROOT/'t03/results'/f'{mode}_n{cfg["grid"]}{suffix}';out.mkdir(parents=True,exist_ok=True)
    import shutil
    for name in ['T03_overview.png','T03_overview.pdf','T03_fields.png','T03_fields.pdf','summary.json','radii_comparison.csv']:shutil.copy2(run/name,out/name)
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()

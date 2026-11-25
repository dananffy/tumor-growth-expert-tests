"""Nonlinear PDE parameter/wavenumber scan with independent fitted growth rates."""
import os,sys,json,csv,hashlib,datetime,time,shutil,importlib.util
from pathlib import Path
os.environ['OPENBLAS_NUM_THREADS']='1';ROOT=Path(__file__).resolve().parent;CODE=ROOT.parent

import numpy as np
from scipy.optimize import brentq
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
LIB=ROOT/'sources/planar_front.py'
spec=importlib.util.spec_from_file_location('frozen_planar',LIB);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    rid=datetime.datetime.now().strftime('%Y%m%d_%H%M%S_%f')+'_parameter_scan';run=ROOT/'runs/T02'/rid;run.mkdir(parents=True)
    snap=run/'source_snapshot';snap.mkdir();shutil.copy2(__file__,snap/'t02_parameter_scan.py');shutil.copy2(LIB,snap/'planar_front.py')
    critical=brentq(lambda v:m.linear_rate(v,1),.5,40,xtol=1e-12)
    configs=[(v,1,49,1) for v in [.5,2.,5.,10.,20.,40.,critical*.9,critical*1.1]]+[(40.,k,49,1) for k in [2,3,4,5,6]]
    configs += [(critical*1.1,1,49,2),(40.,4,49,2),(critical*1.1,1,65,1),(40.,4,65,1)]
    rows=[]
    for idx,(lam,k,n,refine) in enumerate(configs):
        case=f'case_{idx:02d}';out=run/case;out.mkdir();eps=1e-4;theory=m.linear_rate(lam,k);end=min(1.,.5/max(abs(theory),1e-14));dtmax=min(end/200,.0025)/refine
        y=np.arange(n)*2*np.pi/n;f=eps*np.cos(k*y);X=0.;t=0.;hist=[];frames=[];maxres=0.;start=time.perf_counter()
        v0,*_=m.solve_front(np.zeros(n),lam);vexact=100/(np.sqrt(lam)*(1+np.sqrt(lam)));assert abs(v0-vexact).max()<1e-9
        while True:
            v,res,ab,hc=m.solve_front(f,lam);maxres=max(maxres,res);amp=2*abs(np.mean(f*np.exp(-1j*k*y)));hist.append([t,X,amp,eps*np.exp(theory*t),res]);frames.append(f.copy())
            if t>=end-1e-13:break
            dt=min(dtmax,end-t);vp,rp,*_=m.solve_front(f+dt*(v-v.mean()),lam);maxres=max(maxres,rp)
            X+=dt*.5*(v.mean()+vp.mean());f+=dt*.5*((v-v.mean())+(vp-vp.mean()));f-=f.mean();t+=dt
        z=np.array(hist);fit=np.polyfit(z[:,0],np.log(z[:,2]/eps),1)[0];off=m.off_grid(f,lam,ab,hc)
        assert np.isfinite(frames).all() and maxres<1e-8 and off<1e-7
        row={'case_id':case,'lambda':lam,'k':k,'nodes':n,'time_refinement':refine,'epsilon':eps,'end':end,'dt_max':dtmax,'sigma_theory':theory,'sigma_fit':float(fit),'absolute_growth_rate_error':float(abs(fit-theory)),'relative_growth_rate_error':float(abs(fit-theory)/abs(theory)),'amplitude_ratio':float(z[-1,2]/eps),'linear_amplitude_ratio':float(np.exp(theory*end)),'max_relative_amplitude_difference':float(np.max(abs(z[:,2]-z[:,3])/z[:,3])),'collocation_residual':maxres,'off_grid_residual':off,'run_seconds':time.perf_counter()-start}
        np.savez_compressed(out/'solution.npz',y=y,history=z,front_deviations=np.array(frames),nutrient_coefficients=ab,pressure_coefficients=hc)
        with (out/'history.csv').open('w',newline='') as fcsv:w=csv.writer(fcsv);w.writerow(['time','mean_front','amplitude','linear_reference','collocation_residual']);w.writerows(hist)
        (out/'summary.json').write_text(json.dumps(row,indent=2));rows.append(row);(run/'partial_summary.json').write_text(json.dumps(rows,indent=2));print(json.dumps(row),flush=True)
    result={'run_id':rid,'source_sha256':{'driver':sha(Path(__file__)),'frozen_planar_solver':sha(LIB)},'critical_lambda_k1_from_theory':critical,'rows':rows,'scope':'Small amplitude nonlinear periodic planar-interface PDE; transverse Fourier coefficient fit. Different observation windows documented, no general nonlinear stability theorem.'}
    (run/'summary.json').write_text(json.dumps(result,indent=2));(ROOT/'t02_scan_summary.json').write_text(json.dumps(result,indent=2))
    with (run/'growth_rates.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    fig,ax=plt.subplots(1,2,figsize=(10,3.3),layout='constrained');lams=np.linspace(.5,40,350)
    ax[0].plot(lams,[m.linear_rate(v,1) for v in lams],label='Linear theory, k=1');primary=[r for r in rows if r['k']==1 and r['nodes']==49 and r['time_refinement']==1];ax[0].plot([r['lambda'] for r in primary],[r['sigma_fit'] for r in primary],'o',label='Nonlinear PDE fit');ax[0].axhline(0,color='gray',lw=.7);ax[0].set(xlabel='lambda',ylabel='Growth rate');ax[0].legend(fontsize=8)
    primary=[r for r in rows if r['lambda']==40 and r['nodes']==49 and r['time_refinement']==1];ax[1].plot([r['k'] for r in primary],[r['sigma_theory'] for r in primary],label='Linear theory');ax[1].plot([r['k'] for r in primary],[r['sigma_fit'] for r in primary],'o',label='Nonlinear PDE fit');ax[1].axhline(0,color='gray',lw=.7);ax[1].set(xlabel='Wavenumber k (lambda=40)',ylabel='Growth rate');ax[1].legend(fontsize=8)
    for a in ax:a.grid(alpha=.2)
    for ext in ['png','pdf']:fig.savefig(run/f'T02_scan.{ext}',dpi=200)
    plt.close(fig);(run/'checksums.json').write_text(json.dumps({p.relative_to(run).as_posix():sha(p) for p in run.rglob('*') if p.is_file() and p.name!='checksums.json'},indent=2));print('COMPLETED',run)
if __name__=='__main__':main()

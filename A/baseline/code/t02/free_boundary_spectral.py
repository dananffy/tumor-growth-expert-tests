"""Independent 2D in-vivo free-boundary PDE solver for T02.

Modified-Helmholtz Fourier/Bessel expansions inside and outside the boundary
enforce nutrient continuity and normal-flux continuity. A harmonic expansion
enforces p=0, with p=-G0/lambda*c+harmonic inside. Boundary motion uses -dn p.
This is a new spectral collocation implementation, not the upstream ADI code.
"""
import os,sys,json,time,datetime,argparse,hashlib,csv
from pathlib import Path
os.environ['OPENBLAS_NUM_THREADS']='1'
ROOT=Path(__file__).resolve().parents[1]
if sys.platform=='win32':sys.path.insert(0,str(ROOT/'.python_packages'))
import numpy as np
from scipy.special import iv,kv,ivp,kvp
from scipy.integrate import solve_ivp
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt

def solve(r,lam,cB=100.):
    N=len(r);theta=np.arange(N)*2*np.pi/N;K=(N-1)//2;R=r.mean();s=np.sqrt(lam)
    modes=np.r_[0,np.arange(1,K+1),np.arange(1,K+1)]
    trig=np.column_stack([np.ones(N),np.cos(theta[:,None]*np.arange(1,K+1)),np.sin(theta[:,None]*np.arange(1,K+1))])
    trig_d=np.column_stack([np.zeros(N),-np.sin(theta[:,None]*np.arange(1,K+1))*np.arange(1,K+1),np.cos(theta[:,None]*np.arange(1,K+1))*np.arange(1,K+1)])
    kfreq=np.fft.fftfreq(N,1/N);rt=np.fft.ifft(1j*kfreq*np.fft.fft(r)).real
    nr=r/np.hypot(r,rt);nt=-rt/np.hypot(r,rt)
    I=iv(modes,s*r[:,None])/iv(modes,s*R);Ip=s*ivp(modes,s*r[:,None])/iv(modes,s*R)
    O=kv(modes,r[:,None])/kv(modes,R);Op=kvp(modes,r[:,None])/kv(modes,R)
    A=I*trig;B=O*trig
    DA=nr[:,None]*Ip*trig+nt[:,None]/r[:,None]*I*trig_d
    DB=nr[:,None]*Op*trig+nt[:,None]/r[:,None]*O*trig_d
    matrix=np.block([[A,-B],[DA,-DB]])
    target=np.r_[np.full(N,cB),np.zeros(N)]
    coeff=np.linalg.solve(matrix,target);inside=coeff[:N];outside=coeff[N:]
    cb=A@inside;dn_c=DA@inside
    H=(r[:,None]/R)**modes;DH=nr[:,None]*modes/r[:,None]*H*trig+nt[:,None]/r[:,None]*H*trig_d
    H=H*trig;hc=np.linalg.solve(H,cb/lam)
    v=dn_c/lam-DH@hc
    residual=max(np.max(abs(matrix@coeff-target))/cB,np.max(abs(H@hc-cb/lam))/max(cB/lam,1))
    return v/nr,float(residual),cb,v,inside,outside,hc

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--nodes',type=int,default=65);ap.add_argument('--dr',type=float,default=.01);a=ap.parse_args()
    N=a.nodes;assert N%2==1;theta=np.arange(N)*2*np.pi/N
    rid=datetime.datetime.now().astimezone().strftime('%Y%m%d_%H%M%S_%f')+f'_spectral_n{N}'
    run=ROOT/'runs/T02'/rid;run.mkdir(parents=True)
    summaries=[];fig,axes=plt.subplots(1,3,figsize=(14,4.2),layout='constrained')
    for lam in [.5,40.]:
        circle=np.full(N,4.);v,res,*_=solve(circle,lam)
        exact=100*kv(1,4)*iv(1,np.sqrt(lam)*4)/(lam*kv(0,4)*iv(1,np.sqrt(lam)*4)+np.sqrt(lam)*kv(1,4)*iv(0,np.sqrt(lam)*4))
        assert np.max(abs(v-exact))<1e-8
        r=4+.01*np.cos(8*theta);t=0.;rows=[];frames=[];start=time.perf_counter();maxres=res
        for step in range(500):
            speed,res,cb,v,ci,co,hc=solve(r,lam);maxres=max(maxres,res)
            amp=2*abs(np.mean(r*np.exp(-8j*theta)));rows.append([t,r.mean(),amp,amp/r.mean(),res]);frames.append(r.copy())
            if r.mean()>=5:break
            dt=min(a.dr/max(np.max(abs(speed)),1e-12),.02)
            trial=r+dt*speed;speed2,res2,*_=solve(trial,lam);maxres=max(maxres,res2)
            r=r+.5*dt*(speed+speed2);t+=dt
            if maxres>1e-8 or r.min()<3:raise RuntimeError('Collocation residual/geometry failure')
        else:raise RuntimeError('Spectral step limit')
        rows=np.array(rows);frames=np.array(frames)
        ref=solve_ivp(lambda t,y:100*kv(1,y)*iv(1,np.sqrt(lam)*y)/(lam*kv(0,y)*iv(1,np.sqrt(lam)*y)+np.sqrt(lam)*kv(1,y)*iv(0,np.sqrt(lam)*y)),(0,t),[4.],t_eval=rows[:,0],rtol=1e-11,atol=1e-12)
        np.savez_compressed(run/f'lambda{lam:g}.npz',theta=theta,radii=frames,history=rows,boundary_nutrient=cb,boundary_normal_velocity=v,inside_coefficients=ci,outside_coefficients=co,harmonic_coefficients=hc,reference_radius=ref.y[0])
        with (run/f'lambda{lam:g}_history.csv').open('w',newline='') as f:w=csv.writer(f);w.writerow(['time','mean_radius','mode8_absolute','mode8_relative','collocation_residual']);w.writerows(rows)
        axes[0].plot(rows[:,0],rows[:,1],label=f'lambda={lam:g}');axes[1].plot(rows[:,1],rows[:,2]/rows[0,2],label=f'lambda={lam:g}');axes[2].plot(rows[:,1],rows[:,3]/rows[0,3],label=f'lambda={lam:g}')
        summary={'lambda':lam,'nodes':N,'max_radial_increment':a.dr,'steps':len(rows)-1,'end_time':t,'initial_radius':4.,'last_mean_radius':float(r.mean()),'absolute_mode8_initial':float(rows[0,2]),'absolute_mode8_final':float(rows[-1,2]),'relative_mode8_initial':float(rows[0,3]),'relative_mode8_final':float(rows[-1,3]),'absolute_amplitude_ratio':float(rows[-1,2]/rows[0,2]),'relative_amplitude_ratio':float(rows[-1,3]/rows[0,3]),'radial_velocity_reference_max_error':float(np.max(abs(exact-solve(circle,lam)[0]))),'max_boundary_collocation_residual':maxres,'run_seconds':time.perf_counter()-start}
        summaries.append(summary);print(json.dumps(summary),flush=True)
    axes[0].set(title='Mean radius vs time',xlabel='Time',ylabel='Mean radius')
    axes[1].set(title='Absolute mode-8 growth',xlabel='Mean radius',ylabel='Amplitude / initial amplitude')
    axes[2].set(title='Relative mode-8 growth',xlabel='Mean radius',ylabel='Relative amplitude / initial value')
    for ax in axes:ax.legend();ax.grid(alpha=.2)
    fig.suptitle(f'T02 in-vivo free-boundary PDE | Fourier/Bessel collocation, N={N}, RK2')
    for ext in ['png','pdf']:fig.savefig(run/f'T02_spectral.{ext}',dpi=180)
    plt.close(fig)
    s={'run_id':rid,'summaries':summaries,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'model':'2D in-vivo elliptic nutrient inside/outside, p=-c/lambda+harmonic, V=-dn p, G0=1,cB=100; unbounded exterior.','scope':'New nonlinear free-boundary PDE discretization; original ADI reproduction is recorded separately.'}
    (run/'summary.json').write_text(json.dumps(s,indent=2))
    import shutil
    out=ROOT/f't02/results/spectral_n{N}_dr{a.dr:g}';out.mkdir(parents=True,exist_ok=True)
    for name in ['summary.json','T02_spectral.png','T02_spectral.pdf']:shutil.copy2(run/name,out/name)
if __name__=='__main__':main()

"""Nonlinear periodic planar-front PDE collocation (ZAMP Sections 2,4).

c_i uses exp(sqrt(lambda+k^2)*xi), c_o-cB uses
exp(-sqrt(1+k^2)*xi); pressure=-c_i/lambda+harmonic inside.
The mean front location is integrated independently of its periodic shape.
"""
import os,sys,argparse,json,datetime,csv,time,hashlib,shutil
from pathlib import Path
os.environ['OPENBLAS_NUM_THREADS']='1'
ROOT=Path(__file__).resolve().parents[1]
if sys.platform=='win32':sys.path.insert(0,str(ROOT/'.python_packages'))
import numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt

def basis(y,m):
    k=np.arange(1,m+1)
    return np.column_stack([np.ones(len(y)),np.cos(y[:,None]*k),np.sin(y[:,None]*k)]),np.column_stack([np.zeros(len(y)),-np.sin(y[:,None]*k)*k,np.cos(y[:,None]*k)*k])

def solve_front(f,lam):
    n=len(f);m=(n-1)//2;y=np.arange(n)*2*np.pi/n
    modes=np.r_[0,np.arange(1,m+1),np.arange(1,m+1)]
    trig,deriv=basis(y,m)
    fp=np.fft.ifft(1j*np.fft.fftfreq(n,1/n)*np.fft.fft(f)).real
    nx=1/np.hypot(1,fp);ny=-fp*nx
    alpha=np.sqrt(lam+modes*modes);beta=np.sqrt(1+modes*modes)
    ei=np.exp(f[:,None]*alpha);eo=np.exp(-f[:,None]*beta)
    A=ei*trig;B=eo*trig
    DA=ei*(nx[:,None]*alpha*trig+ny[:,None]*deriv)
    DB=eo*(-nx[:,None]*beta*trig+ny[:,None]*deriv)
    matrix=np.block([[A,-B],[DA,-DB]]);target=np.r_[np.full(n,100.),np.zeros(n)]
    ab=np.linalg.solve(matrix,target);ci=A@ab[:n]
    eh=np.exp(f[:,None]*modes);H=eh*trig;DH=eh*(nx[:,None]*modes*trig+ny[:,None]*deriv)
    hc=np.linalg.solve(H,ci/lam)
    velocity=(DA@ab[:n]/lam-DH@hc)/nx
    residual=max(float(abs(matrix@ab-target).max()/100),float(abs(H@hc-ci/lam).max()/max(100/lam,1)))
    return velocity,residual,ab,hc

def linear_rate(lam,k=1):
    s=np.sqrt(lam);A=100/(1+s);V=A/s;alpha=np.sqrt(lam+k*k);beta=np.sqrt(1+k*k)
    q=-100*s/(alpha+beta)
    return A-k*V+q/lam*(alpha-k)

def off_grid(f,lam,ab,hc):
    n=len(f);m=(n-1)//2;yd=np.arange(4*n)*2*np.pi/(4*n);freq=np.fft.fftfreq(n,1/n);weights=np.fft.fft(f)/n
    E=np.exp(1j*yd[:,None]*freq);fd=(E@weights).real;fp=(E@(1j*freq*weights)).real
    nx=1/np.hypot(1,fp);ny=-fp*nx;trig,deriv=basis(yd,m);k=np.r_[0,np.arange(1,m+1),np.arange(1,m+1)]
    alpha=np.sqrt(lam+k*k);beta=np.sqrt(1+k*k)
    ei=np.exp(fd[:,None]*alpha);eo=np.exp(-fd[:,None]*beta)
    ci=(ei*trig)@ab[:n];co=100+(eo*trig)@ab[n:]
    di=(ei*(nx[:,None]*alpha*trig+ny[:,None]*deriv))@ab[:n]
    do=(eo*(-nx[:,None]*beta*trig+ny[:,None]*deriv))@ab[n:]
    p=-ci/lam+(np.exp(fd[:,None]*k)*trig)@hc
    return max(float(abs(ci-co).max()/100),float(abs(di-do).max()/100),float(abs(p).max()/max(100/lam,1)))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--nodes',type=int,default=33);ap.add_argument('--refinement',type=int,default=1);ap.add_argument('--label',default='');a=ap.parse_args();n=a.nodes
    assert n%2==1 and n>=17
    y=np.arange(n)*2*np.pi/n;rid=datetime.datetime.now().astimezone().strftime('%Y%m%d_%H%M%S_%f')+f'_planar_n{n}'
    run=ROOT/'runs/T02'/rid;run.mkdir(parents=True);summary=[]
    (run/'source_snapshot').mkdir();shutil.copy2(__file__,run/'source_snapshot/planar_front.py')
    fig,axes=plt.subplots(2,2,figsize=(11,7),layout='constrained')
    for row,(lam,end,dtmax) in enumerate([(.5,.02,.00025/a.refinement),(40.,1.,.01/a.refinement)]):
        flat=np.zeros(n);flat_velocity,res,*_=solve_front(flat,lam);exact_speed=100/(np.sqrt(lam)*(1+np.sqrt(lam)))
        assert abs(flat_velocity-exact_speed).max()<1e-9
        f=.001*np.cos(y);X=0.;t=0.;rows=[];frames=[];start=time.perf_counter();max_res=res
        while True:
            v,res,ab,hc=solve_front(f,lam);max_res=max(max_res,res)
            amp=2*abs(np.mean(f*np.exp(-1j*y)));linear=.001*np.exp(linear_rate(lam)*t)
            rows.append([t,X,amp,linear,res]);frames.append(f.copy())
            if t>=end-1e-13:break
            step=min(dtmax,end-t)
            vp,resp,*_=solve_front(f+step*(v-v.mean()),lam)
            X+=step*.5*(v.mean()+vp.mean())
            f+=step*.5*((v-v.mean())+(vp-vp.mean()));f-=f.mean();t+=step
            max_res=max(max_res,resp)
        z=np.array(rows);frames=np.array(frames);off=off_grid(f,lam,ab,hc)
        assert np.isfinite(frames).all() and max_res<1e-9 and off<1e-8
        np.savez_compressed(run/f'lambda{lam:g}.npz',y=y,front_deviations=frames,history=z,inside_outside_coefficients=ab,harmonic_coefficients=hc)
        with (run/f'lambda{lam:g}.csv').open('w',newline='') as stream:w=csv.writer(stream);w.writerow(['time','mean_front_position','amplitude_numerical','amplitude_linear_reference','collocation_residual']);w.writerows(rows)
        s={'lambda':lam,'nodes':n,'dt_max':dtmax,'end':end,'initial_amplitude':.001,'final_amplitude':float(z[-1,2]),'amplitude_ratio':float(z[-1,2]/.001),'linear_growth_rate':linear_rate(lam),'linear_amplitude_ratio':float(np.exp(linear_rate(lam)*end)),'max_relative_error_vs_linear_amplitude':float((abs(z[:,2]-z[:,3])/z[:,3]).max()),'flat_front_speed_max_error':float(abs(flat_velocity-exact_speed).max()),'max_collocation_residual':max_res,'off_grid_interface_residual':off,'run_seconds':time.perf_counter()-start}
        summary.append(s);print(json.dumps(s),flush=True)
        axes[row,0].plot(z[:,0],z[:,2]/.001,label='Nonlinear PDE');axes[row,0].plot(z[:,0],z[:,3]/.001,'--',label='Linear reference');axes[row,0].set(title=f'lambda={lam:g}: planar mode-1 amplitude',xlabel='Time',ylabel='Amplitude / initial');axes[row,0].legend(fontsize=8)
        for k in [0,len(frames)//2,len(frames)-1]:axes[row,1].plot(y,frames[k]/.001,label=f't={z[k,0]:.4g}')
        axes[row,1].set(title='Actual front deviation',xlabel='Periodic transverse coordinate y',ylabel='(x - mean front) / initial amplitude');axes[row,1].legend(fontsize=8)
    for axis in axes.flat:axis.grid(alpha=.2)
    fig.suptitle(f'T02 planar traveling-front perturbation | nonlinear elliptic PDE, N={n}')
    for ext in ['png','pdf']:fig.savefig(run/f'T02_planar.{ext}',dpi=180)
    plt.close(fig)
    result={'run_id':rid,'config':vars(a),'summaries':summary,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'model':'Periodic planar in-vivo free boundary; nutrient and pressure solved on the perturbed interface; G0=1,cB=100, transverse wavenumber=1.','reference':'ZAMP traveling-wave stability; rate derived independently from flat-interface flux matching.','scope':'Small periodic perturbation, not a general nonlinear stability proof.'}
    (run/'summary.json').write_text(json.dumps(result,indent=2));(run/'output_checksums.json').write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in run.iterdir() if p.is_file() and p.name!='output_checksums.json'},indent=2))
    suffix=('_'+a.label) if a.label else ''
    out=ROOT/f't02/results/planar_n{n}{suffix}';out.mkdir(parents=True,exist_ok=True)
    for name in ['summary.json','T02_planar.png','T02_planar.pdf']:shutil.copy2(run/name,out/name)
if __name__=='__main__':main()

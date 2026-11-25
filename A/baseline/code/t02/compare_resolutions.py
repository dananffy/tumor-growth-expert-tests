import sys,json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
if sys.platform=='win32':sys.path.insert(0,str(ROOT/'.python_packages'))
import numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import iv,kv,ivp,kvp
def main():
    runs={}
    for p in sorted((ROOT/'runs/T02').glob('*_spectral_n*/summary.json')):
        s=json.loads(p.read_text());runs[s['summaries'][0]['nodes']]=p.parent
    coarse,fine=runs[65],runs[97];records=[]
    fig,axes=plt.subplots(2,2,figsize=(11,7),layout='constrained')
    for row,lam in enumerate([.5,40.]):
        dc=np.load(coarse/f'lambda{lam:g}.npz');df=np.load(fine/f'lambda{lam:g}.npz')
        rc=dc['history'][:,1];rf=df['history'][:,1];common=np.linspace(4,min(rc[-1],rf[-1]),200)
        ac=np.interp(common,rc,dc['history'][:,2]);af=np.interp(common,rf,df['history'][:,2])
        relative=float(np.max(abs(ac-af)/af))
        theta=df['theta'];r=df['radii'][-1];N=len(r);m=(N-1)//2
        td=np.arange(4*N)*2*np.pi/(4*N);freq=np.fft.fftfreq(N,1/N);weights=np.fft.fft(r)/N
        dense=np.exp(1j*td[:,None]*freq)
        rd=(dense@weights).real;rt=(dense@(1j*freq*weights)).real;R=r.mean();nr=rd/np.hypot(rd,rt);nt=-rt/np.hypot(rd,rt)
        modes=np.r_[0,np.arange(1,m+1),np.arange(1,m+1)]
        trig=np.column_stack([np.ones(len(td)),np.cos(td[:,None]*np.arange(1,m+1)),np.sin(td[:,None]*np.arange(1,m+1))])
        trigd=np.column_stack([np.zeros(len(td)),-np.sin(td[:,None]*np.arange(1,m+1))*np.arange(1,m+1),np.cos(td[:,None]*np.arange(1,m+1))*np.arange(1,m+1)])
        s=np.sqrt(lam);I=iv(modes,s*rd[:,None])/iv(modes,s*R);O=kv(modes,rd[:,None])/kv(modes,R)
        Ip=s*ivp(modes,s*rd[:,None])/iv(modes,s*R);Op=kvp(modes,rd[:,None])/kv(modes,R)
        ci=df['inside_coefficients'];co=df['outside_coefficients'];hc=df['harmonic_coefficients']
        cin=(I*trig)@ci;cout=100+(O*trig)@co
        dcin=(nr[:,None]*Ip*trig+nt[:,None]/rd[:,None]*I*trigd)@ci
        dcout=(nr[:,None]*Op*trig+nt[:,None]/rd[:,None]*O*trigd)@co
        pres=-cin/lam+((rd[:,None]/R)**modes*trig)@hc
        offgrid=max(float(abs(cin-cout).max()/100),float(abs(dcin-dcout).max()/100),float(abs(pres).max()/max(100/lam,1)))
        assert relative<1e-3 and offgrid<1e-8
        records.append({'lambda':lam,'coarse_nodes':65,'fine_nodes':97,'comparison_radius_interval':[4,float(common[-1])],'max_relative_mode8_difference_same_radius':relative,'fine_off_collocation_boundary_residual':offgrid})
        axes[row,0].plot(common,ac/.01,label='N=65, dr<=0.01');axes[row,0].plot(common,af/.01,'--',label='N=97, dr<=0.005');axes[row,0].set(title=f'lambda={lam:g}: amplitude vs mean radius',xlabel='Mean radius',ylabel='Mode-8 / initial value');axes[row,0].legend(fontsize=8)
        for k in [0,len(df['radii'])//2,len(df['radii'])-1]:
            snap=df['radii'][k];axes[row,1].plot(theta,(snap-snap.mean())/.01,label=f'R={snap.mean():.3f}')
        axes[row,1].set(title=f'lambda={lam:g}: actual boundary deviation',xlabel='Angle theta',ylabel='(r - mean r) / initial amplitude');axes[row,1].legend(fontsize=8)
    for ax in axes.flat:ax.grid(alpha=.2)
    out=ROOT/'t02/results';fig.suptitle('T02 spatial/time refinement and boundary perturbation profiles')
    for ext in ['png','pdf']:fig.savefig(out/f'T02_resolution_and_shapes.{ext}',dpi=180)
    plt.close(fig);(out/'resolution_validation.json').write_text(json.dumps(records,indent=2));print(json.dumps(records,indent=2))
if __name__=='__main__':main()

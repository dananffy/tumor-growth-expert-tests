"""Independent finite-gamma regression with constant nutrient c=1.

This isolates n_t=(n^gamma)_xx+n from the coupled nutrient problem and compares
the AP prediction/correction algorithm to an independent method-of-lines BDF
integrator. It does not replace validation of the coupled nutrient equations.
"""
import os,sys,json
from pathlib import Path
os.environ['OPENBLAS_NUM_THREADS']='1'
ROOT=Path(__file__).resolve().parents[1]
if sys.platform=='win32':sys.path.insert(0,str(ROOT/'.python_packages'))
sys.path.insert(0,str(ROOT/'t04'))
import run_t04 as ap
import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import diags
def main():
    # Process-local test reduction; main coupled-program sources are unchanged.
    ap.nutrient=lambda n,n2,h,a:np.ones_like(n)
    records=[]
    for gamma in [2,8]:
        errors=[]
        for h in [.08,.04,.02]:
            data,metrics=ap.simulate(gamma,h,.0005,.1,False)
            initial=data['frames'][0.][0]
            def rhs(t,n):
                z=np.maximum(n,0)**gamma
                return (np.r_[z[0],z[:-1]]-2*z+np.r_[z[1:],z[-1]])/h**2+n
            N=len(initial);pattern=diags([np.ones(N-1),np.ones(N),np.ones(N-1)],[-1,0,1],format='csc')
            reference=solve_ivp(rhs,(0,.1),initial,method='BDF',rtol=1e-9,atol=1e-11,jac_sparsity=pattern)
            assert reference.success
            nref=reference.y[:,-1];error=float(np.sum(abs(data['n']-nref))*h);errors.append(error)
            records.append({'gamma':gamma,'dx':h,'AP_dt_max':.0005,'end':.1,'density_L1_difference_vs_BDF':error,'reference_min_density':float(nref.min()),'reference_steps':len(reference.t)})
        assert errors[2]<errors[1]<errors[0],errors
    out=ROOT/'t04/results/finite_gamma_independent_regression.json';out.write_text(json.dumps({'model':'constant nutrient c=1, no autophagy: n_t=(n^gamma)_xx+n','records':records,'checks':'Independent BDF spatial discretization; AP/BDF difference decreases on grid refinement for both gamma values.'},indent=2));print(json.dumps(records,indent=2))
if __name__=='__main__':main()

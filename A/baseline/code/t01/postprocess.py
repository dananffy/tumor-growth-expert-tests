"""Independent ODE reference, data checks, CSV tables, and scientific figures.

Works with a normal scientific Python or the project-local Windows packages.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
local = ROOT / ".python_packages"
if local.exists() and sys.platform == "win32":
    sys.path.insert(0, str(local))
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.special import iv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def spline_radius(points):
    n = len(points)
    theta = np.linspace(0, 2*np.pi, n+1)
    curve = CubicSpline(theta, np.vstack([points, points[0]]), bc_type="periodic")
    # Match upstream: three-point Gauss integration on every spline interval.
    q, w = np.polynomial.legendre.leggauss(3)
    half = np.pi/n
    x = (theta[:-1, None]+theta[1:, None])/2 + half*q
    speed = np.linalg.norm(curve(x, 1), axis=2)
    return float(np.sum(speed*w)*half/(2*np.pi))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", help="Run ID, or defaults to latest_run.txt")
    args = parser.parse_args()
    parent = ROOT / "runs/T01"
    run = parent / (args.run or (parent/"latest_run.txt").read_text().strip())
    meta = json.loads((run/"config_and_environment.json").read_text())
    config = meta["config"]
    p = meta["parameters"]
    if meta["exit_code"] != 0:
        raise RuntimeError("Solver did not complete successfully")
    raw = np.genfromtxt(run/"radius_diagnostics.csv", delimiter=",", names=True)
    times = raw["time"]
    expected = int(round(config["final_time"]/config["dt"]))+1
    if len(times) != expected or not np.allclose(times,np.arange(expected)*config["dt"],rtol=0,atol=1e-13):
        raise RuntimeError("Saved steps do not match configured times")
    if not np.isfinite(np.column_stack([raw[n] for n in raw.dtype.names])).all():
        raise RuntimeError("Non-finite diagnostics")

    def rhs(t,r):
        z = np.sqrt(p["lambda"])*r
        return p["G0"]*p["c_B"]*iv(1,z)/(np.sqrt(p["lambda"])*iv(0,z))
    ref = solve_ivp(rhs,(0,config["final_time"]),[p["R0"]],t_eval=times,method="DOP853",rtol=1e-12,atol=1e-14,dense_output=True)
    tighter = solve_ivp(rhs,(0,config["final_time"]),[p["R0"]],t_eval=times,method="DOP853",rtol=3e-14,atol=1e-15)
    if not ref.success or not tighter.success:
        raise RuntimeError("Reference integration failed")
    error = np.abs(raw["radius_spline"]-ref.y[0])
    radius_crosscheck = []
    files = []
    for i in range(expected):
        for kind in ["ctrl_points","ctrl_nxny","bdry_points"]:
            f = run/f"result/points/{kind}_{i}.txt"
            a = np.loadtxt(f)
            if not np.isfinite(a).all(): raise RuntimeError(f"Non-finite data: {f}")
            files.append(f)
            if kind == "ctrl_points":
                if a.shape != (p["n_ctrl"],2): raise RuntimeError("Wrong boundary shape")
                radius_crosscheck.append(abs(spline_radius(a)-raw["radius_spline"][i]))
        for kind in ["c","p","v"]:
            f = run/f"result/functions/{kind}_{i}.txt"
            a = np.loadtxt(f)
            if not np.isfinite(a).all(): raise RuntimeError(f"Non-finite data: {f}")
            if kind in ("c","p") and a.shape != (config["grid"]+1,config["grid"]+1):
                raise RuntimeError("Wrong field shape")
            files.append(f)
    if max(radius_crosscheck)>1e-10: raise RuntimeError("Independent spline-radius crosscheck failed")
    with (run/"radius_comparison.csv").open("w", newline="") as f:
        writer=csv.writer(f)
        writer.writerow(["step","time","radius_numerical_spline","radius_reference_ode","absolute_error","relative_error","radius_min","radius_max","mean_velocity","reference_velocity"])
        for i,t in enumerate(times):
            writer.writerow([i,t,raw["radius_spline"][i],ref.y[0,i],error[i],error[i]/ref.y[0,i],raw["radius_min"][i],raw["radius_max"][i],raw["velocity_mean"][i],rhs(t,ref.y[0,i])])
    last=expected-1
    points=np.loadtxt(run/f"result/points/ctrl_points_{last}.txt")
    theta=np.linspace(0,2*np.pi,500)
    tdense=np.linspace(0,config["final_time"],300)
    fig,ax=plt.subplots(1,3,figsize=(14,4.2),layout="constrained")
    ax[0].plot(tdense,ref.sol(tdense)[0],label="ODE reference",lw=2)
    ax[0].plot(times,raw["radius_spline"],"o--",ms=4,label="Numerical spline radius")
    ax[0].set(xlabel="Time",ylabel="Radius",title="Radius evolution")
    ax[0].legend(fontsize=8)
    ax[1].plot(times,error,"o-",ms=4)
    ax[1].set(xlabel="Time",ylabel="Absolute radius error",title="Error against independent ODE")
    ax[1].ticklabel_format(axis="y",style="sci",scilimits=(0,0))
    ax[2].plot(ref.y[0,-1]*np.cos(theta),ref.y[0,-1]*np.sin(theta),label="Reference circle")
    ax[2].plot(np.r_[points[:,0],points[0,0]],np.r_[points[:,1],points[0,1]],"o--",ms=3,label="Control points")
    ax[2].set(xlabel="x",ylabel="y",title=f"Boundary at t={times[-1]:g}")
    ax[2].set_aspect("equal"); ax[2].legend(fontsize=8)
    for a in ax: a.grid(alpha=.25)
    fig.suptitle(f"T01 radial benchmark | {config['grid']} x {config['grid']}, dt={config['dt']:g}, Nctrl=32")
    for ext in ["png","pdf"]: fig.savefig(run/f"T01_overview.{ext}",dpi=180)
    plt.close(fig)

    # Fields are stored as [x-index, y-index]. Compare in the common interior
    # of the numerical boundary and the reference circle, excluding a grid-
    # width boundary band to avoid mixing domain-position and field errors.
    grid=config["grid"]
    coords=np.linspace(-1.5,1.5,grid+1)
    xx,yy=np.meshgrid(coords,coords,indexing="ij")
    rr=np.hypot(xx,yy)
    R=float(ref.y[0,-1])
    interior=rr < min(float(raw["radius_min"][-1]),R)-3/grid
    c=np.loadtxt(run/f"result/functions/c_{last}.txt")
    pressure=np.loadtxt(run/f"result/functions/p_{last}.txt")
    cexact=p["c_B"]*iv(0,np.sqrt(p["lambda"])*rr)/iv(0,np.sqrt(p["lambda"])*R)
    pexact=p["G0"]/p["lambda"]*(p["c_B"]-cexact)
    fig,axes=plt.subplots(2,2,figsize=(9,8),layout="constrained")
    mask=rr<=float(raw["radius_min"][-1])
    for a,values,title in zip(axes.flat,[np.where(mask,c,np.nan),np.where(mask,pressure,np.nan),np.where(interior,np.abs(c-cexact),np.nan),np.where(interior,np.abs(pressure-pexact),np.nan)],
                               ["Nutrient c (numerical)","Pressure p (numerical)","|c - reference| (common interior)","|p - reference| (common interior)"]):
        im=a.pcolormesh(xx,yy,values,shading="auto",rasterized=True)
        a.set(title=title,xlabel="x",ylabel="y",xlim=(-.75,.75),ylim=(-.75,.75),aspect="equal")
        fig.colorbar(im,ax=a,shrink=.8)
    fig.suptitle(f"T01 final fields | t={times[-1]:g}, grid={grid}")
    for ext in ["png","pdf"]: fig.savefig(run/f"T01_fields.{ext}",dpi=180)
    plt.close(fig)
    log=(run/"run.log").read_text()
    if "failed to converge" in log or "Invalid velocity" in log or "Could not open file" in log:
        raise RuntimeError("Solver log contains a failure")
    summary={"run_id":meta["run_id"],"grid":grid,"dt":config["dt"],"final_time":float(times[-1]),
             "saved_states":expected,"run_seconds":meta["run_seconds"],"build_seconds":meta["build_seconds"],
             "radius_reference_final":R,"radius_numerical_final":float(raw["radius_spline"][-1]),
             "absolute_radius_error_final":float(error[-1]),"relative_radius_error_final":float(error[-1]/R),
             "max_radius_error_over_time":float(error.max()),
             "ode_tolerance_crosscheck_max_difference":float(np.max(np.abs(ref.y-tighter.y))),
             "spline_radius_independent_crosscheck_max_difference":max(radius_crosscheck),
             "final_radial_spread_max_minus_min":float(raw["radius_max"][-1]-raw["radius_min"][-1]),
             "nutrient_max_abs_error_common_interior":float(np.max(np.abs(c-cexact)[interior])),
             "pressure_max_abs_error_common_interior":float(np.max(np.abs(pressure-pexact)[interior])),
             "all_saved_numeric_arrays_finite":True,"checked_numeric_files":len(files),
             "log_warning_lines":[line for line in log.splitlines() if "warning" in line.lower()],
             "postprocess_versions":{"python":sys.version,"numpy":np.__version__,"scipy":scipy.__version__,"matplotlib":matplotlib.__version__},
             "scope":"Internal radial benchmark reproduction only; no expert acceptance conclusion or separate time/space convergence order."}
    (run/"summary.json").write_text(json.dumps(summary,indent=2))
    with (run/"summary.csv").open("w",newline="") as f:
        writer=csv.writer(f); writer.writerow(["quantity","value"])
        writer.writerows((k,v) for k,v in summary.items() if not isinstance(v,(dict,list)))
    manifest={str(f.relative_to(run)):hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(run.rglob("*")) if f.is_file() and f.name!="output_checksums.json"}
    (run/"output_checksums.json").write_text(json.dumps(manifest,indent=2))
    print(json.dumps(summary,indent=2))


if __name__=="__main__":
    main()

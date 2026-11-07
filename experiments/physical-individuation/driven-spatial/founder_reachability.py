"""Declared zero-L founder-budget assay of the separate open ODE benchmark.

This is not the new spatial experiment. It addresses an explicit missing
reachability test without changing the benchmark's parameters.
"""
from __future__ import annotations

from dataclasses import asdict
from hashlib import sha256
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parent
PRIOR = ROOT.parent/"renewal-turnover"
sys.path.insert(0, str(PRIOR))
from verify_threshold import M, ALPHA, ETA, AFFINITY, vector


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def material_free_energy(r, ell):
    d = M-r-ell
    return 20*(r+ell)+sum(c*(math.log(c)-1) for c in (d, r, ell) if c > 0)


def main():
    protocol = dict(founder_R=[.01, .03, .05, .10, .20, .50, .90], initial_L=0.,
                    M=M, alpha=ALPHA, eta=ETA, fuel_affinity=AFFINITY,
                    arms=[dict(produce=p, protect=b) for p, b in [(True, True), (False, True), (True, False), (False, False)]],
                    horizon=500., sample_dt=.5, terminal_high_tolerance=1e-7,
                    benchmark_source_sha256=digest(PRIOR/"verify_threshold.py"),
                    interval_certificate_sha256=digest(PRIOR/"threshold-interval-certificate.json"),
                    runner_sha256=digest(Path(__file__)),
                    scope="Numerical zero-L formation examples in the supplied open wellmixed benchmark, without founder-budget optimization or spatial-individuality inference.")
    protocol_path=ROOT/"founder-protocol.json"
    encoded=json.dumps(protocol, indent=2)+"\n"
    if protocol_path.exists():
        assert protocol_path.read_text()==encoded
    else:
        protocol_path.write_text(encoded)
    assert not (ROOT/"founder-results.json").exists(), "refusing to replace evidence"
    high_data=json.loads((PRIOR/"threshold-verification.json").read_text())["reversible_equilibria"][2]
    high=np.array([high_data[k] for k in ("R", "L")])
    rows=[]
    times=np.arange(0.,500.1,.5)
    for founder in protocol["founder_R"]:
        for arm in protocol["arms"]:
            def rhs(t, state):
                r, ell = state[:2]
                d=M-r-ell
                chemical_power=AFFINITY*(r*d-ETA*r*r+(r*d-ETA*r*ell if arm["produce"] else 0.))
                return np.r_[vector(np.array([r,ell]), **arm), chemical_power]
            solved=solve_ivp(rhs,(0.,500.),[founder,0.,0.],t_eval=times,rtol=1e-10,atol=1e-12)
            assert solved.success
            material=solved.y[:2]
            assert material.min()>=-1e-11 and material.sum(axis=0).max()<=M+1e-11
            final=material[:,-1]
            delta_g=material_free_energy(*final)-material_free_energy(founder,0.)
            work=float(solved.y[2,-1])
            assert work-delta_g>=-1e-9
            rows.append(dict(founder_R=founder,initial_L=0.,**arm,final_R=float(final[0]),final_L=float(final[1]),
                             reaches_high_state=bool(np.max(abs(final-high))<protocol["terminal_high_tolerance"]),
                             initial_material_free_energy_above_all_D=material_free_energy(founder,0.)-material_free_energy(0.,0.),
                             supplied_fuel_work=work,material_free_energy_change=delta_g,
                             integrated_entropy_production_from_balance=work-delta_g,
                             samples=[dict(time=float(t),R=float(v[0]),L=float(v[1]),chemical_work=float(v[2]))
                                      for t,v in zip(times,solved.y.T)]))
    result=dict(protocol=protocol,protocol_sha256=digest(protocol_path),runs=rows,
                preparations=7,physical_ODE_trajectories=len(rows),
                initialization_note="Templates are supplied. Their free energy difference is counted, but their actual construction protocol/work is not modeled; this is not an abiogenesis derivation.",
                entropy_note="Integrated entropy production is computed through the exact free-energy balance. Direct entropy-rate quadrature at zero L is not claimed.",
                full_goal_achieved=False)
    assert digest(PRIOR/"verify_threshold.py")==protocol["benchmark_source_sha256"]
    (ROOT/"founder-results.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps([{k:v for k,v in row.items() if k!="samples"} for row in rows],indent=2))


if __name__=="__main__":
    main()

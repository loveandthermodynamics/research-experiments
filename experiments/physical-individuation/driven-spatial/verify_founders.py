"""Cross-check founder trajectories with reconstructed equations and DOP853."""
from hashlib import sha256
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

ROOT=Path(__file__).resolve().parent


def g(r,ell,m):
    d=m-r-ell
    return 20*(r+ell)+sum(x*(math.log(x)-1) for x in (r,ell,d) if x>0)


def main():
    source=ROOT/"founder-results.json"
    data=json.loads(source.read_text()); p=data["protocol"]
    assert sha256((ROOT/"founder-protocol.json").read_bytes()).hexdigest()==data["protocol_sha256"]
    assert sha256((ROOT/"founder_reachability.py").read_bytes()).hexdigest()==p["runner_sha256"]
    expected={(r,produce,protect) for r in p["founder_R"] for produce in (True,False) for protect in (True,False)}
    assert {(r["founder_R"],r["produce"],r["protect"]) for r in data["runs"]}==expected
    assert len(data["runs"])==len(expected)==28
    maxstate=maxwork=0.; classifications=[]
    for row in data["runs"]:
        def rhs(t,y):
            r,ell,work=y; d=p["M"]-r-ell
            h=math.exp(-p["alpha"]*ell) if row["protect"] else 1.
            copy=r*d-p["eta"]*r*r
            production=(r*d-p["eta"]*r*ell) if row["produce"] else 0.
            return [copy-h*(r-p["eta"]*d),production-ell+p["eta"]*d,
                    p["fuel_affinity"]*(copy+production)]
        times=np.array([s["time"] for s in row["samples"]])
        sol=solve_ivp(rhs,(0.,p["horizon"]),[row["founder_R"],0.,0.],method="DOP853",
                      t_eval=times,rtol=1e-12,atol=1e-14)
        assert sol.success
        saved=np.array([[s["R"],s["L"],s["chemical_work"]] for s in row["samples"]]).T
        state_error=float(np.max(abs(sol.y[:2]-saved[:2])))
        work_error=float(np.max(abs(sol.y[2]-saved[2])/(1+abs(saved[2]))))
        assert state_error<1e-8 and work_error<1e-8
        maxstate=max(maxstate,state_error); maxwork=max(maxwork,work_error)
        delta=g(row["final_R"],row["final_L"],p["M"])-g(row["founder_R"],0.,p["M"])
        assert abs(delta-row["material_free_energy_change"])<1e-12
        assert abs(row["supplied_fuel_work"]-delta-row["integrated_entropy_production_from_balance"])<1e-9
        assert abs(g(row["founder_R"],0.,p["M"])-g(0.,0.,p["M"])-row["initial_material_free_energy_above_all_D"])<1e-12
        classifications.append(dict(founder_R=row["founder_R"],produce=row["produce"],protect=row["protect"],
                                    reaches_high_state=row["reaches_high_state"]))
    report=dict(checks_passed=True,source_sha256=sha256(source.read_bytes()).hexdigest(),trajectories=28,
                initial_L=0,high_full_founders=[r["founder_R"] for r in data["runs"] if r["produce"] and r["protect"] and r["reaches_high_state"]],
                high_ablations=sum(r["reaches_high_state"] for r in data["runs"] if not(r["produce"] and r["protect"])),
                maximum_state_difference=maxstate,maximum_scaled_work_difference=maxwork,
                method="Directly reconstructed currents; DOP853 versus recorded RK45; all sampled states and work compared",
                scope="Numerical reachability examples, not rigorous threshold-founder optimization, spatial formation or a biodiversity proof.")
    (ROOT/"founder-verification.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()

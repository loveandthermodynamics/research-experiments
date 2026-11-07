"""Deterministic reciprocal-maintenance CRN with a chemical-work ledger.

This is a concentration/large-copy limit, not a four-molecule CTMC. The supplied
reaction grammar is cooperative catalytic maintenance, not a derived membrane,
group reproduction mechanism, or proof of individuality. Fuel and waste are
chemostatted; their net chemical work is explicitly counted.

Run: python3 -W error agent_closure.py
Writes only agent-closure-results.json beside this owned script.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from pathlib import Path
import json

import numpy as np
from scipy.integrate import solve_ivp


@dataclass(frozen=True)
class Parameters:
    material: float = 4.0
    active_energy: float = 4.0
    fuel_energy: float = 6.0
    thermal: float = 0.5
    copy_coefficient: float = 1.0
    fuel: float = 1.0
    waste: float = 1.0
    decay: float = 0.2
    channel_x: float = 1.0
    channel_y: float = 1.0

    @property
    def reverse_coefficient(self) -> float:
        return self.copy_coefficient * np.exp(
            -(self.fuel_energy - self.active_energy) / self.thermal
        )

    @property
    def source(self) -> float:
        return self.decay * np.exp(-self.active_energy / self.thermal)

    @property
    def chemical_gap(self) -> float:
        return self.fuel_energy + self.thermal * np.log(self.fuel / self.waste)


def fluxes(p: Parameters, xy):
    x, y = xy
    d = p.material - x - y
    assert min(x, y, d) > 0
    forward = np.array([
        p.channel_x * p.copy_coefficient * p.fuel * x * y * d,
        p.channel_y * p.copy_coefficient * p.fuel * x * y * d,
        p.decay * x,
        p.decay * y,
    ])
    reverse = np.array([
        p.channel_x * p.reverse_coefficient * p.waste * x * x * y,
        p.channel_y * p.reverse_coefficient * p.waste * x * y * y,
        p.source * d,
        p.source * d,
    ])
    return forward, reverse


def dynamics(p: Parameters, xy):
    forward, reverse = fluxes(p, xy)
    net = forward - reverse
    return np.array([net[0] - net[2], net[1] - net[3]])


def jacobian(p: Parameters, xy):
    x, y = xy
    d = p.material - x - y
    s = p.copy_coefficient * p.fuel
    r = p.reverse_coefficient * p.waste
    return np.array([
        [p.channel_x * (s * y * (d-x) - 2*r*x*y) - p.decay - p.source,
         p.channel_x * (s * x * (d-y) - r*x*x) - p.source],
        [p.channel_y * (s * y * (d-x) - r*y*y) - p.source,
         p.channel_y * (s * x * (d-y) - 2*r*x*y) - p.decay - p.source],
    ])


def free_energy(p: Parameters, xy):
    x, y = xy
    concentrations = np.array([x, y, p.material-x-y])
    assert concentrations.min() > 0
    return float(p.active_energy * (x+y) + p.thermal * np.sum(
        concentrations * (np.log(concentrations)-1)
    ))


def equilibrium_material(p: Parameters):
    ratio = np.exp(-p.active_energy/p.thermal)
    u = ratio * p.material / (1+2*ratio)
    return np.array([u, u])


def physical_rates(p: Parameters, xy):
    x, y = xy
    d = p.material-x-y
    forward, reverse = fluxes(p, xy)
    net = forward-reverse
    gradient = np.array([
        p.active_energy+p.thermal*np.log(x/d),
        p.active_energy+p.thermal*np.log(y/d),
    ])
    dx = np.array([net[0]-net[2], net[1]-net[3]])
    active = forward > 0
    dissipation = float(p.thermal*np.sum(
        net[active]*np.log(forward[active]/reverse[active])
    ))
    work = float(p.chemical_gap*np.sum(net[:2]))
    gdot = float(np.sum(gradient*dx))
    heat = float(p.fuel_energy*np.sum(net[:2])-p.active_energy*np.sum(dx))
    return work, dissipation, gdot, heat


def symmetric_fixed_points(p: Parameters):
    assert p.channel_x == p.channel_y == 1
    s = p.copy_coefficient*p.fuel
    r = p.reverse_coefficient*p.waste
    a = s*p.material
    b = 2*s+r
    roots = np.roots([b, -a, p.decay+2*p.source, -p.source*p.material])
    values = sorted(float(z.real) for z in roots
                    if abs(z.imag)<1e-10 and 0<z.real<p.material/2)
    result = []
    for u in values:
        parallel = 2*a*u-3*b*u*u-p.decay-2*p.source
        transverse = -p.decay-r*u*u
        numerical = np.sort(np.linalg.eigvals(jacobian(p, [u, u])).real)
        assert np.max(np.abs(numerical-np.sort([parallel, transverse]))) < 1e-10
        assert np.max(np.abs(dynamics(p, [u,u]))) < 1e-10
        work, diss, gdot, _ = physical_rates(p, [u,u])
        assert abs(work-diss)<1e-10 and abs(gdot)<1e-10
        result.append({
            "u": u, "x": u, "y": u, "D": p.material-2*u,
            "parallel_eigenvalue": float(parallel),
            "transverse_eigenvalue": float(transverse),
            "classification": "stable" if parallel<0 else "saddle",
            "chemical_power": work,
            "dissipation_power": diss,
            "material_availability": free_energy(p,[u,u])-free_energy(p,equilibrium_material(p)),
        })
    return result


def trajectory(p: Parameters, start, horizon=120.0):
    # Augmented integrals: chemical work, theta*sigma, standard-energy heat,
    # then all four forward fluxes and all four reverse fluxes.
    def ode(_, state):
        xy = state[:2]
        forward, reverse = fluxes(p, xy)
        work, diss, _, heat = physical_rates(p, xy)
        return np.r_[dynamics(p,xy), work, diss, heat, forward, reverse]

    initial = np.r_[start,np.zeros(11)]
    times = np.linspace(0,horizon,1201)
    sol = solve_ivp(ode,[0,horizon],initial,method="DOP853",t_eval=times,
                    rtol=2e-10,atol=2e-12,max_step=.5)
    assert sol.success
    xy = sol.y[:2].T
    d = p.material-np.sum(xy,axis=1)
    assert min(float(xy.min()),float(d.min())) > 0
    g = np.array([free_energy(p,z) for z in xy])
    gref = free_energy(p,equilibrium_material(p))
    assert (g-gref).min()>-1e-10
    work,diss,heat = map(float,sol.y[2:5,-1])
    plus = sol.y[5:9,-1]
    minus = sol.y[9:13,-1]
    net = plus-minus
    final = xy[-1]
    fuel_used = float(np.sum(net[:2]))
    chemical_balance_error = abs(float(g[-1]-g[0])-work+diss)
    energy_balance_error = abs(p.active_energy*float(np.sum(final-start))
                               -p.fuel_energy*fuel_used+heat)
    x_balance_error = abs(float(final[0]-start[0])-float(net[0]-net[2]))
    y_balance_error = abs(float(final[1]-start[1])-float(net[1]-net[3]))
    assert chemical_balance_error<2e-7
    assert energy_balance_error<2e-8
    assert max(x_balance_error,y_balance_error)<2e-9
    work_end,diss_end,gdot_end,_ = physical_rates(p,final)
    return {
        "parameters":asdict(p),"initial":list(map(float,start)),
        "final":final.tolist(),"horizon":horizon,
        "final_derivative_norm":float(np.linalg.norm(dynamics(p,final))),
        "initial_material_availability":float(g[0]-gref),
        "final_material_availability":float(g[-1]-gref),
        "net_fuel_consumed_from_chemostat":fuel_used,
        "integrated_chemical_work":work,"integrated_dissipation":diss,
        "integrated_heat_to_bath":heat,
        "final_chemical_power":work_end,"final_dissipation_power":diss_end,
        "final_free_energy_derivative":gdot_end,
        "G_work_balance_error":chemical_balance_error,
        "standard_energy_heat_balance_error":energy_balance_error,
        "species_flux_balance_error":max(x_balance_error,y_balance_error),
        "minimum_concentration":min(float(xy.min()),float(d.min())),
        "largest_free_energy_increase_on_grid":float(np.max(np.diff(g))),
        "trajectory":{
            "time":times[::30].tolist(),"x":xy[::30,0].tolist(),
            "y":xy[::30,1].tolist(),"D":d[::30].tolist(),
            "material_availability":(g[::30]-gref).tolist(),
        },
    }


def audit_random_states(p: Parameters, count=500):
    rng = np.random.default_rng(3801)
    states = rng.dirichlet([1.5,1.5,1.5],size=count)*p.material
    max_ldb = max_ledger = max_difference = max_jacobian = 0.0
    min_diss = np.inf
    for x,y,d in states:
        forward,reverse = fluxes(p,[x,y])
        mux = p.active_energy+p.thermal*np.log(x)
        muy = p.active_energy+p.thermal*np.log(y)
        mud = p.thermal*np.log(d)
        affinities = np.array([
            mud+p.chemical_gap-mux,mud+p.chemical_gap-muy,mux-mud,muy-mud,
        ])
        max_ldb = max(max_ldb,float(np.max(np.abs(
            np.log(forward/reverse)-affinities/p.thermal))))
        work,diss,gdot,_ = physical_rates(p,[x,y])
        max_ledger = max(max_ledger,abs(gdot-work+diss))
        min_diss = min(min_diss,diss)
        derivative = dynamics(p,[x,y])
        predicted = -(p.decay+p.reverse_coefficient*p.waste*x*y)*(x-y)
        max_difference = max(max_difference,abs(float(derivative[0]-derivative[1])-predicted))
        if min(x,y,d)>1e-4:
            step=1e-6
            numeric=np.column_stack([
                (dynamics(p,np.array([x,y])+step*np.eye(2)[j])
                 -dynamics(p,np.array([x,y])-step*np.eye(2)[j]))/(2*step)
                for j in range(2)
            ])
            max_jacobian=max(max_jacobian,float(np.max(np.abs(numeric-jacobian(p,[x,y])))))
    assert max_ldb<1e-11 and max_ledger<1e-10
    assert max_difference<1e-12 and max_jacobian<2e-8
    assert min_diss>=0
    return {"states":count,"max_local_detailed_balance_log_error":max_ldb,
            "max_G_work_identity_error":max_ledger,
            "max_exact_difference_law_error":max_difference,
            "max_jacobian_finite_difference_error":max_jacobian,
            "minimum_dissipation":float(min_diss)}


def main():
    p=Parameters()
    fixed=symmetric_fixed_points(p)
    assert len(fixed)==3
    low,middle,high=[r["u"] for r in fixed]
    cases={
        "small_symmetric_lesion":(p,np.array([.75*high,.75*high])),
        "asymmetric_lesion":(p,np.array([.5*high,high])),
        "below_symmetric_threshold":(p,np.array([.98*middle]*2)),
        "above_symmetric_threshold":(p,np.array([1.02*middle]*2)),
        "near_zero_initial_activity":(p,np.array([1e-10,1e-10])),
        "x_reciprocal_channel_removed":(replace(p,channel_x=0),np.array([high,high])),
        "y_reciprocal_channel_removed":(replace(p,channel_y=0),np.array([high,high])),
        "both_reciprocal_channels_removed":(replace(p,channel_x=0,channel_y=0),np.array([high,high])),
        "zero_chemical_drive":(replace(p,fuel=p.waste*np.exp(-p.fuel_energy/p.thermal)),np.array([high,high])),
    }
    results={name:trajectory(params,start) for name,(params,start) in cases.items()}
    for name in ("small_symmetric_lesion", "asymmetric_lesion"):
        start = np.array(results[name]["initial"])
        converted = 2*high-float(start.sum())
        results[name]["conditioned_preparation"] = {
            "reference_active_concentrations": [high, high],
            "active_concentration_converted_to_D": converted,
            "standard_energy_released_per_volume": p.active_energy*converted,
            "material_removed": 0,
            "scope": "Conditioned post-lesion initial state; preparation/controller work is not modeled. The subsequent chemical-work ledger starts here.",
        }
    for name in ["small_symmetric_lesion","asymmetric_lesion","above_symmetric_threshold"]:
        assert np.max(np.abs(np.array(results[name]["final"])-high))<2e-7
    for name in ["below_symmetric_threshold","near_zero_initial_activity"]:
        assert np.max(np.abs(np.array(results[name]["final"])-low))<2e-7
    for name in ["x_reciprocal_channel_removed","y_reciprocal_channel_removed",
                 "both_reciprocal_channels_removed","zero_chemical_drive"]:
        assert max(results[name]["final"])<.002
    eq=results["zero_chemical_drive"]
    assert np.max(np.abs(np.array(eq["final"])-equilibrium_material(p)))<2e-7
    assert eq["largest_free_energy_increase_on_grid"]<1e-10
    assert abs(eq["integrated_chemical_work"])<1e-12
    q=p.source/p.decay
    r=p.reverse_coefficient*p.waste
    s=p.copy_coefficient*p.fuel
    zero_drive=cases["zero_chemical_drive"][0]
    test_u=np.linspace(1e-8,p.material/2-1e-8,100)
    factor_error=max(abs(float(dynamics(zero_drive,[u,u])[0])
                         -(r*u*u+p.decay)*(q*(p.material-2*u)-u))
                     for u in test_u)
    assert factor_error<1e-12
    source_zero=np.sort(np.roots([2*s+r,-s*p.material,p.decay]))
    output={
        "status":"Concentration/large-copy driven reciprocal-maintenance model; local detailed balance, explicit chemostat work",
        "parameters":asdict(p),"reverse_coefficient":p.reverse_coefficient,
        "thermal_source_coefficient":p.source,"chemical_work_gap":p.chemical_gap,
        "negligible_source_limiting_nonzero_roots":source_zero.tolist(),
        "negligible_source_saddle_node_condition":{
            "a_squared":float((s*p.material)**2),
            "four_b_delta":float(4*(2*s+r)*p.decay),
        },
        "symmetric_fixed_points":fixed,
        "zero_drive_equilibrium_u":float(equilibrium_material(p)[0]),
        "zero_drive_symmetric_factorization_max_error":float(factor_error),
        "one_channel_ablation_bound":{
            "scope":"stationary concentrations after removing one reciprocal reaction pair",
            "ratio_eta_over_delta":float(q),
            "max_remaining_self_growth_coefficient":float(s*q*p.material**2),
            "decay":p.decay,
            "inactive_channel_species_upper_bound":float(q*p.material),
            "other_species_upper_bound":float(p.source*p.material/(p.decay-s*q*p.material**2)),
        },
        "random_state_audit":audit_random_states(p),"cases":results,
        "interpretation_limits":[
            "N=4 is a concentration inventory, not four individual molecules.",
            "Free energies, chemical work, heat and powers are normalized per-volume quantities.",
            "Bistability is deterministic; a finite positive-rate count-state chain has a unique stationary distribution, with possible metastability.",
            "The supplied cooperative catalytic grammar and well-mixed reactor do not derive a membrane, spatial boundary, reproduction, or individuality.",
            "Removing a reciprocal channel removes both directions; its remaining partner cannot sustain the high state for these parameters.",
            "The source-free threshold is a negligible-source limiting approximation, not exact finite-temperature reversibility.",
            "Fuel and waste are chemostats: continuous high-state maintenance consumes explicitly counted external chemical work.",
        ],"checks_passed":True,
    }
    path=Path(__file__).with_name("agent-closure-results.json")
    path.write_text(json.dumps(output,indent=2)+"\n")
    print(json.dumps({"fixed_points":fixed,"random_audit":output["random_state_audit"],
                      "final_cases":{k:v["final"] for k,v in results.items()},"saved":str(path)},indent=2))


if __name__=="__main__":main()

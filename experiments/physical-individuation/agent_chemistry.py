"""Isothermal, material-balanced alternative regenerative chemistries.

Two equal-volume patches have the same feed, temperature, and reaction grammar.
Only initial chemical identity differs. Diversity is never a dynamical reward.
All activity, stored-energy, and dissipation scores are observations of the ODE.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root

NAMES = ["F", "A0", "B0", "C0", "W0", "A1", "B1", "C1", "W1"]
MASS = np.array([4, 1, 1, 2, 2, 1, 1, 2, 2.], float)
ENERGY = np.array([3, .6, .6, 1, 1, .6, .6, 1, 1.], float)
STOICH = np.zeros((9, 4))
for i, a in enumerate((1, 5)):
    STOICH[[a, a+1, a+2], 2*i] = [-1, -1, 1]
    STOICH[[0, a, a+1, a+3], 2*i+1] = [-1, 1, 1, 1]
RELEASE = -ENERGY @ STOICH
OBS = ["assembly_time_integral", "component_standard_energy_time_integral",
       "fuel_converted", "standard_reaction_energy_release", "energy_outflow",
       "free_energy_time_integral", "reaction_dissipation", "diffusion_dissipation",
       "inlet_mixing_dissipation", "free_energy_outflow", "effective_catalyst_time_integral"]


@dataclass(frozen=True)
class Parameters:
    feed: float = 6.
    dilution: float = 1.
    thermal: float = .05
    shock_start: float = 20.
    shock_duration: float = 6.
    shock_factor: float = 1e-4
    mixing_duration: float = 2.
    mixing_rate: float = 1.
    horizon: float = 120.


def free_energy(x: np.ndarray, p: Parameters) -> float:
    positive = np.maximum(x, 1e-300)
    return float(ENERGY @ x + p.thermal * np.sum(x * (np.log(positive)-1)))


def fluxes(x: np.ndarray, p: Parameters, barriers: np.ndarray):
    plus = np.zeros(4); minus = np.zeros(4)
    for i, a in enumerate((1, 5)):
        aa, bb, cc, ww = x[a:a+4]
        plus[2*i] = aa*bb
        minus[2*i] = np.exp(-.2/p.thermal)*cc
        plus[2*i+1] = barriers[i]*x[0]*cc
        minus[2*i+1] = barriers[i]*np.exp(-.8/p.thermal)*cc*aa*bb*ww
    return plus, minus


def dynamics(state: np.ndarray, p: Parameters, barriers: np.ndarray, mixing: float):
    patches = state[:18].reshape(2, 9)
    dpatches = np.zeros_like(patches)
    aggregate = np.zeros(len(OBS))
    for k, x in enumerate(patches):
        plus, minus = fluxes(x, p, barriers)
        net = plus-minus
        dpatches[k] = STOICH @ net-p.dilution*x
        dpatches[k, 0] += p.dilution*p.feed
        dpatches[k] += mixing*(patches[1-k]-x)
        g = free_energy(x, p)
        mu = ENERGY+p.thermal*np.log(np.maximum(x, 1e-300))
        reaction_diss = float(net @ (-mu @ STOICH))
        inlet_gap = p.thermal*(p.feed*np.log(p.feed/max(x[0], 1e-300))-p.feed+np.sum(x))
        aggregate += .5*np.array([
            x[3]+x[7], .6*(x[1]+x[2]+x[5]+x[6])+x[3]+x[7],
            net[1]+net[3], RELEASE@net, p.dilution*(ENERGY@x), g,
            reaction_diss, 0., p.dilution*inlet_gap, p.dilution*g,
            barriers[0]*x[3]+barriers[1]*x[7],
        ])
    # Equal-volume reciprocal diffusion; no molecule is created or replaced.
    if mixing:
        mu0 = ENERGY+p.thermal*np.log(np.maximum(patches[0], 1e-300))
        mu1 = ENERGY+p.thermal*np.log(np.maximum(patches[1], 1e-300))
        aggregate[7] = .5*mixing*np.dot(patches[0]-patches[1], mu0-mu1)
    return np.r_[dpatches.ravel(), aggregate]


def initial_state(kind: str, p: Parameters):
    # A mature resident of either identity has identical concentrations/energy.
    def resident_rhs(y):
        f, a, b, c, w = y
        ja = a*b-np.exp(-.2/p.thermal)*c
        jg = f*c-np.exp(-.8/p.thermal)*c*a*b*w
        return [p.dilution*(p.feed-f)-jg, jg-ja-p.dilution*a,
                jg-ja-p.dilution*b, ja-p.dilution*c, jg-p.dilution*w]
    sol = root(resident_rhs, [2, 1, 1, 1, 2], tol=1e-12)
    assert sol.success and np.max(np.abs(sol.fun)) < 1e-9 and min(sol.x)>0
    patches = np.zeros((2, 9))
    identities = {"mono0": (0, 0), "mono1": (1, 1), "diverse": (0, 1)}[kind]
    for k, identity in enumerate(identities):
        patches[k, 0] = sol.x[0]
        a = 1+4*identity
        patches[k, a:a+4] = sol.x[1:]
    return np.r_[patches.ravel(), np.zeros(len(OBS))]


def run(kind: str, target: int | str, p: Parameters):
    state = initial_state(kind, p)
    initial = state[:18].reshape(2, 9).copy()
    initial_g = np.mean([free_energy(x,p) for x in initial])
    initial_energy = float(np.mean(initial @ ENERGY))
    shock_end = p.shock_start+p.shock_duration
    phases = [(0, p.shock_start, False, 0.),
              (p.shock_start, shock_end, True, 0.),
              (shock_end, shock_end+p.mixing_duration, False, p.mixing_rate),
              (shock_end+p.mixing_duration, p.horizon, False, 0.)]
    ts=[]; ys=[]
    for start, stop, shock, mixing in phases:
        if stop<=start: continue
        barriers = np.ones(2)
        if shock:
            if target == "common": barriers[:] = p.shock_factor
            elif target != "none": barriers[int(target)] = p.shock_factor
        sol = solve_ivp(lambda t,s: dynamics(s,p,barriers,mixing), (start, stop), state,
                        method="DOP853", rtol=2e-9, atol=2e-11,
                        max_step=.5, dense_output=True)
        assert sol.success, sol.message
        state = sol.y[:,-1]
        sampled_t=np.linspace(start, stop, max(2,round((stop-start)*10)+1))
        ts.extend(sampled_t.tolist());ys.extend(sol.sol(sampled_t).T.tolist())
    history=np.array(ys); times=np.array(ts)
    patches=state[:18].reshape(2,9)
    observations=dict(zip(OBS,map(float,state[18:]),strict=True))
    energy_final=float(np.mean(patches@ENERGY))
    g_final=float(np.mean([free_energy(x,p) for x in patches]))
    gin=p.feed*(3+p.thermal*(np.log(p.feed)-1))
    energy_error=abs(energy_final-initial_energy-(3*p.dilution*p.feed*p.horizon
                      -observations["energy_outflow"]-observations["standard_reaction_energy_release"]))
    g_error=abs(g_final-initial_g-(p.dilution*gin*p.horizon
               -observations["free_energy_outflow"]-observations["reaction_dissipation"]
               -observations["diffusion_dissipation"]-observations["inlet_mixing_dissipation"]))
    assemblies=.5*(history[:,3]+history[:,7]+history[:,12]+history[:,16])
    final_count=float(np.mean(patches[:,3]+patches[:,7]))
    summary={"kind":kind,"target":target,"parameters":asdict(p),
             "initial_material":float(np.mean(initial@MASS)),
             "initial_standard_energy":initial_energy,"initial_free_energy":initial_g,
             "initial_assemblies":float(np.mean(initial[:,3]+initial[:,7])),
             "observations":observations,"final_patches":[dict(zip(NAMES,map(float,x),strict=True)) for x in patches],
             "final_assemblies":final_count,"minimum_assemblies":float(min(assemblies)),
             "postshock_assembly_minimum":float(min(assemblies[times>=shock_end])),
             "maximum_material_error":float(np.max(np.abs(.5*(np.sum(history[:,:9]*MASS,axis=1)
                                                  +np.sum(history[:,9:18]*MASS,axis=1))-4*p.feed))),
             "minimum_numeric_concentration":float(np.min(history[:,:18])),
             "standard_energy_balance_error":energy_error,"mixing_aware_free_energy_balance_error":g_error,
             "final_free_energy":g_final}
    assert summary["maximum_material_error"] < 1e-7
    assert summary["minimum_numeric_concentration"] > -1e-8
    assert energy_error < 1e-6 and g_error < 2e-5, (energy_error,g_error)
    return summary, times, history


def main():
    assert np.allclose(MASS@STOICH,0)
    assert np.allclose(RELEASE,[.2,.8,.2,.8])
    p=Parameters()
    cases={}; trajectories={}
    for kind in ("mono0","mono1","diverse"):
        for target in (0,1,"common","none"):
            summary,t,history=run(kind,target,p)
            key=kind+"_"+str(target)
            cases[key]=summary
            if target==0:
                take=np.unique(np.linspace(0,len(t)-1,241,dtype=int))
                trajectories[key]={"time":t[take].tolist(),
                    "assembled_C":(.5*(history[take,3]+history[take,7]
                                         +history[take,12]+history[take,16])).tolist()}
            print(key, summary["final_assemblies"], summary["observations"]["assembly_time_integral"], flush=True)
    for kind in ("mono0","mono1","diverse"):
        pp=Parameters(mixing_rate=0.)
        cases[kind+"_no_mixing"],_,_=run(kind,0,pp)
    cases["low_feed_diverse_unshocked"],_,_=run("diverse","none",Parameters(feed=4.))
    cases["short_shock_mono0"],_,_=run("mono0",0,Parameters(shock_duration=.5))
    # Equal prior probability of the two initial monotypes is equivalent to an
    # unknown, uniform target identity; no event was chosen using the score.
    mean_mono=.5*(cases["mono0_0"]["observations"]["assembly_time_integral"]
                   +cases["mono1_0"]["observations"]["assembly_time_integral"])
    mean_mono_no_mix=.5*(cases["mono0_no_mixing"]["observations"]["assembly_time_integral"]
                          +cases["mono1_no_mixing"]["observations"]["assembly_time_integral"])
    diverse=cases["diverse_0"]["observations"]["assembly_time_integral"]
    diverse_no_mix=cases["diverse_no_mixing"]["observations"]["assembly_time_integral"]
    assert np.isclose(diverse_no_mix,mean_mono_no_mix,rtol=2e-8)
    assert diverse>mean_mono*1.5
    assert cases["diverse_0"]["final_assemblies"]>2.7
    assert cases["mono0_0"]["final_assemblies"]<1e-8
    assert cases["diverse_common"]["final_assemblies"]<1e-8
    assert cases["low_feed_diverse_unshocked"]["final_assemblies"]<1e-8
    assert cases["short_shock_mono0"]["final_assemblies"]>2.7
    for key in ("initial_material","initial_standard_energy","initial_free_energy","initial_assemblies"):
        assert np.isclose(cases["mono0_0"][key],cases["diverse_0"][key],rtol=1e-14)
    # Analytical free-energy and material identities at arbitrary positive states.
    rng=np.random.default_rng(20261021)
    max_matter=max_energy=max_free=0.;min_rxn=min_diff=min_inlet=np.inf
    for _ in range(250):
        patches=np.exp(rng.uniform(-3,1,(2,9)))
        barriers=np.exp(-rng.uniform(0,10,2));mixing=float(rng.uniform(0,2))
        state=np.r_[patches.ravel(),np.zeros(len(OBS))]
        deriv=dynamics(state,p,barriers,mixing)
        dp=deriv[:18].reshape(2,9);obs=deriv[18:]
        matter=np.mean(patches@MASS);energy=np.mean(patches@ENERGY)
        mus=ENERGY+p.thermal*np.log(patches)
        g=np.mean([free_energy(x,p) for x in patches])
        gin=p.feed*(3+p.thermal*(np.log(p.feed)-1))
        max_matter=max(max_matter,abs(np.mean(dp@MASS)-p.dilution*(4*p.feed-matter)))
        max_energy=max(max_energy,abs(np.mean(dp@ENERGY)
             -(3*p.dilution*p.feed-p.dilution*energy-obs[3])))
        max_free=max(max_free,abs(.5*np.sum(mus*dp)
             -(p.dilution*(gin-g)-obs[6]-obs[7]-obs[8])))
        min_rxn=min(min_rxn,obs[6]);min_diff=min(min_diff,obs[7]);min_inlet=min(min_inlet,obs[8])
    assert max(max_matter,max_energy,max_free)<1e-10
    assert min(min_rxn,min_diff,min_inlet)>=0
    result={"status":"Synthetic reversible chemistry; no diversity reward or external reseeding",
            "species":NAMES,"material_weights":MASS.tolist(),"standard_energies":ENERGY.tolist(),
            "two_patch_volume":[.5,.5],"cases":cases,"trajectories":trajectories,
            "comparison":{"equal_prior_mean_monoculture_integrated_C":mean_mono,
                "diverse_integrated_C":diverse,"relative_gain":diverse/mean_mono-1,
                "no_mixing_mean_mono":mean_mono_no_mix,"no_mixing_diverse":diverse_no_mix,
                "no_mixing_mean_identity_error":abs(diverse_no_mix-mean_mono_no_mix)},
            "random_state_identity_checks":{"states":250,"max_material_error":float(max_matter),
                "max_standard_energy_error":float(max_energy),"max_free_energy_error":float(max_free),
                "minimum_reaction_dissipation":float(min_rxn),"minimum_diffusion_dissipation":float(min_diff),
                "minimum_inlet_mixing_dissipation":float(min_inlet)},
            "interpretation_limits":[
                "Reaction grammar, isomer shock responses, compartment volumes and environmental transport schedule are supplied.",
                "No diversity reward, external reseeding, or replacement of destroyed molecules is used.",
                "Passive diffusion transfers actual donor material and dissipates free energy; continuing feed funds regrowth.",
                "Raw reactor free energy includes unconsumed feed and is not a scalar measure of living organization.",
                "C inventory and catalytic throughput are measured physical outcomes, not objectives driving the ODE.",
                "Global common-mode loss, insufficient feed and indiscriminate mixing can defeat diversity.",
                "The result is finite-horizon insurance with physical recolonization, not a universal biodiversity maximum or proof of antifragility.",
            ],"checks_passed":True}
    Path(__file__).with_name("agent-chemistry-results.json").write_text(json.dumps(result,indent=2)+"\n")
    print("saved",Path(__file__).with_name("agent-chemistry-results.json"))


if __name__=="__main__":main()

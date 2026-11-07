"""Bounded joined-generator software checks; no phase/population experiment."""
from dataclasses import asdict, replace
import gzip
from hashlib import sha256
from itertools import product
import json
import math
from pathlib import Path
import time

import numpy as np

from productive_model import A, B, D, X, Y, COPY, TEMPLATE_PASSIVE, Model, Params, Flags, simulate


ROOT = Path(__file__).resolve().parent


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def hashes():
    return {name:digest(ROOT/name) for name in ("productive_model.py","model-spec.md","software_preflight.py")}


def direct_energy(state,p):
    b = (state.material.reshape(p.side,p.side) == B).astype(float)
    templates = state.templates.reshape(p.side,p.side)
    interactions = np.sum(b*np.roll(b,1,axis=0))+np.sum(b*np.roll(b,1,axis=1))
    return float(p.epsilon_B*b.sum()-p.J_BB*interactions
                 +p.E_X*np.sum(templates == X)+p.E_Y*np.sum(templates == Y)
                 -p.g_X*np.sum((templates == X)*b)-p.g_Y*np.sum((templates == Y)*b))


def main():
    output = ROOT/"software-preflight.json"
    if output.exists(): raise RuntimeError("software preflight already frozen")
    bound = hashes()
    rng = np.random.default_rng(856000)
    pairs,states = 0,0
    max_energy_error,max_log_error,max_drift_error,max_heff_error = 0.,0.,0.,0.
    flags_list = (Flags(),Flags(enable_Y=False),Flags(enable_X=False))
    for mode,mu,kappa,flags in product(("specialist","generalist"),(0.,12.),(0.,.0002),flags_list):
        p = Params(side=4,material_count=5,initial_X=2,initial_Y=1,mu=mu,background_AB_rate=kappa)
        model = Model(p,flags,mode)
        state = model.initial(856001)
        carriers = state.material != 0
        state.material[carriers] = rng.choice((A,B),int(carriers.sum()))
        state.templates[:] = rng.choice((D,X,Y),model.n)
        live = np.flatnonzero(state.templates != D)
        state.instance_tags[:] = 0
        state.instance_tags[live] = np.arange(1,len(live)+1)
        state.next_tag = len(live)+1
        model.validate_state(state)
        states += 1
        before_b = int((state.material == B).sum())
        before_n = int((state.templates != D).sum())
        rates = model.rates(state)[1]
        thermal = 0.
        parent_forward = np.zeros(model.n)
        for k in np.flatnonzero(rates > 0):
            edge = model.events[int(k)]
            changed = model.copy_state(state)
            h0 = direct_energy(changed,p)
            record = model.apply(changed,int(k))
            h1 = direct_energy(changed,p)
            reverse = int(model.reverse_index[k])
            reverse_rate = model.event_rate(changed,reverse)
            log_ratio = math.log(rates[k])-math.log(reverse_rate)
            expected = (record["chemical_work"]-record["delta_H"])/p.theta
            max_energy_error = max(max_energy_error,abs(h1-h0-record["delta_H"]))
            max_log_error = max(max_log_error,abs(log_ratio-expected))
            assert abs(log_ratio-expected) < 2e-13
            delta_b = int((changed.material == B).sum())-before_b
            delta_n = int((changed.templates != D).sum())-before_n
            if edge.kind == COPY:
                assert delta_b == edge.direction*(1 if edge.activity == X else -1)
                assert delta_n == edge.direction
                assert record["chemical_work"] == p.mu*edge.direction
                assert record["catalyst_instance"] == int(state.instance_tags[edge.j])
                if edge.direction == 1:
                    assert record["parent_instance"] == int(state.instance_tags[edge.j])
                    parent_forward[edge.j] += rates[k]
                else:
                    assert record["parent_instance"] == 0
                    assert record["instance_removed"] == int(state.instance_tags[edge.i])
            if edge.kind == TEMPLATE_PASSIVE and edge.before == D:
                thermal += rates[k]
            if kappa == 0 and flags.enable_X != flags.enable_Y:
                sign = 1 if flags.enable_X else -1
                heff_delta = h1-h0-sign*p.mu*delta_b
                max_heff_error = max(max_heff_error,abs(log_ratio+heff_delta/p.theta))
            model.apply(changed,reverse)
            assert np.array_equal(changed.material,state.material)
            assert np.array_equal(changed.templates,state.templates)
            assert np.array_equal(changed.material_tags,state.material_tags)
            model.validate_state(changed)
            # IDs/counters are observer history: a thermal rebirth can have a
            # new ID even when the physical occupancy returns to its start.
            pairs += 1
        assert np.max(parent_forward) <= 4*p.copy_rate+1e-14
        c0 = math.exp(-p.E_X/p.theta)+math.exp(-p.E_Y/p.theta)
        cB = math.exp(-(p.E_X-p.g_X)/p.theta)+math.exp(-(p.E_Y-p.g_Y)/p.theta)
        nD = int((state.templates == D).sum())
        nDB = int(np.sum((state.templates == D)&(state.material == B)))
        expected_thermal = p.loss_rate*(c0*nD+(cB-c0)*nDB)
        max_drift_error = max(max_drift_error,abs(thermal-expected_thermal))
        assert abs(thermal-expected_thermal) < 1e-13
    assert max_energy_error < 1e-12 and max_heff_error < 1e-12
    # Explicit local full cycle at target0 with two distinct NN parents1 and4.
    p = Params(side=4,material_count=1,initial_X=1,initial_Y=1)
    model = Model(p)
    state = model.initial(856002)
    state.material[:] = 0; state.material[0] = A
    state.material_tags[:] = 0; state.material_tags[0] = 1
    state.templates[:] = D; state.templates[1] = X; state.templates[4] = Y
    state.instance_tags[:] = 0; state.instance_tags[1] = 1; state.instance_tags[4] = 2
    state.next_tag = 3
    original = model.copy_state(state)
    cycle = []
    for activity,parent in ((X,1),(Y,4)):
        copy = next(k for k,e in enumerate(model.events) if e.kind == COPY and e.i == 0
                    and e.j == parent and e.activity == activity and e.child == activity and e.direction == 1)
        loss = next(k for k,e in enumerate(model.events) if e.kind == TEMPLATE_PASSIVE
                    and e.i == 0 and e.before == activity and e.after == D)
        cycle.extend((copy,loss))
    log_affinity = 0.
    for k in cycle:
        rate = model.event_rate(state,k)
        model.apply(state,k)
        reverse = model.event_rate(state,int(model.reverse_index[k]))
        log_affinity += math.log(rate)-math.log(reverse)
    assert np.array_equal(state.material,original.material)
    assert np.array_equal(state.templates,original.templates)
    assert state.cumulative_work == 2*p.mu
    assert abs(log_affinity-2*p.mu/p.theta) < 1e-13
    # Short software traces only, with the same microscopic initialization and
    # proposal stream across flags/mode. No formation outcome is tested.
    directory = ROOT/"software-preflight"
    directory.mkdir(exist_ok=True)
    traces = []
    small = Params(side=6,material_count=11,initial_X=2,initial_Y=0)
    for mode,flags,mu in (("specialist",Flags(),12.),("generalist",Flags(),12.),
                          ("specialist",Flags(enable_Y=False),12.),("specialist",Flags(enable_X=False),12.),
                          ("specialist",Flags(),0.)):
        params = replace(small,mu=mu)
        seed = 856003
        start = time.monotonic()
        run = simulate(params,seed,1.,flags,mode,sample_dt=5.)
        elapsed = time.monotonic()-start
        name = f"side6-{mode}-x{int(flags.enable_X)}y{int(flags.enable_Y)}-mu{mu:g}.json.gz"
        path = directory/name
        if path.exists(): raise RuntimeError("software trace already exists")
        with gzip.open(path,"wt",encoding="utf-8") as stream:
            json.dump(run,stream,separators=(",",":"),allow_nan=False)
        assert all(v != B for v in run["initial_state"]["material"])
        assert all(v != Y for v in run["initial_state"]["templates"])
        traces.append(dict(file=str(path.relative_to(ROOT)),sha256=digest(path),seed=seed,params=asdict(params),
                           flags=asdict(flags),activity_mode=mode,horizon=1.,elapsed_seconds=elapsed,
                           accepted_events=len(run["ledger"]),proposals=run["proposal_statistics"]["proposals"],
                           first_law_error=run["first_law_error"]))
    if hashes() != bound: raise RuntimeError("source changed during software preflight")
    result = dict(checks_passed=True,source_hashes=bound,random_fixture_states=states,
                  enabled_reverse_pairs=pairs,maximum_full_energy_delta_error=max_energy_error,
                  maximum_log_LDB_error=max_log_error,maximum_thermal_intensity_error=max_drift_error,
                  maximum_one_family_effective_H_DB_error=max_heff_error,
                  full_cycle_log_affinity=log_affinity,full_cycle_expected_affinity=2*p.mu/p.theta,
                  per_parent_forward_hazard_cap_checked=True,traces=traces,
                  scope="Bounded software and structural-identity checks only. No spatial formation, high population, coexistence, lineage establishment, phase size, generalist superiority or individuality is tested. Exclude these traces from any prospective physical outcome study.")
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(json.dumps(result,indent=2))


if __name__ == "__main__":
    main()

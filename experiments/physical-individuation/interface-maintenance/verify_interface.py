"""Independent physical audit of interface_model.py and saved SSA trajectories.

This verifier rebuilds lattice geometry, Hamiltonian, changes, propensities and
instance bookkeeping. Trajectory replay never calls Model.apply/delta_energy.
Outputs JSON to stdout and does not modify model, archives or result files.

Run with a NumPy-equipped Python, optionally followed by result JSON/JSON.GZ,
summary files containing trajectory_file references, directories or globs.
"""
from __future__ import annotations

import argparse
import copy
from dataclasses import asdict, replace
import glob
import gzip
from hashlib import sha256
import itertools
import json
import math
from pathlib import Path
import random
import sys

sys.dont_write_bytecode = True
import numpy as np
from interface_model import Flags, Model, Params


NAMES = ("passive", "copy_R", "produce_L", "material_swap", "fuel_swap")
R = (1, 2)


def close(a, b, tolerance=2e-8):
    assert math.isfinite(float(a)) and math.isfinite(float(b)), (a, b)
    assert abs(float(a)-float(b)) <= tolerance, (a, b, abs(float(a)-float(b)))


class Reference:
    """Scalar, independently indexed implementation of the specified physics."""

    def __init__(self, params, flags=None, environment=None):
        self.p = params
        self.flags = flags or Flags()
        self.env = [1., 1.] if environment is None else list(environment)
        s = params.side
        self.n = s*s
        self.neighbors = []
        self.edges = []
        for i in range(self.n):
            row, col = divmod(i, s)
            self.neighbors.append((((row-1) % s)*s+col, ((row+1) % s)*s+col,
                                   row*s+(col-1) % s, row*s+(col+1) % s))
            self.edges.extend(((i, ((row+1) % s)*s+col), (i, row*s+(col+1) % s)))
        self.incident = [set() for _ in range(self.n)]
        for k, (i, j) in enumerate(self.edges):
            self.incident[i].add(k)
            self.incident[j].add(k)
        self.events = []
        for i in range(self.n):
            for species in (1, 2, 3):
                self.events.extend(((0, i, -1, 0, species), (0, i, -1, species, 0)))
            for j in self.neighbors[i]:
                for species in R:
                    self.events.extend(((1, i, j, 0, species), (1, i, j, species, 0)))
                self.events.extend(((2, i, j, 0, 3), (2, i, j, 3, 0)))
        for i, j in self.edges:
            self.events.extend(((3, i, j, -1, -1), (4, i, j, -1, -1)))
        self.lookup = {e: k for k, e in enumerate(self.events)}

    def onsite(self, species):
        return self.p.epsilon_R if species in R else self.p.epsilon_L if species == 3 else 0.

    def bond(self, a, b):
        if a == b == 3:
            return -self.p.J_LL
        if (a in R and b == 3) or (a == 3 and b in R):
            return -self.p.J_RL
        return -self.p.J_RR if a in R and b in R else 0.

    def energy(self, state):
        m, f = state["material"], state["fuel"]
        return (sum(self.onsite(a) for a in m) + self.p.E_F*sum(f)
                + sum(self.bond(m[i], m[j]) for i, j in self.edges))

    def changes(self, state, event):
        kind, i, j, before, after = event
        m, f = state["material"], state["fuel"]
        if kind == 3:
            return {i: m[j], j: m[i]}, {}
        if kind == 4:
            return {}, {i: f[j], j: f[i]}
        return {i: after}, ({i: 0 if before == 0 else 1} if kind in (1, 2) else {})

    def delta(self, state, event):
        material_changes, fuel_changes = self.changes(state, event)
        m, f = state["material"], state["fuel"]
        value = sum(self.onsite(a)-self.onsite(m[i]) for i, a in material_changes.items())
        value += self.p.E_F*sum(a-f[i] for i, a in fuel_changes.items())
        affected = set()
        for i in material_changes:
            affected.update(self.incident[i])
        for k in affected:
            i, j = self.edges[k]
            value += (self.bond(material_changes.get(i, m[i]), material_changes.get(j, m[j]))
                      - self.bond(m[i], m[j]))
        return value

    def rate(self, state, event):
        kind, i, j, before, after = event
        m, f, p = state["material"], state["fuel"], self.p
        if kind <= 2:
            if m[i] != before:
                return 0.
            product = max(before, after)
            if kind == 0:
                if product == 3:
                    mobility = p.passive_L
                else:
                    allele = product-1
                    mobility = p.passive_R*p.allele_passive_multipliers[allele]*self.env[allele]
                    if self.flags.protect_R:
                        mobility *= math.exp(-p.protection_strength*sum(m[k] == 3 for k in self.neighbors[i]))
            else:
                if m[j] not in R or f[i] != (1 if before == 0 else 0):
                    return 0.
                mobility = (p.copy_rate*p.mutation[m[j]-1][product-1]*self.flags.copy_R
                            if kind == 1 else p.production_rate*self.flags.produce_L)
        elif kind == 3:
            if m[i] == m[j]:
                return 0.
            mobility = p.material_mobility
            if self.flags.protect_R and (m[i] in R or m[j] in R):
                union = set((i, j, *self.neighbors[i], *self.neighbors[j]))
                mobility *= math.exp(-p.transport_protection_strength*sum(m[k] == 3 for k in union))
        else:
            return p.fuel_mobility if f[i] != f[j] else 0.
        return mobility*math.exp(-max(self.delta(state, event), 0.)/p.theta)

    def validate(self, state):
        m, f, tags = state["material"], state["fuel"], state["instance_tags"]
        assert len(m) == len(f) == len(tags) == self.n
        assert all(a in (0, 1, 2, 3) for a in m)
        assert all(a in (0, 1) for a in f)
        assert all((a == 0) == (tag == 0) for a, tag in zip(m, tags))
        live = [tag for tag in tags if tag]
        assert len(live) == len(set(live))
        assert state["next_tag"] > max(live, default=0)

    def apply(self, state, event_index):
        """Independent replay, including observational identity bookkeeping."""
        e = self.events[event_index]
        kind, i, j, before, after = e
        assert self.rate(state, e) > 0, ("disabled/unphysical event", e)
        delta = self.delta(state, e)
        changes, fuel_changes = self.changes(state, e)
        old_fuel = sum(state["fuel"])
        tags = state["instance_tags"]
        created = removed = parent = catalyst = 0
        if kind <= 2:
            if before:
                removed, tags[i] = tags[i], 0
            if after:
                created = state["next_tag"]
                state["next_tag"] += 1
                tags[i] = created
                if kind in (1, 2):
                    catalyst = tags[j]
                    assert catalyst > 0
                    if kind == 1:
                        parent = catalyst
        elif kind == 3:
            tags[i], tags[j] = tags[j], tags[i]
        for k, value in changes.items():
            state["material"][k] = value
        for k, value in fuel_changes.items():
            state["fuel"][k] = value
        state["cumulative_heat"] += -delta
        state["event_count"] += 1
        fuel_delta = self.p.E_F*(sum(state["fuel"])-old_fuel)
        return {"time": state["time"], "event_index": event_index, "kind": NAMES[kind],
                "i": i, "j": j, "before": before, "after": after,
                "delta_H": delta, "heat_to_bath": -delta, "external_work": 0.,
                "fuel_before": old_fuel, "fuel_after": sum(state["fuel"]),
                "fuel_standard_energy_change": fuel_delta, "material_energy_change": delta-fuel_delta,
                "instance_created": created, "instance_removed": removed,
                "parent_instance": parent, "catalyst_instance": catalyst}

    def summary(self, state):
        m = state["material"]
        rl = sum((m[i] in R and m[j] == 3) or (m[j] in R and m[i] == 3) for i, j in self.edges)
        ll = sum(m[i] == m[j] == 3 for i, j in self.edges)
        return {"time": state["time"], "D": m.count(0), "R0": m.count(1), "R1": m.count(2),
                "L": m.count(3), "fuel": sum(state["fuel"]), "H": self.energy(state),
                "heat_to_bath": state["cumulative_heat"], "external_work": state["external_work"],
                "RL_contacts": rl, "LL_contacts": ll, "events": state["event_count"]}


def compare_record(expected, observed, max_errors=None):
    for key, value in expected.items():
        assert key in observed, ("missing ledger field", key)
        if key in ("time", "delta_H", "heat_to_bath", "external_work", "H",
                   "fuel_standard_energy_change", "material_energy_change"):
            error = abs(float(value)-float(observed[key]))
            if max_errors is not None:
                max_errors[key] = max(max_errors.get(key, 0.), error)
            close(value, observed[key])
        else:
            assert value == observed[key], (key, value, observed[key])


def compare_state(expected, observed):
    for key in ("material", "fuel", "instance_tags", "next_tag", "event_count"):
        assert expected[key] == observed[key], ("final state mismatch", key)
    for key in ("time", "cumulative_heat", "external_work"):
        close(expected[key], observed[key])


def verify_simulation(result):
    """Replay one complete public simulate() record, without Model.apply."""
    p = Params(**result["params"])
    flags = Flags(**result["flags"])
    ref = Reference(p, flags, result.get("environment"))
    state = copy.deepcopy(result["initial_state"])
    ref.validate(state)
    initial_h, initial_q, initial_w = ref.energy(state), state["cumulative_heat"], state["external_work"]
    initial_fuel = sum(state["fuel"])
    initial_count = state["event_count"]
    direct_h = initial_h
    errors = {}
    samples = result.get("samples", [])
    sample_index = 0
    counts = dict.fromkeys(NAMES, 0)
    forward_fuel = reverse_fuel = 0

    def check_samples(until):
        nonlocal sample_index
        while sample_index < len(samples) and samples[sample_index]["time"] <= until:
            expected = ref.summary(state)
            expected["time"] = samples[sample_index]["time"]
            compare_record(expected, samples[sample_index], errors)
            sample_index += 1

    for k, row in enumerate(result["ledger"]):
        assert state["time"] <= row["time"] <= result["horizon"]
        check_samples(row["time"])
        state["time"] = row["time"]
        index = int(row["event_index"])
        assert 0 <= index < len(ref.events)
        event = ref.events[index]
        expected = ref.apply(state, index)
        compare_record(expected, row, errors)
        counts[NAMES[event[0]]] += 1
        if event[0] in (1, 2):
            if event[3] == 0:
                forward_fuel += 1
            else:
                reverse_fuel += 1
        assert sum(state["fuel"]) == initial_fuel-forward_fuel+reverse_fuel
        direct_h += expected["delta_H"]
        close(initial_h+state["external_work"]-initial_w-direct_h-state["cumulative_heat"]+initial_q, 0.)
        if k % 97 == 0 or k == len(result["ledger"])-1:
            close(ref.energy(state), direct_h)
            ref.validate(state)
    check_samples(result["horizon"])
    assert sample_index == len(samples)
    state["time"] = float(result["horizon"])
    compare_state(state, result["final_state"])
    compare_record(ref.summary(state), result["final_summary"], errors)
    assert state["event_count"]-initial_count == len(result["ledger"])
    close(ref.energy(state), direct_h)
    return {"seed": result.get("seed"), "flags": asdict(flags), "events_replayed": len(result["ledger"]),
            "samples_checked": len(samples), "event_counts": counts,
            "gross_fuel_conversions": forward_fuel, "gross_reverse_fuel_conversions": reverse_fuel,
            "net_fuel_consumed": forward_fuel-reverse_fuel, "initial_H": initial_h,
            "final_H": ref.energy(state), "maximum_record_errors": errors}


def verify_model(configurations=8, seed=19873):
    rng = random.Random(seed)
    tested_pairs = state_cases = rate_comparisons = 0
    maxima = {"energy": 0., "delta_H": 0., "rate": 0., "LDB_log": 0.}
    settings = [Params(side=3, initial_templates=4),
                replace(Params(), side=4, initial_templates=5, J_RR=.8,
                        mutation=((.7, .3), (.2, .8))),
                replace(Params(), side=4, initial_templates=5, theta=.7,
                        epsilon_R=3.2, epsilon_L=2.1, E_F=5.2, J_LL=-.2, J_RL=1.3, J_RR=.3)]
    for p in settings:
        all_flags = [Flags(*values) for values in itertools.product((False, True), repeat=3)]
        initial_states = [Model(p, flags).state_to_dict(Model(p, flags).initial(71)) for flags in all_flags]
        assert all(x == initial_states[0] for x in initial_states)
        for configuration in range(configurations):
            base_model = Model(p)
            state = base_model.state_to_dict(base_model.initial(100+configuration))
            state["material"] = [rng.randrange(4) for _ in range(p.side*p.side)]
            state["fuel"] = [rng.randrange(2) for _ in state["material"]]
            state["instance_tags"] = [i+1 if a else 0 for i, a in enumerate(state["material"])]
            state["next_tag"] = len(state["material"])+1
            environment = (0.37, 2.1) if configuration % 2 else (1., 1.)
            outputs = {}
            for flags in all_flags:
                model, ref = Model(p, flags), Reference(p, flags, environment)
                observed_state = Model.state_from_dict(state)
                observed_events, rates = model.rates(observed_state, environment)
                assert [tuple(asdict(e).values()) for e in observed_events] == ref.events
                assert model.events == observed_events
                expected = np.array([ref.rate(state, e) for e in ref.events])
                rate_error = float(np.max(np.abs(expected-rates)))
                maxima["rate"] = max(maxima["rate"], rate_error)
                assert np.allclose(expected, rates, rtol=2e-12, atol=2e-13)
                rate_comparisons += len(rates)
                outputs[flags] = rates
                energy_error = abs(ref.energy(state)-model.energy(observed_state))
                maxima["energy"] = max(maxima["energy"], energy_error)
                close(ref.energy(state), model.energy(observed_state), 2e-10)
                ref.validate(state)
                selected = []
                for kind in range(5):
                    available = [i for i, e in enumerate(ref.events) if e[0] == kind and rates[i] > 0]
                    selected.extend(rng.sample(available, min(3, len(available))))
                for index in selected:
                    e = ref.events[index]
                    changed = copy.deepcopy(state)
                    predicted_record = ref.apply(changed, index)
                    delta = ref.energy(changed)-ref.energy(state)
                    close(delta, predicted_record["delta_H"], 2e-10)
                    close(delta, model.delta_energy(observed_state, index), 2e-10)
                    maxima["delta_H"] = max(maxima["delta_H"], abs(delta-model.delta_energy(observed_state, index)))
                    actual_state = Model.state_from_dict(state)
                    actual_record = model.apply(actual_state, index)
                    compare_record(predicted_record, actual_record)
                    compare_state(changed, Model.state_to_dict(actual_state))
                    inverse = (e[0], e[1], e[2], e[4], e[3]) if e[0] <= 2 else e
                    reverse_index = ref.lookup[inverse]
                    assert int(model.reverse_index[index]) == reverse_index
                    assert tuple(asdict(model.reverse(index)).values()) == inverse
                    reverse_rates = model.rates(Model.state_from_dict(changed), environment)[1]
                    assert reverse_rates[reverse_index] > 0
                    error = abs(math.log(rates[index]/reverse_rates[reverse_index])+delta/p.theta)
                    maxima["LDB_log"] = max(maxima["LDB_log"], error)
                    assert error < 2e-10
                    tested_pairs += 1
                # Metadata, clock, alleged ancestry and accumulated costs have no
                # effect on the Hamiltonian or any propensity.
                altered = Model.state_from_dict(state)
                altered.instance_tags[:] = np.arange(model.n)+10000
                altered.next_tag = 999999
                altered.time, altered.event_count = 1e9, 987654
                altered.cumulative_heat, altered.external_work = -1e7, 1e8
                assert np.array_equal(model.rates(altered, environment)[1], rates)
                close(model.energy(altered), model.energy(observed_state), 0.)
                state_cases += 1
            # Same microscopic state and potentials across factorial conditions.
            for flags, rates in outputs.items():
                for index, e in enumerate(Reference(p).events):
                    if not flags.produce_L and e[0] == 2:
                        assert rates[index] == 0
                    if not flags.copy_R and e[0] == 1:
                        assert rates[index] == 0
                    if e[0] not in (0, 3):
                        other = outputs[replace(flags, protect_R=not flags.protect_R)]
                        close(rates[index], other[index], 0.)
                    if e[0] != 2:
                        other = outputs[replace(flags, produce_L=not flags.produce_L)]
                        close(rates[index], other[index], 0.)
                    if e[0] != 1:
                        other = outputs[replace(flags, copy_R=not flags.copy_R)]
                        close(rates[index], other[index], 0.)
    return {"parameter_settings": len(settings), "random_state_flag_cases": state_cases,
            "individual_propensities_compared": rate_comparisons, "enabled_reverse_pairs": tested_pairs,
            "maximum_errors": maxima, "factorial_initial_states_equal": True,
            "instance_tags_and_age_metadata_do_not_affect_rates_or_H": True,
            "notes": ["Protection changes passive stabilization and transport retention together.",
                      "Attraction remains in H in every flag condition.",
                      "This verifies local physics, not an enclosing membrane or functional information."]}


def discover(path, seen):
    """Find embedded simulations or relative trajectory_file references."""
    path = Path(path).resolve()
    if path in seen:
        return
    seen.add(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as stream:
        data = json.load(stream)

    def walk(value, label="$"):
        if isinstance(value, dict):
            if {"params", "flags", "initial_state", "final_state", "ledger"} <= value.keys():
                yield str(path)+":"+label, value
                return
            if isinstance(value.get("trajectory_file"), str):
                yield from discover(path.parent/value["trajectory_file"], seen)
            for key, item in value.items():
                yield from walk(item, label+"."+key)
        elif isinstance(value, list):
            for index, item in enumerate(value):
                yield from walk(item, label+"["+str(index)+"]")
    yield from walk(data)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="*")
    parser.add_argument("--configurations", type=int, default=8)
    args = parser.parse_args()
    assert args.configurations > 0
    model_source = Path(__file__).with_name("interface_model.py")
    model_hash = sha256(model_source.read_bytes()).hexdigest()
    model_checks = verify_model(args.configurations)
    replays, seen = [], set()
    for pattern in args.paths:
        matches = sorted(glob.glob(pattern, recursive=True))
        assert matches, ("input path/glob not found", pattern)
        for match in matches:
            path = Path(match)
            files = sorted(set(path.rglob("*.json")) | set(path.rglob("*.json.gz"))) if path.is_dir() else [path]
            for file in files:
                for label, result in discover(file, seen):
                    replays.append({"record": label, **verify_simulation(result)})
    assert sha256(model_source.read_bytes()).hexdigest() == model_hash, "model source changed during audit"
    print(json.dumps({"checks_passed": True,
                      "model_sha256": model_hash,
                      "verifier_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
                      "model_checks": model_checks, "trajectory_records_replayed": len(replays),
                      "total_events_replayed": sum(row["events_replayed"] for row in replays),
                      "replays": replays,
                      "scope": "Independent local-energy/rate/replay audit; not a demonstration of collective individuality, aging, agency, membrane integrity or useful information."}, indent=2))


if __name__ == "__main__":
    main()

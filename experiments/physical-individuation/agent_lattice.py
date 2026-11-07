"""Finite-fuel local copying and binding on a two-dimensional lattice.

Each material slot is D, R0, or R1; an independent fuel slot is F or W.
All transport and reactions obey local detailed balance for one Hamiltonian.
Only local rules change the lattice. Connected components are observations.
The autonomous environmental gate changes reaction barriers and has explicit
well energies. No bulk clearing, replacement seed, or component split exists.

This is a finite coarse lattice chemistry, not an identified molecular system.
The simulation does not certify fission, collective inheritance, or exact
distribution-level availability for the large lattice state space.
Run: python3 agent_lattice.py  (requires numpy)
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from pathlib import Path
import json
import math

import numpy as np


@dataclass(frozen=True)
class Parameters:
    side: int = 20
    thermal: float = 1.0
    binding: float = 2.4
    template_energy: float = 13.0
    fuel_energy: float = 16.0
    copy_rate: float = 0.02
    mutation: float = 0.05
    passive_rate: float = 0.004
    bad_passive_rate: float = 0.18
    material_transport: float = 0.1
    fuel_transport: float = 0.2
    gate_return: float = 0.25
    gate_bad_weight: float = 0.06
    initial_radius: float = 3.0
    initial_layout: str = "compact"
    horizon: float = 120.0
    observations: int = 61


class LocalModel:
    def __init__(self, p: Parameters):
        assert p.side >= 4 and p.thermal > 0 and p.binding >= 0
        assert p.template_energy > 4*p.binding
        assert p.fuel_energy > p.template_energy
        assert 0 <= p.mutation <= 1 and 0 < p.gate_bad_weight <= 1
        self.p = p
        self.n = p.side*p.side
        cells = np.arange(self.n).reshape(p.side, p.side)
        self.neighbors = np.stack([np.roll(cells, step, axis).ravel()
                                   for axis, step in ((0, 1), (0, -1), (1, 1), (1, -1))], axis=1)
        self.left = np.repeat(np.arange(self.n), 2)
        self.right = self.neighbors[:, (1, 3)].ravel()
        self.edge_count = len(self.left)
        self.offset_passive = 2*self.edge_count
        self.offset_copy = self.offset_passive + 2*self.n
        self.size = self.offset_copy + 8*self.n
        self.mu = np.array([[1-p.mutation, p.mutation],
                            [p.mutation, 1-p.mutation]])
        self.gate_energy = -p.thermal*math.log(p.gate_bad_weight)

    def initial(self):
        rows, cols = np.indices((self.p.side, self.p.side))
        centered = (rows-self.p.side//2)**2+(cols-self.p.side//2)**2
        material = np.where(centered.ravel() <= self.p.initial_radius**2, 0, -1)
        assert self.p.initial_layout in ("compact", "dispersed")
        if self.p.initial_layout == "dispersed":
            number = int(np.sum(material == 0))
            material[:] = -1
            selected = np.random.default_rng(771).choice(self.n, number, replace=False)
            material[selected] = 0
        fuel = np.ones(self.n, dtype=np.int8)
        return material.astype(np.int8), fuel

    def energy(self, material, fuel, environment):
        active = material >= 0
        return float(self.p.template_energy*active.sum()
                     + self.p.fuel_energy*fuel.sum()
                     - self.p.binding*np.sum(active[self.left]*active[self.right])
                     + (self.gate_energy if environment else 0))

    def rates(self, material, fuel, environment):
        p = self.p
        active = (material >= 0).astype(np.int8)
        neighbor_count = active[self.neighbors].sum(axis=1)
        ai, aj = active[self.left], active[self.right]
        si, sj = neighbor_count[self.left]-aj, neighbor_count[self.right]-ai
        swap_delta = -p.binding*(aj-ai)*(si-sj)
        material_rates = p.material_transport*np.exp(-np.maximum(swap_delta, 0)/p.thermal)
        material_rates *= material[self.left] != material[self.right]
        fuel_rates = p.fuel_transport*(fuel[self.left] != fuel[self.right])
        activation = p.template_energy-p.binding*neighbor_count
        passive = np.empty((self.n, 2))
        damage = np.full(2, p.passive_rate)
        if environment:
            damage[environment-1] = p.bad_passive_rate
        for allele in range(2):
            passive[:, allele] = damage[allele]*np.where(material == allele, 1.,
                np.where(material == -1, np.exp(-activation/p.thermal), 0.))
        helpers = material[self.neighbors]
        copies = np.zeros((self.n, 4, 2))
        for allele in range(2):
            mobility = p.copy_rate*self.mu[np.maximum(helpers, 0), allele]*(helpers >= 0)
            plus = ((material == -1) & (fuel == 1))[:, None]
            minus = ((material == allele) & (fuel == 0))[:, None]
            copies[:, :, allele] = mobility*(plus+minus*np.exp(
                (activation-p.fuel_energy)/p.thermal)[:, None])
        rates = np.concatenate((material_rates, fuel_rates, passive.ravel(), copies.ravel()))
        assert rates.shape == (self.size,) and np.all(rates >= 0)
        # Row-stochastic mutation changes product identity, not total forward rate.
        intended_forward = p.copy_rate*np.sum((helpers >= 0)
            * (((material == -1) & (fuel == 1))[:, None]))
        observed_forward = copies[((material == -1) & (fuel == 1))].sum()
        assert abs(intended_forward-observed_forward) < 1e-10
        return rates, neighbor_count, swap_delta

    def apply(self, event, material, fuel, neighbor_count, swap_delta):
        """Apply one molecular event; return delta H, label, and copied allele."""
        p = self.p
        if event < self.edge_count:
            i, j = self.left[event], self.right[event]
            material[i], material[j] = material[j], material[i]
            return float(swap_delta[event]), "material_swaps", None
        if event < self.offset_passive:
            edge = event-self.edge_count
            i, j = self.left[edge], self.right[edge]
            fuel[i], fuel[j] = fuel[j], fuel[i]
            return 0., "fuel_swaps", None
        if event < self.offset_copy:
            cell, allele = divmod(event-self.offset_passive, 2)
            activation = p.template_energy-p.binding*neighbor_count[cell]
            if material[cell] == -1:
                material[cell] = allele
                return float(activation), "passive_births", allele
            assert material[cell] == allele
            material[cell] = -1
            return float(-activation), "passive_deaths", None
        cell, within = divmod(event-self.offset_copy, 8)
        helper_slot, allele = divmod(within, 2)
        helper = material[self.neighbors[cell, helper_slot]]
        assert helper >= 0
        forward_delta = p.template_energy-p.binding*neighbor_count[cell]-p.fuel_energy
        if material[cell] == -1:
            assert fuel[cell] == 1
            material[cell], fuel[cell] = allele, 0
            return float(forward_delta), ("erroneous_forward_copies" if helper != allele
                                         else "correct_forward_copies"), allele
        assert material[cell] == allele and fuel[cell] == 0
        material[cell], fuel[cell] = -1, 1
        return float(-forward_delta), "reverse_copies", None


def gate_schedule(model, seed):
    """Pre-sampling an independent autonomous gate equals the joint CTMC law."""
    rng = np.random.default_rng(seed)
    t, state = 0., 0
    schedule = []
    while True:
        rate = (2*model.p.gate_return*model.p.gate_bad_weight if state == 0
                else model.p.gate_return)
        t += rng.exponential(1/rate)
        if t > model.p.horizon:
            return schedule
        state = int(rng.integers(1, 3)) if state == 0 else 0
        schedule.append((float(t), state))


def components(material, neighbors):
    """Physical nearest-neighbor contact components; readout only."""
    remaining = set(map(int, np.flatnonzero(material >= 0)))
    sizes, compositions = [], []
    while remaining:
        origin = remaining.pop()
        stack, found = [origin], [origin]
        while stack:
            cell = stack.pop()
            for neighbor in neighbors[cell]:
                neighbor = int(neighbor)
                if neighbor in remaining:
                    remaining.remove(neighbor)
                    stack.append(neighbor)
                    found.append(neighbor)
        sizes.append(len(found))
        compositions.append(int(np.sum(material[found] == 1)))
    order = np.argsort(sizes)[::-1]
    return ([sizes[i] for i in order], [compositions[i] for i in order])


def thermodynamic_checks(seed=1701, configurations=16, checks_per_configuration=12):
    rng = np.random.default_rng(seed)
    model = LocalModel(replace(Parameters(), side=8))
    maximum_ratio_error = maximum_energy_error = 0.
    checks = 0
    for _ in range(configurations):
        material = rng.integers(-1, 2, size=model.n, dtype=np.int8)
        fuel = rng.integers(0, 2, size=model.n, dtype=np.int8)
        environment = int(rng.integers(0, 3))
        rates, neighbor_count, swap_delta = model.rates(material, fuel, environment)
        candidates = np.flatnonzero(rates > 0)
        for event in rng.choice(candidates, size=checks_per_configuration, replace=True):
            changed, changed_fuel = material.copy(), fuel.copy()
            before = model.energy(material, fuel, environment)
            delta, _, _ = model.apply(int(event), changed, changed_fuel, neighbor_count, swap_delta)
            after = model.energy(changed, changed_fuel, environment)
            reverse = model.rates(changed, changed_fuel, environment)[0][event]
            assert reverse > 0
            ratio_error = abs(math.log(rates[event]/reverse)+delta/model.p.thermal)
            energy_error = abs(after-before-delta)
            maximum_ratio_error = max(maximum_ratio_error, ratio_error)
            maximum_energy_error = max(maximum_energy_error, energy_error)
            checks += 1
    gate_error = abs(math.log(model.p.gate_bad_weight)+model.gate_energy/model.p.thermal)
    assert maximum_ratio_error < 1e-11 and maximum_energy_error < 1e-10
    return {"random_transition_pairs": checks,
            "maximum_local_detailed_balance_log_error": maximum_ratio_error,
            "maximum_direct_Hamiltonian_change_error": maximum_energy_error,
            "gate_detailed_balance_log_error": gate_error,
            "forward_copy_total_propensity_identity_checked": True,
            "method": "independent direct Hamiltonian differences and same-event reverse propensities"}


def simulate(p, seed, environment_seed=20261021, retain_snapshots=False):
    model = LocalModel(p)
    rng = np.random.default_rng(seed)
    material, fuel = model.initial()
    environment, time, gate_index = 0, 0., 0
    gates = gate_schedule(model, environment_seed)
    observation_times = np.linspace(0, p.horizon, p.observations)
    initial_h = model.energy(material, fuel, environment)
    initial_fuel = int(fuel.sum())
    fuel_count = initial_fuel
    tags = np.where(material >= 0, np.arange(model.n), -1)
    next_tag = model.n
    event_ledger = []
    counters = {name: 0 for name in ("material_swaps", "fuel_swaps", "passive_births",
        "passive_deaths", "correct_forward_copies", "erroneous_forward_copies", "reverse_copies",
        "gate_entries", "gate_returns")}
    heat, max_ledger_error = 0., 0.
    first_variant = None
    trace, snapshots = [], []
    total_events = 0

    def observe(observation_time):
        nonlocal max_ledger_error
        sizes, novel_counts = components(material, model.neighbors)
        stored = model.energy(material, fuel, environment)
        ledger_error = abs(initial_h-stored-heat)
        max_ledger_error = max(max_ledger_error, ledger_error)
        forward = counters["correct_forward_copies"]+counters["erroneous_forward_copies"]
        assert int(fuel.sum()) == initial_fuel-forward+counters["reverse_copies"]
        trace.append({"time": float(observation_time), "active_templates": int(np.sum(material >= 0)),
            "allele1_templates": int(np.sum(material == 1)), "contact_components": len(sizes),
            "largest_component": max(sizes, default=0), "component_sizes": sizes,
            "component_allele1_counts": novel_counts, "fuel_remaining": int(fuel.sum()),
            "environment": environment, "stored_internal_energy": stored,
            "net_heat_to_bath": heat, "forward_copies": forward,
            "reverse_copies": counters["reverse_copies"]})
        if retain_snapshots and len(trace)-1 in (0, (p.observations-1)//2, p.observations-1):
            snapshots.append({"time": float(observation_time),
                              "material": material.reshape(p.side, p.side).tolist(),
                              "fuel": fuel.reshape(p.side, p.side).tolist(),
                              "template_instance_tags": tags.reshape(p.side, p.side).tolist()})

    observe(0.)
    observation_index = 1
    while time < p.horizon:
        rates, neighbor_count, swap_delta = model.rates(material, fuel, environment)
        total = float(rates.sum())
        reaction_time = time+rng.exponential(1/total) if total else float("inf")
        gate_time = gates[gate_index][0] if gate_index < len(gates) else float("inf")
        next_time = min(reaction_time, gate_time, p.horizon)
        while observation_index < len(observation_times) and observation_times[observation_index] <= next_time:
            observe(observation_times[observation_index])
            observation_index += 1
        time = next_time
        if time == p.horizon:
            break
        if gate_time <= reaction_time:
            new_environment = gates[gate_index][1]
            delta = model.gate_energy*((new_environment != 0)-(environment != 0))
            counters["gate_entries" if new_environment else "gate_returns"] += 1
            environment = new_environment
            gate_index += 1
            event = -new_environment if new_environment else -3
            product_tag, parent_tag = -1, -1
        else:
            event = int(np.searchsorted(np.cumsum(rates), rng.random()*total, side="right"))
            product_tag, parent_tag = -1, -1
            if event < model.edge_count:
                i, j = model.left[event], model.right[event]
                tags[i], tags[j] = tags[j], tags[i]
            elif model.offset_passive <= event < model.offset_copy:
                cell, _ = divmod(event-model.offset_passive, 2)
                if material[cell] == -1:
                    product_tag = next_tag
                    next_tag += 1
                    tags[cell] = product_tag
                else:
                    tags[cell] = -1
            elif event >= model.offset_copy:
                cell, within = divmod(event-model.offset_copy, 8)
                helper_slot, _ = divmod(within, 2)
                if material[cell] == -1:
                    parent_tag = int(tags[model.neighbors[cell, helper_slot]])
                    assert parent_tag >= 0
                    product_tag = next_tag
                    next_tag += 1
                    tags[cell] = product_tag
                else:
                    tags[cell] = -1
            delta, label, allele = model.apply(event, material, fuel, neighbor_count, swap_delta)
            counters[label] += 1
            if "forward_copies" in label:
                fuel_count -= 1
            elif label == "reverse_copies":
                fuel_count += 1
            if first_variant is None and allele == 1:
                first_variant = {"time": time,
                    "source": "paid_copy" if "forward_copies" in label else "passive_thermal"}
        heat -= delta
        event_ledger.append((time, event, delta, heat, fuel_count, environment,
                             product_tag, parent_tag))
        total_events += 1
    assert len(trace) == len(observation_times)
    assert max_ledger_error < 1e-8
    final = trace[-1]
    times = np.array([x["time"] for x in trace])
    active = np.array([x["active_templates"] for x in trace])
    area = float(np.sum((active[1:]+active[:-1])*np.diff(times)/2))
    coexist = any(0 < row["allele1_templates"] < row["active_templates"] for row in trace)
    return {"parameters": asdict(p), "reaction_seed": seed, "environment_seed": environment_seed,
        "initial_templates": trace[0]["active_templates"], "initial_fuel": initial_fuel,
        "initial_stored_internal_energy": initial_h, "events": counters,
        "total_molecular_and_gate_events": total_events,
        "first_new_allele": first_variant, "observed_two_alleles_simultaneously": coexist,
        "final_active_templates": final["active_templates"],
        "final_allele1_templates": final["allele1_templates"],
        "final_contact_components": final["contact_components"],
        "final_largest_component": final["largest_component"],
        "integrated_active_templates_trapezoid": area,
        "net_fuel_consumed": initial_fuel-final["fuel_remaining"],
        "net_heat_to_bath": heat, "final_stored_internal_energy": final["stored_internal_energy"],
        "maximum_first_law_error_on_observation_grid": max_ledger_error,
        "gate_schedule": gates, "trajectory": trace, "snapshots": snapshots,
        "_event_ledger": np.asarray(event_ledger, dtype=float),
        "interpretation": "Contact components only; no fission, collective heredity, or exact availability claim."}


def main():
    base = Parameters()
    runs = {}
    ledgers = {}
    for layout in ("compact", "dispersed"):
        for mutation in (0., base.mutation):
            label = layout+("_no_mutation" if mutation == 0 else "_paid_copying_errors")
            result = simulate(replace(base, mutation=mutation, initial_layout=layout), 101,
                              retain_snapshots=True)
            ledgers[label] = result.pop("_event_ledger")
            result["event_ledger_key"] = label
            runs[label] = [result]
            print(label, result["final_active_templates"],
                  result["net_fuel_consumed"], result["first_new_allele"], flush=True)
    for label, p in (("no_binding_errors", replace(base, binding=0)),
                     ("high_loss_errors", replace(base, passive_rate=.15))):
        result = simulate(p, 101, retain_snapshots=True)
        ledgers[label] = result.pop("_event_ledger")
        result["event_ledger_key"] = label
        runs[label] = [result]
        print(label, result["final_active_templates"], result["net_fuel_consumed"], flush=True)
    critical = 2*math.log(1+math.sqrt(2))
    if base.binding/base.thermal > critical:
        magnetization = (1-math.sinh(base.binding/(2*base.thermal))**-4)**.125
        dilute, dense = (1-magnetization)/2, (1+magnetization)/2
    else:
        dilute = dense = .5
    result = {"status": "Local finite-fuel stochastic lattice prototype; no component operations in dynamics",
        "thermodynamic_checks": thermodynamic_checks(),
        "phase_turnover_diagnostic": {
            "zero_field_coexistence_reference_dilute_fraction": dilute,
            "zero_field_coexistence_reference_dense_fraction": dense,
            "early_fuel_occupancy": 1., "dense_phase_mean_field_copy_rate_per_template": 4*base.copy_rate*(1-dense),
            "dilute_phase_mean_field_copy_rate_per_template": 4*base.copy_rate*(1-dilute),
            "healthy_passive_loss_rate": base.passive_rate,
            "warning": "Zero-field coexistence reference only; template energy and allele degeneracy shift the field. Actual reacting-lattice phases differ. Uncorrelated source estimates at these reference densities are diagnostics, not a phase/fission theorem."},
        "runs": runs, "checks_passed": True,
        "event_ledger_file": "agent-lattice-event-ledgers.npz",
        "event_ledger_columns": ["time", "event_index", "delta_H", "cumulative_heat_to_bath",
                                 "fuel_remaining", "environment", "product_template_tag", "parent_template_tag"],
        "negative_event_codes": {"-1": "enter_bad0", "-2": "enter_bad1", "-3": "return_healthy"},
        "scope": ["Same local reversible generator supplies binding, transport, copying, mutation, and local damage.",
            "All material slots and fuel/waste slots are conserved. No new fuel enters.",
            "Each event's heat is minus its Hamiltonian change, including environmental gate well energy.",
            "The gate is an independently pre-sampled part of an autonomous detailed-balance joint process.",
            "A finite allele grammar, attractive contact energy, and catalytic copying are supplied chemistry.",
            "The compact founder is a supplied monomorphic 29-template seed; dispersed controls have the same inventory and a different binding energy.",
            "Template-instance and copying-ancestry tags are passive readouts, not conserved atom or material identities; decay erases a template instance.",
            "Connected components are readouts and never trigger copying, splitting, relocation, or selection.",
            "A local bath-induced new allele remains possible; first appearance records its source.",
            "No claim of successful division, collective inheritance, antifragility, or open-ended evolution.",
            "Full state-distribution KL is not estimated from one spatial trajectory."]}
    destination = Path(__file__).with_name("agent-lattice-results.json")
    np.savez_compressed(destination.with_name("agent-lattice-event-ledgers.npz"), **ledgers)
    destination.write_text(json.dumps(result, indent=2)+"\n")
    print("saved", destination, flush=True)


if __name__ == "__main__":
    main()

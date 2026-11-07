"""One bounded causal wound-recovery assay of the existing local lattice.

Factorial comparison: intact sham versus four-site internal lesion, crossed
with local catalytic copying versus paired copying-channel ablation. All
other propensities and the Hamiltonian are inherited from agent_lattice.py.
The gate stays in its healthy state. This is an assay of one partial integrity
criterion, not a model of collective fission, interdependent alleles or selves.

The wounded initial state is conditioned on four allowed passive decay events.
Their energy release is recorded; no matter or fuel is erased. Selecting or
preparing the conditioned state is not modeled as a costless autonomous act.
Each post-preparation trajectory is autonomous and closed to material/fuel.
Run: python3 agent_integrity.py
"""
from __future__ import annotations

from dataclasses import asdict, replace
from pathlib import Path
from hashlib import sha256
import json
import math
import sys

sys.dont_write_bytecode = True
import numpy as np
from agent_lattice import Parameters, LocalModel, thermodynamic_checks


REPLICATES = 64
HORIZON = 40.0


def lesion_cells(model):
    center = model.p.side//2
    return np.array([(center+di)*model.p.side+center+dj
                     for di, dj in ((-1, -1), (-1, 0), (0, -1), (0, 0))])


def prepare(model, wounded):
    material, fuel = model.initial()
    pristine_energy = model.energy(material, fuel, 0)
    preparation = []
    if wounded:
        for cell in lesion_cells(model):
            allele = int(material[cell])
            assert allele >= 0
            rates, counts, swaps = model.rates(material, fuel, 0)
            event = model.offset_passive+2*int(cell)+allele
            assert rates[event] > 0
            delta, label, _ = model.apply(event, material, fuel, counts, swaps)
            assert label == "passive_deaths"
            preparation.append({"cell": int(cell), "allowed_decay_event_index": event,
                                "delta_H": delta, "heat_to_bath": -delta,
                                "R_to_D_units": 1, "material_removed": 0, "fuel_consumed": 0})
    prepared_energy = model.energy(material, fuel, 0)
    heat = sum(row["heat_to_bath"] for row in preparation)
    assert abs(prepared_energy+heat-pristine_energy) < 1e-10
    return material, fuel, {"pristine_energy": pristine_energy,
        "prepared_energy": prepared_energy, "passive_path_heat_to_bath": heat,
        "events": preparation, "scope": "Conditioned starting-state assay; no autonomous wound mechanism or controller-cost claim."}


def void_count(material, model):
    """Inactive cells disconnected from a fixed far-exterior reference region."""
    cells = np.arange(model.n).reshape(model.p.side, model.p.side)
    row, col = np.indices(cells.shape)
    far = (row-model.p.side//2)**2+(col-model.p.side//2)**2 >= 4.5**2
    unseen = set(map(int, np.flatnonzero(material < 0)))
    stack = [x for x in np.flatnonzero(far.ravel()) if int(x) in unseen]
    for x in stack:
        unseen.discard(int(x))
    while stack:
        cell = int(stack.pop())
        for neighbor in model.neighbors[cell]:
            neighbor = int(neighbor)
            if neighbor in unseen:
                unseen.remove(neighbor)
                stack.append(neighbor)
    return len(unseen)


def run_one(p, wounded, seed, keep_snapshots=False):
    model = LocalModel(p)
    material, fuel, preparation = prepare(model, wounded)
    initial_material, initial_fuel = material.copy(), fuel.copy()
    tags = np.where(material >= 0, np.arange(model.n), -1)
    tag_sources = {int(tag): "prepared_template" for tag in tags if tag >= 0}
    next_tag = model.n
    rng = np.random.default_rng(seed)
    roi = lesion_cells(model)
    initial_energy = model.energy(material, fuel, 0)
    time, heat, work = 0., 0., 0.
    forward = reverse = passive_births = passive_deaths = 0
    integral_active = integral_roi = integral_void = integral_perimeter = 0.
    observations = np.linspace(0, p.horizon, 21)
    observation_index = 0
    trajectory, snapshots, ledger = [], [], []
    maximum_first_law_error = maximum_direct_delta_error = 0.
    spatial_catalysis = {}

    def observables():
        active = material >= 0
        perimeter = int(np.sum(active[model.left] != active[model.right]))
        return int(active.sum()), int(active[roi].sum()), void_count(material, model), perimeter

    def observe(t):
        active, occupied, void, perimeter = observables()
        living_tags = tags[tags >= 0]
        original = sum(tag_sources[int(tag)] == "prepared_template" for tag in living_tags)
        copied = sum(tag_sources[int(tag)] == "paid_copy" for tag in living_tags)
        roi_copied = int(sum(tag >= 0 and tag_sources[int(tag)] == "paid_copy" for tag in tags[roi]))
        trajectory.append({"time": float(t), "active_templates": active,
            "lesion_ROI_occupancy": occupied, "enclosed_inactive_cells": void,
            "active_inactive_contact_perimeter": perimeter,
            "retained_prepared_template_instances": original,
            "living_paid_copy_instances": copied,
            "lesion_ROI_paid_copy_instances": roi_copied,
            "fuel_remaining": int(fuel.sum()), "net_heat_to_bath": heat,
            "stored_energy": model.energy(material, fuel, 0)})
        if keep_snapshots and len(trajectory)-1 in (0, 10, 20):
            snapshots.append({"time": float(t), "material": material.reshape(p.side, p.side).tolist(),
                "fuel": fuel.reshape(p.side, p.side).tolist(),
                "template_instance_tags": tags.reshape(p.side, p.side).tolist()})

    while time < p.horizon:
        rates, counts, swaps = model.rates(material, fuel, 0)
        total = float(rates.sum())
        event_time = time+rng.exponential(1/total) if total else float("inf")
        stop = min(event_time, p.horizon)
        active, occupied, void, perimeter = observables()
        duration = stop-time
        integral_active += active*duration
        integral_roi += occupied*duration
        integral_void += void*duration
        integral_perimeter += perimeter*duration
        while observation_index < len(observations) and observations[observation_index] <= stop:
            observe(observations[observation_index])
            observation_index += 1
        time = stop
        if stop == p.horizon:
            break
        event = int(np.searchsorted(np.cumsum(rates), rng.random()*total, side="right"))
        product_tag = parent_tag = -1
        if event < model.edge_count:
            i, j = model.left[event], model.right[event]
            tags[i], tags[j] = tags[j], tags[i]
        elif model.offset_passive <= event < model.offset_copy:
            cell, _ = divmod(event-model.offset_passive, 2)
            if material[cell] == -1:
                product_tag = next_tag
                next_tag += 1
                tags[cell] = product_tag
                tag_sources[product_tag] = "passive_thermal"
            else:
                tags[cell] = -1
        elif event >= model.offset_copy:
            cell, within = divmod(event-model.offset_copy, 8)
            helper_slot, _ = divmod(within, 2)
            helper_cell = int(model.neighbors[cell, helper_slot])
            if material[cell] == -1:
                product_tag, parent_tag = next_tag, int(tags[helper_cell])
                assert parent_tag >= 0
                next_tag += 1
                tags[cell] = product_tag
                tag_sources[product_tag] = "paid_copy"
                key = (helper_cell, int(cell))
                spatial_catalysis[key] = spatial_catalysis.get(key, 0)+1
            else:
                tags[cell] = -1
        before = model.energy(material, fuel, 0)
        delta, label, _ = model.apply(event, material, fuel, counts, swaps)
        after = model.energy(material, fuel, 0)
        maximum_direct_delta_error = max(maximum_direct_delta_error, abs(after-before-delta))
        if "forward_copies" in label:
            forward += 1
        elif label == "reverse_copies":
            reverse += 1
        elif label == "passive_births":
            passive_births += 1
        elif label == "passive_deaths":
            passive_deaths += 1
        heat -= delta
        maximum_first_law_error = max(maximum_first_law_error, abs(initial_energy+work-after-heat))
        assert int(fuel.sum()) == int(initial_fuel.sum())-forward+reverse
        assert int(np.sum(material >= 0)) == int(np.sum(initial_material >= 0))+forward-reverse+passive_births-passive_deaths
        ledger.append((time, event, delta, heat, int(fuel.sum()), product_tag, parent_tag))
    assert observation_index == len(observations)
    assert maximum_direct_delta_error < 1e-10 and maximum_first_law_error < 1e-8
    final = trajectory[-1]
    reciprocal = sum(min(count, spatial_catalysis.get((target, helper), 0))
                     for (helper, target), count in spatial_catalysis.items() if helper < target)
    return {"parameters": asdict(p), "wounded": wounded, "seed": seed,
        "preparation": preparation, "initial_active_templates": int(np.sum(initial_material >= 0)),
        "initial_F_inventory": int(initial_fuel.sum()), "initial_D_inventory": int(np.sum(initial_material == -1)),
        "conserved_material_slots": model.n, "conserved_fuel_plus_waste_slots": model.n,
        "final_active_templates": final["active_templates"],
        "final_lesion_ROI_occupancy": final["lesion_ROI_occupancy"],
        "final_enclosed_inactive_cells": final["enclosed_inactive_cells"],
        "final_perimeter": final["active_inactive_contact_perimeter"],
        "final_living_paid_copy_instances": final["living_paid_copy_instances"],
        "final_retained_prepared_template_instances": final["retained_prepared_template_instances"],
        "final_lesion_ROI_paid_copy_instances": final["lesion_ROI_paid_copy_instances"],
        "integral_active_templates": integral_active,
        "integral_lesion_ROI_occupancy": integral_roi,
        "integral_enclosed_inactive_cells": integral_void,
        "integral_perimeter": integral_perimeter,
        "forward_copies": forward, "reverse_copies": reverse,
        "passive_thermal_births": passive_births, "passive_deaths": passive_deaths,
        "net_fuel_consumed": forward-reverse, "net_heat_to_bath": heat,
        "external_work_post_preparation": work,
        "maximum_event_first_law_error": maximum_first_law_error,
        "maximum_event_Hamiltonian_delta_error": maximum_direct_delta_error,
        "reciprocal_spatial_catalytic_event_pairs": reciprocal,
        "trajectory": trajectory, "snapshots": snapshots,
        "_ledger": np.asarray(ledger, dtype=float)}


def mean_se(values):
    values = np.asarray(values, dtype=float)
    mean = float(values.mean())
    se = float(values.std(ddof=1)/math.sqrt(len(values))) if len(values) > 1 else 0.
    return {"mean": mean, "monte_carlo_standard_error": se,
            "approximate_95_percent_normal_interval": [mean-1.96*se, mean+1.96*se]}


def immediate_response(p):
    """Exact generator-level feedback before any subsequent stochastic event."""
    model = LocalModel(p)
    rates = {}
    for wounded in (False, True):
        material, fuel, _ = prepare(model, wounded)
        q = model.rates(material, fuel, 0)[0]
        copy_rates = q[model.offset_copy:].reshape(model.n, 4, 2)
        rates["wounded" if wounded else "sham"] = {
            "forward_copy_rate_in_lesion_ROI": float(copy_rates[lesion_cells(model)].sum()),
            "total_forward_copy_rate": float(copy_rates.sum())}
    extra = rates["wounded"]["forward_copy_rate_in_lesion_ROI"]-rates["sham"]["forward_copy_rate_in_lesion_ROI"]
    assert abs(extra-8*p.copy_rate) < 1e-12
    return {"rates": rates, "extra_ROI_paid_renewal_propensity": extra,
        "derivation": "Four newly inactive targets each touch two surviving catalysts; the eight local channels have total rate 8*k because each mutation row sums to one.",
        "scope": "An exact initial local repair response, not a macroscopic homeostatic or individuality theorem."}


def main():
    p = replace(Parameters(), side=12, horizon=HORIZON)
    cases, ledgers = {}, {}
    for copying in (False, True):
        for wounded in (False, True):
            name = ("copying" if copying else "paired_copying_ablation")+("_wound" if wounded else "_sham")
            case_p = p if copying else replace(p, copy_rate=0.)
            cases[name] = []
            for replicate in range(REPLICATES):
                seed = 9017+replicate
                result = run_one(case_p, wounded, seed, keep_snapshots=(replicate == 0))
                key = name+"_"+str(replicate)
                ledgers[key] = result.pop("_ledger")
                result["event_ledger_key"] = key
                cases[name].append(result)
            print(name, "mean final active", np.mean([x["final_active_templates"] for x in cases[name]]), flush=True)
    metrics = ["final_active_templates", "final_lesion_ROI_occupancy", "final_enclosed_inactive_cells",
        "final_perimeter", "integral_active_templates", "integral_lesion_ROI_occupancy",
        "integral_enclosed_inactive_cells", "net_fuel_consumed", "final_living_paid_copy_instances",
        "final_retained_prepared_template_instances", "final_lesion_ROI_paid_copy_instances"]
    summaries = {name: {metric: mean_se([r[metric] for r in rows]) for metric in metrics}
                 for name, rows in cases.items()}
    # Difference in differences: whether copying specifically reduces the wound
    # deficit, after its generic effect on the corresponding intact sham.
    contrasts = {}
    for metric in ("final_active_templates", "final_lesion_ROI_occupancy", "integral_active_templates", "integral_lesion_ROI_occupancy"):
        deficits_off = np.array([a[metric]-b[metric] for a, b in zip(cases["paired_copying_ablation_sham"], cases["paired_copying_ablation_wound"])])
        deficits_on = np.array([a[metric]-b[metric] for a, b in zip(cases["copying_sham"], cases["copying_wound"])])
        contrasts[metric] = {"sham_minus_wound_without_copying": mean_se(deficits_off),
            "sham_minus_wound_with_copying": mean_se(deficits_on),
            "copying_reduction_of_wound_deficit": mean_se(deficits_off-deficits_on)}
    void_off = np.array([b["integral_enclosed_inactive_cells"]-a["integral_enclosed_inactive_cells"]
                        for a, b in zip(cases["paired_copying_ablation_sham"], cases["paired_copying_ablation_wound"])])
    void_on = np.array([b["integral_enclosed_inactive_cells"]-a["integral_enclosed_inactive_cells"]
                       for a, b in zip(cases["copying_sham"], cases["copying_wound"])])
    contrasts["copying_reduction_of_wound_excess_void_time"] = mean_se(void_off-void_on)
    maximum_error = max(r["maximum_event_first_law_error"] for rows in cases.values() for r in rows)
    source = Path(__file__).with_name("agent_lattice.py")
    output = {"status": "One matched factorial causal recovery assay; partial integrity criterion only",
        "base_source_sha256": sha256(source.read_bytes()).hexdigest(),
        "replicates_per_case": REPLICATES, "horizon": HORIZON,
        "thermodynamic_checks": thermodynamic_checks(),
        "exact_initial_repair_response": immediate_response(p),
        "maximum_event_first_law_error_all_runs": maximum_error,
        "sham_definition": "Identical prepared seed, without the four allowed R-to-D lesion events; no extra perturbation.",
        "matched_comparison": "Copying versus ablation has exactly the same H, initial microscopic state, inventory, and Gibbs reference within wounded or sham conditions.",
        "unmatched_injury_difference": "The lesion intentionally changes 4 R into 4 D and releases 23.2 energy along the stated passive path; wounded and sham do not have equal initial stored energy H.",
        "coupling_ablation": "Set only the common mobility of forward and reverse neighbor-catalyzed copying to zero; leave binding, transport, loss, thermal activation, and H unchanged.",
        "statistics": "Paired seed contrasts; iid seeds across replicates. Approximate Monte Carlo intervals, not rigorous probability bounds.",
        "summaries": summaries, "contrasts": contrasts, "cases": cases,
        "event_ledger_file": "agent-integrity-event-ledgers.npz",
        "event_ledger_columns": ["time", "event_index", "delta_H", "cumulative_heat_to_bath", "fuel_remaining", "product_template_tag", "parent_template_tag"],
        "scope": ["The same local kinetics operate after a conditioned lesion or intact sham.",
            "Copying is metabolic template renewal; disabling it is not a test of organism-level sterility.",
            "Zero copying mobility preserves the formal Gibbs stationary reference but disconnects fuel-conversion sectors; no reachable-equilibrium or availability equality is asserted.",
            "An inactive unit remains in the lattice as D; no atoms, fuel, or waste are erased.",
            "Void regions, contact perimeter, fixed lesion ROI and lineage tags only measure outcomes; none changes a rate.",
            "The copying-mobility ablation tests renewal in addition to unchanged transport; occupant tags identify birth provenance but do not isolate in-situ repair from paid copies transported into the lesion ROI.",
            "Existing alleles remain individually autocatalytic; repair cannot establish obligatory interdependence of differentiated units.",
            "No maintained Markov blanket, recursive hierarchy, homeostatic set point, collective reproduction, or subjectivity is established.",
            "The experiment does not assert costless preparation or estimate total distribution-level availability."]}
    destination = Path(__file__).with_name("agent-integrity-results.json")
    np.savez_compressed(destination.with_name(output["event_ledger_file"]), **ledgers)
    destination.write_text(json.dumps(output, indent=2)+"\n")
    print(json.dumps({"contrasts": contrasts, "maximum_energy_error": maximum_error}, indent=2), flush=True)
    print("saved", destination, flush=True)


if __name__ == "__main__":
    main()

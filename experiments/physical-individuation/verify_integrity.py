"""Independently replay the prepared-wound assay and recompute contrasts."""
from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from agent_lattice import LocalModel, Parameters


def hamiltonian(model, material, fuel):
    active = (material >= 0).reshape(model.p.side, model.p.side)
    bonds = sum(int((active * np.roll(active, 1, axis=axis)).sum()) for axis in (0, 1))
    return float(model.p.template_energy * active.sum() + model.p.fuel_energy * fuel.sum()
                 - model.p.binding * bonds)


def replay(record, events):
    model = LocalModel(Parameters(**record["parameters"]))
    material, fuel = model.initial()
    center = model.p.side // 2
    roi = np.array([(center + a) * model.p.side + center + b
                    for a, b in ((-1, -1), (-1, 0), (0, -1), (0, 0))])
    pristine = hamiltonian(model, material, fuel)
    if record["wounded"]:
        material[roi] = -1
    initial = hamiltonian(model, material, fuel)
    assert abs(pristine - initial - record["preparation"]["passive_path_heat_to_bath"]) < 1e-9
    assert abs(initial - record["preparation"]["prepared_energy"]) < 1e-9
    tags = np.where(material >= 0, np.arange(model.n), -1)
    time = heat = area = roi_area = 0.
    initial_fuel = int(fuel.sum())
    forward = reverse = births = deaths = 0
    maximum_delta_error = maximum_first_law_error = 0.
    for entry in events:
        t, event, recorded_delta, recorded_heat, recorded_fuel, child, parent = entry
        event, child, parent = int(event), int(child), int(parent)
        assert t >= time
        area += int((material >= 0).sum()) * (t - time)
        roi_area += int((material[roi] >= 0).sum()) * (t - time)
        time = t
        assert model.rates(material, fuel, 0)[0][event] > 0
        before = hamiltonian(model, material, fuel)
        if event < model.edge_count:
            i, j = model.left[event], model.right[event]
            material[i], material[j] = material[j], material[i]
            tags[i], tags[j] = tags[j], tags[i]
        elif event < model.offset_passive:
            i, j = model.left[event - model.edge_count], model.right[event - model.edge_count]
            fuel[i], fuel[j] = fuel[j], fuel[i]
        elif event < model.offset_copy:
            cell, allele = divmod(event - model.offset_passive, 2)
            if material[cell] == -1:
                assert child >= model.n and parent == -1
                material[cell], tags[cell] = allele, child
                births += 1
            else:
                assert material[cell] == allele and child == parent == -1
                material[cell], tags[cell] = -1, -1
                deaths += 1
        else:
            cell, rem = divmod(event - model.offset_copy, 8)
            helper, allele = divmod(rem, 2)
            if material[cell] == -1:
                assert parent == tags[model.neighbors[cell, helper]] and parent >= 0
                assert child >= model.n and fuel[cell] == 1
                material[cell], fuel[cell], tags[cell] = allele, 0, child
                forward += 1
            else:
                assert material[cell] == allele and fuel[cell] == 0 and child == parent == -1
                material[cell], fuel[cell], tags[cell] = -1, 1, -1
                reverse += 1
        after = hamiltonian(model, material, fuel)
        maximum_delta_error = max(maximum_delta_error, abs(after - before - recorded_delta))
        heat -= recorded_delta
        maximum_first_law_error = max(maximum_first_law_error, abs(initial - after - heat))
        assert abs(heat - recorded_heat) < 1e-8 and int(fuel.sum()) == int(recorded_fuel)
        assert int(fuel.sum()) == initial_fuel - forward + reverse
    area += int((material >= 0).sum()) * (model.p.horizon - time)
    roi_area += int((material[roi] >= 0).sum()) * (model.p.horizon - time)
    assert forward == record["forward_copies"] and reverse == record["reverse_copies"]
    assert births == record["passive_thermal_births"] and deaths == record["passive_deaths"]
    assert int((material >= 0).sum()) == record["final_active_templates"]
    assert int((material[roi] >= 0).sum()) == record["final_lesion_ROI_occupancy"]
    assert abs(area - record["integral_active_templates"]) < 1e-7
    assert abs(roi_area - record["integral_lesion_ROI_occupancy"]) < 1e-7
    if record["snapshots"]:
        final = record["snapshots"][-1]
        assert np.array_equal(material.reshape(model.p.side, model.p.side), final["material"])
        assert np.array_equal(fuel.reshape(model.p.side, model.p.side), final["fuel"])
        assert np.array_equal(tags.reshape(model.p.side, model.p.side), final["template_instance_tags"])
    assert maximum_delta_error < 1e-8 and maximum_first_law_error < 1e-8
    return len(events), maximum_delta_error, maximum_first_law_error


def main():
    root = Path(__file__).resolve().parent
    data = json.loads((root / "agent-integrity-results.json").read_text())
    with np.load(root / "agent-integrity-event-ledgers.npz") as logs:
        results = [replay(record, logs[record["event_ledger_key"]])
                   for records in data["cases"].values() for record in records]
    contrasts = {}
    for metric in ("final_active_templates", "final_lesion_ROI_occupancy"):
        values = {}
        for prefix in ("copying", "paired_copying_ablation"):
            sham = data["cases"][prefix + "_sham"]
            wound = data["cases"][prefix + "_wound"]
            values[prefix] = np.array([a[metric] - b[metric] for a, b in zip(sham, wound)])
        difference = values["paired_copying_ablation"] - values["copying"]
        contrasts[metric] = {"difference_in_differences": float(difference.mean()),
                             "monte_carlo_SE": float(difference.std(ddof=1) / np.sqrt(len(difference)))}
    output = {"status": "passed", "trajectories": len(results),
              "events_replayed": sum(item[0] for item in results),
              "maximum_independent_energy_delta_error": max(item[1] for item in results),
              "maximum_independent_first_law_error": max(item[2] for item in results),
              "preparation_material_energy_fuel_and_final_snapshots": "passed",
              "exact_activity_and_ROI_integrals": "passed", "contrasts": contrasts,
              "sha256": {name: sha256((root / name).read_bytes()).hexdigest()
                         for name in ("agent_lattice.py", "agent_integrity.py", "agent-integrity-results.json", "agent-integrity-event-ledgers.npz")},
              "scope": "Partial local repair assay; not a demonstration of a membrane, nested individuality or selective boundary maintenance."}
    (root / "integrity-verification.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()

"""Independent replay of lattice event logs, with exact activity integrals."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from agent_lattice import LocalModel, Parameters


def direct_energy(model, material, fuel, environment):
    p = model.p
    grid = (material >= 0).reshape(p.side, p.side)
    bonds = (grid * np.roll(grid, 1, axis=0)).sum()
    bonds += (grid * np.roll(grid, 1, axis=1)).sum()
    return float(p.template_energy * grid.sum() + p.fuel_energy * fuel.sum()
                 - p.binding * bonds + model.gate_energy * bool(environment))


def replay(record, ledger):
    model = LocalModel(Parameters(**record["parameters"]))
    material, fuel = model.initial()
    tags = np.where(material >= 0, np.arange(model.n), -1)
    environment, time, exact_area, heat = 0, 0., 0., 0.
    initial_h = direct_energy(model, material, fuel, environment)
    maximum_delta_error = maximum_first_law_error = 0.
    counters = {key: 0 for key in record["events"]}
    observation_index = 0
    snapshots = {row["time"]: row for row in record["snapshots"]}

    def check_observation(row):
        assert row["active_templates"] == int((material >= 0).sum())
        assert row["allele1_templates"] == int((material == 1).sum())
        assert row["fuel_remaining"] == int(fuel.sum())
        assert row["environment"] == environment
        assert abs(row["stored_internal_energy"] - direct_energy(model, material, fuel, environment)) < 1e-8
        if row["time"] in snapshots:
            snapshot = snapshots[row["time"]]
            assert np.array_equal(material.reshape(model.p.side, model.p.side), snapshot["material"])
            assert np.array_equal(fuel.reshape(model.p.side, model.p.side), snapshot["fuel"])
            tag_grid = snapshot.get("template_instance_tags", snapshot.get("material_tags"))
            assert np.array_equal(tags.reshape(model.p.side, model.p.side), tag_grid)

    for entry in ledger:
        event_time, event, recorded_delta, cumulative_heat, fuel_remaining, env_after, product_tag, parent_tag = entry
        event, env_after = int(event), int(env_after)
        product_tag, parent_tag = int(product_tag), int(parent_tag)
        assert event_time >= time
        while observation_index < len(record["trajectory"]) and record["trajectory"][observation_index]["time"] <= event_time:
            check_observation(record["trajectory"][observation_index])
            observation_index += 1
        exact_area += float((material >= 0).sum()) * (event_time - time)
        time = event_time
        before = direct_energy(model, material, fuel, environment)
        if event < 0:
            next_environment = {-1: 1, -2: 2, -3: 0}[event]
            assert (environment == 0) != (next_environment == 0)
            environment = next_environment
            counters["gate_entries" if environment else "gate_returns"] += 1
            assert product_tag == parent_tag == -1
        else:
            assert model.rates(material, fuel, environment)[0][event] > 0
            if event < model.edge_count:
                i, j = model.left[event], model.right[event]
                material[i], material[j] = material[j], material[i]
                tags[i], tags[j] = tags[j], tags[i]
                label = "material_swaps"
            elif event < model.offset_passive:
                edge = event - model.edge_count
                i, j = model.left[edge], model.right[edge]
                fuel[i], fuel[j] = fuel[j], fuel[i]
                label = "fuel_swaps"
            elif event < model.offset_copy:
                cell, allele = divmod(event - model.offset_passive, 2)
                if material[cell] == -1:
                    assert product_tag >= model.n and parent_tag == -1
                    material[cell], tags[cell] = allele, product_tag
                    label = "passive_births"
                else:
                    assert material[cell] == allele and product_tag == parent_tag == -1
                    material[cell], tags[cell] = -1, -1
                    label = "passive_deaths"
            else:
                cell, within = divmod(event - model.offset_copy, 8)
                helper_slot, allele = divmod(within, 2)
                helper_cell = model.neighbors[cell, helper_slot]
                helper = material[helper_cell]
                assert helper >= 0
                if material[cell] == -1:
                    assert fuel[cell] == 1 and parent_tag == tags[helper_cell]
                    assert product_tag >= model.n
                    material[cell], fuel[cell], tags[cell] = allele, 0, product_tag
                    label = "erroneous_forward_copies" if helper != allele else "correct_forward_copies"
                else:
                    assert material[cell] == allele and fuel[cell] == 0
                    assert product_tag == parent_tag == -1
                    material[cell], fuel[cell], tags[cell] = -1, 1, -1
                    label = "reverse_copies"
            counters[label] += 1
        after = direct_energy(model, material, fuel, environment)
        delta = after - before
        maximum_delta_error = max(maximum_delta_error, abs(delta - recorded_delta))
        heat -= recorded_delta
        maximum_first_law_error = max(maximum_first_law_error, abs(initial_h - after - heat))
        assert abs(heat - cumulative_heat) < 1e-8
        assert int(fuel.sum()) == int(fuel_remaining) and environment == env_after
        assert np.all(tags[material < 0] == -1)
        assert len(np.unique(tags[material >= 0])) == int((material >= 0).sum())
    exact_area += float((material >= 0).sum()) * (model.p.horizon - time)
    while observation_index < len(record["trajectory"]):
        check_observation(record["trajectory"][observation_index])
        observation_index += 1
    assert counters == record["events"]
    assert maximum_delta_error < 1e-8 and maximum_first_law_error < 1e-8
    return {
        "events_replayed": len(ledger),
        "maximum_direct_energy_delta_error": maximum_delta_error,
        "maximum_first_law_error_all_events": maximum_first_law_error,
        "positive_recorded_reaction_propensities": True,
        "all_observations_and_snapshots_match": True,
        "template_instance_and_copying_ancestry_valid": True,
        "exact_event_integrated_active_templates": exact_area,
        "observation_grid_trapezoid_difference": record["integrated_active_templates_trapezoid"] - exact_area,
        "final_active_templates": record["final_active_templates"],
        "net_fuel_consumed": record["net_fuel_consumed"],
    }


def main():
    root = Path(__file__).resolve().parent
    data = json.loads((root / "agent-lattice-results.json").read_text())
    ledgers = np.load(root / data["event_ledger_file"])
    results = {}
    for label, records in data["runs"].items():
        assert len(records) == 1
        record = records[0]
        results[label] = replay(record, ledgers[record["event_ledger_key"]])
    hashes = {name: hashlib.sha256((root / name).read_bytes()).hexdigest()
              for name in ("agent_lattice.py", "agent-lattice-results.json", data["event_ledger_file"])}
    output = {"status": "passed", "source_and_data_sha256": hashes, "runs": results,
              "scope": "Exact replay of six illustrative trajectories; no population-level advantage, exergy or collective-reproduction inference."}
    (root / "lattice-replay-verification.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()

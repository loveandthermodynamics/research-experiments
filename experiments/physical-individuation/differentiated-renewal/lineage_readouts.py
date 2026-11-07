"""Offline ancestry and finite-time measurements. Never called by rates."""
from __future__ import annotations

from collections import Counter
from typing import Any


def base_state(state: dict[str, Any]) -> dict[str, Any]:
    return state.get("base_state", state)


def measure(run: dict[str, Any], appearance_deadline: float = 20., establishment_size: int = 4) -> dict[str, Any]:
    initial, final = base_state(run["initial_state"]), base_state(run["final_state"])
    living = {int(tag): int(species) for tag, species in zip(initial["instance_tags"], initial["material"]) if tag}
    founders = {tag for tag, species in living.items() if species in (1, 2)}
    # A molecule's observational ancestry follows copying only. Neither reused
    # precursor sites nor unrelated thermal births inherit a removed tag.
    founder_ancestry = set(founders)
    paid_r, paid_l, thermal_r = set(), set(), set()
    focal_ancestry: set[int] = set()
    focal_live: set[int] = set()
    focal = focal_time = focal_extinction = establishment_time = None
    focal_parent_founder_ancestry = None
    first_r1_time = 0. if any(s == 2 for s in living.values()) else None
    first_r1_origin = "initial" if first_r1_time is not None else None
    focal_r1_first_loss = None
    focal_r1_restorations = 0
    founder_clearance = None if founders else 0.
    births = Counter()
    integrals = Counter()
    previous = 0.
    first_inactivity = None if founders else 0.
    max_focal = 0
    work = heat = 0.
    gate_events = 0
    after_clearance_copies = 0
    initial_gate = run["initial_state"].get("gate")
    gate = initial_gate
    gate_times = Counter()
    event_counts = Counter()

    def accumulate(dt: float) -> None:
        assert dt >= -1e-12
        dt = max(0., dt)
        counts = Counter(living.values())
        r = counts[1] + counts[2]
        integrals["R"] += r * dt
        integrals["R0"] += counts[1] * dt
        integrals["R1"] += counts[2] * dt
        integrals["L"] += counts[3] * dt
        integrals["both_types_present_time"] += int(counts[1] > 0 and counts[2] > 0) * dt
        integrals["founder_descendant_R1"] += sum(s == 2 and tag in founder_ancestry for tag, s in living.items()) * dt
        integrals["focal_descendants"] += len(focal_live) * dt
        integrals["focal_R1"] += sum(living.get(tag) == 2 for tag in focal_live) * dt
        integrals["founder_instances_remaining"] += len(founders & living.keys()) * dt
        if gate is not None:
            gate_times[str(gate)] += dt

    for event in run["ledger"]:
        time = float(event["time"])
        accumulate(time - previous)
        previous = time
        kind = event["kind"]
        event_counts[kind] += 1
        work += float(event.get("external_work", 0.))
        heat += float(event.get("heat_to_bath", 0.))
        if kind == "gate_flip":
            gate_events += 1
            gate = event.get("gate_after", event.get("after"))
            assert gate in (-1, 1)
        created, removed = int(event.get("instance_created", 0)), int(event.get("instance_removed", 0))
        parent = int(event.get("parent_instance", 0))
        before, after = int(event.get("before", -1)), int(event.get("after", -1))
        previous_focal_r1 = sum(living.get(tag) == 2 for tag in focal_live)
        if removed:
            assert removed in living
            living.pop(removed)
            focal_live.discard(removed)
        if created:
            assert created not in living and after in (1, 2, 3)
            living[created] = after
            if after == 2 and first_r1_time is None:
                first_r1_time, first_r1_origin = time, kind
            if kind == "copy_R":
                assert parent > 0 and after in (1, 2)
                paid_r.add(created)
                births["paid_R"] += 1
                helper = int(event["helper_species"])
                births[f"paid_R{helper-1}_to_R{after-1}"] += 1
                if parent in founder_ancestry:
                    founder_ancestry.add(created)
                if focal is None and helper == 1 and after == 2 and time <= appearance_deadline:
                    focal, focal_time = created, time
                    focal_parent_founder_ancestry = parent in founder_ancestry
                    focal_ancestry.add(created)
                if created in focal_ancestry or parent in focal_ancestry:
                    focal_ancestry.add(created)
                    focal_live.add(created)
                    births["focal_R0_to_R1" if helper == 1 and after == 2 else "focal_other_copy"] += 1
                elif helper == 1 and after == 2:
                    births["outside_focal_R0_to_R1"] += 1
                if founder_clearance is not None:
                    after_clearance_copies += 1
            elif kind == "produce_L":
                assert after == 3
                paid_l.add(created)
                births["paid_L"] += 1
            elif kind == "passive":
                assert parent == 0
                births["thermal_R" if after in (1, 2) else "thermal_L"] += 1
                if after in (1, 2):
                    thermal_r.add(created)
            else:
                raise AssertionError(f"unexpected birth channel {kind}")
        if focal is not None:
            current_focal_r1 = sum(living.get(tag) == 2 for tag in focal_live)
            if current_focal_r1 == 0 and focal_r1_first_loss is None:
                focal_r1_first_loss = time
            if previous_focal_r1 == 0 and current_focal_r1 > 0 and time > focal_time:
                focal_r1_restorations += 1
            max_focal = max(max_focal, len(focal_live))
            if establishment_time is None and len(focal_live) >= establishment_size:
                establishment_time = time
            if focal_extinction is None and not focal_live:
                focal_extinction = time
        if founder_clearance is None and not founders.intersection(living):
            founder_clearance = time
        if first_inactivity is None and not any(s in (1, 2) for s in living.values()):
            first_inactivity = time

    accumulate(float(run["horizon"]) - previous)
    expected_living = {int(t): int(s) for t, s in zip(final["instance_tags"], final["material"]) if t}
    assert living == expected_living
    assert focal_extinction is None or not focal_live
    final_counts = Counter(living.values())
    initial_h = run["samples"][0]["H"]
    final_h = run["final_summary"]["H"]
    first_law = abs(final_h - initial_h - work + heat)
    assert first_law < 3e-7
    assert abs(final["external_work"] - initial["external_work"] - work) < 3e-7
    assert abs(final["cumulative_heat"] - initial["cumulative_heat"] - heat) < 3e-7
    return {
        "integral_R": integrals["R"], "integral_R0": integrals["R0"],
        "integral_R1": integrals["R1"], "integral_L": integrals["L"],
        "both_types_present_time": integrals["both_types_present_time"],
        "integral_founder_descendant_R1": integrals["founder_descendant_R1"],
        "integral_focal_descendants": integrals["focal_descendants"],
        "integral_focal_R1": integrals["focal_R1"],
        "focal_paid_mutant_appeared_by_deadline": focal is not None,
        "first_R1_appearance_time": first_r1_time, "first_R1_appearance_origin": first_r1_origin,
        "focal_tag": focal, "focal_appearance_time": focal_time,
        "focal_parent_has_initial_founder_ancestry": focal_parent_founder_ancestry,
        "focal_established": establishment_time is not None,
        "establishment_population": "copying descendants over both current alleles",
        "focal_establishment_time": establishment_time,
        "focal_extinction_time": focal_extinction,
        "focal_survives_to_horizon": bool(focal_live), "maximum_focal_descendants": max_focal,
        "focal_R1_first_loss_time": focal_r1_first_loss,
        "focal_R1_uninterrupted_to_horizon": focal is not None and focal_r1_first_loss is None,
        "focal_R1_restorations_after_loss": focal_r1_restorations,
        "final_focal_descendants": len(focal_live),
        "final_focal_R1": sum(living[t] == 2 for t in focal_live),
        "founder_instance_clearance_time": founder_clearance,
        "paid_R_births_after_founder_instance_clearance": after_clearance_copies,
        "final_R": final_counts[1] + final_counts[2], "final_R1": final_counts[2],
        "final_L": final_counts[3], "first_R_inactivity": first_inactivity,
        "first_episode_duration": float(run["horizon"]) if first_inactivity is None else first_inactivity,
        "chemical_work": work, "heat_to_bath": heat, "first_law_error": first_law,
        "gate_events": gate_events, "gate_time": dict(gate_times),
        "birth_counts": dict(births), "event_counts": dict(event_counts),
        "scope": "Offline finite-horizon genealogy. Descendants follow copying through backmutation, never thermal births or recycled precursor. Presence, establishment and functional benefit are distinct; no rate consults these labels."
    }

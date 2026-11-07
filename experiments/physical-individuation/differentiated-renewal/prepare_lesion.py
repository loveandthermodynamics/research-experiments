"""Pure preparation of external fixed-site production-mobility interventions.

No simulation, physical state reset, genealogy-dependent rate or experiment is
performed here. Matching concerns instantaneous forward+reverse propensity,
not net productive current, future hazard, apparatus cost or ownership.
"""
from __future__ import annotations

from dataclasses import asdict
from hashlib import sha256
import math
from pathlib import Path
from typing import Any

import numpy as np

from differentiated_model import DifferentiatedModel, JointState
from interface_model import Event, PRODUCE_L


def _squares(side: int) -> list[dict[str, Any]]:
    result = []
    for row in range(side):
        for col in range(side):
            cells = tuple(sorted({((row + dr) % side) * side + (col + dc) % side
                                  for dr in (0, 1) for dc in (0, 1)}))
            result.append({"row": row, "col": col, "cells": cells,
                           "center": (row + 0.5, col + 0.5)})
    return result


def _distance(first, second, side):
    offsets = [min(abs(a - b), side - abs(a - b)) for a, b in zip(first["center"], second["center"])]
    return float(sum(offsets))


def prepare_lesion(model: DifferentiatedModel, state: JointState,
                   preparation_time: float | None = None, min_center_distance: float = 3.0) -> dict[str, Any]:
    """Return fixed near/remote producer-site mobilities and matching metadata.

    q_near/q_remote below mean fractions REMOVED. Actual retained mobilities
    equal the baseline array multiplied by (1-q) at the selected producer sites.
    Event j is the producer/helper; its target i need not be inside the square.
    Call before re-zeroing the preparation state's observer clock.
    """
    when = state.time if preparation_time is None else float(preparation_time)
    if not math.isfinite(when) or when < 0 or abs(state.time - when) > 1e-9:
        raise ValueError("state must be at the declared preparation time before bookkeeping reset")
    if not math.isfinite(min_center_distance) or min_center_distance < 0:
        raise ValueError("min_center_distance must be finite and nonnegative")
    before = model.state_to_dict(state)
    initial_h = model.energy(state)
    material, side = state.base_state.material, model.params.side
    candidates = _squares(side)
    _, rates = model.rates(state)
    forward_terms = [[] for _ in range(model.n)]
    reverse_terms = [[] for _ in range(model.n)]
    for event, rate in zip(model.events, rates):
        if isinstance(event, Event) and event.kind == PRODUCE_L:
            (forward_terms if event.before == 0 else reverse_terms)[event.j].append(float(rate))
    forward = [math.fsum(values) for values in forward_terms]
    reverse = [math.fsum(values) for values in reverse_terms]

    def description(square):
        cells = square["cells"]
        species = material[list(cells)]
        f, r = math.fsum(forward[i] for i in cells), math.fsum(reverse[i] for i in cells)
        return {"top_left": [square["row"], square["col"]], "sites": list(cells),
                "center": list(square["center"]), "R": int(np.isin(species, (1, 2)).sum()),
                "R0": int((species == 1).sum()), "R1": int((species == 2).sum()),
                "L": int((species == 3).sum()), "D": int((species == 0).sum()),
                "raw_paid_L_forward_propensity": f, "raw_paid_L_reverse_propensity": r,
                "raw_paid_L_total_propensity": math.fsum((f, r)), "raw_paid_L_net_current": f - r}

    descriptions = {(square["row"], square["col"]): description(square) for square in candidates}
    near = min(candidates, key=lambda sq: (-descriptions[(sq["row"], sq["col"])]["R"],
                                          -descriptions[(sq["row"], sq["col"])]["L"], sq["row"], sq["col"]))
    near_info = descriptions[(near["row"], near["col"])]
    near_sites = set(near["cells"])
    eligible = [square for square in candidates if not near_sites.intersection(square["cells"])
                and _distance(near, square, side) >= min_center_distance]
    h_near = near_info["raw_paid_L_total_propensity"]
    remote = min(eligible, key=lambda sq: (
        abs(descriptions[(sq["row"], sq["col"])]["raw_paid_L_total_propensity"] - h_near),
        -_distance(near, sq, side), sq["row"], sq["col"])) if eligible else None
    remote_info = descriptions[(remote["row"], remote["col"])] if remote else None
    h_remote = remote_info["raw_paid_L_total_propensity"] if remote_info else 0.0
    reasons = []
    if near_info["R"] == 0:
        reasons.append("near_square_has_no_templates")
    if remote is None:
        reasons.append("no_disjoint_remote_square_at_declared_distance")
    if h_near == 0:
        reasons.append("zero_near_initial_paid_L_propensity")
    if h_remote == 0:
        reasons.append("zero_remote_initial_paid_L_propensity")
    matched = not reasons
    q_near = min(1.0, h_remote / h_near) if matched else 0.0
    q_remote = min(1.0, h_near / h_remote) if matched else 0.0
    baseline = np.asarray(model.production_site_mobility)
    near_mask, remote_mask = baseline.copy(), baseline.copy()
    near_mask[list(near["cells"])] *= 1 - q_near
    if remote:
        remote_mask[list(remote["cells"])] *= 1 - q_remote

    def removed(info, fraction):
        if info is None:
            return {"forward": 0.0, "reverse": 0.0, "total": 0.0, "net_current": 0.0}
        f = fraction * info["raw_paid_L_forward_propensity"]
        r = fraction * info["raw_paid_L_reverse_propensity"]
        return {"forward": f, "reverse": r, "total": f + r, "net_current": f - r}

    removed_near, removed_remote = removed(near_info, q_near), removed(remote_info, q_remote)
    assert abs(removed_near["total"] - removed_remote["total"]) <= 1e-12 * max(1, h_near, h_remote)
    neighbors = model.driven.base.neighbors
    near_targets = {int(i) for j in near["cells"] for i in neighbors[j]}
    remote_targets = {int(i) for j in remote["cells"] for i in neighbors[j]} if remote else set()
    assert model.state_to_dict(state) == before, "lesion preparation must not mutate the physical snapshot"
    assert model.energy(state) == initial_h
    metadata = {
        "preparation_time": when, "gate": state.gate, "initial_stored_energy_H": initial_h,
        "initial_material_units": model.n, "initial_fuel": int(state.base_state.fuel.sum()),
        "params": asdict(model.params), "flags": asdict(model.flags), "response": asdict(model.response),
        "response_mode": model.response_mode, "transport_mode": model.transport_mode,
        "near": near_info, "remote": remote_info, "minimum_periodic_center_distance": min_center_distance,
        "periodic_center_distance_metric": "Manhattan (L1)",
        "actual_periodic_center_distance": _distance(near, remote, side) if remote else None,
        "eligible_remote_squares": len(eligible), "selection_hazard_difference": abs(h_remote - h_near) if remote else None,
        "near_selection": "Maximum R count, then maximum L count, then lexicographic top-left row/column.",
        "remote_selection": "Disjoint 2x2 squares at distance>=threshold: closest absolute total paid-L propensity, then farther periodic center, then lexicographic row/column; no rounded score ties.",
        "matched": matched, "unmatched_reasons": reasons,
        "near_removed_mobility_fraction": q_near, "remote_removed_mobility_fraction": q_remote,
        "near_removed_initial_propensity": removed_near, "remote_removed_initial_propensity": removed_remote,
        "removed_total_matching_error": removed_near["total"] - removed_remote["total"],
        "possible_L_target_neighborhood_overlap": sorted(near_targets & remote_targets),
        "site_paid_L_forward_propensities": forward, "site_paid_L_reverse_propensities": reverse,
        "state_energy_jump": 0.0, "intervention_apparatus_work_modeled": False,
        "fixture_qualification": "Fixed external kinetic fixture. State/H are unchanged, but switching/construction work is excluded and is not asserted to be zero.",
        "scope": "Initial total propensity matching is not matching net production, future flux, accepted transport, genealogical ownership or individual boundaries. Producer j lies inside the square; target i may be outside. No-op unmatched preparations remain part of the assay.",
        "helper_source_sha256": sha256(Path(__file__).read_bytes()).hexdigest()}
    return {"baseline_production_site_mobility": tuple(float(v) for v in baseline),
            "near_production_site_mobility": tuple(float(v) for v in near_mask),
            "remote_production_site_mobility": tuple(float(v) for v in remote_mask), "metadata": metadata}

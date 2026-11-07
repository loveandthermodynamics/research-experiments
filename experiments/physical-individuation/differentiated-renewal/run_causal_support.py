"""Prospective fixed-location support assay; execution requires separate review.

The protocol is frozen before this study's outcomes. Observer-selected squares
only set a fixed external kinetic intervention and readouts; no group or lineage
label enters the natural model's rates. No prior physical record is modified.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import gzip
from hashlib import sha256
import json
import math
from pathlib import Path
import time
from typing import Any

import numpy as np
from scipy.stats import t

from differentiated_model import (
    DifferentiatedModel, DriveParams, Flags, Params, ResponseParams,
    simulate_differentiated,
)
from cohort_support_readouts import support_readouts
from prepare_lesion import prepare_lesion


ROOT = Path(__file__).resolve().parent
PROTOCOL = ROOT / "causal-support-protocol.json"
RESULTS = ROOT / "causal-support-results.json"
MODES = ("opposed", "aligned")
LESIONS = ("sham", "near", "remote")
PRIMARY = ("opposed_sham_minus_near", "opposed_remote_minus_near")
SECONDARY = ("remote_minus_near_response_interaction", "sham_opposed_minus_aligned")
OUTCOME = "outcome_integral_local_R"
READOUT_FIELDS = (
    OUTCOME, "integral_global_R", "integral_local_L", "integral_global_cohort_R",
    "integral_local_cohort_R", "integral_global_cohort_R1",
    "integral_local_R_after_original_cohort_clearance",
    "integral_cohort_R_after_original_cohort_clearance", "initial_cohort_size",
    "final_cohort_size", "final_cohort_R1", "final_local_R", "final_global_R",
    "cohort_forward_L", "cohort_reverse_L", "cohort_net_L", "near_site_forward_L",
    "near_site_reverse_L", "cohort_copies", "cohort_copies_after_original_clearance",
    "cohort_outward_crossings", "cohort_inward_crossings", "chemical_work",
    "heat_to_bath", "total_preparation_and_continuation_work",
    "total_preparation_and_continuation_heat",
)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def write_once(path: Path, value: Any) -> None:
    encoded = json.dumps(value, indent=2, allow_nan=False) + "\n"
    if path.exists():
        if path.read_text() != encoded:
            raise RuntimeError(f"frozen file differs: {path}")
    else:
        path.write_text(encoded)


def source_hashes() -> dict[str, str]:
    files = [ROOT / name for name in (
        "run_causal_support.py", "causal-support-design.md", "cohort_support_readouts.py",
        "prepare_lesion.py", "differentiated_model.py", "model-spec.md", "lineage_readouts.py",
    )]
    files.extend((ROOT.parent / "driven-spatial/driven_model.py",
                  ROOT.parent / "interface-maintenance/interface_model.py"))
    return {str(path.relative_to(ROOT)) if path.is_relative_to(ROOT)
            else "../" + str(path.relative_to(ROOT.parent)): digest(path) for path in files}


def freeze_protocol() -> None:
    selected_path = ROOT / "confirmation-protocol.json"
    selected = load(selected_path)
    if selected["settings"] != [{"asymmetry": .8, "flip_rate": .025, "initial_gate": None}]:
        raise RuntimeError("selected response setting differs from the authorized fixed setting")
    protocol = dict(
        study="external_fixed_location_causal_support", stage="fixed_confirmation",
        source_hashes=source_hashes(), selected_parameter_protocol_file=selected_path.name,
        selected_parameter_protocol_sha256=digest(selected_path),
        params=selected["params"], drive=selected["drive"], response=selected["settings"][0],
        flags=asdict(Flags(produce_L=True, protect_R=True, copy_R=True)),
        seeds=list(range(1070001, 1070065)), preparation_horizon=100.,
        continuation_horizon=80., post_seed_offset=3000000, sample_dt=2.,
        max_events=5000000, max_gate_events=1000000,
        preparation_response_mode="opposed", transport_mode="local",
        response_modes=list(MODES), lesions=list(LESIONS),
        selection=dict(
            near="Periodic 2x2 square: maximum R count, then L count, then lexicographic top-left.",
            remote="Disjoint 2x2 squares at periodic Manhattan center distance >=3: minimum absolute difference in initial total paid-L propensity; then larger distance; then lexicographic top-left.",
            minimum_periodic_manhattan_center_distance=3.,
            removed_fraction_near="min(1, Hremote/Hnear)",
            removed_fraction_remote="min(1, Hnear/Hremote)",
            matching="Forward plus reverse initial total propensity only; not net current, strength, future flux or outcome.",
            unmatched="Keep every preparation; absent templates, zero propensity or no eligible remote square gives two no-op lesions and an explicit flag.",
        ),
        initial_state_policy="Same realized material, fuel, instance tags, next_tag and gate in all six branches. Reset only time and cumulative heat/work/event counters; never reset physical gate or select successful preparations.",
        gate_pairing="Same post seed origin+3000000 across all six arms; the independent autonomous gate clock has identical realized gate history.",
        primary_outcome=OUTCOME,
        primary_region="Near 2x2 square plus its periodic nearest-neighbor shell, frozen once per preparation.",
        primary_contrasts=list(PRIMARY),
        primary_inference=dict(
            n=64, degrees_freedom=63, interval="Two-sided t63 95% confidence interval for paired seed means.",
            joint_claim="Both primary lower endpoints must be strictly positive. No outcome-dependent choice, filtering, pooling, extra scan or rescue endpoint.",
        ),
        primary_bootstrap_sensitivity=dict(
            replicates=10000, seed=1080002, confidence=.95,
            resampling="Paired seeds; one shared index matrix for both primary contrasts; numpy percentile quantiles with default linear interpolation.",
            rule="Sensitivity only; cannot replace the declared two-sided t63 decision.",
        ),
        named_secondary_contrasts=list(SECONDARY),
        secondary_lineage="Prepared near-square R instances and their paid-copy descendants, through backmutation; unrelated thermal births do not join. Track local/global counts, crossings, original-instance clearance and post-clearance copying. Ancestry does not enter rates.",
        cost_policy="Retain the complete preparation and each continuation ledger. Report initial supplied inventory and stored H, preparation and branch chemical work/heat, and their totals. Initial resources are supplied; stored H is not distribution-level exergy. Shared reservoir parameters do not imply equal realized work.",
        fixture_scope="Response-map and fixed-site paired-mobility changes are external counterfactual interventions. State/H remain unchanged, but apparatus construction, measurement and switching cost are excluded and not asserted to be zero.",
        supported_inference="A positive joint result supports continuing production at these fixed nearby locations during this finite protocol, conditional on prepared states and supplied chemistry. It does not establish ownership, a self-bounded individual, organizational memory, aging, reproduction or ethical participant status.",
        expected_preparations=64, expected_continuations=384,
        history_directory="causal-support", result_file=RESULTS.name,
    )
    write_once(PROTOCOL, protocol)
    print(json.dumps({"protocol_file": PROTOCOL.name, "protocol_sha256": digest(PROTOCOL),
                      "runner_sha256": digest(Path(__file__).resolve()),
                      "ensemble_executed": False}, indent=2))


def estimate(values: list[float]) -> dict[str, Any]:
    data = np.asarray(values, dtype=float)
    if data.ndim != 1 or len(data) < 2 or not np.all(np.isfinite(data)):
        raise ValueError("finite paired values from at least two preparations required")
    n = len(data)
    mean = float(data.mean())
    se = float(data.std(ddof=1) / np.sqrt(n))
    critical = float(t.ppf(.975, n - 1))
    return dict(n=n, mean=mean, monte_carlo_standard_error=se, degrees_freedom=n - 1,
                t_95_interval=[mean - critical * se, mean + critical * se])


def physical_accounting(run: dict[str, Any]) -> dict[str, Any]:
    model = DifferentiatedModel(
        Params(**run["params"]), Flags(**run["flags"]), DriveParams(**run["drive"]),
        ResponseParams(**run["response"]), run["response_mode"], run["transport_mode"],
        run["production_site_mobility"],
    )
    initial = model.state_from_dict(run["initial_state"])
    final = model.state_from_dict(run["final_state"])
    initial_h, final_h = model.energy(initial), model.energy(final)
    work = math.fsum(float(event.get("external_work", 0.)) for event in run["ledger"])
    heat = math.fsum(float(event["heat_to_bath"]) for event in run["ledger"])
    first_law_error = abs(final_h - initial_h - work + heat)
    if first_law_error >= 3e-7:
        raise RuntimeError("trajectory first law failed")
    if abs(final.base_state.external_work - work) >= 3e-7 or abs(final.base_state.cumulative_heat - heat) >= 3e-7:
        raise RuntimeError("trajectory cumulative counters disagree with complete ledger")
    material = initial.base_state.material
    return dict(
        initial_stored_energy_H=initial_h, final_stored_energy_H=final_h,
        initial_material_units=model.n, initial_fuel=int(initial.base_state.fuel.sum()),
        final_fuel=int(final.base_state.fuel.sum()),
        initial_species_counts={name: int((material == species).sum())
                                for name, species in (("D", 0), ("R0", 1), ("R1", 2), ("L", 3))},
        initial_gate=initial.gate, final_gate=final.gate,
        chemical_work=work, heat_to_bath=heat, first_law_error=first_law_error,
        material_reservoir_transport_events=len(run["ledger"]) - final.gate_event_count,
        gate_events=final.gate_event_count, ledger_events=len(run["ledger"]),
    )


def reset_bookkeeping(model: DifferentiatedModel, endpoint: dict[str, Any]):
    state = model.state_from_dict(endpoint)
    before = model.state_to_dict(state)
    state.time = 0.
    state.base_state.cumulative_heat = 0.
    state.base_state.external_work = 0.
    state.base_state.event_count = 0
    state.gate_event_count = 0
    after = model.state_to_dict(state)
    expected = json.loads(json.dumps(before))
    expected["base_state"].update(time=0., cumulative_heat=0., external_work=0., event_count=0)
    expected["gate_event_count"] = 0
    if after != expected:
        raise RuntimeError("bookkeeping reset changed physical state")
    return state


def save_record(path: Path, run: dict[str, Any]) -> dict[str, str]:
    if path.exists():
        raise RuntimeError(f"physical record already exists: {path}")
    with gzip.open(path, "wt", encoding="utf-8") as stream:
        json.dump(run, stream, separators=(",", ":"), allow_nan=False)
    return {"trajectory_file": str(path.relative_to(ROOT)), "trajectory_sha256": digest(path)}


def gate_trace(run: dict[str, Any]):
    return [(event["time"], event["gate_after"]) for event in run["ledger"]
            if event["kind"] == "gate_flip"]


def summarize(rows: list[dict[str, Any]], protocol: dict[str, Any]):
    groups = []
    paired = {}
    for row in rows:
        key = (row["response_mode"], row["lesion"])
        target = paired.setdefault(row["origin_seed"], {})
        if key in target:
            raise RuntimeError("duplicate continuation arm")
        target[key] = row
    expected = {(mode, lesion) for mode in MODES for lesion in LESIONS}
    if sorted(paired) != protocol["seeds"] or any(set(arms) != expected for arms in paired.values()):
        raise RuntimeError("missing or extra seed/arm cells")
    for mode in MODES:
        for lesion in LESIONS:
            arm = [paired[seed][(mode, lesion)] for seed in protocol["seeds"]]
            groups.append(dict(response_mode=mode, lesion=lesion,
                               statistics={field: estimate([row[field] for row in arm])
                                           for field in READOUT_FIELDS},
                               original_cohort_present=sum(row["original_cohort_present"] for row in arm),
                               original_cohort_cleared_after_presence=sum(row["original_cohort_cleared_after_presence"] for row in arm),
                               final_zero_global_R=sum(row["final_global_R"] == 0 for row in arm)))
    differences = {name: [] for name in PRIMARY + SECONDARY}
    seed_values = []
    for seed in protocol["seeds"]:
        arm = paired[seed]
        y = {key: row[OUTCOME] for key, row in arm.items()}
        values = (
            y[("opposed", "sham")] - y[("opposed", "near")],
            y[("opposed", "remote")] - y[("opposed", "near")],
            y[("opposed", "remote")] - y[("opposed", "near")]
            - y[("aligned", "remote")] + y[("aligned", "near")],
            y[("opposed", "sham")] - y[("aligned", "sham")],
        )
        seed_values.append({"origin_seed": seed, **dict(zip(PRIMARY + SECONDARY, values))})
        for name, value in zip(PRIMARY + SECONDARY, values):
            differences[name].append(value)
    contrasts = [dict(contrast=name, field=OUTCOME, primary=name in PRIMARY, **estimate(values))
                 for name, values in differences.items()]
    spec = protocol["primary_bootstrap_sensitivity"]
    vectors = np.asarray([differences[name] for name in PRIMARY], dtype=float).T
    rng = np.random.default_rng(spec["seed"])
    indices = rng.integers(0, len(vectors), size=(spec["replicates"], len(vectors)))
    bootmeans = vectors[indices].mean(axis=1)
    bootstrap = [dict(contrast=name, replicates=spec["replicates"], seed=spec["seed"],
                      percentile_95_interval=np.quantile(bootmeans[:, i], [.025, .975]).tolist(),
                      sensitivity_only=True) for i, name in enumerate(PRIMARY)]
    return groups, contrasts, bootstrap, seed_values


def execute() -> None:
    if RESULTS.exists():
        raise RuntimeError("causal support results are already frozen")
    protocol = load(PROTOCOL)
    if protocol["source_hashes"] != source_hashes():
        raise RuntimeError("bound source changed after protocol freeze")
    selected_path = ROOT / protocol["selected_parameter_protocol_file"]
    if digest(selected_path) != protocol["selected_parameter_protocol_sha256"]:
        raise RuntimeError("selected primary parameter protocol changed")
    started = time.monotonic()
    directory = ROOT / protocol["history_directory"]
    directory.mkdir(exist_ok=True)
    params, drive = Params(**protocol["params"]), DriveParams(**protocol["drive"])
    response, flags = ResponseParams(**protocol["response"]), Flags(**protocol["flags"])
    prep_model = DifferentiatedModel(params, flags, drive, response,
                                    protocol["preparation_response_mode"], protocol["transport_mode"])
    preparations, rows = [], []
    for seed in protocol["seeds"]:
        prep = simulate_differentiated(
            params, seed, protocol["preparation_horizon"], flags=flags, drive=drive,
            response=response, response_mode=protocol["preparation_response_mode"],
            transport_mode=protocol["transport_mode"], sample_dt=protocol["sample_dt"],
            max_events=protocol["max_events"], max_gate_events=protocol["max_gate_events"],
        )
        prep_physics = physical_accounting(prep)
        endpoint = prep_model.state_from_dict(prep["final_state"])
        plan = prepare_lesion(
            prep_model, endpoint, preparation_time=protocol["preparation_horizon"],
            min_center_distance=protocol["selection"]["minimum_periodic_manhattan_center_distance"],
        )
        prep_ref = save_record(directory / f"s{seed}-preparation.json.gz", prep)
        prepared = reset_bookkeeping(prep_model, prep["final_state"])
        expected_initial = prep_model.state_to_dict(prepared)
        prep_entry = dict(origin_seed=seed, **prep_ref, physical_accounting=prep_physics,
                          branch_initial_state=expected_initial, intervention_plan=plan,
                          final_zero_R=plan["metadata"]["near"]["R"] == 0)
        preparations.append(prep_entry)
        mobility = {"sham": plan["baseline_production_site_mobility"],
                    "near": plan["near_production_site_mobility"],
                    "remote": plan["remote_production_site_mobility"]}
        common_gates = None
        for mode in protocol["response_modes"]:
            for lesion in protocol["lesions"]:
                run = simulate_differentiated(
                    params, seed + protocol["post_seed_offset"], protocol["continuation_horizon"],
                    flags=flags, drive=drive, response=response, response_mode=mode,
                    transport_mode=protocol["transport_mode"], production_site_mobility=mobility[lesion],
                    initial_state=prepared, sample_dt=protocol["sample_dt"],
                    max_events=protocol["max_events"], max_gate_events=protocol["max_gate_events"],
                )
                if run["initial_state"] != expected_initial:
                    raise RuntimeError("branch physical initial states differ")
                gates = gate_trace(run)
                if common_gates is None:
                    common_gates = gates
                elif common_gates != gates:
                    raise RuntimeError("paired branches have different independent gate histories")
                physics = physical_accounting(run)
                observations = support_readouts(run, plan["metadata"]["near"]["sites"])
                # Immutable annotations for independent readout replay; none
                # of these observer-defined fields were available to rates.
                run["near_sites"] = plan["metadata"]["near"]["sites"]
                run["lesion_metadata"] = plan["metadata"]
                run["causal_support_context"] = dict(
                    origin_seed=seed, lesion=lesion,
                    preparation_trajectory_file=prep_ref["trajectory_file"],
                    preparation_trajectory_sha256=prep_ref["trajectory_sha256"],
                    protocol_file=PROTOCOL.name, protocol_sha256=digest(PROTOCOL),
                    observer_annotations_only=True,
                )
                total_work = prep_physics["chemical_work"] + physics["chemical_work"]
                total_heat = prep_physics["heat_to_bath"] + physics["heat_to_bath"]
                total_error = abs(physics["final_stored_energy_H"] - prep_physics["initial_stored_energy_H"]
                                  - total_work + total_heat)
                if total_error >= 5e-7:
                    raise RuntimeError("preparation plus branch first law failed")
                ref = save_record(directory / f"s{seed}-{mode}-{lesion}.json.gz", run)
                rows.append(dict(
                    origin_seed=seed, post_seed=seed + protocol["post_seed_offset"],
                    response_mode=mode, lesion=lesion, matched=plan["metadata"]["matched"],
                    preparation_trajectory_file=prep_ref["trajectory_file"],
                    preparation_trajectory_sha256=prep_ref["trajectory_sha256"], **ref,
                    **observations, **physics,
                    preparation_chemical_work=prep_physics["chemical_work"],
                    preparation_heat_to_bath=prep_physics["heat_to_bath"],
                    original_initial_stored_energy_H=prep_physics["initial_stored_energy_H"],
                    total_preparation_and_continuation_work=total_work,
                    total_preparation_and_continuation_heat=total_heat,
                    total_preparation_and_continuation_first_law_error=total_error,
                ))
        print(f"causal support: seed={seed} preparation and six branches complete", flush=True)
    groups, contrasts, bootstrap, seed_values = summarize(rows, protocol)
    primary = [contrast for contrast in contrasts if contrast["primary"]]
    if len(preparations) != protocol["expected_preparations"] or len(rows) != protocol["expected_continuations"]:
        raise RuntimeError("physical record count differs from predeclared protocol")
    if source_hashes() != protocol["source_hashes"]:
        raise RuntimeError("bound source changed during execution")
    result = dict(
        study=protocol["study"], source_hashes=source_hashes(), protocol_file=PROTOCOL.name,
        protocol_sha256=digest(PROTOCOL), selected_parameter_protocol_sha256=digest(selected_path),
        all_preparations_included=True, preparations=preparations, runs=rows,
        groups=groups, paired_contrasts=contrasts, paired_seed_differences=seed_values,
        primary_results=primary, primary_bootstrap_sensitivity=bootstrap,
        declared_joint_claim_supported=all(contrast["t_95_interval"][0] > 0 for contrast in primary),
        matching_summary=dict(
            matched=sum(prep["intervention_plan"]["metadata"]["matched"] for prep in preparations),
            unmatched=sum(not prep["intervention_plan"]["metadata"]["matched"] for prep in preparations),
            final_zero_R_preparations=sum(prep["final_zero_R"] for prep in preparations),
            maximum_removed_total_matching_error=max(abs(prep["intervention_plan"]["metadata"]["removed_total_matching_error"]) for prep in preparations),
        ),
        physical_totals=dict(
            unique_preparation_records=len(preparations), continuation_records=len(rows),
            unique_ledger_events=sum(prep["physical_accounting"]["ledger_events"] for prep in preparations)
                                 + sum(row["ledger_events"] for row in rows),
            preparation_work_counted_once=sum(prep["physical_accounting"]["chemical_work"] for prep in preparations),
            continuation_work_all_counterfactual_arms=sum(row["chemical_work"] for row in rows),
            maximum_trajectory_first_law_error=max([prep["physical_accounting"]["first_law_error"] for prep in preparations]
                                                  + [row["first_law_error"] for row in rows]),
            maximum_preparation_plus_branch_first_law_error=max(row["total_preparation_and_continuation_first_law_error"] for row in rows),
            scope="Preparation resources are shared across counterfactual branches, not physically cloned in a single realization. Summed work across alternative arms is only a simulation-data total. Bath heat alone is not total irreversible entropy production.",
        ),
        elapsed_seconds=time.monotonic() - started,
        scope=protocol["supported_inference"], fixture_scope=protocol["fixture_scope"],
    )
    write_once(RESULTS, result)
    print(json.dumps({"results_file": RESULTS.name, "results_sha256": digest(RESULTS),
                      "primary_results": primary,
                      "declared_joint_claim_supported": result["declared_joint_claim_supported"]}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("freeze-protocol", "run"))
    action = parser.parse_args().action
    if action == "freeze-protocol":
        freeze_protocol()
    else:
        execute()

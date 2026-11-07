"""Bounded software/physics checks, not a formation or size experiment."""
from dataclasses import asdict, replace
import gzip
from hashlib import sha256
import json
import math
from pathlib import Path
import time

import numpy as np

from active_phase_model import A, B, S, CATALYZED, PASSIVE, Model, Params, simulate


ROOT = Path(__file__).resolve().parent


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def hashes():
    return {name: digest(ROOT / name) for name in ("active_phase_model.py", "model-spec.md", "software_smoke.py")}


def direct_energy(material, enzyme, params):
    """Independent full-state energy expression, not a local delta shortcut."""
    side = params.side
    b = (material.reshape(side, side) == B).astype(float)
    bonds = np.sum(b * np.roll(b, -1, axis=0)) + np.sum(b * np.roll(b, -1, axis=1))
    return float(params.epsilon_B * b.sum() - params.J_BB * bonds
                 - params.J_EB * np.sum(enzyme.reshape(side, side) * b))


def main():
    target = ROOT / "software-smoke.json"
    if target.exists():
        raise RuntimeError("software smoke already frozen; do not rerun or overwrite")
    bound = hashes()
    rng = np.random.default_rng(745000)
    pairs = 0
    max_delta_error, max_log_error, max_rate_error = 0., 0., 0.
    for side in (4, 6):
        for mu in (0., 12.):
            p = Params(side=side, material_count=side * side // 3, enzyme_count=2, mu=mu)
            model = Model(p)
            for trial in range(12):
                state = model.initial(1000 + trial)
                occupied = np.flatnonzero(state.material != S)
                state.material[occupied] = rng.choice((A, B), len(occupied))
                events, rates = model.rates(state)
                assert model.energy(state) == direct_energy(state.material, state.enzyme, p)
                for k in np.flatnonzero(rates > 0):
                    scalar = model.event_rate(state, int(k))
                    max_rate_error = max(max_rate_error, abs(scalar - rates[k]))
                    candidate = model.copy_state(state)
                    h0 = direct_energy(candidate.material, candidate.enzyme, p)
                    record = model.apply(candidate, int(k))
                    h1 = direct_energy(candidate.material, candidate.enzyme, p)
                    max_delta_error = max(max_delta_error, abs(h1 - h0 - record["delta_H"]))
                    reverse = int(model.reverse_index[k])
                    q_reverse = model.event_rate(candidate, reverse)
                    log_ratio = math.log(scalar) - math.log(q_reverse)
                    expected = (record["chemical_work"] - record["delta_H"]) / p.theta
                    max_log_error = max(max_log_error, abs(log_ratio - expected))
                    assert abs(log_ratio - expected) < 3e-13
                    model.apply(candidate, reverse)
                    assert np.array_equal(candidate.material, state.material)
                    assert np.array_equal(candidate.enzyme, state.enzyme)
                    assert np.array_equal(candidate.material_tags, state.material_tags)
                    model.validate_state(candidate)
                    pairs += 1
    assert max_delta_error < 1e-12 and max_rate_error < 1e-13
    cycle_model = Model(Params(side=4, material_count=5, enzyme_count=2))
    cycle_state = cycle_model.initial(745004)
    # Construct a software fixture with colocated A and E, preserving inventories.
    site = int(np.flatnonzero(cycle_state.material == A)[0])
    if cycle_state.enzyme[site] == 0:
        other = int(np.flatnonzero(cycle_state.enzyme == 1)[0])
        cycle_state.enzyme[site], cycle_state.enzyme[other] = 1, 0
    original = cycle_model.copy_state(cycle_state)
    passive = next(k for k, e in enumerate(cycle_model.events)
                   if e.kind == PASSIVE and e.i == site and e.before == A)
    catalytic = next(k for k, e in enumerate(cycle_model.events)
                     if e.kind == CATALYZED and e.i == site and e.before == B)
    cycle_log_affinity = 0.
    for k in (passive, catalytic):
        q = cycle_model.event_rate(cycle_state, k)
        cycle_model.apply(cycle_state, k)
        qr = cycle_model.event_rate(cycle_state, int(cycle_model.reverse_index[k]))
        cycle_log_affinity += math.log(q) - math.log(qr)
    assert np.array_equal(cycle_state.material, original.material)
    assert np.array_equal(cycle_state.enzyme, original.enzyme)
    assert np.array_equal(cycle_state.material_tags, original.material_tags)
    assert cycle_state.net_fuel_to_waste == 1
    assert abs(cycle_log_affinity - cycle_model.params.mu / cycle_model.params.theta) < 1e-13
    directory = ROOT / "software-smoke"
    directory.mkdir(exist_ok=True)
    traces = []
    for p, seed in ((Params(side=6, material_count=11, enzyme_count=2), 745001),
                    (Params(side=6, material_count=11, enzyme_count=2, mu=0.), 745002),
                    (Params(), 745003)):
        start = time.monotonic()
        run = simulate(p, seed, 1., sample_dt=5.)
        elapsed = time.monotonic() - start
        assert all(m != B for m in run["initial_state"]["material"])
        path = directory / f"side{p.side}-mu{p.mu:g}-s{seed}.json.gz"
        if path.exists():
            raise RuntimeError("software trace already exists")
        with gzip.open(path, "wt", encoding="utf-8") as stream:
            json.dump(run, stream, separators=(",", ":"), allow_nan=False)
        traces.append(dict(params=asdict(p), seed=seed, horizon=1., elapsed_seconds=elapsed,
                           trajectory_file=str(path.relative_to(ROOT)), trajectory_sha256=digest(path),
                           events=len(run["ledger"]), proposals=run["proposal_statistics"]["proposals"],
                           first_law_error=run["first_law_error"],
                           initial_zero_B=True, final_summary=run["final_summary"]))
    if hashes() != bound:
        raise RuntimeError("source changed during software checks")
    result = dict(
        checks_passed=True, source_hashes=bound, numpy_version=np.__version__,
        random_software_states=48, enabled_reverse_pairs=pairs,
        maximum_full_energy_delta_error=max_delta_error,
        maximum_log_detailed_balance_error=max_log_error,
        maximum_vector_scalar_rate_error=max_rate_error,
        chemical_cycle_log_affinity=cycle_log_affinity,
        chemical_cycle_expected_affinity=cycle_model.params.mu / cycle_model.params.theta,
        traces=traces,
        scope="Bounded software checks only. Random microscopic fixtures and three one-unit traces do not test phase formation, selected size, coarsening arrest, sustained turnover or individuation; exclude them from prospective morphology-study totals.",
    )
    target.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

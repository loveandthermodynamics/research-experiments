# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy", "scipy", "matplotlib"]
# ///
"""Finite chemical-resource model with reversible repair and assembling contacts.

Three labeled molecular sites retain distinct response templates in both active
and inactive states. They are not biological individuals. All contacts are
ordinary reversible conformation/binding states. No group/diversity score enters
the generator. Reactions are thermodynamically realizable coarse transitions,
not an atomistically identified chemistry or an origin of response templates.

The system is closed to matter, coupled to a thermal bath, and starts with a
finite fuel inventory. A degenerate reversible environmental gate changes only
transition-state prefactors, never well energies or forward/reverse affinities.
It is a joint detailed-balance Markov model, not prescribed external work.
The initial fuel disequilibrium is the source of transient maintenance/repair.

Run: uv run --no-project finite_resource_groups.py
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import expm_multiply
from scipy.special import gammaln, logsumexp, xlogy


N = 3
EDGES = ((0, 1), (0, 2), (1, 2))
NCONFIG = 2**N
NBOND = 2**len(EDGES)


@dataclass(frozen=True)
class Parameters:
    fuel_packets: int = 12
    thermal: float = 0.4
    active_energy: float = 4.0
    fuel_energy: float = 7.0
    contact_energy: float = 0.5
    active_binding: float = 0.8
    contact_rate: float = 0.5
    repair_rate: float = 0.08
    contagion_rate: float = 0.3
    permeability: float = 1.0
    damage_good: float = 0.01
    damage_bad: float = 1.0
    environment_rate: float = 0.2
    correlated_shock_rate: float = 0.0


def active_energy(config: int, bonds: int, p: Parameters) -> float:
    bits = [(config >> j) & 1 for j in range(N)]
    bound = [(bonds >> k) & 1 for k in range(len(EDGES))]
    return (p.active_energy * sum(bits) + p.contact_energy * sum(bound)
            - p.active_binding * sum(on * bits[i] * bits[j]
                                     for on, (i, j) in zip(bound, EDGES)))


def largest_component(bonds: int) -> int:
    parts = [{i} for i in range(N)]
    for k, (i, j) in enumerate(EDGES):
        if (bonds >> k) & 1:
            a = next(s for s in parts if i in s)
            b = next(s for s in parts if j in s)
            if a is not b:
                a.update(b)
                parts.remove(b)
    return max(map(len, parts))


def make_model(p: Parameters, responses: tuple[int, ...]):
    assert len(responses) == N and all(s in (-1, 0, 1) for s in responses)
    kfuel = p.fuel_packets
    size = NCONFIG * NBOND * (kfuel + 1) * 2
    def index(a, b, f, e):
        return (((a * NBOND + b) * (kfuel + 1) + f) * 2 + e)

    h_a = np.array([[active_energy(a, b, p) for b in range(NBOND)]
                    for a in range(NCONFIG)])
    log_pi_a = -h_a.reshape(-1) / p.thermal
    log_pi_a -= logsumexp(log_pi_a)
    fuel_count = np.arange(kfuel + 1)
    log_comb = (gammaln(kfuel + 1) - gammaln(fuel_count + 1)
                - gammaln(kfuel - fuel_count + 1))
    h_f = p.fuel_energy * fuel_count - p.thermal * log_comb
    log_pi_f = -h_f / p.thermal
    log_pi_f -= logsumexp(log_pi_f)
    log_pi = (log_pi_a[:, None, None] + log_pi_f[None, :, None]
              + np.full((1, 1, 2), -np.log(2))).reshape(-1)

    rows, cols, rates = [], [], []
    out = np.zeros(size)
    repair_forward = np.zeros(size)
    repair_reverse = np.zeros(size)
    max_balance_error = 0.0

    # Add paired transitions together; the pairing proves local detailed balance
    # and makes a quantitative log-ratio audit possible without searching edges.
    def pair(i, j, qij, qji):
        nonlocal max_balance_error
        if qij == 0 and qji == 0:
            return
        assert qij > 0 and qji > 0
        max_balance_error = max(max_balance_error,
                                abs(np.log(qij / qji) - (log_pi[j] - log_pi[i])))
        rows.extend((i, j)); cols.extend((j, i)); rates.extend((qij, qji))
        out[i] += qij; out[j] += qji

    for a in range(NCONFIG):
        bits = np.array([(a >> j) & 1 for j in range(N)])
        for b in range(NBOND):
            neighbors = [[] for _ in range(N)]
            for edge, (i, j) in enumerate(EDGES):
                if (b >> edge) & 1:
                    neighbors[i].append(j); neighbors[j].append(i)
            for f in range(kfuel + 1):
                for e in range(2):
                    i = index(a, b, f, e)
                    if e == 0:
                        pair(i, index(a, b, f, 1), p.environment_rate,
                             p.environment_rate)
                    # A shared vulnerability is a collective denaturation
                    # channel. Each vulnerability class has the same shock rate;
                    # hence every active site has the same marginal exposure.
                    # Its rare thermal reverse is retained, not discarded.
                    if p.correlated_shock_rate:
                        for response_class in set(responses):
                            mask = sum(1 << j for j, r in enumerate(responses)
                                       if r == response_class)
                            aa = a & ~mask
                            if aa != a:
                                drop = h_a[a, b] - h_a[aa, b]
                                assert drop > 0
                                pair(i, index(aa, b, f, e),
                                     p.correlated_shock_rate,
                                     p.correlated_shock_rate * np.exp(-drop / p.thermal))
                    for edge in range(len(EDGES)):
                        if not (b >> edge) & 1:
                            bb = b | (1 << edge)
                            dh = h_a[a, bb] - h_a[a, b]
                            pair(i, index(a, bb, f, e),
                                 p.contact_rate * np.exp(-dh / p.thermal),
                                 p.contact_rate)
                    for target in range(N):
                        if bits[target]:
                            continue
                        aa = a | (1 << target)
                        dh = h_a[aa, b] - h_a[a, b]
                        live = sum(bits[j] for j in neighbors[target]) / (N - 1)
                        dead = sum(1 - bits[j] for j in neighbors[target]) / (N - 1)
                        # Response 0 is a constant, same-mean generalist control.
                        if responses[target] == 0:
                            damage = (p.damage_good + p.damage_bad) / 2
                        else:
                            bad = responses[target] == (1 if e == 0 else -1)
                            damage = p.damage_bad if bad else p.damage_good
                        passive = damage + p.permeability * p.contagion_rate * dead
                        pair(i, index(aa, b, f, e),
                             passive * np.exp(-dh / p.thermal), passive)
                        if f and live:
                            j = index(aa, b, f - 1, e)
                            forward = p.permeability * p.repair_rate * f * live
                            reverse = (p.permeability * p.repair_rate
                                       * np.exp(-(p.fuel_energy - dh) / p.thermal)
                                       * (kfuel - f + 1) * live)
                            pair(i, j, forward, reverse)
                            repair_forward[i] += forward
                            repair_reverse[j] += reverse

    rows.extend(range(size)); cols.extend(range(size)); rates.extend(-out)
    q = coo_matrix((rates, (rows, cols)), shape=(size, size)).tocsr()
    assert max_balance_error < 8e-13
    assert np.max(np.abs(q.sum(axis=1))) < 2e-12
    pi = np.exp(log_pi)
    assert np.max(np.abs(q.T @ pi)) < 2e-14
    f_observable = np.tile(np.repeat(fuel_count, 2), NCONFIG * NBOND)
    assert np.max(np.abs(q @ f_observable + repair_forward - repair_reverse)) < 2e-11
    initial = np.zeros(size)
    initial[index(NCONFIG - 1, 0, kfuel, 0)] = 0.5
    initial[index(NCONFIG - 1, 0, kfuel, 1)] = 0.5
    return q, initial, log_pi, log_pi_a, log_pi_f, repair_forward, repair_reverse, max_balance_error


def area(values, times):
    return float(np.sum((values[1:] + values[:-1]) * np.diff(times) / 2))


def run(p: Parameters, responses: tuple[int, ...], horizon=30.0, samples=121):
    q, initial, lp, la, lf, rf, rr, ldb_error = make_model(p, responses)
    times = np.linspace(0, horizon, samples)
    path = expm_multiply(q.T, initial, start=0, stop=horizon, num=samples)
    assert path.min() > -1e-12
    path = np.maximum(path, 0)
    assert np.max(np.abs(path.sum(axis=1) - 1)) < 2e-10
    shaped = path.reshape(samples, NCONFIG * NBOND, p.fuel_packets + 1, 2)
    pa = shaped.sum(axis=(2, 3))
    pf = shaped.sum(axis=(1, 3))
    pe = shaped.sum(axis=(1, 2))
    total_b = p.thermal * (np.sum(xlogy(path, path), axis=1) - path @ lp)
    active_b = p.thermal * (np.sum(xlogy(pa, pa), axis=1) - pa @ la)
    fuel_b = p.thermal * (np.sum(xlogy(pf, pf), axis=1) - pf @ lf)
    env_b = p.thermal * (np.sum(xlogy(pe, pe), axis=1) + np.log(2))
    mutual = total_b - active_b - fuel_b - env_b
    assert mutual.min() > -2e-10
    assert np.diff(total_b).max() < 1e-9
    active_n = np.repeat([int(a).bit_count() for a in range(NCONFIG)], NBOND)
    live = pa @ active_n
    fuel = pf @ np.arange(p.fuel_packets + 1)
    bonds = pa @ np.tile([int(b).bit_count() for b in range(NBOND)], NCONFIG)
    group_size = pa @ np.tile([largest_component(b) for b in range(NBOND)], NCONFIG)
    forward = path @ rf
    reverse = path @ rr
    net_integral = area(forward - reverse, times)
    # This numerical quadrature check is separate from exact generator identities.
    assert abs((fuel[0] - fuel[-1]) - net_integral) < 0.03
    values = {
        "responses": list(responses), "horizon": horizon,
        "initial_total_available_free_energy": float(total_b[0]),
        "initial_active_available_free_energy": float(active_b[0]),
        "active_exergy_time_integral": area(active_b, times),
        "active_site_time_integral": area(live, times),
        "fuel_consumed_net": float(fuel[0] - fuel[-1]),
        "total_free_energy_dissipated": float(total_b[0] - total_b[-1]),
        "final_active_sites": float(live[-1]),
        "final_active_exergy": float(active_b[-1]),
        "final_unused_fuel": float(fuel[-1]),
        "average_bonds": area(bonds, times) / horizon,
        "average_largest_contact_component": area(group_size, times) / horizon,
        "maximum_local_detailed_balance_log_error": ldb_error,
        "net_fuel_balance_quadrature_error": abs((fuel[0] - fuel[-1]) - net_integral),
        "initial_zero_contacts": True,
        "no_ecological_reward": True,
        "state_count": q.shape[0],
    }
    tracks = {"time": times, "active_exergy": active_b, "active_sites": live,
              "fuel_remaining": fuel, "total_exergy": total_b,
              "contacts": bonds, "largest_component": group_size,
              "forward_repair_flux": forward, "reverse_repair_flux": reverse}
    return values, tracks


def main():
    p = Parameters()
    results = []
    tracks = {}
    for label, response in (("monoculture", (1, 1, 1)),
                            ("response_diverse", (1, 1, -1)),
                            ("same_mean_generalist", (0, 0, 0))):
        val, trajectory = run(p, response)
        val["label"] = label
        results.append(val); tracks[label] = trajectory
        print(json.dumps(val), flush=True)
    assert np.ptp([r["initial_total_available_free_energy"] for r in results]) < 1e-12
    assert np.ptp([r["initial_active_available_free_energy"] for r in results]) < 1e-12
    scan = []
    for permeability in (0.0, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0):
        pp = replace(p, permeability=permeability)
        for label, response in (("monoculture", (1, 1, 1)),
                                ("response_diverse", (1, 1, -1))):
            val, _ = run(pp, response, samples=241)
            val.update(label=label, permeability=permeability)
            scan.append(val)
        print(json.dumps({"permeability": permeability,
                          "active_exergy_integrals": [v["active_exergy_time_integral"]
                                                     for v in scan[-2:]]}), flush=True)
    output = {"status": "finite-resource reversible coarse chemical model",
              "parameters": asdict(p), "baseline": results, "permeability_scan": scan,
              "limitations": [
                  "Fixed site templates and response types; no new hereditary types.",
                  "Thermal revival is present; zero activity is not absorbing extinction.",
                  "Only a finite fuel resource; no eternal growth or cosmic theorem.",
                  "Contact components are assemblies, not demonstrated collective individuals.",
                  "Optimal fixed permeability is not evolution of permeability.",
                  "Barrier gating is a reversible kinetic idealization, not supplied external work.",
              ], "thermodynamic_checks_passed": True}
    folder = Path(__file__).parent
    (folder / "finite-resource-results.json").write_text(json.dumps(output, indent=2) + "\n")
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.1), constrained_layout=True)
    for label, tr in tracks.items():
        axes[0].plot(tr["time"], tr["active_exergy"], label=label.replace("_", " "))
        axes[1].plot(tr["time"], tr["fuel_remaining"])
    axes[0].set(xlabel="Time", ylabel="Available free energy of active subsystem",
                title="Finite fuel pays for repair")
    axes[0].legend(fontsize=8)
    axes[1].set(xlabel="Time", ylabel="Fuel packets remaining", title="No resource reset")
    for label in ("monoculture", "response_diverse"):
        rows = [v for v in scan if v["label"] == label]
        axes[2].plot([v["permeability"] for v in rows],
                     [v["active_exergy_time_integral"] for v in rows], "o-",
                     label=label.replace("_", " "))
    axes[2].set(xscale="symlog", xlabel="Exchange / contagion permeability",
                ylabel="Integral of active available free energy",
                title="Grouping has physical consequences")
    axes[2].legend(fontsize=8)
    fig.savefig(folder / "finite-resource-dynamics.png", dpi=170)


if __name__ == "__main__":
    main()

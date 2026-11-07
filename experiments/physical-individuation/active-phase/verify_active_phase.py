"""Independent scalar physical reference for active-phase lattice chemistry.

No phase labels, clusters, preferred sizes or observer masks enter this
implementation. Model comparisons import the production module only locally;
the scalar Hamiltonian, rates, work and state updates are independent.
"""
from __future__ import annotations

import copy
from bisect import bisect_right
from dataclasses import asdict
import gzip
from hashlib import sha256
import itertools
import json
import math
from pathlib import Path
import random
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    path = Path(path)
    if path.suffix == ".gz":
        with gzip.open(path, "rt") as stream:
            return json.load(stream)
    return json.loads(path.read_text())


def close(a, b, tolerance=2e-8):
    assert math.isfinite(float(a)) and math.isfinite(float(b)), (a, b)
    assert abs(float(a)-float(b)) <= tolerance, (a, b, abs(float(a)-float(b)))


def compare(expected, observed, path=""):
    if isinstance(expected, dict):
        for key, value in expected.items():
            compare(value, observed[key], path+"."+key)
    elif isinstance(expected, (list, tuple)):
        assert len(expected) == len(observed), path
        for i, (a, b) in enumerate(zip(expected, observed)):
            compare(a, b, path+f"[{i}]")
    elif expected is None or isinstance(expected, (str, bool)):
        assert expected == observed, (path, expected, observed)
    else:
        close(expected, observed)


class Reference:
    def __init__(self, params):
        self.p = dict(params)
        self.side, self.n = self.p["side"], self.p["side"]**2
        self.neighbors, self.edges = [], []
        for i in range(self.n):
            row, col = divmod(i, self.side)
            self.neighbors.append((((row-1) % self.side)*self.side+col,
                                   ((row+1) % self.side)*self.side+col,
                                   row*self.side+(col-1) % self.side,
                                   row*self.side+(col+1) % self.side))
            self.edges.extend(((i, row*self.side+(col+1) % self.side),
                               (i, ((row+1) % self.side)*self.side+col)))
        self.incident = [set() for _ in range(self.n)]
        for index, (i, j) in enumerate(self.edges):
            self.incident[i].add(index)
            self.incident[j].add(index)
        self.events = [("material_swap", i, j, -1, -1) for i, j in self.edges]
        self.events += [("enzyme_swap", i, j, -1, -1) for i, j in self.edges]
        self.events += [("passive", i, -1, before, after) for i in range(self.n) for before, after in ((1, 2), (2, 1))]
        self.events += [("catalyzed", i, -1, before, after) for i in range(self.n) for before, after in ((2, 1), (1, 2))]
        self.lookup = {event: index for index, event in enumerate(self.events)}
        self.reverse_index = [index if event[0].endswith("swap") else
                              self.lookup[(event[0], event[1], event[2], event[4], event[3])]
                              for index, event in enumerate(self.events)]

    def energy(self, state):
        b = [int(a == 2) for a in state["material"]]
        return (self.p["epsilon_B"]*sum(b)
                -self.p["J_BB"]*sum(b[i]*b[j] for i, j in self.edges)
                -self.p["J_EB"]*sum(e*a for e, a in zip(state["enzyme"], b)))

    def changes(self, state, event):
        kind, i, j, before, after = event
        if kind == "material_swap":
            return {i: state["material"][j], j: state["material"][i]}, {}
        if kind == "enzyme_swap":
            return {}, {i: state["enzyme"][j], j: state["enzyme"][i]}
        return {i: after}, {}

    def delta(self, state, event):
        changes, enzyme_changes = self.changes(state, event)
        m, e = state["material"], state["enzyme"]
        b = lambda i: int(m[i] == 2)
        changed_b = lambda i: int(changes.get(i, m[i]) == 2)
        delta = self.p["epsilon_B"]*sum(changed_b(i)-b(i) for i in changes)
        affected = set()
        for i in changes:
            affected.update(self.incident[i])
        delta -= self.p["J_BB"]*sum(changed_b(i)*changed_b(j)-b(i)*b(j)
                                    for index in affected for i, j in (self.edges[index],))
        delta -= self.p["J_EB"]*sum(enzyme_changes.get(i, e[i])*changed_b(i)-e[i]*b(i)
                                    for i in set(changes) | set(enzyme_changes))
        return delta

    def work(self, event):
        return self.p["mu"]*(1 if event[3] == 2 else -1) if event[0] == "catalyzed" else 0.

    def fuel_count(self, event):
        return (1 if event[3] == 2 else -1) if event[0] == "catalyzed" else 0

    def mobility(self, event):
        return self.p[{"material_swap": "material_mobility", "enzyme_swap": "enzyme_mobility",
                       "passive": "passive_rate", "catalyzed": "catalytic_rate"}[event[0]]]

    def rate(self, state, event):
        kind, i, j, before, after = event
        if kind == "material_swap":
            if state["material"][i] == state["material"][j]:
                return 0.
        elif kind == "enzyme_swap":
            if state["enzyme"][i] == state["enzyme"][j]:
                return 0.
        elif state["material"][i] != before or (kind == "catalyzed" and state["enzyme"][i] != 1):
            return 0.
        exponent = (self.work(event)-self.delta(state, event))/self.p["theta"]
        return self.mobility(event)*math.exp(min(0., exponent))

    def validate(self, state):
        m, e, tags = state["material"], state["enzyme"], state["material_tags"]
        assert len(m) == len(e) == len(tags) == self.n
        assert all(a in (0, 1, 2) for a in m) and all(a in (0, 1) for a in e)
        assert sum(a != 0 for a in m) == self.p["material_count"]
        assert sum(e) == self.p["enzyme_count"]
        assert all((a == 0) == (tag == 0) for a, tag in zip(m, tags))
        live = [tag for tag in tags if tag]
        assert len(live) == len(set(live)) == self.p["material_count"]

    def apply_physical(self, state, index):
        """State and thermodynamic increment only, independent of record naming."""
        event = self.events[index]
        assert self.rate(state, event) > 0, ("disabled event", event)
        kind, i, j, before, after = event
        delta, work, fuel = self.delta(state, event), self.work(event), self.fuel_count(event)
        material_changes, enzyme_changes = self.changes(state, event)
        if kind == "material_swap":
            state["material_tags"][i], state["material_tags"][j] = state["material_tags"][j], state["material_tags"][i]
        for site, value in material_changes.items():
            state["material"][site] = value
        for site, value in enzyme_changes.items():
            state["enzyme"][site] = value
        state["cumulative_work"] += work
        state["cumulative_heat"] += work-delta
        state["net_fuel_to_waste"] += fuel
        state["event_count"] += 1
        return dict(delta_H=delta, chemical_work=work, heat_to_bath=work-delta, fuel_to_waste=fuel)

    def summary(self, state, include_state=False):
        m, e = state["material"], state["enzyme"]
        b = [a == 2 for a in m]
        out = dict(time=state["time"], H=self.energy(state), S=m.count(0), A=m.count(1), B=m.count(2),
                   enzymes=sum(e), enzymes_on_B=sum(x*int(y) for x, y in zip(e, b)),
                   BB_contacts=sum(b[i] and b[j] for i, j in self.edges),
                   chemical_work=state["cumulative_work"], heat_to_bath=state["cumulative_heat"],
                   net_fuel_to_waste=state["net_fuel_to_waste"], events=state["event_count"])
        if include_state:
            out.update(material=m[:], enzyme=e[:], material_tags=state["material_tags"][:])
        return out

    def apply(self, state, index):
        event = self.events[index]
        kind, i, j, _, _ = event
        sites = [i] if j == -1 else [i, j]

        def snapshot():
            return {key: [state[key][site] for site in sites]
                    for key in ("material", "enzyme", "material_tags")}

        before = snapshot()
        increments = self.apply_physical(state, index)
        return dict(time=state["time"], event_index=index, kind=kind, i=i, j=j,
                    before=before, after=snapshot(), **increments,
                    cumulative_work=state["cumulative_work"], cumulative_heat=state["cumulative_heat"],
                    net_fuel_to_waste=state["net_fuel_to_waste"], event_count=state["event_count"])


def verify_trajectory(run, full_proposals=True):
    """Replay every accepted event; optionally every rejected proposal too."""
    import numpy as np
    ref, state = Reference(run["params"]), copy.deepcopy(run["initial_state"])
    ref.validate(state)
    initial_h = ref.energy(state)
    assert all(state[key] == 0 for key in ("time", "cumulative_work", "cumulative_heat", "net_fuel_to_waste", "event_count"))
    original_tags = set(state["material_tags"])-{0}
    samples, sample_index, counts = run["samples"], 0, {}
    direct_h, maximum_first_law_error = initial_h, 0.

    def samples_until(until):
        nonlocal sample_index
        while sample_index < len(samples) and samples[sample_index]["time"] <= until:
            observed = samples[sample_index]
            expected = ref.summary(state, True)
            expected["time"] = observed["time"]
            compare(expected, observed)
            sample_index += 1

    def accepted(saved, proposal_number=None, rejected_since_last=None):
        nonlocal direct_h, maximum_first_law_error
        assert state["time"] <= saved["time"] <= run["horizon"]
        samples_until(saved["time"])
        state["time"] = saved["time"]
        index = saved["event_index"]
        assert 0 <= index < len(ref.events)
        before_rate = ref.rate(state, ref.events[index])
        expected = ref.apply(state, index)
        compare(expected, saved)
        if proposal_number is not None:
            assert saved["proposal_number"] == proposal_number
            assert saved["rejected_since_last_accepted"] == rejected_since_last
        reverse_rate = ref.rate(state, ref.events[ref.reverse_index[index]])
        close(math.log(before_rate/reverse_rate), expected["heat_to_bath"]/ref.p["theta"], 2e-11)
        direct_h += expected["delta_H"]
        error = abs(direct_h-initial_h-state["cumulative_work"]+state["cumulative_heat"])
        maximum_first_law_error = max(maximum_first_law_error, error)
        assert error < 1e-6
        close(state["cumulative_work"], ref.p["mu"]*state["net_fuel_to_waste"], 2e-8)
        counts[expected["kind"]] = counts.get(expected["kind"], 0)+1
        if state["event_count"] % 97 == 0:
            close(ref.energy(state), direct_h, 2e-8)
            ref.validate(state)

    stats = run["proposal_statistics"]
    bounds = [ref.p[key] for key in ("material_mobility", "enzyme_mobility", "passive_rate", "catalytic_rate")]
    cumulative = list(itertools.accumulate(bounds))
    total_bound = 2*ref.n*cumulative[-1]
    close(total_bound, stats["total_bound"])
    if full_proposals:
        rng = np.random.default_rng(np.random.SeedSequence([run["seed"], 74503]))
        proposal = accepted_index = rejected = 0
        proposed_by_family, accepted_by_family = [0]*4, [0]*4
        clock, censored = 0., None
        while clock < run["horizon"]:
            next_time = clock+float(rng.exponential(1./total_bound)) if total_bound else math.inf
            if next_time > run["horizon"]:
                censored = next_time if math.isfinite(next_time) else None
                break
            clock = next_time
            proposal += 1
            family = bisect_right(cumulative, rng.random()*cumulative[-1])
            offset = int(rng.integers(2*ref.n))
            index = family*2*ref.n+offset
            proposed_by_family[family] += 1
            probability = ref.rate(state, ref.events[index])/bounds[family]
            if rng.random() < probability:
                assert accepted_index < len(run["ledger"]), "Missing accepted proposal"
                saved = run["ledger"][accepted_index]
                assert index == saved["event_index"]
                close(clock, saved["time"], 1e-12)
                accepted(saved, proposal, rejected)
                accepted_by_family[family] += 1
                accepted_index += 1
                rejected = 0
            else:
                rejected += 1
        assert accepted_index == len(run["ledger"])
        families = ("material_swap", "enzyme_swap", "passive", "catalyzed")
        compare(dict(total_bound=total_bound, proposals=proposal, accepted=accepted_index,
            rejected=proposal-accepted_index, terminal_rejected_since_last_accepted=rejected,
            proposed_by_family=dict(zip(families, proposed_by_family)),
            accepted_by_family=dict(zip(families, accepted_by_family)), terminal_censored_proposal_time=censored), stats)
    else:
        previous_proposal, rejected_sum = 0, 0
        for saved in run["ledger"]:
            assert saved["proposal_number"] == previous_proposal+saved["rejected_since_last_accepted"]+1
            rejected_sum += saved["rejected_since_last_accepted"]
            previous_proposal = saved["proposal_number"]
            accepted(saved)
        assert stats["accepted"] == len(run["ledger"])
        assert stats["rejected"] == rejected_sum+stats["terminal_rejected_since_last_accepted"]
        assert stats["proposals"] == stats["accepted"]+stats["rejected"]
    samples_until(run["horizon"])
    assert sample_index == len(samples)
    state["time"] = run["horizon"]
    compare(state, run["final_state"])
    compare(ref.summary(state), run["final_summary"])
    close(ref.energy(state), direct_h, 2e-8)
    ref.validate(state)
    assert original_tags == set(state["material_tags"])-{0}
    final_error = abs(ref.energy(state)-initial_h-state["cumulative_work"]+state["cumulative_heat"])
    close(final_error, run["first_law_error"])
    return dict(seed=run["seed"], accepted_events_replayed=len(run["ledger"]),
        proposals=stats["proposals"], rejected_proposals=stats["rejected"], exact_all_proposals_regenerated=full_proposals,
        samples_checked=len(samples), event_counts=counts, initial_H=initial_h, final_H=ref.energy(state),
        chemical_work=state["cumulative_work"], heat_to_bath=state["cumulative_heat"],
        net_fuel_to_waste=state["net_fuel_to_waste"], maximum_running_first_law_error=maximum_first_law_error,
        final_first_law_error=final_error, conserved_material_tags=True,
        scope="Complete accepted physical ledger replay; exact rejected-clock replay is stated separately. Signed fuel counts refer to supplied ideal chemostats, not a finite closed fuel inventory.")


def check_model(params, configurations=3, seed=91037):
    """Production imports are limited to comparison, never physical replay."""
    from active_phase_model import Model, Params
    rng = random.Random(seed)
    maxima = dict(H=0., delta_H=0., rate=0., LDB_log=0., cycle_affinity_log=0.)
    comparisons = pairs = cases = 0
    # Both supplied thermodynamic controls retain the same material inventories.
    parameter_sets = [dict(params), {**params, "mu": 0.}, {**params, "J_EB": 0.}]
    for p in parameter_sets:
        ref, model = Reference(p), Model(Params(**p))
        actual_events = [tuple(asdict(event).values()) for event in model.events]
        assert actual_events == ref.events
        assert list(model.reverse_index) == ref.reverse_index
        for _ in range(configurations):
            # Arbitrary random physical states with exactly conserved inventories.
            sites = list(range(ref.n)); rng.shuffle(sites)
            m = [0]*ref.n
            for i in sites[:p["material_count"]]:
                m[i] = rng.choice((1, 2))
            rng.shuffle(sites)
            e = [int(i in set(sites[:p["enzyme_count"]])) for i in range(ref.n)]
            tags, tag = [], 0
            for a in m:
                tag += int(a != 0)
                tags.append(tag if a else 0)
            state = dict(material=m, enzyme=e, material_tags=tags, time=0., cumulative_work=0.,
                         cumulative_heat=0., net_fuel_to_waste=0, event_count=0)
            converted = model.state_from_dict(state)
            close(model.energy(converted), ref.energy(state))
            maxima["H"] = max(maxima["H"], abs(model.energy(converted)-ref.energy(state)))
            _, rates = model.rates(converted)
            for index, event in enumerate(ref.events):
                expected = ref.rate(state, event)
                close(expected, rates[index], 2e-12)
                close(expected, model.event_rate(converted, index), 2e-12)
                maxima["rate"] = max(maxima["rate"], abs(expected-float(rates[index])))
                comparisons += 1
                if expected <= 0:
                    continue
                delta = ref.delta(state, event)
                close(delta, model.delta_energy(converted, index))
                maxima["delta_H"] = max(maxima["delta_H"], abs(delta-model.delta_energy(converted, index)))
                after = copy.deepcopy(state)
                increment = ref.apply_physical(after, index)
                close(ref.energy(after)-ref.energy(state), delta)
                reverse = ref.events[ref.reverse_index[index]]
                ratio = math.log(expected/ref.rate(after, reverse))
                error = abs(ratio-(ref.work(event)-delta)/p["theta"])
                assert error < 2e-11
                maxima["LDB_log"] = max(maxima["LDB_log"], error)
                ref.validate(after)
                pairs += 1
                # Explicit two-channel physical cycle, not merely a biased well.
                if event[0] == "passive" and event[3:] == (1, 2) and e[event[1]]:
                    catalyst = ("catalyzed", event[1], -1, 2, 1)
                    passive_reverse = ("passive", event[1], -1, 2, 1)
                    catalyst_reverse = ("catalyzed", event[1], -1, 1, 2)
                    affinity = math.log(expected/ref.rate(after, passive_reverse)) + math.log(ref.rate(after, catalyst)/ref.rate(state, catalyst_reverse))
                    error = abs(affinity-p["mu"]/p["theta"])
                    assert error < 2e-11
                    maxima["cycle_affinity_log"] = max(maxima["cycle_affinity_log"], error)
            altered = model.copy_state(converted)
            altered.material_tags *= 19
            assert (model.rates(altered)[1] == rates).all(), "Conserved tags must not affect physical rates"
            cases += 1
    return dict(cases=cases, scalar_rates_checked=comparisons, enabled_reverse_pairs_checked=pairs,
                maximum_errors=maxima, zero_mu_control=True, no_enzyme_affinity_control=True,
                tags_nonfunctional=True,
                scope="Conditional specified-generator thermodynamics. Heat/theta is medium entropy flow; system Shannon entropy is not supplied by an individual-path heat ledger. No phase formation, size selection or individuation is proved by this check.")


def observe_independent(material, enzyme):
    """Union-find clusters and half-spectrum observables, not production BFS."""
    import numpy as np
    array, enzymes = np.asarray(material), np.asarray(enzyme)
    if array.ndim == 1:
        side = math.isqrt(len(array))
        assert side*side == len(array)
        array, enzymes = array.reshape(side, side), enzymes.reshape(side, side)
    assert array.ndim == 2 and array.shape[0] == array.shape[1] and enzymes.shape == array.shape
    side, n = array.shape[0], array.size
    b = (array == 2).ravel()
    count, ne = int(b.sum()), int(enzymes.sum())
    parent = list(range(n))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    bonds = 0
    for i in range(n):
        if not b[i]:
            continue
        row, col = divmod(i, side)
        for j in (row*side+(col+1) % side, ((row+1) % side)*side+col):
            if b[j]:
                bonds += 1
                a, c = root(i), root(j)
                if a != c:
                    parent[c] = a
    sizes = {}
    for i in range(n):
        if b[i]:
            sizes[root(i)] = sizes.get(root(i), 0)+1
    sizes = sorted(sizes.values(), reverse=True)
    ein = int(sum(int(e) for e, sticky in zip(enzymes.ravel(), b) if sticky))
    eout = ne-ein
    pin, pout = (ein/count if count else None), (eout/(n-count) if count < n else None)
    centered = b.reshape(side, side).astype(float)-count/n
    transform = np.fft.rfft2(centered)
    multiplicity = np.full(transform.shape[1], 2.)
    multiplicity[0] = 1.
    if side % 2 == 0:
        multiplicity[-1] = 1.
    powers = np.real(transform*np.conj(transform))*multiplicity[None, :]/n
    krow, kcol = np.fft.fftfreq(side)*2*math.pi, np.fft.rfftfreq(side)*2*math.pi
    radius = np.hypot(krow[:, None], kcol[None, :])
    positive = radius > 0
    low = positive & (radius <= math.pi/2+1e-12)
    mass = float(powers[positive].sum())
    lowmass = float(powers[low].sum())
    firstmoment = float((powers[low]*radius[low]).sum())
    neighbor = 2*bonds/(4*count) if count else None
    return dict(side=side, n_B=count, B_fraction=count/n, n_E=ne, cluster_sizes=sizes,
                largest_cluster=sizes[0] if sizes else 0,
                largest_B_fraction=sizes[0]/count if count else None,
                mass_weighted_cluster_size=sum(size**2 for size in sizes)/count if count else None,
                clusters_ge_9=sum(size >= 9 for size in sizes), B_mass_in_clusters_ge_9=sum(size for size in sizes if size >= 9),
                neighbor_B_fraction=neighbor,
                excess_neighbor_B_fraction=neighbor-(count-1)/(n-1) if count else None,
                E_on_B=ein, E_on_nonB=eout, E_density_B=pin, E_density_nonB=pout,
                E_enrichment=pin/pout if pin is not None and pout else None,
                E_on_B_excess=ein-ne*count/n,
                E_on_B_excess_fraction=ein/ne-count/n if ne else 0.,
                spectral_wavelength_low_k=2*math.pi*lowmass/firstmoment if firstmoment else None,
                low_k_power_fraction=lowmass/mass if mass else None,
                spectral_parseval_residual=float(powers.sum())-count*(1-count/n))


def check_morphology(seed=319011):
    import numpy as np
    from morphology import observe
    rng = np.random.default_rng(seed)
    cases, max_error = 0, 0.
    for side in (4, 8, 16, 32):
        for case in range(6):
            material = rng.integers(0, 3, (side, side))
            enzyme = rng.integers(0, 2, (side, side))
            if case == 0:
                material[:] = 0
            elif case == 1:
                material[:] = 2
            elif case == 2:
                material[:] = 0
                material[0, 0] = material[-1, 0] = 2
            elif case == 3:
                material[:] = 0
                material[:3, :3] = 2
                enzyme[:] = 0
                enzyme[0, 0] = 1  # Perfect enzyme partition; enrichment is undefined, excess is positive.
            expected = observe_independent(material, enzyme)
            observed = observe(material, enzyme)
            compare(expected, observed)
            for key, value in expected.items():
                if isinstance(value, float) and value is not None:
                    max_error = max(max_error, abs(value-observed[key]))
            cases += 1
    return dict(cases=cases, cluster_algorithm="independent union-find",
                spectrum_algorithm="independent real half-spectrum with explicit conjugate multiplicities",
                maximum_absolute_error=max_error,
                scope="Software observables only; the low-k cutoff and cluster>=9 threshold are declared analysis choices, not generated physical regulation.")

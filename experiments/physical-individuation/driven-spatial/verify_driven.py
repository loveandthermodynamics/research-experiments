"""Independent scalar energy, rate, work and trajectory audit of driven lattice.

Uses the frozen independent lattice geometry and Hamiltonian implementation;
never uses the driven production apply/delta implementation for ledger replay.
"""
from __future__ import annotations

import copy
import argparse
from dataclasses import asdict, replace
import gzip
from hashlib import sha256
import itertools
import json
import math
from pathlib import Path
import random
import statistics
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
BASE = ROOT.parent / "interface-maintenance"
sys.path.insert(0, str(BASE))
from verify_interface import Reference, Params, Flags, close, compare_record, compare_state
sys.path.insert(0, str(ROOT))
from driven_model import DrivenModel, DriveParams

REFERENCE_SHA256 = "ccbaac8114d5f4ad48cf1bd655de7c186ddca0f4e99197322ed3fc4e43aef0b6"
BASE_MODEL_SHA256 = "08c5d9948e84d62ddc77a0265e5e91e805b2badf47bc516542a36fa918201927"
DRIVEN_MODEL_SHA256 = "3dd7b498f0f1b31739aa706966c220320e10aaecc31a68523ce4c1a0866b228f"
PILOT_PROTOCOL_SHA256 = "e5ebd316d7532dea68f4a662139ffdc725c2f43a7d5ac8511ddfb16ec980f6b8"
CONFIRMATION_PROTOCOL_SHA256 = "4bf7630fcdcb17b4ab91288a5d7f33ac1cc5d847add3fa77b264ea97ebaef30a"


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    path = Path(path)
    if path.suffix == ".gz":
        with gzip.open(path, "rt") as stream:
            return json.load(stream)
    return json.loads(path.read_text())


class DrivenReference(Reference):
    def __init__(self, params, flags=None, drive=None, transport_mode="local"):
        super().__init__(params, flags)
        self.mode = transport_mode
        self.drive = drive or DriveParams()
        if self.drive.mu is None:
            self.probability = self.drive.fuel_probability
            self.mu = params.E_F+params.theta*math.log(self.probability/(1-self.probability))
        else:
            self.mu = self.drive.mu
            x = (self.mu-params.E_F)/params.theta
            self.probability = 1/(1+math.exp(-x)) if x >= 0 else math.exp(x)/(1+math.exp(x))
        self.base_count = len(self.events)
        self.events += [("reservoir_exchange", i, -1, before, 1-before)
                        for i in range(self.n) for before in (0, 1)]
        self.events += [("mixed_L_swap", i, j, -1, -1)
                        for i, j in itertools.combinations(range(self.n), 2)]
        lookup = {e: k for k, e in enumerate(self.events)}
        self.reverse_index = []
        for k, (kind, i, j, before, after) in enumerate(self.events):
            if kind in (0, 1, 2, "reservoir_exchange"):
                self.reverse_index.append(lookup[(kind, i, j, after, before)])
            else:
                self.reverse_index.append(k)

    def delta(self, state, event):
        kind, i, j, before, after = event
        if kind == "reservoir_exchange":
            return self.p.E_F*(after-before)
        if kind == "mixed_L_swap":
            return super().delta(state, (3, i, j, -1, -1))
        return super().delta(state, event)

    def rate(self, state, event):
        kind, i, j, before, after = event
        if kind == "reservoir_exchange":
            if state["fuel"][i] != before:
                return 0.
            return self.drive.exchange_rate*(self.probability if before == 0 else 1-self.probability)
        if kind == "mixed_L_swap":
            if self.mode != "mixed_L" or 3 not in (state["material"][i], state["material"][j]):
                return 0.
            return super().rate(state, (3, i, j, -1, -1))*4/(self.n-1)
        if kind == 3 and self.mode == "mixed_L" and 3 in (state["material"][i], state["material"][j]):
            return 0.
        return super().rate(state, event)

    def work(self, event):
        return self.mu*(event[4]-event[3]) if event[0] == "reservoir_exchange" else 0.

    def apply(self, state, index):
        event = self.events[index]
        kind, i, j, before, after = event
        if index < self.base_count:
            record = super().apply(state, index)
            record.update(chemical_work=0., subsystem="internal")
            return record
        assert self.rate(state, event) > 0
        old_fuel = sum(state["fuel"])
        delta, work = self.delta(state, event), self.work(event)
        if kind == "reservoir_exchange":
            state["fuel"][i] = after
        else:
            state["material"][i], state["material"][j] = state["material"][j], state["material"][i]
            state["instance_tags"][i], state["instance_tags"][j] = state["instance_tags"][j], state["instance_tags"][i]
        heat = work-delta
        state["cumulative_heat"] += heat
        state["external_work"] += work
        state["event_count"] += 1
        new_fuel = sum(state["fuel"])
        fuel_delta = self.p.E_F*(new_fuel-old_fuel)
        return dict(time=state["time"], event_index=index, kind=kind, i=i, j=j,
                    before=before, after=after, delta_H=delta, heat_to_bath=heat,
                    external_work=work, chemical_work=work, fuel_before=old_fuel,
                    fuel_after=new_fuel, fuel_standard_energy_change=fuel_delta,
                    material_energy_change=delta-fuel_delta, instance_created=0,
                    instance_removed=0, parent_instance=0, catalyst_instance=0,
                    subsystem="reservoir" if kind == "reservoir_exchange" else "transport")

    def summary(self, state, include_state=False):
        out = super().summary(state)
        out["chemical_work"] = state["external_work"]
        if include_state:
            out["material"] = state["material"][:]
            out["instance_tags"] = state["instance_tags"][:]
        return out


def verify_trajectory(run):
    ref = DrivenReference(Params(**run["params"]), Flags(**run["flags"]),
                          DriveParams(**run["drive"]), run["transport_mode"])
    close(ref.mu, run["effective_mu"])
    close(ref.probability, run["effective_reservoir_fuel_probability"])
    state = copy.deepcopy(run["initial_state"])
    ref.validate(state)
    initial_h = ref.energy(state)
    initial_q, initial_w = state["cumulative_heat"], state["external_work"]
    initial_count, initial_fuel = state["event_count"], sum(state["fuel"])
    direct_h, sample_index = initial_h, 0
    samples, errors = run["samples"], {}
    counts, reservoir_balance, internal_fuel = {}, 0, 0

    def samples_until(until):
        nonlocal sample_index
        while sample_index < len(samples) and samples[sample_index]["time"] <= until:
            observed = samples[sample_index]
            expected = ref.summary(state, True)
            expected["time"] = observed["time"]
            compare_record(expected, observed, errors)
            sample_index += 1

    for n, saved in enumerate(run["ledger"]):
        assert state["time"] <= saved["time"] <= run["horizon"]
        samples_until(saved["time"])
        state["time"] = saved["time"]
        index = int(saved["event_index"])
        assert 0 <= index < len(ref.events)
        expected = ref.apply(state, index)
        compare_record(expected, saved, errors)
        counts[saved["kind"]] = counts.get(saved["kind"], 0)+1
        if saved["kind"] == "reservoir_exchange":
            reservoir_balance += saved["after"]-saved["before"]
        elif saved["kind"] in ("copy_R", "produce_L"):
            internal_fuel += -1 if saved["before"] == 0 else 1
        assert sum(state["fuel"]) == initial_fuel+reservoir_balance+internal_fuel
        close(state["external_work"]-initial_w, ref.mu*reservoir_balance)
        direct_h += expected["delta_H"]
        close(initial_h+state["external_work"]-initial_w-direct_h-state["cumulative_heat"]+initial_q, 0.)
        if n % 97 == 0 or n == len(run["ledger"])-1:
            close(ref.energy(state), direct_h)
            ref.validate(state)
    samples_until(run["horizon"])
    assert sample_index == len(samples)
    state["time"] = run["horizon"]
    compare_state(state, run["final_state"])
    compare_record(ref.summary(state), run["final_summary"], errors)
    assert state["event_count"]-initial_count == len(run["ledger"])
    close(ref.energy(state), direct_h)
    return dict(seed=run["seed"], flags=run["flags"], transport_mode=run["transport_mode"],
                events_replayed=len(run["ledger"]), event_counts=counts,
                samples_checked=len(samples), net_reservoir_fuel_exchange=reservoir_balance,
                internal_fuel_change=internal_fuel, initial_H=initial_h, final_H=ref.energy(state),
                chemical_work=state["external_work"]-initial_w,
                heat_to_bath=state["cumulative_heat"]-initial_q, maximum_record_errors=errors)


def check_physics(params, drive, configurations=2, seed=5471901):
    rng = random.Random(seed)
    p, supplied_drive = Params(**params), DriveParams(**drive)
    drives = (supplied_drive, replace(supplied_drive, mu=0.))
    maxima = dict(H=0., delta_H=0., rate=0., local_detailed_balance_log=0., zero_mu_full_H_DB_log=0.)
    comparisons = pairs = nonadjacent = cases = 0
    flags = [Flags(produce_L=a, protect_R=b) for a, b in itertools.product((False, True), repeat=2)]
    initial = [DrivenModel(p, f, supplied_drive, mode).state_to_dict(DrivenModel(p, f, supplied_drive, mode).initial(137))
               for mode in ("local", "mixed_L") for f in flags]
    assert all(state == initial[0] for state in initial)
    for _ in range(configurations):
        state = copy.deepcopy(initial[0])
        state["material"] = [rng.randrange(4) for _ in state["material"]]
        state["fuel"] = [rng.randrange(2) for _ in state["fuel"]]
        state["instance_tags"] = [j+1 if a else 0 for j, a in enumerate(state["material"])]
        state["next_tag"] = len(state["material"])+1
        for dr, mode, f in itertools.product(drives, ("local", "mixed_L"), flags):
            ref, model = DrivenReference(p, f, dr, mode), DrivenModel(p, f, dr, mode)
            converted = model.state_from_dict(state)
            events, rates = model.rates(converted)
            assert [tuple(asdict(event).values()) for event in events] == ref.events
            assert list(model.reverse_index) == ref.reverse_index
            close(model.energy(converted), ref.energy(state))
            maxima["H"] = max(maxima["H"], abs(model.energy(converted)-ref.energy(state)))
            enabled = []
            for index, event in enumerate(ref.events):
                expected = ref.rate(state, event)
                close(expected, rates[index], 2e-12)
                maxima["rate"] = max(maxima["rate"], abs(expected-float(rates[index])))
                comparisons += 1
                if expected > 0:
                    enabled.append(index)
            # Bias the challenge toward the new nonadjacent swap implementation.
            rng.shuffle(enabled)
            selected = enabled[:60]
            extra = [k for k in enabled if ref.events[k][0] == "mixed_L_swap"
                     and ref.events[k][2] not in ref.neighbors[ref.events[k][1]]]
            selected += extra[:40]
            for index in set(selected):
                event = ref.events[index]
                before_h = ref.energy(state)
                delta = ref.delta(state, event)
                actual_delta = model.delta_energy(converted, index)
                close(delta, actual_delta)
                maxima["delta_H"] = max(maxima["delta_H"], abs(delta-actual_delta))
                after = copy.deepcopy(state)
                ref.apply(after, index)
                close(ref.energy(after)-before_h, delta)
                reverse = ref.events[ref.reverse_index[index]]
                ratio = math.log(ref.rate(state, event)/ref.rate(after, reverse))
                discrepancy = abs(ratio-(ref.work(event)-delta)/p.theta)
                assert discrepancy < 2e-11
                maxima["local_detailed_balance_log"] = max(maxima["local_detailed_balance_log"], discrepancy)
                if dr.mu == 0:
                    full_h_error = abs(ratio+delta/p.theta)
                    assert full_h_error < 2e-11
                    maxima["zero_mu_full_H_DB_log"] = max(maxima["zero_mu_full_H_DB_log"], full_h_error)
                nonadjacent += int(event[0] == "mixed_L_swap" and event[2] not in ref.neighbors[event[1]])
                pairs += 1
            metadata = copy.deepcopy(state)
            metadata["time"], metadata["event_count"] = 1e6, 4567
            metadata["external_work"], metadata["cumulative_heat"] = 9999., -123.
            metadata["instance_tags"] = [t+9000 if t else 0 for t in metadata["instance_tags"]]
            metadata["next_tag"] += 18000
            _, unchanged_rates = model.rates(model.state_from_dict(metadata))
            assert all(a == b for a, b in zip(rates, unchanged_rates))
            close(model.energy(model.state_from_dict(metadata)), ref.energy(state))
            cases += 1
    assert nonadjacent > 0
    return dict(copy_rate=p.copy_rate, parameter_values=params, state_drive_mode_flag_cases=cases,
                individual_propensities_compared=comparisons, reverse_pairs=pairs,
                nonadjacent_swap_reverse_pairs=nonadjacent, maximum_errors=maxima,
                all_initial_physical_states_identical=True, observer_metadata_not_in_rates_or_H=True,
                zero_mu_recovers_full_H_detailed_balance=True,
                mixed_per_L_unaccepted_proposal_budget=4*p.material_mobility,
                mixing_changes_transport_graph_and_accepted_flux=True)


def spatial(material, ref):
    nr, nl = sum(a in (1, 2) for a in material), material.count(3)
    rl = sum((material[i] in (1, 2) and material[j] == 3)
             or (material[j] in (1, 2) and material[i] == 3) for i, j in ref.edges)
    unseen = {i for i, a in enumerate(material) if a}
    sizes, mixed = [], 0
    while unseen:
        start = unseen.pop()
        queue, members = [start], {start}
        while queue:
            for j in ref.neighbors[queue.pop()]:
                if j in unseen:
                    unseen.remove(j)
                    members.add(j)
                    queue.append(j)
        sizes.append(len(members))
        kinds = {material[j] for j in members}
        mixed += int(3 in kinds and bool(kinds & {1, 2}))
    null = 4*nr*nl/(ref.n-1)
    return dict(R0=material.count(1), R1=material.count(2), R=nr, L=nl, D=material.count(0),
                RL_edges=rl, LL_edges=sum(material[i] == material[j] == 3 for i, j in ref.edges),
                RR_edges=sum(material[i] in (1, 2) and material[j] in (1, 2) for i, j in ref.edges),
                R_with_L_neighbor=sum(a in (1, 2) and any(material[j] == 3 for j in ref.neighbors[i])
                                      for i, a in enumerate(material)),
                L_with_R_neighbor=sum(a == 3 and any(material[j] in (1, 2) for j in ref.neighbors[i])
                                      for i, a in enumerate(material)),
                RL_enrichment_over_random_placement=rl/null if null else None,
                contact_components=len(sizes), largest_contact_component=max(sizes, default=0),
                mixed_RL_components=mixed)


def trapezoid(values, times):
    return sum((a+b)*.5*(v-u) for a, b, u, v in zip(values[:-1], values[1:], times[:-1], times[1:]))


def readouts(run, original_founders=None):
    initial, final = run["initial_state"], run["final_state"]
    material = initial["material"]
    nr, nl = sum(a in (1, 2) for a in material), material.count(3)
    integral_r = integral_l = previous = work = heat = 0.
    first_zero = 0. if nr == 0 else None
    counts = dict(paid_R_births=0, reverse_R_copies=0, paid_L_births=0, reverse_L_production=0,
                  passive_R_births=0, passive_R_deaths=0, passive_L_births=0, passive_L_deaths=0,
                  reservoir_fuel_up=0, reservoir_fuel_down=0, mixed_L_swaps=0)
    paid_r_tags, paid_l_tags, removal_times = set(), set(), {}
    first_l_birth = None
    for event in run["ledger"]:
        t = event["time"]
        integral_r += nr*(t-previous)
        integral_l += nl*(t-previous)
        previous = t
        work += event["external_work"]
        heat += event["heat_to_bath"]
        kind, before, after = event["kind"], event["before"], event["after"]
        if event["instance_removed"]:
            removal_times[event["instance_removed"]] = t
        if kind == "reservoir_exchange":
            counts["reservoir_fuel_up" if after == 1 else "reservoir_fuel_down"] += 1
        elif kind == "mixed_L_swap":
            counts["mixed_L_swaps"] += 1
        elif kind in ("passive", "copy_R", "produce_L"):
            nr += int(after in (1, 2))-int(before in (1, 2))
            nl += int(after == 3)-int(before == 3)
            if after == 3 and first_l_birth is None:
                first_l_birth = dict(time=t, channel=kind)
            if kind == "copy_R":
                counts["paid_R_births" if after in (1, 2) else "reverse_R_copies"] += 1
                if after in (1, 2):
                    paid_r_tags.add(event["instance_created"])
            elif kind == "produce_L":
                counts["paid_L_births" if after == 3 else "reverse_L_production"] += 1
                if after == 3:
                    paid_l_tags.add(event["instance_created"])
            else:
                label = "passive_R_" if before in (1, 2) or after in (1, 2) else "passive_L_"
                counts[label+("births" if before == 0 else "deaths")] += 1
        if first_zero is None and nr == 0:
            first_zero = t
    integral_r += nr*(run["horizon"]-previous)
    integral_l += nl*(run["horizon"]-previous)
    assert nr == sum(a in (1, 2) for a in final["material"])
    assert nl == final["material"].count(3)
    ref = DrivenReference(Params(**run["params"]), Flags(**run["flags"]),
                          DriveParams(**run["drive"]), run["transport_mode"])
    initial_r = {tag for a, tag in zip(material, initial["instance_tags"]) if a in (1, 2)}
    initial_l = {tag for a, tag in zip(material, initial["instance_tags"]) if a == 3}
    founders = initial_r if original_founders is None else set(original_founders)
    live_tags = set(initial["instance_tags"])-{0}
    births_after = {}
    for label, group, channel in (("original_founder", founders, "copy_R"),
                                  ("preparation_R", initial_r, "copy_R"),
                                  ("preparation_L", initial_l, "produce_L")):
        group = group & live_tags
        clearance = (max((removal_times[tag] for tag in group), default=0.)
                     if all(tag in removal_times for tag in group) else None)
        births_after[("paid_R_births_after_" if channel == "copy_R" else "paid_L_births_after_")+label+"_clearance"] = sum(
            event["kind"] == channel and event["before"] == 0 and clearance is not None and event["time"] >= clearance
            for event in run["ledger"])
    sample_spatial = [spatial(sample["material"], ref) for sample in run["samples"]]
    times = [sample["time"] for sample in run["samples"]]
    contact = [s["RL_edges"] for s in sample_spatial]
    excess = [s["RL_edges"]-4*s["R"]*s["L"]/(ref.n-1) for s in sample_spatial]
    enrichments = [s["RL_enrichment_over_random_placement"] for s in sample_spatial
                   if s["RL_enrichment_over_random_placement"] is not None]
    powers = []
    side = ref.p.side
    # Direct Fourier sums, independent of the runner's FFT implementation.
    waves = [[complex(math.cos(2*math.pi*(kx*(j//side)+ky*(j % side))/side),
                      math.sin(2*math.pi*(kx*(j//side)+ky*(j % side))/side))
              for j in range(ref.n)] for kx, ky in ((1, 0), (-1, 0), (0, 1), (0, -1))]
    for sample, geometry in zip(run["samples"], sample_spatial):
        rcount = geometry["R"]
        if 0 < rcount < ref.n:
            amplitudes = [sum(weight for a, weight in zip(sample["material"], wave) if a in (1, 2)) for wave in waves]
            powers.append(statistics.mean(abs(z)**2 for z in amplitudes)/(rcount*(ref.n-rcount)/(ref.n-1)))
    final_r_tags = [tag for a, tag in zip(final["material"], final["instance_tags"]) if a in (1, 2)]
    final_l_tags = [tag for a, tag in zip(final["material"], final["instance_tags"]) if a == 3]
    out = dict(**counts, integral_R=integral_r, integral_L=integral_l, final_R=nr, final_L=nl,
               initial_R=len(initial_r), initial_L=len(initial_l), initial_fuel=sum(initial["fuel"]),
               final_fuel=sum(final["fuel"]), initial_H=ref.energy(initial), final_H=ref.energy(final),
               chemical_work=work, heat_to_bath=heat, events=len(run["ledger"]),
               first_R_inactivity=first_zero, first_episode_duration=run["horizon"] if first_zero is None else first_zero,
               final_paid_R_instances=sum(t in paid_r_tags for t in final_r_tags),
               final_paid_L_instances=sum(t in paid_l_tags for t in final_l_tags),
               founder_R_instances_remaining=sum(t in founders for t in final_r_tags),
               sampled_integral_RL_contacts=trapezoid(contact, times),
               sampled_integral_RL_excess_over_fixed_counts=trapezoid(excess, times),
               sampled_mean_RL_enrichment=statistics.mean(enrichments) if enrichments else None,
               sampled_mean_low_wave_Rpower_normalized=statistics.mean(powers) if powers else None,
               final_spatial=spatial(final["material"], ref), **births_after)
    close(out["final_H"]-out["initial_H"]-work+heat, 0.)
    return out, first_l_birth


def compare_values(expected, observed):
    if isinstance(expected, dict):
        for key, value in expected.items():
            compare_values(value, observed[key])
    elif isinstance(expected, (list, tuple)):
        assert len(expected) == len(observed)
        for value, other in zip(expected, observed):
            compare_values(value, other)
    elif expected is None:
        assert observed is None
    elif isinstance(expected, (float, int)):
        close(expected, observed)
    else:
        assert expected == observed


def estimate(values):
    n = len(values)
    mean = statistics.mean(values)
    se = statistics.stdev(values)/math.sqrt(n) if n > 1 else None
    return dict(n=n, mean=mean, monte_carlo_standard_error=se,
                approximate_95_percent_normal_interval=[mean-1.96*se, mean+1.96*se] if se is not None else None)


WEIGHTS = {
    "production_effect_when_protected": {("local", True, True): 1, ("local", False, True): -1},
    "production_protection_interaction": {("local", True, True): 1, ("local", True, False): -1,
                                          ("local", False, True): -1, ("local", False, False): 1},
    "locality_effect_full": {("local", True, True): 1, ("mixed_L", True, True): -1},
}
FULL = dict(produce_L=True, protect_R=True, copy_R=True)


def verify_study(path):
    result = load(path)
    protocol_path = ROOT/result["protocol_file"]
    protocol = load(protocol_path)
    assert digest(protocol_path) == result["protocol_sha256"]
    if result["stage"] == "pilot":
        assert digest(protocol_path) == PILOT_PROTOCOL_SHA256
    elif result["stage"] == "confirmation":
        assert digest(protocol_path) == CONFIRMATION_PROTOCOL_SHA256
    else:
        raise AssertionError("undeclared stage")
    assert protocol["source_hashes"] == result["source_hashes"]
    for name, expected in protocol["source_hashes"].items():
        assert digest(ROOT/name) == expected, ("changed source", name)
    for field in ("stage", "copy_rate_levels", "seeds", "preparation_horizon"):
        assert protocol[field] == result[field]
    assert result["horizon"] == protocol["continuation_horizon"] == 80.
    assert protocol["preparation_horizon"] == 20.
    assert protocol["transport_modes"] == ["local", "mixed_L"]
    assert set(protocol["primary_contrasts"]) == set(WEIGHTS)
    assert protocol["primary_field"] == "integral_R"
    rates, seeds = protocol["copy_rate_levels"], protocol["seeds"]
    cells = {(mode, a, b) for mode, a, b in itertools.product(protocol["transport_modes"], (False, True), (False, True))}
    assert {(a["produce_L"], a["protect_R"]) for a in protocol["arms"]} == {(False, False), (False, True), (True, False), (True, True)}
    assert all(a["copy_R"] for a in protocol["arms"])
    preparations = {(r["copy_rate"], r["origin_seed"]): r for r in result["preparations"]}
    expected_preps = set(itertools.product(rates, seeds))
    assert len(preparations) == len(result["preparations"]) and set(preparations) == expected_preps
    branches = {(r["copy_rate"], r["origin_seed"], r["transport_mode"], r["flags"]["produce_L"], r["flags"]["protect_R"]): r
                for r in result["runs"]}
    assert len(branches) == len(result["runs"])
    assert set(branches) == {(rate, seed, *cell) for rate, seed in expected_preps for cell in cells}
    replays, values, initial_by_seed, formations = [], {}, {}, []
    no_r = no_l = 0
    for rate, seed in sorted(expected_preps):
        prep_row = preparations[(rate, seed)]
        prep_path = ROOT/prep_row["trajectory_file"]
        assert digest(prep_path) == prep_row["trajectory_sha256"]
        prep = load(prep_path)
        params = {**protocol["base_params"], "copy_rate": rate}
        assert prep["params"] == params and prep["drive"] == protocol["drive"]
        assert prep["flags"] == FULL and prep["transport_mode"] == "local"
        assert prep["seed"] == prep_row["seed"] == seed and prep["horizon"] == 20.
        assert prep["copy_rate"] == rate and prep["phase"] == "preparation"
        assert prep["stage"] == result["stage"] and prep["source_hashes"] == protocol["source_hashes"]
        assert prep["protocol_sha256"] == digest(protocol_path)
        initial = prep["initial_state"]
        assert initial["material"].count(1) == 8 and initial["material"].count(2) == initial["material"].count(3) == 0
        assert sum(initial["fuel"]) == params["initial_fuel"] == 58
        if seed in initial_by_seed:
            assert initial == initial_by_seed[seed]
        initial_by_seed[seed] = initial
        replay = verify_trajectory(prep)
        replay.update(record=str(prep_path), phase="preparation", copy_rate=rate, origin_seed=seed)
        replays.append(replay)
        measured_prep, first_l = readouts(prep)
        compare_values(measured_prep, prep_row)
        formations.append(dict(copy_rate=rate, origin_seed=seed, initial_L=0, first_L_birth=first_l,
                               final_R=measured_prep["final_R"], final_L=measured_prep["final_L"],
                               initial_H=measured_prep["initial_H"], chemical_work=measured_prep["chemical_work"]))
        no_r += int(measured_prep["final_R"] == 0)
        no_l += int(measured_prep["final_L"] == 0)
        reset = copy.deepcopy(prep["final_state"])
        for key in ("time", "cumulative_heat", "external_work", "event_count"):
            reset[key] = 0
        founders = {t for a, t in zip(initial["material"], initial["instance_tags"]) if a in (1, 2)}
        for cell in sorted(cells):
            row = branches[(rate, seed, *cell)]
            run_path = ROOT/row["trajectory_file"]
            assert digest(run_path) == row["trajectory_sha256"]
            run = load(run_path)
            mode, a, b = cell
            assert run["params"] == params and run["drive"] == protocol["drive"]
            assert run["flags"] == row["flags"] == dict(produce_L=a, protect_R=b, copy_R=True)
            assert run["transport_mode"] == mode
            assert run["initial_state"] == reset
            assert run["seed"] == row["seed"] == seed+protocol["post_seed_offset"]
            assert run["origin_seed"] == seed and run["copy_rate"] == rate
            assert run["horizon"] == 80. and run["phase"] == "continuation"
            assert run["stage"] == result["stage"] and run["source_hashes"] == protocol["source_hashes"]
            assert run["protocol_sha256"] == digest(protocol_path)
            assert run["preparation_file"] == row["preparation_file"] == prep_row["trajectory_file"]
            assert run["preparation_sha256"] == row["preparation_sha256"] == digest(prep_path)
            replay = verify_trajectory(run)
            replay.update(record=str(run_path), phase="continuation", copy_rate=rate, origin_seed=seed)
            replays.append(replay)
            measured, _ = readouts(run, founders)
            compare_values(measured, row)
            assert row["prepared_R"] == measured_prep["final_R"]
            assert row["prepared_L"] == measured_prep["final_L"]
            assert row["prepared_fuel"] == measured_prep["final_fuel"]
            close(row["supplied_initial_H"], measured_prep["initial_H"])
            close(row["preparation_chemical_work"], measured_prep["chemical_work"])
            close(row["preparation_heat_to_bath"], measured_prep["heat_to_bath"])
            for total, field in (("total_chemical_work", "chemical_work"), ("total_heat_to_bath", "heat_to_bath")):
                measured[total] = measured_prep[field]+measured[field]
                close(measured[total], row[total])
            close(measured_prep["initial_H"]+measured["total_chemical_work"]-measured["final_H"]-measured["total_heat_to_bath"], 0.)
            values[(rate, seed, *cell)] = measured
    for group in result["groups"]:
        cell = group["transport_mode"], group["flags"]["produce_L"], group["flags"]["protect_R"]
        for field, observed in group["statistics"].items():
            compare_values(estimate([values[(group["copy_rate"], seed, *cell)][field] for seed in seeds]), observed)
        for field, observed in group["morphology_statistics_conditioned_on_defined_samples"].items():
            defined = [values[(group["copy_rate"], seed, *cell)][field] for seed in seeds
                       if values[(group["copy_rate"], seed, *cell)][field] is not None]
            compare_values(estimate(defined) if defined else dict(n=0, mean=None), observed)
    contrasts = []
    for saved in result["paired_contrasts"]:
        expected = estimate([sum(weight*values[(saved["copy_rate"], seed, *cell)][saved["field"]]
                                 for cell, weight in WEIGHTS[saved["contrast"]].items()) for seed in seeds])
        compare_values(expected, saved)
        if saved["field"] == "integral_R":
            contrasts.append(dict(copy_rate=saved["copy_rate"], contrast=saved["contrast"], **expected))
    joint = None
    if result["stage"] == "confirmation":
        assert len(contrasts) == 3
        joint = all(c["approximate_95_percent_normal_interval"][0] > 0 for c in contrasts)
        assert result["declared_joint_claim_supported"] == joint
        for expected in contrasts:
            compare_values(expected, next(c for c in result["primary_confirmation_results"] if c["contrast"] == expected["contrast"]))
    return dict(checks_passed=True, stage=result["stage"], result_file=str(path), result_sha256=digest(path),
                protocol_sha256=digest(protocol_path), preparation_trajectories=len(preparations),
                continuation_trajectories=len(branches), total_events_replayed=sum(r["events_replayed"] for r in replays),
                samples_checked=sum(r["samples_checked"] for r in replays), all_seed_rate_mode_flag_cells_present=True,
                all_eight_branch_physical_states_identical=True, no_R_preparations_retained=no_r,
                no_L_preparations_retained=no_l, supplied_initial_H_and_all_preparation_work_counted=True,
                zero_L_formation_histories=formations, all_primary_and_secondary_statistics_recomputed=True,
                primary_R_time_contrasts=contrasts, declared_joint_claim_supported=joint, replays=replays)


def check_selection():
    pilot_protocol_path = ROOT/"pilot-protocol.json"
    assert digest(pilot_protocol_path) == PILOT_PROTOCOL_SHA256
    pilot_protocol, pilot = load(pilot_protocol_path), load(ROOT/"pilot-results.json")
    protocol_path = ROOT/"confirmation-protocol.json"
    protocol = load(protocol_path)
    assert digest(ROOT/"pilot-results.json") == protocol["pilot_results_sha256"]
    assert digest(pilot_protocol_path) == protocol["pilot_protocol_sha256"]
    for field in ("base_params", "drive", "source_hashes", "preparation_horizon", "continuation_horizon",
                  "arms", "transport_modes", "post_seed_offset", "selection_rule", "primary_contrasts", "confirmation_rule"):
        assert protocol[field] == pilot_protocol[field]
    assert protocol["seeds"] == pilot_protocol["confirmation_seeds"]
    assert not set(protocol["seeds"]) & set(pilot_protocol["seeds"])
    scores = []
    for rate in pilot_protocol["copy_rate_levels"]:
        index = {(r["origin_seed"], r["transport_mode"], r["flags"]["produce_L"], r["flags"]["protect_R"]): r["integral_R"]
                 for r in pilot["runs"] if r["copy_rate"] == rate}
        stats = [estimate([sum(weight*index[(seed, *cell)] for cell, weight in WEIGHTS[name].items())
                           for seed in pilot_protocol["seeds"]]) for name in pilot_protocol["primary_contrasts"]]
        zs = [s["mean"]/s["monte_carlo_standard_error"] if s["monte_carlo_standard_error"] else 0. for s in stats]
        scores.append(dict(copy_rate=rate, selection_score=min(zs), primary_z_scores=zs))
    selected = max(scores, key=lambda s: (s["selection_score"], -s["copy_rate"]))["copy_rate"]
    assert protocol["copy_rate_levels"] == [selected] and protocol["selected_copy_rate"] == selected
    for expected, observed in zip(scores, protocol["selection_scores"]):
        compare_values(expected, observed)
    return dict(checks_passed=True, selected_copy_rate=selected, scores_recomputed=scores,
                fresh_confirmation_seeds=True, confirmation_protocol_sha256=digest(protocol_path))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", nargs="+", type=Path)
    args = parser.parse_args()
    source_hash = digest(Path(__file__))
    assert digest(BASE/"verify_interface.py") == REFERENCE_SHA256
    assert digest(BASE/"interface_model.py") == BASE_MODEL_SHA256
    assert digest(ROOT/"driven_model.py") == DRIVEN_MODEL_SHA256
    studies = [verify_study(path.resolve()) for path in args.results]
    records = [r["record"] for s in studies for r in s["replays"]]
    assert len(records) == len(set(records))
    selection = check_selection() if any(s["stage"] == "confirmation" for s in studies) else None
    actual = {}
    for path in args.results:
        result = load(path)
        protocol = load(ROOT/result["protocol_file"])
        for rate in protocol["copy_rate_levels"]:
            actual[rate] = ({**protocol["base_params"], "copy_rate": rate}, protocol["drive"])
    model_checks = [check_physics(params, drive) for _, (params, drive) in sorted(actual.items())]
    stress = asdict(Params(side=4, initial_templates=4, J_RR=.7, J_LL=1.2, J_RL=.37))
    model_checks.append(check_physics(stress, asdict(DriveParams()), configurations=1))
    assert digest(Path(__file__)) == source_hash and digest(ROOT/"driven_model.py") == DRIVEN_MODEL_SHA256
    assert digest(BASE/"interface_model.py") == BASE_MODEL_SHA256 and digest(BASE/"verify_interface.py") == REFERENCE_SHA256
    print(json.dumps(dict(checks_passed=True, verifier_sha256=source_hash, model_sha256=DRIVEN_MODEL_SHA256,
                          independent_reference_sha256=REFERENCE_SHA256, base_model_sha256=BASE_MODEL_SHA256,
                          trajectory_records_replayed=len(records),
                          preparation_trajectories=sum(s["preparation_trajectories"] for s in studies),
                          continuation_trajectories=sum(s["continuation_trajectories"] for s in studies),
                          total_events_replayed=sum(s["total_events_replayed"] for s in studies),
                          model_checks=model_checks, selection_audit=selection, studies=studies,
                          scope="Matched reservoir potential and exchange capacity, supplied initial founders/fuel/H and transport fixtures. Actual chemical work can differ by arm. Full preparation+continuation first law and conservative transport verified. Zero-L local formation and postformation transport comparisons are distinct. Bath heat is not total entropy production; no boundary, individuality, biodiversity or closed cosmic theorem is established."), indent=2))


if __name__ == "__main__":
    main()

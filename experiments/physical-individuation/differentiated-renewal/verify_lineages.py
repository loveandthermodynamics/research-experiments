"""Adversarial unit fixtures for offline ancestry readouts, not physical runs.

The synthetic ledgers carry zero placeholder energies. They test identity,
time integrals and interpretation only; physical event feasibility and energy
accounting belong to the separate stochastic-model verifier.
"""
from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
from lineage_readouts import measure
import verify_study as independent

ROOT = Path(__file__).resolve().parent


def birth(time, tag, species, parent=0, helper=0, kind="copy_R"):
    return dict(time=time, kind=kind, before=0, after=species,
                instance_created=tag, instance_removed=0, parent_instance=parent,
                helper_species=helper, external_work=0., heat_to_bath=0.)


def death(time, tag, species):
    return dict(time=time, kind="passive", before=species, after=0,
                instance_created=0, instance_removed=tag, parent_instance=0,
                external_work=0., heat_to_bath=0.)


def gate(time, before, after):
    return dict(time=time, kind="gate_flip", before=before, after=after,
                gate_before=before, gate_after=after, instance_created=0,
                instance_removed=0, parent_instance=0, external_work=0., heat_to_bath=0.)


def fixture(events, horizon=10.):
    living = {1: 1}
    initial = dict(material=[1]+[0]*15, instance_tags=[1]+[0]*15,
                   external_work=0., cumulative_heat=0.)
    g, gate_count = -1, 0
    previous = 0.
    for event in events:
        assert previous <= event["time"] <= horizon
        previous = event["time"]
        if event["instance_removed"]:
            assert living.pop(event["instance_removed"]) == event["before"]
        if event["instance_created"]:
            if event["kind"] == "copy_R":
                assert living[event["parent_instance"]] == event["helper_species"]
            assert event["instance_created"] not in living
            living[event["instance_created"]] = event["after"]
        if event["kind"] == "gate_flip":
            assert event["before"] == g
            g, gate_count = event["after"], gate_count+1
    tags = sorted(living)
    final = dict(material=[living[t] for t in tags]+[0]*(16-len(tags)),
                 instance_tags=tags+[0]*(16-len(tags)), external_work=0., cumulative_heat=0.)
    return dict(initial_state=dict(base_state=initial, gate=-1, gate_event_count=0),
                final_state=dict(base_state=final, gate=g, gate_event_count=gate_count),
                samples=[dict(H=0.)], final_summary=dict(H=0.), ledger=events, horizon=horizon)


def check(expected, observed, path=""):
    if isinstance(expected, dict):
        for key, value in expected.items():
            check(value, observed[key], path+"."+key)
    elif isinstance(expected, float):
        assert abs(expected-observed) < 1e-10, (path, expected, observed)
    else:
        assert expected == observed, (path, expected, observed)


def main():
    source = ROOT/"lineage_readouts.py"
    source_hash = sha256(source.read_bytes()).hexdigest()
    tests = []

    class SyntheticZeroEnergy:
        def energy(self, state):
            return 0.

    # Synthetic fixtures intentionally do not carry a Hamiltonian/valid rates.
    # Their graph identities/integrals are cross-checked without pretending to
    # replay physical energies; the separate model verifier handles that.
    independent.make_reference = lambda run: SyntheticZeroEnergy()

    def test(name, events, expected, horizon=10.):
        run = fixture(events, horizon)
        observed = measure(run)
        check(expected, observed, name)
        check(independent.graph_readouts(run), observed, name+".independent_lifetime_graph")
        assert observed["establishment_population"] == "copying descendants over both current alleles"
        assert abs(observed["integral_R"]-observed["integral_R0"]-observed["integral_R1"]) < 1e-10
        assert observed["integral_focal_R1"] <= observed["integral_focal_descendants"]
        tests.append(dict(name=name, checks_passed=True, expected_fields=expected))

    test("late_thermal_immigration_is_not_a_focal_mutant",
         [death(5, 1, 1), birth(25, 2, 2, kind="passive")],
         dict(focal_tag=None, focal_paid_mutant_appeared_by_deadline=False,
              first_R1_appearance_time=25., first_R1_appearance_origin="passive",
              integral_R=80., integral_R1=75., integral_founder_descendant_R1=0.,
              final_R1=1, first_episode_duration=5., focal_established=False), horizon=100.)
    test("late_paid_mutant_does_not_retroactively_meet_deadline",
         [birth(21, 2, 2, 1, 1)],
         dict(focal_tag=None, focal_paid_mutant_appeared_by_deadline=False,
              first_R1_appearance_time=21., first_R1_appearance_origin="copy_R",
              integral_R=179., integral_founder_descendant_R1=79.,
              birth_counts=dict(outside_focal_R0_to_R1=1)), horizon=100.)
    test("first_any_origin_can_precede_first_paid_mutant",
         [birth(1, 2, 2, kind="passive"), birth(2, 3, 2, 1, 1)],
         dict(first_R1_appearance_time=1., first_R1_appearance_origin="passive",
              focal_tag=3, focal_appearance_time=2., integral_R1=17.,
              integral_focal_R1=8., integral_founder_descendant_R1=8.))
    test("backmutation_preserves_subtree_but_R1_presence_can_break",
         [birth(1, 2, 2, 1, 1), birth(2, 3, 1, 2, 2), death(3, 2, 2),
          birth(4, 4, 2, kind="passive"), birth(5, 5, 2, 1, 1),
          birth(6, 6, 2, 3, 1), death(7, 6, 2), birth(8, 7, 2, 3, 1)],
         dict(focal_tag=2, focal_extinction_time=None, focal_survives_to_horizon=True,
              focal_R1_first_loss_time=3., focal_R1_restorations_after_loss=2,
              focal_R1_uninterrupted_to_horizon=False, final_focal_descendants=2,
              final_focal_R1=1, integral_focal_descendants=13., integral_focal_R1=5.,
              integral_founder_descendant_R1=10., both_types_present_time=8.,
              integral_R0=18., integral_R1=16., integral_R=34.,
              birth_counts=dict(outside_focal_R0_to_R1=1)))
    test("subtree_extinction_is_not_reversed_by_new_immigration",
         [birth(1, 2, 2, 1, 1), death(2, 2, 2),
          birth(3, 3, 2, kind="passive"), birth(4, 4, 2, 1, 1)],
         dict(focal_tag=2, focal_extinction_time=2., focal_survives_to_horizon=False,
              final_focal_descendants=0, integral_focal_descendants=1., integral_focal_R1=1.,
              focal_R1_restorations_after_loss=0, final_R1=2,
              integral_founder_descendant_R1=7.))
    test("four_descendants_can_be_established_without_four_R1",
         [birth(1, 2, 2, 1, 1), birth(2, 3, 1, 2, 2),
          birth(3, 4, 1, 3, 1), birth(4, 5, 1, 3, 1), death(5, 2, 2)],
         dict(focal_established=True, focal_establishment_time=4., maximum_focal_descendants=4,
              final_focal_descendants=3, final_focal_R1=0, final_R1=0,
              focal_R1_first_loss_time=5., focal_R1_uninterrupted_to_horizon=False,
              integral_focal_descendants=25., integral_focal_R1=4.))
    test("paid_mutant_can_descend_from_a_thermal_R0",
         [death(1, 1, 1), birth(2, 2, 1, kind="passive"), birth(3, 3, 2, 2, 1)],
         dict(focal_tag=3, focal_appearance_time=3., integral_focal_R1=7.,
              focal_parent_has_initial_founder_ancestry=False,
              integral_founder_descendant_R1=0., founder_instance_clearance_time=1.,
              paid_R_births_after_founder_instance_clearance=1, first_episode_duration=1.))
    test("first_focal_failure_is_not_replaced_by_a_later_successful_mutant",
         [birth(1, 2, 2, 1, 1), death(2, 2, 2), birth(3, 3, 2, 1, 1),
          birth(4, 4, 2, 3, 2), birth(5, 5, 2, 3, 2), birth(6, 6, 2, 3, 2)],
         dict(focal_tag=2, focal_extinction_time=2., focal_established=False,
              maximum_focal_descendants=1, final_focal_descendants=0,
              final_R1=4, focal_R1_restorations_after_loss=0))
    test("catalyzed_L_is_not_a_copying_descendant",
         [birth(1, 2, 2, 1, 1), birth(2, 3, 3, helper=2, kind="produce_L"), death(3, 2, 2)],
         dict(focal_extinction_time=3., final_focal_descendants=0, final_L=1,
              integral_focal_descendants=2., integral_focal_R1=2.))
    test("gate_residence_time_does_not_create_or_relabel_descendants",
         [gate(2, -1, 1), gate(5, 1, -1)],
         dict(gate_events=2, gate_time={"-1": 7., "1": 3.},
              integral_R=10., final_R=1, focal_tag=None))
    test("copying_after_founder_instance_turnover_retains_founder_ancestry",
         [birth(1, 2, 1, 1, 1), death(2, 1, 1), birth(3, 3, 2, 2, 1)],
         dict(founder_instance_clearance_time=2., paid_R_births_after_founder_instance_clearance=1,
              focal_parent_has_initial_founder_ancestry=True,
              focal_tag=3, integral_founder_descendant_R1=7., integral_focal_R1=7.))
    assert sha256(source.read_bytes()).hexdigest() == source_hash
    print(json.dumps(dict(checks_passed=True, source_sha256=source_hash,
                          verifier_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                          adversarial_synthetic_cases=len(tests), tests=tests,
                          independent_lifetime_graph_matches_all_synthetic_cases=True,
                          interpretation_limits=[
                              "Focal establishment means subtree size, including backmutated R0 descendants.",
                              "A paid focal mutant need not have initial-founder ancestry; its parent can be a thermal R0.",
                              "Neither subtree expansion nor persistent R1 occupancy establishes positive selection, benefit or individuality."
                          ], scope="Synthetic offline-readout fixtures only. Not physical trajectories, stochastic experiments or thermodynamic validation; do not add these to physical replay totals."), indent=2))


if __name__ == "__main__":
    main()

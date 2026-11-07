"""Independent fixed-region/cohort reconstruction and adversarial fixtures.

No rates are called and no natural trajectory is simulated. Synthetic ledgers
test observer bookkeeping, not chemical feasibility or energy conservation.
Public API: reconstruct_support(run, near_sites), verify_run(run, near_sites).
Importing this module neither imports the readout under test nor writes files.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import gzip
import json
import math
from pathlib import Path
import sys
from typing import Any, Sequence

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
R_TYPES = (1, 2)
CHEMICAL = {"passive", "copy_R", "produce_L"}
TRANSPORT = {"material_swap", "mixed_L_swap"}
OTHER = {"fuel_swap", "reservoir_exchange", "gate_flip"}


def _base(state: dict) -> dict:
    return state["base_state"] if "base_state" in state else state


def _slots(state: dict, size: int) -> list[tuple[int, int]]:
    material, tags = state["material"], state["instance_tags"]
    assert len(material) == len(tags) == size
    slots = [(int(s), int(t)) for s, t in zip(material, tags)]
    assert all(s in (0, 1, 2, 3) and (t > 0) == (s > 0) for s, t in slots)
    present = [t for _, t in slots if t]
    assert len(present) == len(set(present)), "duplicate living instance"
    return slots


def _region(side: int, core: set[int]) -> set[int]:
    result = set(core)
    for site in core:
        row, col = divmod(site, side)
        for dr, dc in ((0, 1), (0, -1), (1, 0), (-1, 0)):
            result.add(((row + dr) % side) * side + (col + dc) % side)
    return result


def reconstruct_support(run: dict[str, Any], near_sites: Sequence[int]) -> dict[str, Any]:
    """Reconstruct using per-founder roots, life endings and interval summaries."""
    side = int(run["params"]["side"])
    assert side >= 3
    size = side * side
    core = set(map(int, near_sites))
    assert all(0 <= i < size for i in core)
    region = _region(side, core)
    slots = _slots(_base(run["initial_state"]), size)
    originals = {tag for site, (species, tag) in enumerate(slots)
                 if site in core and species in R_TYPES}
    # Store an individual founder root, not merely membership in one union.
    roots: dict[int, int | None] = {tag: tag if tag in originals else None
                                    for species, tag in slots if tag}
    live = {tag: site for site, (species, tag) in enumerate(slots) if tag}
    original_deaths: dict[int, tuple[int, float]] = {}
    copies: list[tuple[int, float]] = []
    counts: Counter = Counter()
    transport: Counter = Counter()
    chemical: Counter = Counter()
    segments: list[tuple[float, float, tuple[int, ...]]] = []
    horizon = float(run["horizon"])
    assert math.isfinite(horizon) and horizon >= 0
    previous = 0.0

    def snapshot() -> tuple[int, ...]:
        local_r = global_r = local_l = cohort_r = local_cohort = cohort_r1 = 0
        for site, (species, tag) in enumerate(slots):
            if species == 3 and site in region:
                local_l += 1
            if species not in R_TYPES:
                continue
            global_r += 1
            local_r += site in region
            if roots[tag] is not None:
                cohort_r += 1
                local_cohort += site in region
                cohort_r1 += species == 2
        return local_r, global_r, local_l, cohort_r, local_cohort, cohort_r1

    initial_counts = snapshot()
    for index, event in enumerate(run["ledger"]):
        now = float(event["time"])
        assert math.isfinite(now) and previous <= now <= horizon + 1e-10
        segments.append((previous, now, snapshot()))
        previous = now
        kind = event["kind"]
        i, j = int(event.get("i", -1)), int(event.get("j", -1))
        removed = int(event.get("instance_removed", 0))
        created = int(event.get("instance_created", 0))
        parent = int(event.get("parent_instance", 0))
        if kind in TRANSPORT:
            assert 0 <= i < size and 0 <= j < size and i != j
            assert not created and not removed and not parent
            departing = [(i, j, slots[i]), (j, i, slots[j])]
            slots[i], slots[j] = slots[j], slots[i]
            for source, target, (species, tag) in departing:
                if tag:
                    live[tag] = target
                if (source in region) == (target in region):
                    continue
                direction = "out" if source in region else "in"
                transport[(species, direction)] += 1
                if species in R_TYPES and roots[tag] is not None:
                    counts[f"cohort_{direction}"] += 1
        elif kind in CHEMICAL:
            assert 0 <= i < size
            before, after = int(event["before"]), int(event["after"])
            assert before in (0, 1, 2, 3) and after in (0, 1, 2, 3)
            assert slots[i][0] == before and (before == 0) != (after == 0)
            old_tag = slots[i][1]
            assert removed == old_tag
            if after:
                assert created > 0 and created not in roots
            else:
                assert created == 0
            helper_tag = 0
            if kind in {"copy_R", "produce_L"}:
                assert 0 <= j < size and i != j and slots[j][0] in R_TYPES
                helper_tag = slots[j][1]
                if "helper_species" in event:
                    assert int(event["helper_species"]) == slots[j][0]
            if kind == "copy_R":
                assert (before or after) in R_TYPES
                assert parent == (helper_tag if after else 0)
            else:
                assert parent == 0
            if kind == "produce_L":
                assert (before or after) == 3
                orientation = "forward" if after == 3 else "reverse"
                # Deliberately ignore the record's optional catalyst field:
                # reverse records may store zero. The helper exists PRE-jump.
                if roots[helper_tag] is not None:
                    counts[f"cohort_{orientation}_L"] += 1
                if j in core:
                    counts[f"near_{orientation}_L"] += 1
            old_cohort = before in R_TYPES and roots[old_tag] is not None
            new_root = roots[parent] if kind == "copy_R" and created else None
            if removed:
                del live[removed]
                if removed in originals:
                    original_deaths[removed] = (index, now)
            if created:
                roots[created] = new_root
                live[created] = i
            if before in R_TYPES:
                chemical["R_death"] += 1
                chemical["local_R_death"] += i in region
            if after in R_TYPES:
                chemical["R_birth"] += 1
                chemical["local_R_birth"] += i in region
            if before == 3:
                chemical["local_L_death"] += i in region
            if after == 3:
                chemical["local_L_birth"] += i in region
            if old_cohort:
                chemical["cohort_death"] += 1
                chemical["local_cohort_death"] += i in region
            if created and new_root is not None:
                copies.append((index, now))
                chemical["local_cohort_birth"] += i in region
            slots[i] = (after, created)
        else:
            assert kind in OTHER, f"unrecognized nonmaterial event: {kind}"
            assert not created and not removed and not parent

    segments.append((previous, horizon, snapshot()))
    assert slots == _slots(_base(run["final_state"]), size), "final physical replay mismatch"
    assert live == {tag: site for site, (_, tag) in enumerate(slots) if tag}
    if not originals:
        cleared, clearance_index = 0.0, -1
    elif originals == set(original_deaths):
        clearance_index, cleared = max(original_deaths.values())
    else:
        cleared = clearance_index = None
    total = [math.fsum((end-start)*values[k] for start, end, values in segments)
             for k in range(6)]
    after_local = after_cohort = 0.0
    if originals and cleared is not None:
        after_local = math.fsum(max(0.0, end-max(start, cleared))*values[0]
                               for start, end, values in segments)
        after_cohort = math.fsum(max(0.0, end-max(start, cleared))*values[3]
                                for start, end, values in segments)
    final_counts = snapshot()
    cohort_sites = [i for i, (s, tag) in enumerate(slots)
                    if s in R_TYPES and roots[tag] is not None]
    final_root_counts = Counter(roots[slots[i][1]] for i in cohort_sites)
    after_copies = (sum(index > clearance_index for index, _ in copies)
                    if originals and clearance_index is not None else 0)
    results = dict(
        outcome_integral_local_R=total[0], integral_global_R=total[1], integral_local_L=total[2],
        integral_global_cohort_R=total[3], integral_local_cohort_R=total[4],
        integral_global_cohort_R1=total[5],
        integral_local_R_after_original_cohort_clearance=after_local,
        integral_cohort_R_after_original_cohort_clearance=after_cohort,
        original_cohort_clearance_time=cleared, initial_cohort_size=len(originals),
        original_cohort_present=bool(originals),
        original_cohort_cleared_after_presence=bool(originals) and cleared is not None,
        final_cohort_size=len(cohort_sites),
        final_cohort_R1=sum(slots[i][0] == 2 for i in cohort_sites),
        final_local_R=final_counts[0], final_global_R=final_counts[1],
        cohort_forward_L=counts["cohort_forward_L"], cohort_reverse_L=counts["cohort_reverse_L"],
        cohort_net_L=counts["cohort_forward_L"]-counts["cohort_reverse_L"],
        near_site_forward_L=counts["near_forward_L"], near_site_reverse_L=counts["near_reverse_L"],
        cohort_copies=len(copies), cohort_copies_after_original_clearance=after_copies,
        cohort_outward_crossings=counts["cohort_out"], cohort_inward_crossings=counts["cohort_in"],
        fixed_core_sites=sorted(core), fixed_readout_region=sorted(region))
    r_in = sum(transport[(s, "in")] for s in R_TYPES)
    r_out = sum(transport[(s, "out")] for s in R_TYPES)
    balances = dict(
        global_cohort=results["final_cohort_size"]-(len(originals)+len(copies)-chemical["cohort_death"]),
        local_cohort=final_counts[4]-(initial_counts[4]+chemical["local_cohort_birth"]
                       -chemical["local_cohort_death"]+counts["cohort_in"]-counts["cohort_out"]),
        global_R=final_counts[1]-(initial_counts[1]+chemical["R_birth"]-chemical["R_death"]),
        local_R=final_counts[0]-(initial_counts[0]+chemical["local_R_birth"]
                               -chemical["local_R_death"]+r_in-r_out),
        local_L=sum(s == 3 for i, (s, _) in enumerate(slots) if i in region)
            -(sum(s == 3 for i, s in enumerate(_base(run["initial_state"])["material"]) if i in region)
              +chemical["local_L_birth"]-chemical["local_L_death"]
              +transport[(3, "in")]-transport[(3, "out")]))
    assert all(value == 0 for value in balances.values()), balances
    return dict(measurements=results, diagnostics=dict(
        ledger_events=len(run["ledger"]), balance_errors=balances,
        per_original_final_R_descendants={str(t): final_root_counts[t] for t in sorted(originals)},
        original_instance_deaths={str(t): time for t, (_, time) in original_deaths.items()},
        cohort_R_losses=chemical["cohort_death"], local_cohort_R_births=chemical["local_cohort_birth"],
        local_cohort_R_losses=chemical["local_cohort_death"], all_R_inward_crossings=r_in,
        all_R_outward_crossings=r_out, L_inward_crossings=transport[(3, "in")],
        L_outward_crossings=transport[(3, "out")],
        near_site_net_L=counts["near_forward_L"]-counts["near_reverse_L"]))


def _compare(expected: dict, observed: dict, context: str) -> None:
    for key, want in expected.items():
        assert key in observed, (context, "missing field", key)
        got = observed[key]
        if isinstance(want, float):
            assert math.isclose(want, got, rel_tol=2e-10, abs_tol=2e-9), (context, key, want, got)
        else:
            assert want == got, (context, key, want, got)


def verify_run(run: dict[str, Any], near_sites: Sequence[int]) -> dict[str, Any]:
    independent = reconstruct_support(run, near_sites)
    from cohort_support_readouts import support_readouts
    _compare(independent["measurements"], support_readouts(run, near_sites), "readout function")
    for name in ("support_readouts", "cohort_support_readouts"):
        if name in run:
            _compare(independent["measurements"], run[name], f"recorded {name}")
    return dict(checks_passed=True, **independent)


def _event(time, kind, i=-1, j=-1, before=-1, after=-1, removed=0, created=0, parent=0):
    return dict(time=float(time), kind=kind, i=i, j=j, before=before, after=after,
                instance_removed=removed, instance_created=created, parent_instance=parent,
                catalyst_instance=parent if kind == "copy_R" and created else 0)


def _fixture(initial: dict[int, tuple[int, int]], final: dict[int, tuple[int, int]],
             events: list[dict], horizon=8.0, side=5) -> dict:
    def state(entries):
        return dict(material=[entries.get(i, (0, 0))[0] for i in range(side*side)],
                    instance_tags=[entries.get(i, (0, 0))[1] for i in range(side*side)])
    return dict(params=dict(side=side), initial_state=dict(base_state=state(initial), gate=-1),
                final_state=dict(base_state=state(final), gate=-1), ledger=events, horizon=horizon)


def fixtures() -> list[dict]:
    results = []
    def test(name, run, core, expected):
        reconstructed = verify_run(run, core)
        _compare(expected, reconstructed["measurements"], name)
        results.append(dict(name=name, explicit_expected=expected, **reconstructed))

    swap = lambda t, i, j, kind="material_swap": _event(t, kind, i, j)
    kill = lambda t, i, s, tag: _event(t, "passive", i, before=s, after=0, removed=tag)
    thermal = lambda t, i, s, tag: _event(t, "passive", i, before=0, after=s, created=tag)
    copy = lambda t, i, j, s, tag, parent: _event(t, "copy_R", i, j, 0, s, created=tag, parent=parent)
    prod = lambda t, i, j, tag: _event(t, "produce_L", i, j, 0, 3, created=tag)
    revl = lambda t, i, j, tag: _event(t, "produce_L", i, j, 3, 0, removed=tag)

    test("leaving_region_is_not_lineage_extinction",
         _fixture({0:(1,1)}, {3:(2,2)},
                  [swap(1,0,1),swap(2,1,2),copy(3,3,2,2,2,1),kill(4,2,1,1)],6), [0],
         dict(outcome_integral_local_R=2.,integral_global_R=7.,integral_global_cohort_R=7.,
              integral_local_cohort_R=2.,integral_global_cohort_R1=3.,
              original_cohort_clearance_time=4.,final_cohort_size=1,final_cohort_R1=1,
              final_local_R=0,final_global_R=1,cohort_outward_crossings=1,cohort_inward_crossings=0,
              integral_local_R_after_original_cohort_clearance=0.,
              integral_cohort_R_after_original_cohort_clearance=2.))
    test("thermal_replacement_does_not_restore_descent",
         _fixture({0:(2,1)}, {0:(2,2)},[kill(1,0,2,1),thermal(2,0,2,2)],5), [0],
         dict(outcome_integral_local_R=4.,integral_global_cohort_R=1.,final_cohort_size=0,
              final_global_R=1,original_cohort_clearance_time=1.,
              integral_local_R_after_original_cohort_clearance=3.,
              integral_cohort_R_after_original_cohort_clearance=0.))
    test("moving_helper_and_pre_event_reverse_production",
         _fixture({0:(1,1),1:(3,10),2:(2,2),3:(3,11)},
                  {0:(3,13),1:(1,1),2:(2,2),5:(3,12)},
                  [revl(1,1,0,10),swap(2,0,5),prod(3,6,5,12),
                   swap(4,5,6,"mixed_L_swap"),revl(5,3,2,11),swap(6,6,1),prod(7,0,1,13)],8), [0],
         dict(outcome_integral_local_R=6.,integral_global_R=16.,integral_local_L=6.,
              integral_global_cohort_R=8.,integral_local_cohort_R=6.,
              cohort_forward_L=2,cohort_reverse_L=1,cohort_net_L=1,
              near_site_forward_L=0,near_site_reverse_L=1,
              cohort_outward_crossings=1,cohort_inward_crossings=1,
              original_cohort_clearance_time=None,final_cohort_size=1))
    test("gross_counterflow_is_not_net_flow",
         _fixture({0:(1,1),1:(2,2)},{2:(1,1),3:(2,2)},
                  [swap(1,0,1),swap(2,1,2),swap(3,2,3),swap(4,0,1),swap(5,1,2),swap(6,2,3)],8), [0,1],
         dict(integral_global_R=16.,outcome_integral_local_R=11.,integral_local_cohort_R=11.,
              cohort_outward_crossings=2,cohort_inward_crossings=1,final_cohort_size=2,
              integral_global_cohort_R1=8.,final_local_R=1))
    test("backmutation_turnover_and_birth_outside_region",
         _fixture({0:(2,1)}, {1:(1,2),2:(2,3)},
                  [copy(1,1,0,1,2,1),kill(2,0,2,1),copy(3,2,1,2,3,2)],5), [0],
         dict(integral_global_R=8.,outcome_integral_local_R=6.,integral_global_cohort_R1=4.,
              original_cohort_clearance_time=2.,cohort_copies=2,cohort_copies_after_original_clearance=1,
              integral_local_R_after_original_cohort_clearance=3.,
              integral_cohort_R_after_original_cohort_clearance=5.,
              cohort_outward_crossings=0,cohort_inward_crossings=0,final_cohort_size=2,final_cohort_R1=1))
    test("periodic_boundary_shell",
         _fixture({0:(1,1)},{15:(1,1)},[swap(1,0,20),swap(2,20,15)],4), [0],
         dict(outcome_integral_local_R=2.,integral_global_cohort_R=4.,
              fixed_readout_region=[0,1,4,5,20],cohort_outward_crossings=1,final_cohort_size=1))
    test("empty_initial_cohort_has_no_turnover_event",
         _fixture({1:(1,7)},{0:(2,8),1:(1,7)},[thermal(1,0,2,8)],3), [0],
         dict(initial_cohort_size=0,original_cohort_present=False,
              original_cohort_cleared_after_presence=False,original_cohort_clearance_time=0.,
              outcome_integral_local_R=5.,integral_global_cohort_R=0.,
              integral_local_R_after_original_cohort_clearance=0.,
              integral_cohort_R_after_original_cohort_clearance=0.,final_cohort_size=0))
    reverse_copy = _event(1,"copy_R",0,1,1,0,removed=1)
    test("reverse_copy_loss_and_unrelated_sources",
         _fixture({0:(1,1),1:(2,2)},{0:(1,3),1:(2,2),2:(2,4)},
                  [reverse_copy,thermal(2,0,1,3),copy(3,2,1,2,4,2)],5), [0],
         dict(outcome_integral_local_R=9.,integral_global_R=11.,integral_global_cohort_R=1.,
              original_cohort_clearance_time=1.,final_cohort_size=0,final_global_R=3,
              cohort_copies=0,integral_local_R_after_original_cohort_clearance=7.))
    return results


def _load(path: Path) -> dict:
    if path.suffix == ".gz":
        with gzip.open(path,"rt") as stream:
            return json.load(stream)
    return json.loads(path.read_text())


def _near(data: dict, explicit: Sequence[int] | None) -> list[int]:
    if explicit is not None:
        return list(map(int,explicit))
    for name in ("near_sites", "fixed_core_sites"):
        if name in data:
            return list(map(int,data[name]))
    for name in ("lesion_metadata", "preparation_metadata", "metadata", "lesion"):
        value = data.get(name,{})
        if isinstance(value,dict):
            if "near" in value and "sites" in value["near"]:
                return list(map(int,value["near"]["sites"]))
            if "metadata" in value and "near" in value["metadata"]:
                return list(map(int,value["metadata"]["near"]["sites"]))
    for name in ("support_readouts", "cohort_support_readouts"):
        if name in data and "fixed_core_sites" in data[name]:
            return list(map(int,data[name]["fixed_core_sites"]))
    raise ValueError("near sites must be supplied explicitly or in immutable raw-run metadata")


def main() -> None:
    parser=argparse.ArgumentParser()
    parser.add_argument("--inputs",nargs="*",default=[])
    parser.add_argument("--near-sites",nargs="+",type=int)
    parser.add_argument("--output",default=str(ROOT/"support-readout-verification.json"))
    args=parser.parse_args()
    sources=[ROOT/"cohort_support_readouts.py",ROOT/"causal-support-design.md",ROOT/"lineage_readouts.py"]
    hashes={str(p.name):sha256(p.read_bytes()).hexdigest() for p in sources}
    tests=fixtures()
    runs=[]
    paths=[]
    for item in args.inputs:
        path=Path(item)
        paths.extend(sorted(path.glob("*.json.gz"))+sorted(path.glob("*.json")) if path.is_dir() else [path])
    for path in paths:
        data=_load(path)
        if "ledger" not in data or "initial_state" not in data:
            continue
        checked=verify_run(data,_near(data,args.near_sites))
        runs.append(dict(path=str(path.resolve()),sha256=sha256(path.read_bytes()).hexdigest(),**checked))
    assert hashes=={str(p.name):sha256(p.read_bytes()).hexdigest() for p in sources}, "source changed during verification"
    result=dict(checks_passed=True,source_hashes=hashes,
                verifier_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                synthetic_fixtures=tests,physical_raw_runs_checked=len(runs),raw_runs=runs,
                scope="Independent offline readout reconstruction. Synthetic histories are bookkeeping fixtures, not physical simulations. No rate, thermodynamic-feasibility, apparatus-cost or individual-boundary claim is made by this verifier.")
    output=Path(args.output)
    assert output.parent.resolve()==ROOT.resolve(), "write only the owned verification output in this directory"
    output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(dict(checks_passed=True,synthetic_fixtures=len(tests),physical_raw_runs_checked=len(runs),
                          output=str(output.resolve()),source_hashes=hashes),indent=2))


if __name__=="__main__":
    main()

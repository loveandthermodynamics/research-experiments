"""Fixed-region and copying-cohort readouts for the external lesion assay."""
from __future__ import annotations

from collections import Counter
from typing import Any, Sequence

from lineage_readouts import base_state


def support_readouts(run: dict[str, Any], near_sites: Sequence[int]) -> dict[str, Any]:
    side = run["params"]["side"]
    core = set(map(int, near_sites))
    region = set(core)
    for i in core:
        a, b = divmod(i, side)
        region.update((((a+1) % side)*side+b, ((a-1) % side)*side+b,
                       a*side+(b+1) % side, a*side+(b-1) % side))
    initial, final = base_state(run["initial_state"]), base_state(run["final_state"])
    material, tags = initial["material"][:], initial["instance_tags"][:]
    original = {tags[i] for i in core if material[i] in (1, 2)}
    ancestry = set(original)
    integrals, counts = Counter(), Counter()
    original_cleared = 0. if not original else None
    previous = 0.

    def accumulate(dt):
        assert dt >= -1e-12
        dt = max(0., dt)
        local_r = sum(material[i] in (1, 2) for i in region)
        global_r = sum(s in (1, 2) for s in material)
        cohort = [i for i, tag in enumerate(tags) if tag in ancestry and material[i] in (1, 2)]
        integrals["local_R"] += dt*local_r
        integrals["global_R"] += dt*global_r
        integrals["local_L"] += dt*sum(material[i] == 3 for i in region)
        integrals["global_cohort_R"] += dt*len(cohort)
        integrals["local_cohort_R"] += dt*sum(i in region for i in cohort)
        integrals["global_cohort_R1"] += dt*sum(material[i] == 2 for i in cohort)
        if original and original_cleared is not None:
            integrals["local_R_after_original_cohort_clearance"] += dt*local_r
            integrals["cohort_R_after_original_cohort_clearance"] += dt*len(cohort)

    for event in run["ledger"]:
        now = float(event["time"])
        accumulate(now-previous)
        previous = now
        kind = event["kind"]
        i, j = int(event.get("i", -1)), int(event.get("j", -1))
        created, removed = int(event.get("instance_created", 0)), int(event.get("instance_removed", 0))
        parent = int(event.get("parent_instance", 0))
        before, after = int(event.get("before", -1)), int(event.get("after", -1))
        if kind == "produce_L":
            # Attribute reverse production using the pre-event helper. The
            # inherited base ledger stores catalyst_instance=0 on reverse.
            if tags[j] in ancestry:
                counts["cohort_forward_L" if after == 3 else "cohort_reverse_L"] += 1
            if j in core:
                counts["near_site_forward_L" if after == 3 else "near_site_reverse_L"] += 1
        if kind in ("material_swap", "mixed_L_swap"):
            for src, dst in ((i, j), (j, i)):
                if tags[src] in ancestry and material[src] in (1, 2):
                    if src in region and dst not in region:
                        counts["cohort_outward_crossings"] += 1
                    if src not in region and dst in region:
                        counts["cohort_inward_crossings"] += 1
            material[i], material[j] = material[j], material[i]
            tags[i], tags[j] = tags[j], tags[i]
        elif kind in ("passive", "copy_R", "produce_L"):
            assert material[i] == before
            if removed:
                assert tags[i] == removed
            if kind == "copy_R" and created and parent in ancestry:
                ancestry.add(created)
                counts["cohort_copies"] += 1
                if original_cleared is not None:
                    counts["cohort_copies_after_original_clearance"] += 1
            material[i], tags[i] = after, created
        if original_cleared is None and not original.intersection(tags):
            original_cleared = now
    accumulate(float(run["horizon"])-previous)
    assert material == final["material"] and tags == final["instance_tags"]
    cohort = [i for i, tag in enumerate(tags) if tag in ancestry and material[i] in (1, 2)]
    return dict(
        outcome_integral_local_R=integrals["local_R"], integral_global_R=integrals["global_R"],
        integral_local_L=integrals["local_L"], integral_global_cohort_R=integrals["global_cohort_R"],
        integral_local_cohort_R=integrals["local_cohort_R"], integral_global_cohort_R1=integrals["global_cohort_R1"],
        integral_local_R_after_original_cohort_clearance=integrals["local_R_after_original_cohort_clearance"],
        integral_cohort_R_after_original_cohort_clearance=integrals["cohort_R_after_original_cohort_clearance"],
        original_cohort_clearance_time=original_cleared, initial_cohort_size=len(original),
        original_cohort_present=bool(original),
        original_cohort_cleared_after_presence=bool(original) and original_cleared is not None,
        final_cohort_size=len(cohort), final_cohort_R1=sum(material[i] == 2 for i in cohort),
        final_local_R=sum(material[i] in (1, 2) for i in region),
        final_global_R=sum(s in (1, 2) for s in material),
        cohort_forward_L=counts["cohort_forward_L"], cohort_reverse_L=counts["cohort_reverse_L"],
        cohort_net_L=counts["cohort_forward_L"]-counts["cohort_reverse_L"],
        near_site_forward_L=counts["near_site_forward_L"], near_site_reverse_L=counts["near_site_reverse_L"],
        cohort_copies=counts["cohort_copies"],
        cohort_copies_after_original_clearance=counts["cohort_copies_after_original_clearance"],
        cohort_outward_crossings=counts["cohort_outward_crossings"],
        cohort_inward_crossings=counts["cohort_inward_crossings"],
        fixed_core_sites=sorted(core), fixed_readout_region=sorted(region),
        scope="External fixed-site lesion, frozen region plus NN shell. Global copying descendants and crossings distinguish lineage loss from transport. This does not establish molecular ownership, a self-produced boundary or a higher individual.")

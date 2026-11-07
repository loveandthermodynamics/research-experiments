"""Offline chemical-state residence, origins and copying genealogy.

No model import, rate call, simulation or output occurs on import. Conserved
carrier B spells are NOT molecular aging. Copying instance genealogy is distinct
from current phenotype and from all later thermal/mutant arrivals.
"""
from __future__ import annotations

from collections import Counter
import math
from typing import Any


ARRAYS = ("material", "templates", "material_tags", "instance_tags")


def _require(test: bool, message: str) -> None:
    if not test:
        raise ValueError(message)


def measure(run: dict[str, Any], appearance_deadline: float = 20.,
            establishment_size: int = 4, *, include_details: bool = False) -> dict[str, Any]:
    """Measure one complete accepted-event history from productive_model.simulate.

    Focal=first X-parent/Y-offspring paid birth by deadline whose parent descends
    by copying from an INITIAL X. Its subtree includes that focal instance and
    all copying descendants over both phenotypes. Samples exactly at a jump time
    are PRE-jump, as in the model's sample loop. Final state is POST all events.
    """
    _require(math.isfinite(appearance_deadline) and appearance_deadline >= 0, "invalid appearance deadline")
    _require(isinstance(establishment_size, int) and establishment_size >= 1, "invalid establishment threshold")
    initial = run["initial_state"]
    horizon = float(run["horizon"])
    start = float(initial.get("time", 0.))
    _require(start == 0. and math.isfinite(horizon) and horizon >= start, "reader requires reset initial time0 and complete finite horizon")
    arrays = {key: list(map(int, initial[key])) for key in ARRAYS}
    size = len(arrays["material"])
    _require(all(len(arrays[key]) == size for key in ARRAYS), "initial array size mismatch")
    _require(all(v in (0,1,2) for v in arrays["material"]+arrays["templates"]), "invalid physical species")
    material_sites: dict[int,int] = {}
    living: dict[int, dict[str, Any]] = {}
    history: dict[int, dict[str, Any]] = {}
    for site in range(size):
        material, tag = arrays["material"][site], arrays["material_tags"][site]
        species, instance = arrays["templates"][site], arrays["instance_tags"][site]
        _require((material == 0) == (tag == 0), "carrier tag/occupancy mismatch")
        _require((species == 0) == (instance == 0), "instance tag/occupancy mismatch")
        if tag:
            _require(tag > 0 and tag not in material_sites, "duplicate carrier tag")
            material_sites[tag] = site
        if instance:
            _require(instance > 0 and instance not in history, "duplicate initial instance")
            data = dict(tag=instance, species=species, birth_time=start, parent=0,
                        origin="initial", initial_X_ancestry=species == 1,
                        focal_ancestry=False, death_time=None, death_kind=None)
            history[instance] = data
            living[instance] = data
    initial_carriers = set(material_sites)
    original_X = {tag for tag, data in living.items() if data["species"] == 1}
    originals_live = set(original_X)
    X_ancestry_live = set(original_X)
    focal_live: set[int] = set()
    next_tag = int(initial["next_tag"])
    _require(next_tag > max(history, default=0), "initial next_tag must exceed all instances")
    counts = Counter(arrays["templates"])
    nB = arrays["material"].count(2)
    nX, nY = counts[1], counts[2]
    n_X_lineage_Y = n_focal_Y = 0
    initial_B, initial_N = nB, nX+nY
    integral = Counter()
    origins, losses, copy_counts, phase_counts = Counter(), Counter(), Counter(), Counter()
    birth_species, death_species = Counter(), Counter()
    episodes: list[dict[str, Any]] = []
    opened: dict[int, tuple[float,bool,str]] = {
        tag:(start,True,"initial_B") for tag,site in material_sites.items()
        if arrays["material"][site] == 2
    }
    episode_entries, episode_exits = Counter(), Counter()
    focal = focal_time = focal_extinction = established_time = None
    focal_Y_first_loss = None
    focal_Y_restorations = max_focal = 0
    founder_clearance = None
    post_clear_copies = post_clear_descendant_copies = post_clear_removed_descendants = 0
    first_post_clear_descendant_copy = None
    first_new_Y = None
    sample_index = 0
    samples = run.get("samples", [])
    sample_times = [float(sample["time"]) for sample in samples]
    _require(all(math.isfinite(t) and start <= t <= horizon for t in sample_times), "invalid sample time")
    _require(all(a <= b for a,b in zip(sample_times,sample_times[1:])), "samples out of order")
    previous = start

    def check_samples(until: float) -> None:
        nonlocal sample_index
        while sample_index < len(samples) and sample_times[sample_index] <= until:
            sample = samples[sample_index]
            for key in ARRAYS:
                if key in sample:
                    _require(list(map(int,sample[key])) == arrays[key], f"sample {sample_index} {key} inconsistent with event history")
            expected = {"B":nB,"X":nX,"Y":nY,"D":size-nX-nY}
            for key,value in expected.items():
                if key in sample:
                    _require(int(sample[key]) == value, f"sample {sample_index} count {key} inconsistent")
            sample_index += 1

    def advance(time: float) -> None:
        nonlocal previous
        _require(math.isfinite(time) and previous <= time <= horizon, "ledger times outside complete ordered horizon")
        dt = time-previous
        integral["B"] += nB*dt
        integral["X"] += nX*dt
        integral["Y"] += nY*dt
        integral["templates"] += (nX+nY)*dt
        integral["initial_X_originals"] += len(originals_live)*dt
        integral["initial_X_lineage"] += len(X_ancestry_live)*dt
        integral["initial_X_lineage_Y"] += n_X_lineage_Y*dt
        integral["focal"] += len(focal_live)*dt
        integral["focal_Y"] += n_focal_Y*dt
        if founder_clearance is not None:
            integral["post_clearance_templates"] += (nX+nY)*dt
            integral["post_clearance_initial_X_lineage"] += len(X_ancestry_live)*dt
        previous = time

    def close_episode(tag: int, time: float, cause: str, censored: bool) -> None:
        _require(tag in opened, "B exit without open carrier episode")
        entry,left,entry_cause = opened.pop(tag)
        episodes.append(dict(carrier=tag,start=entry,end=time,duration=time-entry,
                             left_censored=left,right_censored=censored,
                             entry_cause=entry_cause,exit_cause=cause))
        if not censored:
            episode_exits[cause] += 1

    for event_number,event in enumerate(run["ledger"]):
        time = float(event["time"])
        check_samples(time)
        advance(time)
        kind = event["kind"]
        i,j = int(event["i"]),int(event.get("j",-1))
        sites = [i] if j == -1 else [i,j]
        _require(len(set(sites)) == len(sites) and all(0 <= site < size for site in sites), "invalid local sites")
        before,after = event["before"],event["after"]
        for key in ARRAYS:
            _require(len(before[key]) == len(sites) and len(after[key]) == len(sites), "local array size mismatch")
            _require(list(map(int,before[key])) == [arrays[key][site] for site in sites], f"event {event_number} pre-state {key} mismatch")
        created,removed = int(event.get("instance_created",0)),int(event.get("instance_removed",0))
        parent = int(event.get("parent_instance",0))
        catalyst = int(event.get("catalyst_instance",0))
        origin = event.get("birth_origin")
        activity,direction = int(event.get("activity",0)),int(event.get("direction",0))
        helper,child = int(event.get("helper_species",0)),int(event.get("offspring_species",0))
        _require(not (created and removed), "one event cannot both create and remove an instance")
        expected = {key:list(map(int,before[key])) for key in ARRAYS}
        cause = kind
        is_birth = is_loss = False
        if kind in ("material_swap","template_swap"):
            _require(len(sites)==2 and not any((created,removed,parent,catalyst)), "swap instance metadata invalid")
            keys = ("material","material_tags") if kind=="material_swap" else ("templates","instance_tags")
            for key in keys:
                expected[key].reverse()
        elif kind == "phase_passive":
            _require(len(sites)==1 and before["material"][0] in (1,2), "invalid passive phase substrate")
            _require(not any((created,removed,parent,catalyst)), "phase channel changes instance metadata")
            expected["material"][0] = 3-int(before["material"][0])
            cause = f"phase_passive_{before['material'][0]}_to_{expected['material'][0]}"
            phase_counts[cause] += 1
        elif kind == "template_passive":
            _require(len(sites)==1 and not any((parent,catalyst)), "invalid passive template metadata")
            is_birth = int(before["templates"][0]) == 0
            if is_birth:
                _require(created > 0 and removed == 0 and origin == "thermal", "thermal birth metadata invalid")
                new_species = int(after["templates"][0])
                _require(new_species in (1,2), "thermal birth needs physical species")
                expected["templates"][0] = new_species
                expected["instance_tags"][0] = created
            else:
                is_loss = True
                _require(created == 0 and removed == int(before["instance_tags"][0]) and origin is None, "passive loss removes wrong target")
                expected["templates"][0] = 0
                expected["instance_tags"][0] = 0
        elif kind == "copy":
            _require(len(sites)==2 and activity in (1,2) and direction in (-1,1), "invalid copy activity/direction")
            _require(child in (1,2) and helper == int(before["templates"][1]) and helper in (1,2), "copy helper/offspring species mismatch")
            _require(catalyst == int(before["instance_tags"][1]) and catalyst in living, "copy catalyst must be current live helper")
            substrate,product = (1,2) if activity==1 else (2,1)
            expected["material"][0] = product if direction==1 else substrate
            _require(int(before["material"][0]) == (substrate if direction==1 else product), "copy material substrate mismatch")
            if direction == 1:
                is_birth = True
                _require(int(before["templates"][0]) == 0 and created > 0 and removed == 0,
                         "forward copy target must be empty and create new instance")
                _require(parent == catalyst and origin == "paid_copy", "paid-copy parent/origin mismatch")
                expected["templates"][0] = child
                expected["instance_tags"][0] = created
            else:
                is_loss = True
                _require(int(before["templates"][0]) == child and removed == int(before["instance_tags"][0]), "reverse copy consumes wrong current target")
                _require(created == 0 and parent == 0 and origin is None, "reverse copy is not a birth or historical-parent operation")
                expected["templates"][0] = 0
                expected["instance_tags"][0] = 0
            copy_counts[f"activity_{activity}_{'forward' if direction==1 else 'reverse'}"] += 1
            cause = f"copy_activity_{activity}_{'forward' if direction==1 else 'reverse'}"
        else:
            raise ValueError(f"unknown accepted-event kind: {kind}")
        _require(all(expected[key] == list(map(int,after[key])) for key in ARRAYS), "event post-state differs from declared physical channel")
        if not is_birth and not is_loss:
            _require(origin is None, "non-birth channel has birth origin")

        prior_focal_Y = n_focal_Y
        focal_preexisting = focal is not None
        cleared_before_event = founder_clearance is not None
        if is_loss:
            _require(removed in living and removed == arrays["instance_tags"][i], "loss target is not live")
            data = living.pop(removed)
            losses[kind] += 1
            death_species[str(data["species"])] += 1
            nX -= data["species"] == 1
            nY -= data["species"] == 2
            data["death_time"],data["death_kind"] = time,kind
            if cleared_before_event and kind=="copy" and removed in X_ancestry_live:
                post_clear_removed_descendants += 1
            if removed in X_ancestry_live and data["species"]==2:
                n_X_lineage_Y -= 1
            if removed in focal_live and data["species"]==2:
                n_focal_Y -= 1
            originals_live.discard(removed)
            X_ancestry_live.discard(removed)
            focal_live.discard(removed)
        if is_birth:
            _require(created == next_tag and created not in history, "instance IDs must be fresh monotonically allocated")
            next_tag += 1
            new_species = int(after["templates"][0])
            from_X = origin=="paid_copy" and history[parent]["initial_X_ancestry"]
            from_focal = origin=="paid_copy" and history[parent]["focal_ancestry"]
            eligible = (origin=="paid_copy" and helper==1 and new_species==2 and
                        from_X and time<=appearance_deadline)
            if focal is None and eligible:
                focal,focal_time = created,time
                from_focal = True
            if from_focal:
                _require(focal_extinction is None, "extinct focal subtree cannot resurrect")
            data = dict(tag=created,species=new_species,birth_time=time,parent=parent,
                        origin=origin,initial_X_ancestry=bool(from_X),
                        focal_ancestry=bool(from_focal),death_time=None,death_kind=None)
            history[created] = living[created] = data
            if from_X:
                X_ancestry_live.add(created)
                n_X_lineage_Y += new_species==2
            if from_focal:
                focal_live.add(created)
                n_focal_Y += new_species==2
            nX += new_species==1
            nY += new_species==2
            origins[origin] += 1
            birth_species[f"{origin}_{new_species}"] += 1
            if origin == "paid_copy":
                birth_species[f"paid_helper_{helper}_offspring_{new_species}"] += 1
                if cleared_before_event:
                    post_clear_copies += 1
                    if from_X:
                        post_clear_descendant_copies += 1
                        if first_post_clear_descendant_copy is None:
                            first_post_clear_descendant_copy = time
            if new_species==2 and first_new_Y is None:
                first_new_Y = dict(tag=created,time=time,origin=origin,parent=parent,
                                   parent_initial_X_ancestry=bool(from_X))
        if original_X and not originals_live and founder_clearance is None:
            founder_clearance = time

        # Follow conserved material identities across the entire local change;
        # a B/A swap moves two carriers rather than ending/starting two spells.
        before_by_tag = {int(tag):int(mat) for tag,mat in zip(before["material_tags"],before["material"]) if tag}
        after_by_tag = {int(tag):int(mat) for tag,mat in zip(after["material_tags"],after["material"]) if tag}
        _require(set(before_by_tag) == set(after_by_tag), "event creates or destroys a material carrier")
        for tag,old in before_by_tag.items():
            new = after_by_tag[tag]
            if old==2 and new!=2:
                close_episode(tag,time,cause,False)
                nB -= 1
            elif old!=2 and new==2:
                _require(tag not in opened, "duplicate open B episode")
                opened[tag] = (time,False,cause)
                episode_entries[cause] += 1
                nB += 1
        # Remove local carrier locations before adding their moved locations.
        for tag in before_by_tag:
            material_sites.pop(tag)
        for local,site in enumerate(sites):
            for key in ARRAYS:
                arrays[key][site] = int(after[key][local])
            tag = arrays["material_tags"][site]
            if tag:
                _require(tag not in material_sites, "material tag duplicated outside local event")
                material_sites[tag] = site
        if focal is not None:
            current_Y = n_focal_Y
            if current_Y==0 and focal_Y_first_loss is None:
                focal_Y_first_loss = time
            if focal_preexisting and prior_focal_Y==0 and current_Y>0:
                focal_Y_restorations += 1
            max_focal = max(max_focal,len(focal_live))
            if established_time is None and len(focal_live)>=establishment_size:
                established_time = time
            if not focal_live and focal_extinction is None:
                focal_extinction = time

    check_samples(horizon)
    advance(horizon)
    _require(sample_index == len(samples), "sample history incomplete")
    final = run["final_state"]
    for key in ARRAYS:
        _require(arrays[key] == list(map(int,final[key])), f"final {key} differs from complete event history")
    _require(int(final["next_tag"]) == next_tag, "final instance counter mismatch")
    _require(set(material_sites) == initial_carriers, "global material tags changed")
    _require(len(living)==nX+nY and set(living)==set(tag for tag in arrays["instance_tags"] if tag), "final global living-set mismatch")
    _require(n_X_lineage_Y == sum(living[t]["species"]==2 for t in X_ancestry_live), "initial-X phenotype counter mismatch")
    _require(n_focal_Y == sum(living[t]["species"]==2 for t in focal_live), "focal phenotype counter mismatch")
    for tag in sorted(opened):
        close_episode(tag,horizon,"observation_end",True)
    durations = [episode["duration"] for episode in episodes]
    duration_sum = math.fsum(durations)
    _require(math.isclose(duration_sum,integral["B"],rel_tol=2e-11,abs_tol=1e-8), "truncated B episodes do not exhaust event-integrated B-time")
    netX = copy_counts["activity_1_forward"]-copy_counts["activity_1_reverse"]
    netY = copy_counts["activity_2_forward"]-copy_counts["activity_2_reverse"]
    net0 = phase_counts["phase_passive_1_to_2"]-phase_counts["phase_passive_2_to_1"]
    _require(nB-initial_B == netX-netY+net0, "signed material-current balance fails")
    _require(nX+nY-initial_N == netX+netY+origins["thermal"]-losses["template_passive"], "signed template-renewal balance fails")
    fully_observed = [episode["duration"] for episode in episodes if not episode["left_censored"] and not episode["right_censored"]]
    summary = dict(episodes=len(episodes),left_censored=sum(e["left_censored"] for e in episodes),
                   right_censored=sum(e["right_censored"] for e in episodes),
                   fully_observed_episodes=len(fully_observed),truncated_duration_sum=duration_sum,
                   truncated_duration_mean=duration_sum/len(episodes) if episodes else None,
                   fully_observed_duration_mean=math.fsum(fully_observed)/len(fully_observed) if fully_observed else None,
                   truncated_duration_min=min(durations) if durations else None,
                   truncated_duration_max=max(durations) if durations else None,
                   entry_causes=dict(episode_entries),exit_causes=dict(episode_exits),
                   scope="Carrier chemical-state residence with explicit observation censoring; neither molecular aging nor an uncensored stationary lifetime estimate.")
    result = dict(
        horizon=horizon,events=len(run["ledger"]),samples_checked=len(samples),
        integral_B=integral["B"],integral_X=integral["X"],integral_Y=integral["Y"],
        integral_templates=integral["templates"],B_episode_summary=summary,
        births=dict(origins),birth_species=dict(birth_species),losses=dict(losses),death_species=dict(death_species),
        copy_channel_counts=dict(copy_counts),phase_channel_counts=dict(phase_counts),
        net_X_activity_copy_count=netX,net_Y_activity_copy_count=netY,net_background_A_to_B_count=net0,
        first_new_Y_beyond_initial=first_new_Y,
        appearance_deadline=appearance_deadline,establishment_size=establishment_size,
        focal_paid_Y_from_initial_X_appeared_by_deadline=focal is not None,
        focal_tag=focal,focal_appearance_time=focal_time,
        focal_established=established_time is not None,focal_establishment_time=established_time,
        maximum_live_focal_subtree=max_focal,
        focal_population_definition="Focal instance plus all paid copying descendants, over both current phenotypes; no unrelated thermal/mutant arrivals.",
        focal_extinction_time=focal_extinction,focal_survives_to_horizon=bool(focal_live),
        final_focal_subtree=len(focal_live),final_focal_Y=n_focal_Y,
        integral_focal_subtree=integral["focal"],integral_focal_Y=integral["focal_Y"],
        focal_Y_first_loss_time=focal_Y_first_loss,focal_Y_restorations_after_loss=focal_Y_restorations,
        focal_Y_uninterrupted_to_horizon=focal is not None and focal_Y_first_loss is None,
        initial_X_founders=len(original_X),initial_X_founder_cohort_present=bool(original_X),
        original_X_founders_cleared_after_presence=founder_clearance is not None,
        original_X_founder_clearance_time=founder_clearance,final_original_X_founders=len(originals_live),
        final_initial_X_lineage=len(X_ancestry_live),
        integral_original_X_founders=integral["initial_X_originals"],
        integral_initial_X_lineage=integral["initial_X_lineage"],
        integral_initial_X_lineage_Y=integral["initial_X_lineage_Y"],
        integral_templates_after_original_X_clearance=integral["post_clearance_templates"],
        integral_initial_X_lineage_after_original_X_clearance=integral["post_clearance_initial_X_lineage"],
        forward_paid_copies_after_original_X_clearance=post_clear_copies,
        forward_paid_copies_by_initial_X_descendants_after_clearance=post_clear_descendant_copies,
        reverse_consumed_initial_X_descendants_after_clearance=post_clear_removed_descendants,
        first_paid_copy_by_initial_X_descendant_after_clearance=first_post_clear_descendant_copy,
        final_X=nX,final_Y=nY,final_B=nB,
        checks=dict(complete_local_replay=True,carrier_conservation=True,fresh_instance_ids=True,
                    signed_material_current_balance=True,signed_template_current_balance=True,
                    truncated_episode_integral_identity=True,sample_arrays_and_counts=True),
        scope="Analysis-only complete-event chemistry and genealogy. Occupancy, paid origin, ancestry, phenotype continuity, finite-horizon establishment and causal functional benefit are distinct. Founder clearance windows are descriptive, not an aging or efficiency comparison. No readout enters a physical rate.",
    )
    if include_details:
        result["B_episode_records"] = episodes
        result["instance_histories"] = [history[tag] for tag in sorted(history)]
    return result


analyze_renewal = measure

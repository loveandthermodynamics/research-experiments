"""Independent audit of frozen joined study; replay only, never simulation."""
from __future__ import annotations

from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.stats import t

sys.dont_write_bytecode=True
from verify_productive import ROOT,Reference,compare,close,digest,load,verify_trajectory,observe_spatial_independent

KEYS=("material","templates","material_tags","instance_tags")
METRICS=("B_fraction","large_B_mass_fraction","excess_neighbor_B_fraction","spectral_wavelength_low_k",
         "largest_B_fraction","mass_weighted_cluster_size","template_density","template_density_near_dense_B",
         "n_X","n_Y","Y_on_B_excess_fraction","clusters_ge_9")


def renew_from_lifetimes(run,deadline=20.,threshold=4):
    """Offline birth/death intervals and parent DAG, rather than live counters.

    Material spell endpoints follow carrier IDs. A separate physical replay
    verifies every local array and event before these observational calculations.
    """
    initial,H=run["initial_state"],run["horizon"]
    lives={tag:dict(species=r,parent=0,origin="initial",birth=0.,death=H,born=-1,removed=None)
           for tag,r in zip(initial["instance_tags"],initial["templates"]) if tag}
    originals={tag for tag,data in lives.items() if data["species"]==1}
    lineage={tag for tag in originals}
    focal=None; focal_time=None; cohort=set()
    thermal=Counter();birthspecies=Counter();deathspecies=Counter();losses=Counter();copies=Counter();phases=Counter()
    open_b={tag:(0.,True,"initial_B") for tag,mat in zip(initial["material_tags"],initial["material"]) if tag and mat==2}
    episodes=[]; entries=Counter();exits=Counter();first_y=None
    birth_events=[];reverse_events=[]
    for order,e in enumerate(run["ledger"]):
        now=e["time"];kind=e["kind"]
        cause=kind
        if kind=="copy":
            label=f"activity_{e['activity']}_{'forward' if e['direction']==1 else 'reverse'}"
            copies[label]+=1;cause="copy_"+label
        elif kind=="phase_passive":
            cause=f"phase_passive_{e['before']['material'][0]}_to_{e['after']['material'][0]}";phases[cause]+=1
        created,removed=e["instance_created"],e["instance_removed"]
        if created:
            r=e["after"]["templates"][0];parent=e["parent_instance"];origin=e["birth_origin"]
            assert created not in lives and (parent==0 or parent in lives)
            lives[created]=dict(species=r,parent=parent,origin=origin,birth=now,death=H,born=order,removed=None)
            from_x=parent in lineage
            if from_x:lineage.add(created)
            if focal is None and origin=="paid_copy" and e["helper_species"]==1 and r==2 and from_x and now<=deadline:
                focal,focal_time=created,now;cohort.add(created)
            elif parent in cohort:cohort.add(created)
            thermal[origin]+=1;birthspecies[f"{origin}_{r}"]+=1
            if origin=="paid_copy":birthspecies[f"paid_helper_{e['helper_species']}_offspring_{r}"]+=1
            birth_events.append((now,order,created,parent,origin))
            if r==2 and first_y is None:
                first_y=dict(tag=created,time=now,origin=origin,parent=parent,parent_initial_X_ancestry=from_x)
        if removed:
            assert lives[removed]["removed"] is None
            lives[removed]["death"]=now;lives[removed]["removed"]=order
            deathspecies[str(lives[removed]["species"])]+=1;losses[kind]+=1
            if kind=="copy":reverse_events.append((now,order,removed))
        before={tag:mat for tag,mat in zip(e["before"]["material_tags"],e["before"]["material"]) if tag}
        after={tag:mat for tag,mat in zip(e["after"]["material_tags"],e["after"]["material"]) if tag}
        assert before.keys()==after.keys()
        for tag,old in before.items():
            new=after[tag]
            if old==2 and new!=2:
                entry,left,entrycause=open_b.pop(tag)
                episodes.append(dict(carrier=tag,start=entry,end=now,duration=now-entry,left_censored=left,right_censored=False,
                                     entry_cause=entrycause,exit_cause=cause));exits[cause]+=1
            elif old!=2 and new==2:
                assert tag not in open_b;open_b[tag]=(now,False,cause);entries[cause]+=1
    for tag,(start,left,cause) in sorted(open_b.items()):
        episodes.append(dict(carrier=tag,start=start,end=H,duration=H-start,left_censored=left,right_censored=True,
                             entry_cause=cause,exit_cause="observation_end"))
    # Monotone allocated IDs let this ancestry set include every backmutation.
    duration=lambda tag:lives[tag]["death"]-lives[tag]["birth"]
    integral=lambda tags:math.fsum(duration(tag) for tag in tags)
    living=lambda tag:lives[tag]["removed"] is None
    final_x=sum(living(tag) and d["species"]==1 for tag,d in lives.items())
    final_y=sum(living(tag) and d["species"]==2 for tag,d in lives.items())
    clearance_order=max((lives[tag]["removed"] for tag in originals),default=None) if originals and all(not living(tag) for tag in originals) else None
    clearance=run["ledger"][clearance_order]["time"] if clearance_order is not None else None
    post_interval=lambda tag:max(0.,lives[tag]["death"]-max(lives[tag]["birth"],clearance)) if clearance is not None else 0.
    peak=0;established=extinction=first_loss=None;restorations=0
    live_count=y_count=0
    # Events retain original order, so birth/death endpoints at equal times are
    # treated like the source's post-event observations rather than averaged.
    endpoints=[]
    for tag in cohort:
        d=lives[tag];endpoints.append((d["born"],d["birth"],1,d["species"]))
        if d["removed"] is not None:endpoints.append((d["removed"],d["death"],-1,d["species"]))
    for order,now,sign,species in sorted(endpoints):
        old_y=y_count;live_count+=sign;y_count+=sign*(species==2);peak=max(peak,live_count)
        if established is None and live_count>=threshold:established=now
        if first_loss is None and y_count==0:first_loss=now
        if old_y==0 and y_count>0 and order!=lives[focal]["born"]:restorations+=1
        if live_count==0 and extinction is None:extinction=now
    durations=[d["duration"] for d in episodes]
    fully=[d["duration"] for d in episodes if not d["left_censored"] and not d["right_censored"]]
    cx=copies["activity_1_forward"]-copies["activity_1_reverse"]
    cy=copies["activity_2_forward"]-copies["activity_2_reverse"]
    c0=phases["phase_passive_1_to_2"]-phases["phase_passive_2_to_1"]
    post_paid=[row for row in birth_events if clearance_order is not None and row[1]>clearance_order and row[4]=="paid_copy"]
    post_desc=[row for row in post_paid if row[2] in lineage]
    return dict(
        integral_B=math.fsum(durations),integral_X=integral(tag for tag,d in lives.items() if d["species"]==1),
        integral_Y=integral(tag for tag,d in lives.items() if d["species"]==2),integral_templates=integral(lives),
        births=dict(thermal),birth_species=dict(birthspecies),losses=dict(losses),death_species=dict(deathspecies),
        copy_channel_counts=dict(copies),phase_channel_counts=dict(phases),net_X_activity_copy_count=cx,
        net_Y_activity_copy_count=cy,net_background_A_to_B_count=c0,first_new_Y_beyond_initial=first_y,
        appearance_deadline=deadline,establishment_size=threshold,focal_paid_Y_from_initial_X_appeared_by_deadline=focal is not None,
        focal_tag=focal,focal_appearance_time=focal_time,focal_established=established is not None,focal_establishment_time=established,
        maximum_live_focal_subtree=peak,focal_extinction_time=extinction,focal_survives_to_horizon=any(living(tag) for tag in cohort),
        final_focal_subtree=sum(living(tag) for tag in cohort),final_focal_Y=sum(living(tag) and lives[tag]["species"]==2 for tag in cohort),
        integral_focal_subtree=integral(cohort),integral_focal_Y=integral(tag for tag in cohort if lives[tag]["species"]==2),
        focal_Y_first_loss_time=first_loss,focal_Y_restorations_after_loss=restorations,
        focal_Y_uninterrupted_to_horizon=focal is not None and first_loss is None,
        initial_X_founders=len(originals),initial_X_founder_cohort_present=bool(originals),
        original_X_founders_cleared_after_presence=clearance is not None,original_X_founder_clearance_time=clearance,
        final_original_X_founders=sum(living(tag) for tag in originals),final_initial_X_lineage=sum(living(tag) for tag in lineage),
        integral_original_X_founders=integral(originals),integral_initial_X_lineage=integral(lineage),
        integral_initial_X_lineage_Y=integral(tag for tag in lineage if lives[tag]["species"]==2),
        integral_templates_after_original_X_clearance=math.fsum(post_interval(tag) for tag in lives),
        integral_initial_X_lineage_after_original_X_clearance=math.fsum(post_interval(tag) for tag in lineage),
        forward_paid_copies_after_original_X_clearance=len(post_paid),forward_paid_copies_by_initial_X_descendants_after_clearance=len(post_desc),
        reverse_consumed_initial_X_descendants_after_clearance=sum(clearance_order is not None and order>clearance_order and tag in lineage for now,order,tag in reverse_events),
        first_paid_copy_by_initial_X_descendant_after_clearance=post_desc[0][0] if post_desc else None,
        final_X=final_x,final_Y=final_y,final_B=initial["material"].count(2)+cx-cy+c0,
        B_episode_summary=dict(episodes=len(episodes),left_censored=sum(e["left_censored"] for e in episodes),
            right_censored=sum(e["right_censored"] for e in episodes),fully_observed_episodes=len(fully),
            truncated_duration_sum=math.fsum(durations),truncated_duration_mean=math.fsum(durations)/len(durations) if durations else None,
            fully_observed_duration_mean=math.fsum(fully)/len(fully) if fully else None,truncated_duration_min=min(durations) if durations else None,
            truncated_duration_max=max(durations) if durations else None,entry_causes=dict(entries),exit_causes=dict(exits)),
        B_episode_records=sorted(episodes,key=lambda e:(e["carrier"],e["start"],e["end"])))


def check_lifetimes():
    from renewal_fixtures import histories
    from renewal_readouts import measure
    for name,run,expected in histories():
        actual=measure(run,include_details=True);expected=renew_from_lifetimes(run)
        actual["B_episode_records"]=sorted(actual["B_episode_records"],key=lambda e:(e["carrier"],e["start"],e["end"]))
        compare(expected,actual,name)


def verify_initial(run):
    p=run["params"];seed=run["seed"];n=p["side"]**2
    arrays={key:[0]*n for key in KEYS}
    for label,stream,values,count in (("material",85601,(1,),p["material_count"]),
                                     ("templates",85602,(1,2),p["initial_X"]+p["initial_Y"])):
        rng=np.random.default_rng(np.random.SeedSequence([seed,stream]))
        sites=rng.choice(n,count,replace=False)
        tags="material_tags" if label=="material" else "instance_tags"
        for index,site in enumerate(sites):
            arrays[label][int(site)]=1 if label=="material" or index<p["initial_X"] else 2
            arrays[tags][int(site)]=index+1
    expected=dict(**arrays,next_tag=p["initial_X"]+p["initial_Y"]+1,time=0.,cumulative_work=0.,cumulative_heat=0.,net_fuel_to_waste=0,event_count=0)
    compare(expected,run["initial_state"])
    assert 2 not in arrays["material"] and 2 not in arrays["templates"]
    return digest_bytes(expected)


def digest_bytes(value):
    from hashlib import sha256
    return sha256(json.dumps(value,sort_keys=True,separators=(",",":")).encode()).hexdigest()


def finite_bound(p,horizon):
    # Keep the constant-source and depletion pieces separate until the final
    # division; this expression was independently certified before outcomes.
    c0=math.exp(-p["E_X"]/p["theta"])+math.exp(-p["E_Y"]/p["theta"])
    cb=math.exp(-(p["E_X"]-p["g_X"])/p["theta"])+math.exp(-(p["E_Y"]-p["g_Y"])/p["theta"])
    assert p["background_AB_rate"]==0 and p["E_X"]>=max(0,p["g_X"]) and p["E_Y"]>=max(0,p["g_Y"])
    numerator=c0*p["side"]**2+max(0.,cb-c0)*p["material_count"]
    numerator+=(p["material_count"]+p["initial_X"]+p["initial_Y"])/(p["loss_rate"]*horizon)
    return numerator/(1+c0)


def arithmetic_mean(values):return None if any(x is None for x in values) else math.fsum(values)/len(values)
def ratio(a,b):return None if a is None or b is None or b==0 else a/b
def difference(a,b):return None if a is None or b is None else a-b


def analyze_one(run):
    side=run["params"]["side"];n=side**2
    series=[]
    assert [s["time"] for s in run["samples"]]==[float(t) for t in range(0,401,5)]
    for sample in run["samples"]:
        row=observe_spatial_independent(np.asarray(sample["material"]).reshape(side,side),np.asarray(sample["templates"]).reshape(side,side))
        row["time"]=sample["time"];series.append(row)
    windows={}
    for name,low,high in (("early",100.,200.),("late",300.,400.)):
        rows=[s for s in series if low<=s["time"]<=high];assert len(rows)==21
        windows[name]={key:arithmetic_mean([row[key] for row in rows]) for key in METRICS}
    renewal=renew_from_lifetimes(run);renewal.pop("B_episode_records")
    signed={1:0,2:0}
    for event in run["ledger"]:
        if event["kind"]=="copy" and 300.<event["time"]<=400.:signed[event["activity"]]+=event["direction"]
    integrals=dict(template_time=renewal["integral_templates"],B_time=renewal["integral_B"],
                   whole_mean_template_density=renewal["integral_templates"]/run["horizon"]/n,
                   late_net_X_current_density=signed[1]/(100*n),late_net_Y_current_density=signed[2]/(100*n))
    return dict(side=side,arm=run["arm"],seed=run["seed"],windows=windows,timeseries=series,
                growth=ratio(windows["late"]["spectral_wavelength_low_k"],windows["early"]["spectral_wavelength_low_k"]),
                integrals=integrals,renewal=renewal,finite_ablation_bound_density=finite_bound(run["params"],run["horizon"])/n,
                physical_final=run["final_summary"])


def verify_record(row):
    protocol=load(ROOT/"protocol.json");phash=digest(ROOT/"protocol.json")
    path=ROOT/row["path"];assert digest(path)==row["sha256"]
    run=load(path);assert run["protocol_sha256"]==phash
    config=protocol["configurations"][f"{row['side']}-{row['arm']}"]
    compare(config,run)
    assert run["seed"]==row["seed"] and run["arm"]==row["arm"] and run["horizon"]==protocol["horizon"]
    initial_hash=verify_initial(run)
    physical=verify_trajectory(run,True)
    analysis=analyze_one(run)
    assert len(run["ledger"])==row["events"] and len(run["samples"])==row["samples"]
    close(row["proposals"],run["proposal_statistics"]["proposals"])
    report=dict(path=row["path"],sha256=row["sha256"],side=row["side"],arm=row["arm"],seed=row["seed"],
                checks_passed=True,initial_state_sha256=initial_hash,physical=physical,analysis=analysis,
                verifier_sha256=digest(Path(__file__)),reference_sha256=digest(ROOT/"verify_productive.py"))
    target=ROOT/"audit-records"/(path.name.replace(".json.gz","-verification.json"))
    target.parent.mkdir(exist_ok=True)
    target.write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    return dict(path=row["path"],events=physical["events"],proposals=physical["proposals"],checks_passed=True)


def check_bindings():
    protocol=load(ROOT/"protocol.json");pre=load(ROOT/"preflight-verification.json")
    assert pre["checks_passed"] and digest(ROOT/"preflight-verification.json")==protocol["independent_preflight_sha256"]
    for name,value in protocol["source_sha256"].items():assert digest(ROOT/name)==value,name
    for name,value in pre["source_sha256"].items():assert digest(ROOT/name)==value,name
    assert protocol["seeds"]==list(range(730100,730108)) and protocol["sides"]==[16,32]
    assert protocol["arms"]==["specialist","generalist","X_only","Y_only"]
    assert protocol["required_claim_A_gates"]==20 and protocol["required_claim_B_gates"]==2
    assert protocol["horizon"]==400 and protocol["sample_dt"]==5
    return protocol


def replay(workers=2):
    check_bindings();check_lifetimes()
    execution=load(ROOT/("execution.json" if (ROOT/"execution.json").exists() else "execution-progress.json"))
    jobs=[]
    for row in execution["records"]:
        cache=ROOT/"audit-records"/(Path(row["path"]).name.replace(".json.gz","-verification.json"))
        if cache.exists():
            previous=load(cache)
            assert previous["sha256"]==row["sha256"] and previous["checks_passed"]
            assert previous["verifier_sha256"]==digest(Path(__file__))
        else:jobs.append(row)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for future in as_completed([pool.submit(verify_record,row) for row in jobs]):print(json.dumps(future.result()),flush=True)


def estimate(values):
    if any(value is None or not math.isfinite(value) for value in values):
        return dict(estimable=False,n=len(values),values=values,reason="undefined values retained; no seed filtering")
    n=len(values);average=math.fsum(values)/n
    sem=math.sqrt(math.fsum((value-average)**2 for value in values)/(n*(n-1)))
    width=float(t.ppf(.975,n-1))*sem
    indices=np.random.default_rng(841701).integers(n,size=(10000,n))
    boot=np.asarray(values)[indices].mean(axis=1)
    return dict(estimable=True,n=n,values=values,mean=average,sem=sem,t95=[average-width,average+width],
                paired_seed_bootstrap95=[float(x) for x in np.percentile(boot,[2.5,97.5])])


def gate(name,values,relation,threshold):
    result=estimate(values);result.update(name=name,relation=relation,threshold=threshold)
    result["met"]=bool(result["estimable"] and (result["t95"][0]>threshold if relation=="lower_gt" else result["t95"][1]<=threshold))
    return result


def finish():
    protocol=check_bindings();execution=load(ROOT/"execution.json");analysis=load(ROOT/"analysis.json")
    expected={(side,arm,seed) for side in protocol["sides"] for arm in protocol["arms"] for seed in protocol["seeds"]}
    assert len(execution["records"])==64 and {(r["side"],r["arm"],r["seed"]) for r in execution["records"]}==expected
    records=[]
    for row in execution["records"]:
        cached=load(ROOT/"audit-records"/(Path(row["path"]).name.replace(".json.gz","-verification.json")))
        assert cached["verifier_sha256"]==digest(Path(__file__)) and cached["sha256"]==row["sha256"] and cached["checks_passed"]
        assert digest(ROOT/row["path"])==row["sha256"]
        records.append(cached)
    for side in (16,32):
        for seed in protocol["seeds"]:
            paired=[r for r in records if r["side"]==side and r["seed"]==seed]
            assert len(paired)==4 and len({r["initial_state_sha256"] for r in paired})==1
    actual_index={(r["side"],r["arm"],r["seed"]):r for r in analysis["rows"]}
    rows=[r["analysis"] for r in records]
    for row in rows:compare(row,actual_index[(row["side"],row["arm"],row["seed"])])
    index={(r["side"],r["arm"],r["seed"]):r for r in rows};seeds=protocol["seeds"]
    get=lambda side,arm,seed:index[(side,arm,seed)]
    gates_a=[];gates_b=[]
    for side in (16,32):
        for metric,threshold in (("B_fraction",.05),("large_B_mass_fraction",.5),("excess_neighbor_B_fraction",.15)):
            gates_a.append(gate(f"L{side}_late_{metric}",[get(side,"specialist",s)["windows"]["late"][metric] for s in seeds],"lower_gt",threshold))
        gates_a.append(gate(f"L{side}_wavelength_growth",[get(side,"specialist",s)["growth"] for s in seeds],"upper_le",1.25))
        for activity in ("X","Y"):
            name=f"late_net_{activity}_current_density"
            gates_a.append(gate(f"L{side}_{name}",[get(side,"specialist",s)["integrals"][name] for s in seeds],"lower_gt",0.))
        for control in ("X_only","Y_only"):
            gates_a.append(gate(f"L{side}_late_template_density_minus_{control}",
                [difference(get(side,"specialist",s)["windows"]["late"]["template_density"],get(side,control,s)["windows"]["late"]["template_density"]) for s in seeds],"lower_gt",0.))
        gates_a.append(gate(f"L{side}_whole_template_density_minus_finite_ablation_bound",[get(side,"specialist",s)["integrals"]["whole_mean_template_density"]-get(side,"specialist",s)["finite_ablation_bound_density"] for s in seeds],"lower_gt",0.))
        gates_b.append(gate(f"L{side}_specialist_minus_generalist_late_dense_associated_template_density",[difference(get(side,"specialist",s)["windows"]["late"]["template_density_near_dense_B"],get(side,"generalist",s)["windows"]["late"]["template_density_near_dense_B"]) for s in seeds],"lower_gt",0.))
    gates_a.append(gate("specialist_late_wavelength_L32_over_L16",[ratio(get(32,"specialist",s)["windows"]["late"]["spectral_wavelength_low_k"],get(16,"specialist",s)["windows"]["late"]["spectral_wavelength_low_k"]) for s in seeds],"upper_le",1.5))
    gates_a.append(gate("specialist_late_largest_B_fraction_L16_minus_L32",[difference(get(16,"specialist",s)["windows"]["late"]["largest_B_fraction"],get(32,"specialist",s)["windows"]["late"]["largest_B_fraction"]) for s in seeds],"lower_gt",0.))
    assert len(gates_a)==20 and len(gates_b)==2
    compare(gates_a,analysis["claim_A_gates"]);compare(gates_b,analysis["claim_B_gates"])
    assert all(g["met"] for g in gates_a)==analysis["claim_A_met"] and all(g["met"] for g in gates_b)==analysis["claim_B_met"]
    summaries={}
    for side in (16,32):
        for arm in protocol["arms"]:
            selected=[get(side,arm,s) for s in seeds]
            summaries[f"L{side}-{arm}"]=dict(late={m:estimate([r["windows"]["late"][m] for r in selected]) for m in METRICS},
                growth=estimate([r["growth"] for r in selected]),whole_mean_template_density=estimate([r["integrals"]["whole_mean_template_density"] for r in selected]),
                chemical_work_per_area=estimate([r["physical_final"]["chemical_work"]/side**2 for r in selected]),
                heat_to_bath_per_area=estimate([r["physical_final"]["heat_to_bath"]/side**2 for r in selected]),
                net_X_current=estimate([r["integrals"]["late_net_X_current_density"] for r in selected]),net_Y_current=estimate([r["integrals"]["late_net_Y_current_density"] for r in selected]),
                finite_ablation_bound_density=selected[0]["finite_ablation_bound_density"],renewal=[r["renewal"] for r in selected])
    compare(summaries,analysis["summary"])
    totals=dict(unique_records=len(records),independent_seed_rows=8,events=sum(r["physical"]["events"] for r in records),
                proposals=sum(r["physical"]["proposals"] for r in records),samples=sum(r["physical"]["samples"] for r in records))
    assert totals["events"]==execution["events"] and totals["samples"]==execution["samples"]
    output=dict(checks_passed=True,protocol_sha256=digest(ROOT/"protocol.json"),execution_sha256=digest(ROOT/"execution.json"),
                analysis_sha256=digest(ROOT/"analysis.json"),verifier_sha256=digest(Path(__file__)),reference_sha256=digest(ROOT/"verify_productive.py"),
                totals=totals,all_accepted_events_and_rejected_proposals_replayed=True,all_initial_states_regenerated=True,
                physical_arrays_and_resources_matched_across_arms=True,all_source_hashes_unchanged=True,all_histories_retained=True,
                independent_union_find_and_real_half_spectrum=True,independent_birth_death_lifetime_genealogy=True,
                claim_A_met=analysis["claim_A_met"],claim_B_met=analysis["claim_B_met"],claim_A_gates=gates_a,claim_B_gates=gates_b,
                summaries=summaries,records=[{key:value for key,value in r.items() if key!="analysis"} for r in records],
                scope="Complete finite-horizon software/physical/statistical verification of the declared study. Domain co-occurrence is distinct from tracked component renewal or integrity; association is not taxon necessity, efficiency, inherited individuality, informational aging or moral agency. Initial preparation and chemostat construction are supplied and excluded from chemical work. Tiny preflight fixtures are excluded from outcome totals.")
    (ROOT/"independent-verification.json").write_text(json.dumps(output,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print(json.dumps({key:output[key] for key in ("checks_passed","totals","claim_A_met","claim_B_met")}),flush=True)


if __name__=="__main__":
    if sys.argv[1:]==["--test-readouts"]:check_lifetimes();print("independent lifetime reader passed11manual histories")
    elif sys.argv[1:2]==["--replay"]:replay(int(sys.argv[2]) if len(sys.argv)>2 else 2)
    elif sys.argv[1:]==["--finish"]:finish()
    else:raise SystemExit("use --test-readouts, --replay [workers], or --finish")

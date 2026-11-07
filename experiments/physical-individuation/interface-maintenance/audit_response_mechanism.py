"""Read-only reconstruction of frozen response-ensemble mechanism diagnostics.

No model imports, new simulation, output file, or source/result mutation.
Run: python3 audit_response_mechanism.py
Optional: python3 audit_response_mechanism.py --summary /path/to/summary.json
Trajectory paths are resolved relative to the summary's parent directory.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
from pathlib import Path
import statistics


def mean(values):
    values=list(values)
    return statistics.mean(values) if values else None


def replay(path: Path):
    with gzip.open(path,"rt") as handle:
        run=json.load(handle)
    p=run["params"]
    side=p["side"]
    state=run["initial_state"]["base_state"]
    material=list(state["material"])
    fuel=list(state["fuel"])
    gates=list(run["initial_state"]["gates"])
    neighbors=[[
        ((j//side-1)%side)*side+j%side,
        ((j//side+1)%side)*side+j%side,
        (j//side)*side+(j%side-1)%side,
        (j//side)*side+(j%side+1)%side,
    ] for j in range(side*side)]
    first_gate0_bad=None
    first_paid_variant=None
    first_extinction_fuel=None
    paid_R=paid_L=0
    both_seen=False
    lost_with_partner=False
    r1time=c1time=protected_r1time=0.0
    r1_copies_gate10=r1_to_r0_gate10=r1_L_gate10=0
    r0_reintroduced_by_r1=0
    after_gate0_R=after_gate0_L=0
    last=0.0

    def exposures(dt):
        if gates!=[1,0]:
            return 0.0,0.0,0.0
        cells=[j for j,a in enumerate(material) if a==2]
        c=sum(material[i]==0 and fuel[i]==1 for j in cells for i in neighbors[j])
        protected=sum(math.exp(-p["protection_strength"]*sum(material[i]==3 for i in neighbors[j]))
                      for j in cells)
        return dt*len(cells),dt*c,dt*protected

    for event in run["ledger"]:
        dr,dc,dp=exposures(event["time"]-last)
        r1time+=dr;c1time+=dc;protected_r1time+=dp
        last=event["time"]
        kind=event["kind"]
        if kind=="environment_gate":
            if event["i"]==0 and event["after"]==1 and first_gate0_bad is None:
                first_gate0_bad={"time":last,"fuel":sum(fuel),
                                 "R0":material.count(1),"R1":material.count(2),
                                 "L":material.count(3),"paid_R_before":paid_R,
                                 "paid_L_before":paid_L}
            gates=event["gates_after"]
            continue
        i,j=event["i"],event["j"]
        n0,n1=material.count(1),material.count(2)
        if n0 and n1:both_seen=True
        if both_seen and ((n0==0 and n1>0) or (n1==0 and n0>0)):
            lost_with_partner=True
        if kind=="copy_R" and event["before"]==0:
            paid_R+=1
            if first_gate0_bad is not None:after_gate0_R+=1
            if event["after"]==2 and first_paid_variant is None:
                first_paid_variant={"time":last,"fuel":sum(fuel),
                                    "R":n0+n1,"L":material.count(3)}
            if material[j]==2:
                if gates==[1,0]:
                    r1_copies_gate10+=1
                    r1_to_r0_gate10+=event["after"]==1
                if event["after"]==1 and n0==0 and both_seen:
                    r0_reintroduced_by_r1+=1
        if kind=="produce_L" and event["before"]==0:
            paid_L+=1
            if first_gate0_bad is not None:after_gate0_L+=1
            if material[j]==2 and gates==[1,0]:r1_L_gate10+=1
        if kind in ("passive","copy_R","produce_L"):
            material[i]=event["after"]
            if kind in ("copy_R","produce_L"):
                fuel[i]=0 if event["before"]==0 else 1
            if n0+n1 and material.count(1)+material.count(2)==0 and first_extinction_fuel is None:
                first_extinction_fuel=sum(fuel)
        elif kind=="material_swap":
            material[i],material[j]=material[j],material[i]
        elif kind=="fuel_swap":
            fuel[i],fuel[j]=fuel[j],fuel[i]
    dr,dc,dp=exposures(run["horizon"]-last)
    r1time+=dr;c1time+=dc;protected_r1time+=dp
    assert material==run["final_state"]["base_state"]["material"]
    assert fuel==run["final_state"]["base_state"]["fuel"]
    return {
        "params":p,"first_gate0_bad":first_gate0_bad,
        "first_paid_variant":first_paid_variant,
        "fuel_at_first_R_extinction":first_extinction_fuel,
        "R1_time_gate10":r1time,"R1_fueled_D_neighbor_time_gate10":c1time,
        "R1_protected_death_exposure_gate10":protected_r1time,
        "R1_copies_gate10":r1_copies_gate10,"those_creating_R0":r1_to_r0_gate10,
        "R1_L_births_gate10":r1_L_gate10,
        "R0_reintroduced_by_R1":r0_reintroduced_by_r1,
        "lost_class_with_partner_alive":lost_with_partner,
        "paid_R_after_first_gate0_bad":after_gate0_R,
        "paid_L_after_first_gate0_bad":after_gate0_L,
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--summary",type=Path,
                        default=Path(__file__).with_name("response-ensemble-results.json"))
    args=parser.parse_args()
    summary=json.loads(args.summary.read_text())
    root=args.summary.parent
    output={"summary":str(args.summary),
            "scope":"Read-only ledger reconstruction and conditional state-distribution diagnostics; no counterfactual model run.",
            "groups":[]}
    for mode,varied in [("common",True),("differentiated",True),("differentiated",False)]:
        selected=[row for row in summary["runs"] if row["mode"]==mode and row["varied"]==varied]
        records=[replay(root/row["trajectory_file"]) for row in selected]
        first=[z["first_gate0_bad"] for z in records if z["first_gate0_bad"] is not None]
        alive=[z for z in first if z["R0"]+z["R1"]>0]
        variants=[z["first_paid_variant"] for z in records if z["first_paid_variant"] is not None]
        group={
            "mode":mode,"varied":varied,"runs":len(records),
            "actual_mutation_matrix":records[0]["params"]["mutation"],
            "first_gate0_bad_observed":len(first),"R_alive_then":len(alive),
            "both_labels_alive_then":sum(z["R0"]>0 and z["R1"]>0 for z in alive),
            "first_gate0_bad_mean_snapshot":{k:mean(z[k] for z in first) for k in
                ["time","fuel","R0","R1","L","paid_R_before","paid_L_before"]},
            "first_gate0_bad_R_alive_mean_snapshot":{k:mean(z[k] for z in alive) for k in
                ["time","fuel","R0","R1","L"]},
            "first_paid_variant_mean_snapshot":{k:mean(z[k] for z in variants) for k in
                ["time","fuel","R","L"]},
            "mean_fuel_remaining_at_first_R_extinction":mean(z["fuel_at_first_R_extinction"] for z in records
                if z["fuel_at_first_R_extinction"] is not None),
            "mean_paid_R_after_first_gate0_bad":mean(z["paid_R_after_first_gate0_bad"] for z in records),
            "mean_paid_L_after_first_gate0_bad":mean(z["paid_L_after_first_gate0_bad"] for z in records),
            "R1_copies_gate10":sum(z["R1_copies_gate10"] for z in records),
            "those_creating_R0":sum(z["those_creating_R0"] for z in records),
            "R1_L_births_gate10":sum(z["R1_L_births_gate10"] for z in records),
            "R0_reintroduced_by_R1":sum(z["R0_reintroduced_by_R1"] for z in records),
            "runs_with_class_loss_while_partner_alive":sum(z["lost_class_with_partner_alive"] for z in records),
            "observed_reverse_R_copies":sum(z["reverse_R_copies"] for z in selected),
            "observed_reverse_L_production":sum(z["reverse_L_production"] for z in selected),
        }
        if mode=="differentiated" and varied:
            p=records[0]["params"]
            tr=sum(z["R1_time_gate10"] for z in records)
            tc=sum(z["R1_fueled_D_neighbor_time_gate10"] for z in records)
            td=sum(z["R1_protected_death_exposure_gate10"] for z in records)
            c,a=tc/tr,td/tr
            self_mu=p["mutation"][1][1]
            death=p["passive_R"]*p["allele_passive_multipliers"][1]*a
            group["healthy_R1_observed_context_diagnostic"]={
                "R1_time":tr,"fueled_D_neighbor_time":tc,"protected_death_exposure":td,
                "mean_fueled_D_neighbors":c,"mean_protection_factor":a,
                "self_copy_hazard_per_R1":p["copy_rate"]*self_mu*c,
                "death_hazard_per_R1":death,
                "copy_rate_threshold_if_context_were_frozen":death/(self_mu*c),
                "self_copy_minus_death_at_copy_rate_012_frozen_context":.12*self_mu*c-death,
                "scope":"Conditional exposure-weighted averages, not an exact invasion eigenvalue; altered kinetics must regenerate all contexts.",
            }
        output["groups"].append(group)
    print(json.dumps(output,indent=2))


if __name__=="__main__":main()

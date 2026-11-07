"""Exact finite-resource copying with physical errors and reversible shocks.

No diversity score, external mutation source, fixed restored mixture, or survival reward
enters the generator. Every erroneous copy consumes monomer and one fuel unit.
The finite grammar of response alleles and collective shock kernel are supplied.
Copying and passive decay use mass action; shocks are a reversible aggregate
kernel, not an elementary fixed-order dilute chemical reaction.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict, replace
from itertools import product
from math import comb, lgamma
from pathlib import Path
import json

import numpy as np
from scipy.integrate import simpson
from scipy.sparse import coo_matrix, csr_matrix, bmat
from scipy.sparse.linalg import expm_multiply


@dataclass(frozen=True)
class Parameters:
    material: int = 6
    fuel_inventory: int = 12
    alleles: int = 2
    mutation: float = .2
    copy_rate: float = .04
    decay_rate: float = .03
    shock_rate: float = .2
    shock_mode: str = "allele"
    template_energy: float = .8
    fuel_energy: float = 1.2
    thermal: float = .05
    horizon: float = 8.


EVENTS = ["forward_copies", "reverse_copies", "erroneous_forward_copies",
          "passive_destroyed_templates", "passive_thermal_templates",
          "shock_destroyed_templates", "shock_thermal_templates"]


def compositions(total: int, dimensions: int):
    return [r for r in product(range(total+1),repeat=dimensions) if sum(r)<=total]


def kernel(p: Parameters):
    assert p.material>=1 and p.fuel_inventory>=0 and p.alleles>=2
    assert 0<=p.mutation<=1 and p.copy_rate>=0 and p.decay_rate>=0
    assert p.shock_rate>=0 and p.shock_mode in ("allele","common","none")
    assert p.thermal>0 and p.template_energy>0 and p.fuel_energy>p.template_energy
    materials=compositions(p.material,p.alleles)
    states=[(r,f) for r in materials for f in range(p.fuel_inventory+1)]
    index={s:i for i,s in enumerate(states)}
    mi={r:i for i,r in enumerate(materials)}
    size=len(states)
    rows=[];cols=[];values=[]
    counters=np.zeros((size,len(EVENTS)))
    first_rates=np.zeros((size,3))  # first new allele via copying, passive, shock
    mu=np.full((p.alleles,p.alleles),p.mutation/(p.alleles-1))
    np.fill_diagonal(mu,1-p.mutation)
    beta=1/p.thermal
    back_copy=np.exp(-(p.fuel_energy-p.template_energy)*beta)
    back_passive=np.exp(-p.template_energy*beta)

    def edge(source, rnew, fnew, rate, event, amount=1.,first_kind=None):
        if rate==0:return
        assert rate>0 and 0<=fnew<=p.fuel_inventory
        assert all(x>=0 for x in rnew) and sum(rnew)<=p.material
        target=index[(tuple(rnew),fnew)]
        rows.append(source);cols.append(target);values.append(float(rate))
        counters[source,event]+=rate*amount
        old=states[source][0]
        if first_kind is not None and sum(old[1:])==0 and sum(rnew[1:])>0:
            first_rates[source,first_kind]+=rate

    for s,(r,f) in enumerate(states):
        d=p.material-sum(r);w=p.fuel_inventory-f
        # Copy-template i is a catalyst; product j may equal or differ from i.
        for i in range(p.alleles):
            for j in range(p.alleles):
                if d and f and r[i] and mu[i,j]:
                    new=list(r);new[j]+=1
                    rate=p.copy_rate*mu[i,j]*d*f*r[i]
                    edge(s,new,f-1,rate,0,first_kind=0)
                    if i!=j:counters[s,2]+=rate
                available=r[j]-(1 if i==j else 0)
                if w and r[i] and available>0 and mu[i,j]:
                    new=list(r);new[j]-=1
                    rate=p.copy_rate*mu[i,j]*back_copy*w*r[i]*available
                    edge(s,new,f+1,rate,1)
        for i in range(p.alleles):
            if r[i]:
                new=list(r);new[i]-=1
                edge(s,new,f,p.decay_rate*r[i],3)
            if d:
                new=list(r);new[i]+=1
                edge(s,new,f,p.decay_rate*back_passive*d,4,first_kind=1)
        if p.shock_mode=="allele":
            for i in range(p.alleles):
                if r[i]:
                    k=r[i];new=list(r);new[i]=0
                    edge(s,new,f,p.shock_rate,5,amount=k)
                else:
                    for k in range(1,d+1):
                        new=list(r);new[i]=k
                        rate=p.shock_rate*np.exp(-k*p.template_energy*beta)*comb(d,k)
                        edge(s,new,f,rate,6,amount=k,first_kind=2)
        elif p.shock_mode=="common":
            if sum(r):
                edge(s,[0]*p.alleles,f,p.shock_rate,5,amount=sum(r))
            else:
                for new in materials:
                    k=sum(new)
                    if not k:continue
                    multiplicity=np.exp(lgamma(p.material+1)-lgamma(p.material-k+1)
                                        -sum(lgamma(x+1) for x in new))
                    rate=p.shock_rate*np.exp(-k*p.template_energy*beta)*multiplicity
                    edge(s,new,f,rate,6,amount=k,first_kind=2)
    Q=coo_matrix((values,(rows,cols)),shape=(size,size)).tocsr()
    outgoing=np.asarray(Q.sum(axis=1)).ravel()
    Q=Q-csr_matrix((outgoing,(np.arange(size),np.arange(size))),shape=(size,size))
    # Exact Gibbs distribution includes indistinguishable-particle factorials.
    log_pi_material=[]
    for r in materials:
        d=p.material-sum(r)
        log_pi_material.append(lgamma(p.material+1)-lgamma(d+1)
            -sum(lgamma(x+1) for x in r)-beta*p.template_energy*sum(r)
            -p.material*np.log1p(p.alleles*back_passive))
    log_pi_fuel=np.array([lgamma(p.fuel_inventory+1)-lgamma(f+1)
        -lgamma(p.fuel_inventory-f+1)-beta*p.fuel_energy*f
        -p.fuel_inventory*np.log1p(np.exp(-beta*p.fuel_energy))
        for f in range(p.fuel_inventory+1)])
    log_pi_material=np.array(log_pi_material)
    log_pi=(log_pi_material[:,None]+log_pi_fuel[None,:]).ravel()
    coo=Q.tocoo();positive=coo.row!=coo.col
    a=coo.row[positive];b=coo.col[positive];rates=coo.data[positive]
    reverse=np.asarray(Q[b,a]).ravel()
    log_error=np.max(np.abs(log_pi[a]+np.log(rates)-log_pi[b]-np.log(reverse)))
    assert log_error<5e-12
    assert np.max(np.abs(np.asarray(Q.sum(axis=1))))<1e-12
    assert abs(np.exp(log_pi).sum()-1)<1e-12
    return Q,states,materials,index,mi,counters,first_rates,log_pi,log_pi_material,log_pi_fuel,float(log_error)


def kl(probability: np.ndarray, log_reference: np.ndarray):
    valid=probability>0
    return float(np.sum(probability[valid]*(np.log(probability[valid])-log_reference[valid])))


def run(p: Parameters, grid: int=401):
    Q,states,materials,index,mi,counters,first_rates,lp,lpm,lpf,db_error=kernel(p)
    initial=np.zeros(len(states));ancestor=(1,)+(0,)*(p.alleles-1)
    initial[index[(ancestor,p.fuel_inventory)]]=1
    times=np.linspace(0,p.horizon,grid)
    probs=expm_multiply(Q.T,initial,start=0,stop=p.horizon,num=grid)
    assert probs.min()>-1e-12
    assert np.max(np.abs(probs.sum(axis=1)-1))<1e-11
    active=np.array([sum(r) for r,f in states],float)
    fuel=np.array([f for r,f in states],float)
    variants=np.array([sum(x>0 for x in r) for r,f in states],float)
    alive=(active>0).astype(float)
    novel=np.array([sum(r[1:])>0 for r,f in states],float)
    expected_active=np.sum(probs*active,axis=1)
    expected_alive=np.sum(probs*alive,axis=1)
    expected_novel=np.sum(probs*novel,axis=1)
    expected_variants=np.sum(probs*variants,axis=1)
    material_probs=probs.reshape(grid,len(materials),p.fuel_inventory+1).sum(axis=2)
    fuel_probs=probs.reshape(grid,len(materials),p.fuel_inventory+1).sum(axis=1)
    b_material=np.array([p.thermal*kl(row,lpm) for row in material_probs])
    b_fuel=np.array([p.thermal*kl(row,lpf) for row in fuel_probs])
    b_full=np.array([p.thermal*kl(row,lp) for row in probs])
    assert np.max(np.diff(b_full))<2e-10
    assert min(b_material.min(),b_fuel.min(),b_full.min())>-1e-11
    mutual=b_full-b_material-b_fuel
    assert mutual.min()>-1e-10
    # Exact integrated event counts through a block matrix exponential.
    augmented=bmat([[Q.T,None],[csr_matrix(counters.T),csr_matrix((len(EVENTS),len(EVENTS)))]],format="csr")
    outcome=expm_multiply(augmented*p.horizon,np.r_[initial,np.zeros(len(EVENTS))])
    event_counts=dict(zip(EVENTS,map(float,outcome[-len(EVENTS):]),strict=True))
    fuel_used=p.fuel_inventory-float(np.sum(probs[-1]*fuel))
    assert abs(fuel_used-(event_counts["forward_copies"]-event_counts["reverse_copies"]))<1e-10
    final_mean_energy=p.template_energy*expected_active[-1]+p.fuel_energy*(p.fuel_inventory-fuel_used)
    initial_mean_energy=p.template_energy+p.fuel_energy*p.fuel_inventory
    energy_release=(p.fuel_energy-p.template_energy)*fuel_used+p.template_energy*(
        event_counts["passive_destroyed_templates"]-event_counts["passive_thermal_templates"]
        +event_counts["shock_destroyed_templates"]-event_counts["shock_thermal_templates"])
    energy_balance_error=abs(final_mean_energy-initial_mean_energy+energy_release)
    assert energy_balance_error<1e-10
    # First appearance source attribution, not merely final variant occupancy.
    ancestral=np.array([i for i,(r,f) in enumerate(states) if sum(r[1:])==0])
    first_Q=bmat([[Q[ancestral][:,ancestral],csr_matrix(first_rates[ancestral])],
                  [None,csr_matrix((3,3))]],format="csr")
    assert np.max(np.abs(np.asarray(first_Q.sum(axis=1))))<1e-12
    first_initial=np.r_[initial[ancestral],np.zeros(3)]
    first=expm_multiply(first_Q.T*p.horizon,first_initial)
    first_prob=dict(zip(["paid_copy","passive_thermal","shock_thermal"],map(float,first[-3:]),strict=True))
    assert sum(first_prob.values())<=1+1e-10
    summary={"parameters":asdict(p),"states":len(states),"max_detailed_balance_log_error":db_error,
        "initial_material_exergy":float(b_material[0]),"initial_total_exergy":float(b_full[0]),
        "final_material_exergy":float(b_material[-1]),"integrated_material_exergy":float(simpson(b_material,x=times)),
        "final_fuel_exergy":float(b_fuel[-1]),"final_total_exergy":float(b_full[-1]),
        "final_material_fuel_information_exergy":float(mutual[-1]),
        "final_expected_templates":float(expected_active[-1]),
        "integrated_expected_templates":float(simpson(expected_active,x=times)),
        "final_survival_probability":float(expected_alive[-1]),
        "final_novel_allele_probability":float(expected_novel[-1]),
        "final_expected_realized_alleles":float(expected_variants[-1]),
        "fuel_used":fuel_used,"event_counts":event_counts,
        "initial_mean_internal_energy":float(initial_mean_energy),
        "final_mean_internal_energy":float(final_mean_energy),
        "mean_energy_event_ledger_error":float(energy_balance_error),
        "first_novel_allele_source_probability":first_prob,
        "maximum_total_exergy_increase_on_grid":float(np.max(np.diff(b_full))),
        "trajectory":{"time":times[::5].tolist(),"material_exergy":b_material[::5].tolist(),
            "total_exergy":b_full[::5].tolist(),"expected_templates":expected_active[::5].tolist(),
            "survival_probability":expected_alive[::5].tolist()}}
    return summary


def main():
    base=Parameters()
    specs={"no_mutation":replace(base,mutation=0),"copying_errors":base,
        "no_shock_no_mutation":replace(base,mutation=0,shock_mode="none"),
        "no_shock_copying_errors":replace(base,shock_mode="none"),
        "common_shock_no_mutation":replace(base,mutation=0,shock_mode="common"),
        "common_shock_copying_errors":replace(base,shock_mode="common"),
        "one_site_no_mutation":replace(base,mutation=0,material=1),
        "one_site_copying_errors":replace(base,material=1),
        "no_fuel_no_mutation":replace(base,mutation=0,fuel_inventory=0),
        "no_fuel_copying_errors":replace(base,fuel_inventory=0),
        "three_alleles_no_mutation":replace(base,mutation=0,alleles=3),
        "three_alleles_copying_errors":replace(base,alleles=3),
        "four_alleles_no_mutation":replace(base,mutation=0,alleles=4),
        "four_alleles_copying_errors":replace(base,alleles=4)}
    cases={}
    for name,p in specs.items():
        r=run(p);cases[name]=r
        print(name,{k:r[k] for k in ["final_material_exergy","integrated_material_exergy",
              "final_expected_templates","final_survival_probability"]},flush=True)
    assert np.isclose(cases["one_site_no_mutation"]["final_material_exergy"],
                      cases["one_site_copying_errors"]["final_material_exergy"],rtol=1e-10)
    assert np.isclose(cases["no_fuel_no_mutation"]["final_material_exergy"],
                      cases["no_fuel_copying_errors"]["final_material_exergy"],rtol=1e-10)
    r0=cases["no_mutation"];r1=cases["copying_errors"]
    assert r1["final_material_exergy"]>r0["final_material_exergy"]
    assert r1["first_novel_allele_source_probability"]["paid_copy"]>.5
    assert sum(v for k,v in r1["first_novel_allele_source_probability"].items() if k!="paid_copy")<1e-5
    assert cases["common_shock_copying_errors"]["final_material_exergy"]<cases["common_shock_no_mutation"]["final_material_exergy"]
    assert cases["no_shock_copying_errors"]["final_material_exergy"]<cases["no_shock_no_mutation"]["final_material_exergy"]
    refinement={}
    for name,p in (("no_mutation",specs["no_mutation"]),("copying_errors",base)):
        refined=run(p,grid=801)
        difference=abs(refined["integrated_material_exergy"]-cases[name]["integrated_material_exergy"])
        assert difference<4e-6
        refinement[name]={"coarse_grid":401,"fine_grid":801,"integral_difference":difference}
    result={"status":"Exact finite CTMC; supplied allele grammar; every erroneous copy consumes material and fuel",
        "cases":cases,"integration_refinement":refinement,
        "comparison":{"final_material_exergy_ratio":r1["final_material_exergy"]/r0["final_material_exergy"],
            "integrated_material_exergy_ratio":r1["integrated_material_exergy"]/r0["integrated_material_exergy"],
            "survival_probability_ratio":r1["final_survival_probability"]/r0["final_survival_probability"]},
        "interpretation_limits":[
            "Allele repertoire and the class-specific shared catastrophe response are supplied finite chemistry.",
            "Bulk all-of-type clearing is a collective reversible Markov channel, not elementary fixed-order mass action.",
            "The reverse shock rate includes the factorial degeneracy required by local detailed balance.",
            "Rare thermal reverse channels exist and first-appearance attribution reports them separately.",
            "Material-subsystem exergy may rise through fuel conversion while total closed-system exergy decreases.",
            "No external reseeding, restored mixture, diversity reward, or survival objective drives the generator.",
            "This is paid copying with errors, not open-ended sequence invention or demonstrated collective selfhood.",
        ],"checks_passed":True}
    out=Path(__file__).with_name("agent-mutation-results.json")
    out.write_text(json.dumps(result,indent=2)+"\n")
    print("saved",out)


if __name__=="__main__":main()

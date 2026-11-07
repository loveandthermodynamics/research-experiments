"""Manual adversarial readout histories, not physical simulations or rate tests."""
from __future__ import annotations

from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path

from renewal_readouts import measure

KEYS=("material","templates","material_tags","instance_tags")


class Fixture:
    def __init__(self, material=None, templates=None):
        material=list(material) if material is not None else [1]*16
        templates=list(templates) if templates is not None else [1]+[0]*15
        mt,it=[],[]
        mc=ic=0
        for a,t in zip(material,templates):
            if a:mc+=1
            if t:ic+=1
            mt.append(mc if a else 0);it.append(ic if t else 0)
        self.state=dict(material=material,templates=templates,material_tags=mt,
                        instance_tags=it,next_tag=ic+1,time=0.)
        self.initial=deepcopy(self.state)
        self.events=[]
        self.samples=[self._sample(0.)]

    def _sample(self,time):
        out={key:self.state[key].copy() for key in KEYS}
        out.update(time=time,B=out["material"].count(2),X=out["templates"].count(1),
                   Y=out["templates"].count(2),D=out["templates"].count(0))
        return out

    def sample(self,time):
        self.samples.append(self._sample(time))

    def event(self,time,kind,i,j=-1,activity=0,child=0,direction=0,species=0):
        sites=[i] if j==-1 else [i,j]
        before={key:[self.state[key][site] for site in sites] for key in KEYS}
        created=removed=parent=catalyst=helper=0
        origin=None
        if kind in ("material_swap","template_swap"):
            for key in (("material","material_tags") if kind=="material_swap" else ("templates","instance_tags")):
                self.state[key][i],self.state[key][j]=self.state[key][j],self.state[key][i]
        elif kind=="phase_passive":
            self.state["material"][i]=3-self.state["material"][i]
        elif kind=="template_passive":
            if species:
                assert self.state["templates"][i]==0
                created=self.state["next_tag"];self.state["next_tag"]+=1
                self.state["templates"][i]=species;self.state["instance_tags"][i]=created
                origin="thermal"
            else:
                assert self.state["templates"][i]!=0
                removed=self.state["instance_tags"][i]
                self.state["templates"][i]=self.state["instance_tags"][i]=0
        elif kind=="copy":
            helper=self.state["templates"][j];catalyst=self.state["instance_tags"][j]
            assert helper in (1,2)
            source,target=(1,2) if activity==1 else (2,1)
            assert self.state["material"][i]==(source if direction==1 else target)
            if direction==1:
                assert self.state["templates"][i]==0
                created=self.state["next_tag"];self.state["next_tag"]+=1
                self.state["templates"][i]=child;self.state["instance_tags"][i]=created
                parent=catalyst;origin="paid_copy"
            else:
                assert self.state["templates"][i]==child
                removed=self.state["instance_tags"][i]
                self.state["templates"][i]=self.state["instance_tags"][i]=0
            self.state["material"][i]=target if direction==1 else source
        else:
            raise AssertionError(kind)
        assert time>=self.state["time"]
        self.state["time"]=time
        after={key:[self.state[key][site] for site in sites] for key in KEYS}
        self.events.append(dict(time=time,kind=kind,i=i,j=j,activity=activity,
            direction=direction,offspring_species=child,helper_species=helper,
            catalyst_instance=catalyst,parent_instance=parent,instance_created=created,
            instance_removed=removed,birth_origin=origin,before=before,after=after))

    def finish(self,horizon):
        assert horizon>=self.state["time"]
        self.state["time"]=horizon
        self.samples.append(self._sample(horizon))
        return dict(initial_state=deepcopy(self.initial),final_state=deepcopy(self.state),
                    horizon=horizon,ledger=deepcopy(self.events),samples=deepcopy(self.samples))


def independent_carrier_check(record,result):
    """Separate full-array interval replay and per-carrier post-state traces.

    Unlike the reader's local before/after transitions, this uses only initial
    arrays and successive POST arrays to reconstruct carrier-state histories.
    It is an independently implemented algorithm, not an independent author.
    """
    arrays={key:record["initial_state"][key].copy() for key in KEYS}
    traces={tag:[(0.,mat)] for tag,mat in zip(arrays["material_tags"],arrays["material"]) if tag}
    last={tag:mat for tag,mat in zip(arrays["material_tags"],arrays["material"]) if tag}
    integral=previous=0.
    for e in record["ledger"]:
        integral+=(e["time"]-previous)*arrays["material"].count(2);previous=e["time"]
        sites=[e["i"]] if e["j"]==-1 else [e["i"],e["j"]]
        for key in KEYS:
            for site,value in zip(sites,e["after"][key]):arrays[key][site]=value
        for tag,mat in zip(e["after"]["material_tags"],e["after"]["material"]):
            if tag and last[tag]!=mat:
                traces[tag].append((e["time"],mat));last[tag]=mat
    H=record["horizon"]
    integral+=(H-previous)*arrays["material"].count(2)
    expected=[]
    for tag,trace in traces.items():
        for index,(time,state) in enumerate(trace):
            if state==2:
                end=trace[index+1][0] if index+1<len(trace) else H
                expected.append((tag,time,end,index==0,index+1==len(trace)))
    actual=[(e["carrier"],e["start"],e["end"],e["left_censored"],e["right_censored"])
            for e in result["B_episode_records"]]
    assert sorted(expected)==sorted(actual)
    assert math.isclose(integral,result["integral_B"],abs_tol=1e-12)
    return True


def histories():
    cases=[]
    # Carrier 1 moves in B before its X-reverse exit; Y-reverse starts carrier2.
    mat=[2,1,1]+[0]*13;tpl=[2,1,0,1,2]+[0]*11
    f=Fixture(mat,tpl)
    f.event(1.,"material_swap",0,1)
    f.event(2.,"copy",1,3,activity=1,child=1,direction=-1)
    f.event(3.,"copy",0,4,activity=2,child=2,direction=-1)
    f.event(4.,"material_swap",0,2)
    f.event(5.,"phase_passive",2)
    f.event(6.,"copy",0,3,activity=1,child=1,direction=1)
    cases.append(("carrier_moves_and_both_reverse_transformations",f.finish(8.),
                  {"integral_B":6.,"final_B":1,"net_X_activity_copy_count":0,
                   "net_Y_activity_copy_count":-1,"net_background_A_to_B_count":-1,
                   "B_episode_summary.episodes":3,"B_episode_summary.left_censored":1,
                   "B_episode_summary.right_censored":1,"B_episode_summary.fully_observed_episodes":1}))

    # Initial founder moves independently of carriers. Earlier thermal-X copied
    # Y is excluded; focal later backmutates and retains global ancestry.
    mat=[1]*16;mat[4]=mat[6]=2
    f=Fixture(mat)
    f.event(.25,"template_swap",0,10)
    f.event(.5,"template_passive",1,species=1)
    f.event(1.,"copy",2,1,activity=1,child=2,direction=1)
    f.event(2.,"copy",3,10,activity=1,child=2,direction=1)
    f.event(3.,"copy",4,3,activity=2,child=1,direction=1)
    f.event(4.,"copy",5,4,activity=1,child=2,direction=1)
    f.event(5.,"copy",6,3,activity=2,child=2,direction=1)
    f.event(6.,"copy",3,1,activity=1,child=2,direction=-1)
    f.event(7.,"template_passive",5)
    f.event(8.,"template_passive",6)
    f.event(9.,"copy",7,4,activity=1,child=2,direction=1)
    f.event(9.5,"template_passive",10)
    f.event(9.7,"copy",8,4,activity=1,child=1,direction=1)
    cases.append(("paid_ancestry_backmutation_unrelated_reverse_and_founder_clearance",f.finish(10.),
                  {"focal_tag":4,"focal_appearance_time":2.,"maximum_live_focal_subtree":4,
                   "focal_established":True,"focal_establishment_time":5.,"focal_survives_to_horizon":True,
                   "final_focal_subtree":3,"final_focal_Y":1,"focal_Y_first_loss_time":8.,
                   "focal_Y_restorations_after_loss":1,"original_X_founder_clearance_time":9.5,
                   "forward_paid_copies_by_initial_X_descendants_after_clearance":1,
                   "integral_original_X_founders":9.5,"integral_initial_X_lineage":27.8,
                   "integral_focal_subtree":18.3,"integral_focal_Y":11.,
                   "integral_initial_X_lineage_after_original_X_clearance":1.3,
                   "integral_templates_after_original_X_clearance":2.3,
                   "first_new_Y_beyond_initial.tag":3,
                   "first_new_Y_beyond_initial.parent_initial_X_ancestry":False}))

    f=Fixture()
    f.event(1.,"copy",1,0,activity=1,child=2,direction=1)
    f.event(2.,"template_passive",1)
    f.event(3.,"template_passive",1,species=2)
    f.event(4.,"copy",1,0,activity=1,child=2,direction=-1)
    f.event(5.,"copy",1,0,activity=1,child=2,direction=1)
    cases.append(("thermal_and_later_paid_arrivals_do_not_resurrect_focal",f.finish(6.),
                  {"focal_tag":2,"focal_extinction_time":2.,"focal_survives_to_horizon":False,
                   "final_focal_subtree":0,"final_Y":1,"maximum_live_focal_subtree":1,
                   "focal_Y_restorations_after_loss":0,"integral_B":4.}))

    mat=[1]*16
    for site in (2,3,4,5):mat[site]=2
    f=Fixture(mat)
    f.event(1.,"copy",1,0,activity=1,child=2,direction=1)
    for time,site,child in [(2.,2,1),(4.,3,2),(6.,4,1),(8.,5,2)]:
        f.event(time,"copy",site,1,activity=2,child=child,direction=1)
        f.event(time+1.,"template_passive",site)
    cases.append(("cumulative_births_are_not_peak_live_subtree",f.finish(10.),
                  {"births.paid_copy":5,"maximum_live_focal_subtree":2,
                   "focal_established":False,"final_focal_subtree":1}))

    f=Fixture(templates=[0]*16)
    f.event(1.,"template_passive",0,species=1)
    f.event(2.,"copy",1,0,activity=1,child=2,direction=1)
    cases.append(("empty_initial_X_cohort_not_a_clearance_event",f.finish(3.),
                  {"initial_X_founder_cohort_present":False,"original_X_founders_cleared_after_presence":False,
                   "original_X_founder_clearance_time":None,
                   "focal_paid_Y_from_initial_X_appeared_by_deadline":False,
                   "integral_templates":3.,"integral_templates_after_original_X_clearance":0.,
                   "forward_paid_copies_after_original_X_clearance":0}))

    mat=[1]*16;mat[1]=2
    f=Fixture(mat,[2]+[0]*15)
    f.event(1.,"copy",1,0,activity=2,child=1,direction=1)
    f.event(2.,"copy",2,1,activity=1,child=2,direction=1)
    cases.append(("initial_Y_ancestor_is_not_initial_X_ancestry",f.finish(3.),
                  {"initial_X_founders":0,"focal_paid_Y_from_initial_X_appeared_by_deadline":False,
                   "final_initial_X_lineage":0}))

    for time,appeared in [(20.,True),(20.001,False)]:
        mat=[1]*16;mat[1]=2
        f=Fixture(mat)
        f.event(time,"copy",1,0,activity=2,child=2,direction=1)
        cases.append((f"deadline_{time}_generalist_X_parent_using_Y_activity",f.finish(25.),
                      {"focal_paid_Y_from_initial_X_appeared_by_deadline":appeared,
                       "maximum_live_focal_subtree":1 if appeared else 0,
                       "net_Y_activity_copy_count":1}))

    f=Fixture([2]+[1]*15)
    cases.append(("zero_horizon_left_and_right_censored_spell",f.finish(0.),
                  {"integral_B":0.,"B_episode_summary.episodes":1,
                   "B_episode_summary.left_censored":1,"B_episode_summary.right_censored":1,
                   "B_episode_summary.truncated_duration_sum":0.,
                   "original_X_founders_cleared_after_presence":False}))

    f=Fixture()
    f.sample(1.) # Exact-jump sample is the pre-event state.
    f.event(1.,"copy",1,0,activity=1,child=2,direction=1)
    cases.append(("sample_at_jump_is_pre_event",f.finish(3.),{"integral_B":2.,"integral_templates":5.}))
    f=Fixture()
    f.event(1.,"copy",1,0,activity=1,child=2,direction=1)
    record=f.finish(1.)
    record["samples"][-1]={key:record["initial_state"][key].copy() for key in KEYS}
    record["samples"][-1].update(time=1.,B=0,X=1,Y=0,D=15)
    cases.append(("event_at_horizon_has_post_final_state_but_zero_new_residence",record,
                  {"integral_B":0.,"B_episode_summary.episodes":1,
                   "B_episode_summary.truncated_duration_sum":0.,"B_episode_summary.right_censored":1,
                   "focal_survives_to_horizon":True,"final_Y":1}))
    return cases


def nested(data,path):
    for part in path.split("."):data=data[part]
    return data


def main():
    reports=[]
    cases=histories()
    for name,record,expected in cases:
        result=measure(record,include_details=True)
        for path,value in expected.items():
            actual=nested(result,path)
            assert math.isclose(actual,value,abs_tol=2e-12) if isinstance(value,float) else actual==value, (name,path,actual,value)
        assert independent_carrier_check(record,result)
        reports.append(dict(name=name,passed=True,manual_expected=expected,measurements=result))
    invalid=[]
    base=cases[2][1]
    corruptions={
        "forged_forward_parent":lambda r:r["ledger"][0].__setitem__("parent_instance",999),
        "reverse_removes_historical_not_current_target":lambda r:r["ledger"][3].__setitem__("instance_removed",2),
        "reverse_is_mislabelled_as_paid_birth":lambda r:r["ledger"][3].__setitem__("birth_origin","paid_copy"),
        "reused_dead_instance_ID":lambda r:(r["ledger"][4].__setitem__("instance_created",2),
                                            r["ledger"][4]["after"]["instance_tags"].__setitem__(0,2)),
    }
    for name,corrupt in corruptions.items():
        record=deepcopy(base);corrupt(record)
        try:measure(record)
        except ValueError as error:invalid.append(dict(name=name,rejected=True,reason=str(error)))
        else:raise AssertionError(f"corrupt history accepted: {name}")
    record=deepcopy(cases[-2][1])
    record["samples"][1].update({key:record["final_state"][key].copy() for key in KEYS})
    try:measure(record)
    except ValueError as error:invalid.append(dict(name="post_event_sample_at_exact_jump",rejected=True,reason=str(error)))
    else:raise AssertionError("post-event sample at exact jump was accepted")
    root=Path(__file__).resolve().parent
    output=dict(passed=True,manual_histories=len(reports),invalid_histories=len(invalid),
                fixtures=reports,rejections=invalid,
                source_sha256={name:hashlib.sha256((root/name).read_bytes()).hexdigest()
                               for name in ("renewal_readouts.py","renewal_fixtures.py")},
                scope="Manually specified analysis fixtures and a separate carrier-trace algorithm; no physical trajectory, rate check, energy claim or ensemble.")
    target=root/"renewal-fixture-results.json"
    target.write_text(json.dumps(output,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print(json.dumps(dict(passed=True,manual_histories=len(reports),invalid_histories=len(invalid),output=str(target))))


if __name__=="__main__":main()

"""Independent scalar energy/rate/instance reference for joined chemistry."""
from __future__ import annotations

from bisect import bisect_right
import copy
from dataclasses import asdict, replace
import itertools
import json
import math
from pathlib import Path
import random
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent/"active-phase"))
from verify_active_phase import Reference as Geometry, close, compare, digest, load
from verify_active_phase import observe_independent as observe_base_independent


class Reference(Geometry):
    def __init__(self, params, flags=None, activity_mode="specialist"):
        super().__init__(params)
        self.flags = flags or dict(enable_X=True, enable_Y=True)
        self.mode = activity_mode
        self.E = (0., self.p["E_X"], self.p["E_Y"])
        self.g = (0., self.p["g_X"], self.p["g_Y"])
        self.events = [(kind,i,j,-1,-1,0,0,0) for kind in ("material_swap","template_swap") for i,j in self.edges]
        self.events += [("template_passive",i,-1,before,after,0,0,0)
                        for i in range(self.n) for child in (1,2) for before,after in ((0,child),(child,0))]
        self.events += [("phase_passive",i,-1,before,after,0,0,0)
                        for i in range(self.n) for before,after in ((1,2),(2,1))]
        for i in range(self.n):
            row,col=divmod(i,self.side)
            neighbors=(row*self.side+(col+1)%self.side, ((row+1)%self.side)*self.side+col,
                       row*self.side+(col-1)%self.side, ((row-1)%self.side)*self.side+col)
            self.events += [("copy",i,j,-1,-1,activity,child,direction)
                            for j in neighbors for activity in (1,2) for child in (1,2) for direction in (1,-1)]
        self.lookup={edge:k for k,edge in enumerate(self.events)}
        self.reverse_index=[k if k<4*self.n else k^1 for k in range(len(self.events))]

    def energy(self,state):
        m,t=state["material"],state["templates"]
        return (self.p["epsilon_B"]*m.count(2)-self.p["J_BB"]*sum(m[i]==m[j]==2 for i,j in self.edges)
                +sum(self.E[r]-self.g[r]*int(a==2) for a,r in zip(m,t)))

    def changes(self,state,event):
        kind,i,j,before,after,activity,child,direction=event
        if kind=="material_swap": return {i:state["material"][j],j:state["material"][i]},{}
        if kind=="template_swap": return {},{i:state["templates"][j],j:state["templates"][i]}
        if kind=="template_passive": return {},{i:after}
        if kind=="phase_passive": return {i:after},{}
        substrate,product=(1,2) if activity==1 else (2,1)
        return {i:product if direction==1 else substrate},{i:child if direction==1 else 0}

    def delta(self,state,event):
        mc,tc=self.changes(state,event)
        m,t=state["material"],state["templates"]
        changed_b=lambda i:int(mc.get(i,m[i])==2)
        b=lambda i:int(m[i]==2)
        value=self.p["epsilon_B"]*sum(changed_b(i)-b(i) for i in mc)
        value+=sum(self.E[r]-self.E[t[i]] for i,r in tc.items())
        affected=set()
        for i in mc: affected.update(self.incident[i])
        value-=self.p["J_BB"]*sum(changed_b(i)*changed_b(j)-b(i)*b(j) for k in affected for i,j in (self.edges[k],))
        value-=sum(self.g[tc.get(i,t[i])]*changed_b(i)-self.g[t[i]]*b(i) for i in set(mc)|set(tc))
        return value

    def work(self,event): return self.p["mu"]*event[7] if event[0]=="copy" else 0.

    def rate(self,state,event):
        kind,i,j,before,after,activity,child,direction=event
        m,t,p=state["material"],state["templates"],self.p
        if kind=="material_swap":
            if m[i]==m[j]: return 0.
            mobility=p["material_mobility"]
        elif kind=="template_swap":
            if t[i]==t[j]: return 0.
            mobility=p["template_mobility"]
        elif kind=="template_passive":
            if t[i]!=before: return 0.
            mobility=p["loss_rate"]
        elif kind=="phase_passive":
            if m[i]!=before: return 0.
            mobility=p["background_AB_rate"]
        else:
            if not self.flags["enable_X" if activity==1 else "enable_Y"]: return 0.
            if t[j]==0 or (self.mode=="specialist" and t[j]!=activity): return 0.
            substrate,product=(1,2) if activity==1 else (2,1)
            if t[i]!=(0 if direction==1 else child) or m[i]!=(substrate if direction==1 else product): return 0.
            mobility=p["copy_rate"]*((1-p["mutation_error"]) if t[j]==child else p["mutation_error"])
        return mobility*math.exp(min(0.,(self.work(event)-self.delta(state,event))/p["theta"])) if mobility else 0.

    def validate(self,state):
        m,t,mt,it=(state[key] for key in ("material","templates","material_tags","instance_tags"))
        assert all(len(a)==self.n for a in (m,t,mt,it))
        assert all(a in (0,1,2) for a in m+t)
        assert sum(a!=0 for a in m)==self.p["material_count"]
        for values,tags in ((m,mt),(t,it)):
            assert all((a==0)==(tag==0) for a,tag in zip(values,tags))
            live=[tag for tag in tags if tag]
            assert len(live)==len(set(live)) and all(tag>0 for tag in live)
        assert state["next_tag"]>max(it,default=0)

    def apply(self,state,index):
        edge=self.events[index]
        assert self.rate(state,edge)>0,("disabled",edge)
        kind,i,j,before,after,activity,child,direction=edge
        sites=[i] if j==-1 else [i,j]
        snapshot=lambda:{key:[state[key][site] for site in sites] for key in ("material","templates","material_tags","instance_tags")}
        old=snapshot()
        delta,work=self.delta(state,edge),self.work(edge)
        helper=state["templates"][j] if kind=="copy" else 0
        catalyst=state["instance_tags"][j] if kind=="copy" else 0
        created=removed=parent=0
        origin=None
        mc,tc=self.changes(state,edge)
        if kind=="material_swap":
            state["material_tags"][i],state["material_tags"][j]=state["material_tags"][j],state["material_tags"][i]
        elif kind=="template_swap":
            state["instance_tags"][i],state["instance_tags"][j]=state["instance_tags"][j],state["instance_tags"][i]
        elif kind in ("template_passive","copy"):
            if tc[i]:
                created=state["next_tag"]
                state["next_tag"]+=1
                state["instance_tags"][i]=created
                origin="paid_copy" if kind=="copy" else "thermal"
                parent=catalyst if kind=="copy" else 0
            else:
                removed=state["instance_tags"][i]
                state["instance_tags"][i]=0
        for site,value in mc.items(): state["material"][site]=value
        for site,value in tc.items(): state["templates"][site]=value
        fuel=direction if kind=="copy" else 0
        state["cumulative_work"]+=work
        state["cumulative_heat"]+=work-delta
        state["net_fuel_to_waste"]+=fuel
        state["event_count"]+=1
        return dict(time=state["time"],event_index=index,kind=kind,i=i,j=j,activity=activity,direction=direction,
            offspring_species=child,helper_species=helper,catalyst_instance=catalyst,parent_instance=parent,
            instance_created=created,instance_removed=removed,birth_origin=origin,before=old,after=snapshot(),
            delta_H=delta,chemical_work=work,heat_to_bath=work-delta,fuel_to_waste=fuel,
            cumulative_work=state["cumulative_work"],cumulative_heat=state["cumulative_heat"],
            net_fuel_to_waste=state["net_fuel_to_waste"],event_count=state["event_count"])

    def summary(self,state,include_state=False):
        m,t=state["material"],state["templates"]
        out=dict(time=state["time"],H=self.energy(state),S=m.count(0),A=m.count(1),B=m.count(2),D=t.count(0),X=t.count(1),Y=t.count(2),
                 X_on_B=sum(a==2 and r==1 for a,r in zip(m,t)),Y_on_B=sum(a==2 and r==2 for a,r in zip(m,t)),
                 BB_contacts=sum(m[i]==m[j]==2 for i,j in self.edges),chemical_work=state["cumulative_work"],
                 heat_to_bath=state["cumulative_heat"],net_fuel_to_waste=state["net_fuel_to_waste"],
                 events=state["event_count"],next_tag=state["next_tag"])
        if include_state: out.update({key:state[key][:] for key in ("material","templates","material_tags","instance_tags")})
        return out


def check_model(seed=33701):
    from productive_model import Model,Params,Flags
    rng=random.Random(seed)
    default=Params(side=4,material_count=5,initial_X=2,initial_Y=1)
    parameters=[default,replace(default,mu=0.),replace(default,background_AB_rate=.0002),
                replace(default,E_X=1.,g_X=2.5,E_Y=3.,g_Y=-.5,theta=.7),replace(default,mutation_error=0.),replace(default,mutation_error=.5)]
    comparisons=pairs=cases=0
    maximum=dict(H=0.,delta_H=0.,rate=0.,LDB_log=0.,effective_equilibrium_log=0.)
    for p in parameters:
        for configuration in range(2):
            carrier=rng.sample(range(16),5)
            material=[rng.choice((1,2)) if i in carrier else 0 for i in range(16)]
            templates=[rng.randrange(3) for _ in range(16)]
            mt=[i+1 if a else 0 for i,a in enumerate(material)]
            it=[i+1 if r else 0 for i,r in enumerate(templates)]
            state=dict(material=material,templates=templates,material_tags=mt,instance_tags=it,next_tag=17,time=0.,
                       cumulative_work=0.,cumulative_heat=0.,net_fuel_to_waste=0,event_count=0)
            for mode,a,b in itertools.product(("specialist","generalist"),(False,True),(False,True)):
                flags=dict(enable_X=a,enable_Y=b)
                model,ref=Model(p,Flags(**flags),mode),Reference(asdict(p),flags,mode)
                converted=model.state_from_dict(state)
                assert [tuple(asdict(event).values()) for event in model.events]==ref.events
                assert list(model.reverse_index)==ref.reverse_index
                close(model.energy(converted),ref.energy(state))
                maximum["H"]=max(maximum["H"],abs(model.energy(converted)-ref.energy(state)))
                _,actual=model.rates(converted)
                for index,edge in enumerate(ref.events):
                    rate=ref.rate(state,edge)
                    close(rate,actual[index],2e-12)
                    maximum["rate"]=max(maximum["rate"],abs(rate-float(actual[index])))
                    comparisons+=1
                    if rate<=0: continue
                    delta=ref.delta(state,edge)
                    close(delta,model.delta_energy(converted,index))
                    maximum["delta_H"]=max(maximum["delta_H"],abs(delta-model.delta_energy(converted,index)))
                    following=copy.deepcopy(state)
                    expected=ref.apply(following,index)
                    production=model.copy_state(converted)
                    observed=model.apply(production,index)
                    compare(expected,observed)
                    compare(following,model.state_to_dict(production))
                    close(ref.energy(following)-ref.energy(state),delta)
                    back=ref.rate(following,ref.events[ref.reverse_index[index]])
                    assert back>0
                    error=abs(math.log(rate/back)-(ref.work(edge)-delta)/p.theta)
                    assert error<2e-11
                    maximum["LDB_log"]=max(maximum["LDB_log"],error)
                    if p.background_AB_rate==0 and (a!=b):
                        sign=1 if a else -1
                        db=following["material"].count(2)-state["material"].count(2)
                        error=abs(math.log(rate/back)+(delta-sign*p.mu*db)/p.theta)
                        assert error<2e-11
                        maximum["effective_equilibrium_log"]=max(maximum["effective_equilibrium_log"],error)
                    ref.validate(following)
                    pairs+=1
                altered=model.copy_state(converted)
                altered.instance_tags*=19; altered.material_tags*=23; altered.next_tag*=19
                assert (model.rates(altered)[1]==actual).all()
                cases+=1
    return dict(cases=cases,scalar_rates_checked=comparisons,enabled_reverse_pairs_checked=pairs,maximum_errors=maximum,
                controls="both activity modes/all enable flags, mu0, positive background, m0/.5 and binding exceeding template energy",
                tags_nonfunctional=True,source_law_condition_not_imposed_outside_its_domain=True)


def check_constructed_instances():
    from productive_model import Model,Params
    p=Params(side=4,material_count=3,initial_X=2,initial_Y=1)
    ref,model=Reference(asdict(p)),Model(p)
    state=dict(material=[0]*16,templates=[0]*16,material_tags=[0]*16,instance_tags=[0]*16,next_tag=4,
               time=0.,cumulative_work=0.,cumulative_heat=0.,net_fuel_to_waste=0,event_count=0)
    for tag,i in enumerate((1,5,6),1): state["material"][i]=1 if i==5 else 2; state["material_tags"][i]=tag
    for tag,i,r in ((1,4,1),(2,6,2),(3,1,1)): state["templates"][i]=r; state["instance_tags"][i]=tag
    initial=copy.deepcopy(state)
    actual=model.state_from_dict(state)
    records=[]
    # A different current X helper removes the child; historical-parent equality
    # must not be an eligibility condition for physical reverse chemistry.
    indices=[ref.lookup[("copy",5,4,-1,-1,1,1,1)],ref.lookup[("copy",5,1,-1,-1,1,1,-1)],
             ref.lookup[("template_passive",5,-1,0,2,0,0,0)],ref.lookup[("template_passive",5,-1,2,0,0,0,0)]]
    for index in indices:
        expected=ref.apply(state,index)
        compare(expected,model.apply(actual,index))
        compare(state,model.state_to_dict(actual))
        records.append(expected)
    assert records[0]["parent_instance"]==1 and records[1]["catalyst_instance"]==3
    assert records[1]["parent_instance"]==0 and records[1]["instance_removed"]==records[0]["instance_created"]
    assert records[2]["birth_origin"]=="thermal" and records[2]["parent_instance"]==0
    assert state["material_tags"]==initial["material_tags"] and state["material"]==initial["material"]
    assert state["templates"]==initial["templates"] and state["instance_tags"]==initial["instance_tags"]
    assert state["next_tag"]==initial["next_tag"]+2
    # Full four-step cycle, with the actual implementation and independent H.
    state=copy.deepcopy(initial); actual=model.state_from_dict(state)
    indices=[ref.lookup[("copy",5,4,-1,-1,1,1,1)],ref.lookup[("template_passive",5,-1,1,0,0,0,0)],
             ref.lookup[("copy",5,6,-1,-1,2,2,1)],ref.lookup[("template_passive",5,-1,2,0,0,0,0)]]
    affinity=0.
    for index in indices:
        rate=ref.rate(state,ref.events[index])
        expected=ref.apply(state,index)
        compare(expected,model.apply(actual,index))
        affinity+=math.log(rate/ref.rate(state,ref.events[ref.reverse_index[index]]))
    close(affinity,24.,2e-12)
    assert state["material"]==initial["material"] and state["templates"]==initial["templates"]
    close(state["cumulative_work"],24.)
    close(ref.energy(state),ref.energy(initial))
    return dict(constructed_events=8,unrelated_reverse_catalyst_allowed=True,thermal_birth_parent_zero=True,
                physical_four_step_cycle_affinity=affinity,observer_next_tag_does_not_return=True,
                scope="Constructed physical transition fixtures, not simulated outcomes or phase evidence.")


def verify_trajectory(run,full_proposals=True):
    import numpy as np
    ref=Reference(run["params"],run["flags"],run["activity_mode"])
    state=copy.deepcopy(run["initial_state"]); ref.validate(state)
    assert all(state[key]==0 for key in ("time","cumulative_work","cumulative_heat","net_fuel_to_waste","event_count"))
    h0=ref.energy(state); direct_h=h0
    b0=state["material"].count(2); n0=sum(r!=0 for r in state["templates"])
    bcount,ncount=b0,n0
    initial_carriers=set(state["material_tags"])-{0}
    cx=cy=c0=births=losses=0
    samples,si=run["samples"],0
    event_counts={}; maximum_error=0.

    def samples_until(until):
        nonlocal si
        while si<len(samples) and samples[si]["time"]<=until:
            expected=ref.summary(state,True); expected["time"]=samples[si]["time"]
            compare(expected,samples[si]); si+=1

    def accept(saved,proposal=None,rejected=None):
        nonlocal direct_h,maximum_error,cx,cy,c0,births,losses,bcount,ncount
        assert state["time"]<=saved["time"]<=run["horizon"]
        samples_until(saved["time"]); state["time"]=saved["time"]
        index=saved["event_index"]
        edge=ref.events[index]; rate=ref.rate(state,edge)
        expected=ref.apply(state,index); compare(expected,saved)
        if proposal is not None:
            assert saved["proposal_number"]==proposal and saved["rejected_since_last_accepted"]==rejected
        back=ref.rate(state,ref.events[ref.reverse_index[index]])
        close(math.log(rate/back),expected["heat_to_bath"]/ref.p["theta"],2e-11)
        if edge[0]=="copy":
            if edge[5]==1: cx+=edge[7]; bcount+=edge[7]
            else: cy+=edge[7]; bcount-=edge[7]
            ncount+=edge[7]
        elif edge[0]=="phase_passive":
            change=1 if edge[4]==2 else -1; c0+=change; bcount+=change
        elif edge[0]=="template_passive":
            if edge[3]==0: births+=1; ncount+=1
            else: losses+=1; ncount-=1
        assert bcount-b0==cx-cy+c0 and ncount-n0==cx+cy+births-losses
        assert state["net_fuel_to_waste"]==cx+cy
        close(state["cumulative_work"],ref.p["mu"]*(cx+cy))
        direct_h+=expected["delta_H"]
        error=abs(direct_h-h0-state["cumulative_work"]+state["cumulative_heat"])
        maximum_error=max(maximum_error,error); assert error<2e-6
        event_counts[edge[0]]=event_counts.get(edge[0],0)+1
        if state["event_count"]%97==0:
            close(ref.energy(state),direct_h)
            assert bcount==state["material"].count(2) and ncount==sum(r!=0 for r in state["templates"])
            ref.validate(state)

    p=ref.p
    counts=[2*ref.n,2*ref.n,4*ref.n,2*ref.n,32*ref.n]
    starts=[0,2*ref.n,4*ref.n,8*ref.n,10*ref.n]
    bounds=[p["material_mobility"],p["template_mobility"],p["loss_rate"],p["background_AB_rate"],p["copy_rate"]*(1-p["mutation_error"])]
    cumulative=list(itertools.accumulate(n*m for n,m in zip(counts,bounds))); total=cumulative[-1]
    stats=run["proposal_statistics"]; close(total,stats["total_bound"])
    if full_proposals:
        rng=np.random.default_rng(np.random.SeedSequence([run["seed"],85603]))
        clock=0.; proposal=accepted=rejected=0; censored=None
        proposed_counts=[0]*5; accepted_counts=[0]*5
        while clock<run["horizon"]:
            next_time=clock+float(rng.exponential(1./total)) if total else math.inf
            if next_time>run["horizon"]:
                censored=next_time if math.isfinite(next_time) else None; break
            clock=next_time; proposal+=1
            family=bisect_right(cumulative,rng.random()*total)
            index=starts[family]+int(rng.integers(counts[family])); proposed_counts[family]+=1
            probability=ref.rate(state,ref.events[index])/bounds[family]
            assert probability<=1.+1e-14
            if rng.random()<probability:
                saved=run["ledger"][accepted]; assert saved["event_index"]==index
                close(saved["time"],clock,1e-12); accept(saved,proposal,rejected)
                accepted_counts[family]+=1; accepted+=1; rejected=0
            else: rejected+=1
        assert accepted==len(run["ledger"])
        families=("material_swap","template_swap","template_passive","phase_passive","copy")
        compare(dict(total_bound=total,proposals=proposal,accepted=accepted,rejected=proposal-accepted,
            terminal_rejected_since_last_accepted=rejected,proposed_by_family=dict(zip(families,proposed_counts)),
            accepted_by_family=dict(zip(families,accepted_counts)),terminal_censored_proposal_time=censored),stats)
    else:
        previous_proposal=0
        for saved in run["ledger"]:
            assert saved["proposal_number"]==previous_proposal+saved["rejected_since_last_accepted"]+1
            previous_proposal=saved["proposal_number"]; accept(saved)
        assert stats["accepted"]==len(run["ledger"])
    samples_until(run["horizon"]); assert si==len(samples)
    state["time"]=run["horizon"]; compare(state,run["final_state"]); compare(ref.summary(state),run["final_summary"])
    close(ref.energy(state),direct_h); ref.validate(state)
    assert initial_carriers==set(state["material_tags"])-{0}
    assert bcount==state["material"].count(2) and ncount==sum(r!=0 for r in state["templates"])
    close(abs(ref.energy(state)-h0-state["cumulative_work"]+state["cumulative_heat"]),run["first_law_error"])
    return dict(seed=run["seed"],activity_mode=run["activity_mode"],flags=run["flags"],events=len(run["ledger"]),
        proposals=stats["proposals"],samples=len(samples),every_proposal_regenerated=full_proposals,event_counts=event_counts,
        signed_X_current_count=cx,signed_Y_current_count=cy,signed_background_current_count=c0,
        thermal_births=births,passive_losses=losses,initial_B=b0,final_B=bcount,initial_templates=n0,final_templates=ncount,
        initial_H=h0,final_H=ref.energy(state),chemical_work=state["cumulative_work"],heat_to_bath=state["cumulative_heat"],
        maximum_running_first_law_error=maximum_error,conserved_working_material_tags=True,
        exact_B_template_and_fuel_balances=True,
        scope="Full accepted physical and observational-ID replay; no claim of stationary sampling, phase formation or inherited individual.")


def observe_spatial_independent(material,templates):
    """Union-find membership plus a set-valued shell, not production BFS/roll."""
    import numpy as np
    material,templates=np.asarray(material),np.asarray(templates)
    if material.ndim==1:
        side=math.isqrt(material.size)
        material,templates=material.reshape(side,side),templates.reshape(side,side)
    side=material.shape[0]; n=material.size
    expected=observe_base_independent(material,templates==2)
    expected={"Y_"+key[2:] if key.startswith("E_") else key:value for key,value in expected.items()}
    expected["n_Y"]=expected.pop("n_E")
    occupied=set(int(i) for i in np.flatnonzero(material.ravel()==2))
    parent={i:i for i in occupied}
    def find(i):
        if parent[i]!=i: parent[i]=find(parent[i])
        return parent[i]
    neighbors=lambda i:(((i//side-1)%side)*side+i%side,((i//side+1)%side)*side+i%side,
                        (i//side)*side+(i%side-1)%side,(i//side)*side+(i%side+1)%side)
    for i in occupied:
        for j in neighbors(i):
            if j in occupied: parent[find(j)]=find(i)
    components={}
    for i in occupied: components.setdefault(find(i),set()).add(i)
    core=set().union(*(group for group in components.values() if len(group)>=9)) if components else set()
    near=core|{j for i in core for j in neighbors(i)}
    t=templates.ravel()
    nx,ny=int(np.count_nonzero(t==1)),int(np.count_nonzero(t==2))
    xnear=sum(t[i]==1 for i in near); ynear=sum(t[i]==2 for i in near)
    expected.update(n_X=nx,n_templates=nx+ny,template_density=(nx+ny)/n,
        templates_near_dense_B=int(xnear+ynear),template_density_near_dense_B=(xnear+ynear)/n,
        X_near_dense_B=int(xnear),Y_near_dense_B=int(ynear),dense_neighborhood_sites=len(near),
        large_B_mass_fraction=len(core)/len(occupied) if occupied else 0.)
    return expected


def check_spatial(seed=98134):
    import numpy as np
    from spatial_readouts import observe
    rng=np.random.default_rng(seed)
    cases=0
    for side in (4,8,16,32):
        for fixture in range(8):
            m=rng.integers(0,3,(side,side)); t=rng.integers(0,3,(side,side))
            if fixture==0: m[:]=0
            if fixture==1: m[:]=2
            if fixture==2:
                m[:]=0; m[:3,:3]=2
            if fixture==3:
                m[:]=0; m[[side-1,0,1],:3]=2
            if fixture==4:
                m[:]=0; m[:3,:3]=2; m[-3:,-3:]=2
            compare(observe_spatial_independent(m,t),observe(m,t)); cases+=1
    return dict(cases=cases,independent_membership="union-find",independent_shell="set union of explicit periodic neighbors",
                denominator="entire lattice area, not occupied neighborhood",nullable_morphology_preserved=True)


def reconstruct_renewal_independent(run,deadline=20.,threshold=4):
    """Global post-array census and parent-pointer ancestry; no DUT reader import.

    This deliberately rescans whole physical arrays after each event. It is a
    slower independent oracle for fixtures and optional full study readouts.
    """
    a={key:run["initial_state"][key][:] for key in ("material","templates","material_tags","instance_tags")}
    origins={tag:(0,"initial",r) for tag,r in zip(a["instance_tags"],a["templates"]) if tag}
    initial_x={tag for tag,r in zip(a["instance_tags"],a["templates"]) if r==1}
    initial_b=a["material"].count(2); initial_n=sum(bool(r) for r in a["templates"])
    def ancestor(tag,ancestors):
        while tag:
            if tag in ancestors:return True
            tag=origins[tag][0]
        return False
    def census():
        live={tag:r for tag,r in zip(a["instance_tags"],a["templates"]) if tag}
        from_x={tag for tag in live if ancestor(tag,initial_x)}
        from_focal={tag for tag in live if focal and ancestor(tag,{focal})}
        return live,from_x,from_focal
    focal=focal_time=extinction=established=first_y_loss=clearance=first_post=None
    first_y=None; peak=restorations=postcopies=postdesc=postreverse=0
    integral={key:0. for key in ("B","X","Y","templates","originals","lineage","lineage_Y","focal","focal_Y","post_templates","post_lineage")}
    counts={"thermal":0,"paid_copy":0}; losses={}; copycounts={}; phases={}; births={}; deaths={}
    entry={tag:(0.,True) for tag,mat in zip(a["material_tags"],a["material"]) if tag and mat==2}
    episodes=[]; entrycauses={}; exitcauses={}
    previous=0.
    def advance(now):
        nonlocal previous
        live,from_x,from_focal=census(); dt=now-previous
        values=dict(B=a["material"].count(2),X=sum(r==1 for r in live.values()),Y=sum(r==2 for r in live.values()),
                    templates=len(live),originals=len(set(live)&initial_x),lineage=len(from_x),
                    lineage_Y=sum(live[tag]==2 for tag in from_x),focal=len(from_focal),focal_Y=sum(live[tag]==2 for tag in from_focal),
                    post_templates=len(live) if clearance is not None else 0,post_lineage=len(from_x) if clearance is not None else 0)
        for key in integral:integral[key]+=dt*values[key]
        previous=now
    for e in run["ledger"]:
        time=e["time"]; advance(time)
        live,from_x,from_focal=census(); prior_focal_y=sum(live[tag]==2 for tag in from_focal)
        had_focal=focal is not None; cleared=clearance is not None
        kind=e["kind"]; created=e.get("instance_created",0); removed=e.get("instance_removed",0)
        cause=kind
        if kind=="copy":
            cause=f"copy_activity_{e['activity']}_{'forward' if e['direction']==1 else 'reverse'}"
            key=f"activity_{e['activity']}_{'forward' if e['direction']==1 else 'reverse'}"
            copycounts[key]=copycounts.get(key,0)+1
        elif kind=="phase_passive":
            cause=f"phase_passive_{e['before']['material'][0]}_to_{e['after']['material'][0]}"
            phases[cause]=phases.get(cause,0)+1
        if removed:
            losses[kind]=losses.get(kind,0)+1
            deathkey=str(live[removed]); deaths[deathkey]=deaths.get(deathkey,0)+1
            if cleared and kind=="copy" and removed in from_x:postreverse+=1
        if created:
            species=e["after"]["templates"][0]; parent=e.get("parent_instance",0)
            origin=e["birth_origin"]; origins[created]=(parent,origin,species)
            is_from_x=ancestor(created,initial_x)
            if focal is None and origin=="paid_copy" and e["helper_species"]==1 and species==2 and is_from_x and time<=deadline:
                focal,focal_time=created,time
            counts[origin]+=1
            key=f"{origin}_{species}";births[key]=births.get(key,0)+1
            if origin=="paid_copy":
                key=f"paid_helper_{e['helper_species']}_offspring_{species}";births[key]=births.get(key,0)+1
                if cleared:
                    postcopies+=1
                    if is_from_x:
                        postdesc+=1
                        if first_post is None:first_post=time
            if species==2 and first_y is None:
                first_y=dict(tag=created,time=time,origin=origin,parent=parent,parent_initial_X_ancestry=is_from_x)
        old={tag:mat for tag,mat in zip(a["material_tags"],a["material"]) if tag}
        sites=[e["i"]] if e["j"]==-1 else [e["i"],e["j"]]
        for key in a:
            for site,value in zip(sites,e["after"][key]):a[key][site]=value
        new={tag:mat for tag,mat in zip(a["material_tags"],a["material"]) if tag}
        for tag,mat in new.items():
            if mat==2 and old[tag]!=2:
                entry[tag]=(time,False);entrycauses[cause]=entrycauses.get(cause,0)+1
            if mat!=2 and old[tag]==2:
                start,left=entry.pop(tag);episodes.append((tag,start,time,left,False))
                exitcauses[cause]=exitcauses.get(cause,0)+1
        live,from_x,from_focal=census()
        if initial_x and not set(live)&initial_x and clearance is None:clearance=time
        focal_y=sum(live[tag]==2 for tag in from_focal)
        if focal is not None:
            peak=max(peak,len(from_focal))
            if established is None and len(from_focal)>=threshold:established=time
            if first_y_loss is None and focal_y==0:first_y_loss=time
            if had_focal and prior_focal_y==0 and focal_y>0:restorations+=1
            if not from_focal and extinction is None:extinction=time
    advance(run["horizon"])
    for tag,(start,left) in entry.items():episodes.append((tag,start,run["horizon"],left,True))
    live,from_x,from_focal=census()
    cx=copycounts.get("activity_1_forward",0)-copycounts.get("activity_1_reverse",0)
    cy=copycounts.get("activity_2_forward",0)-copycounts.get("activity_2_reverse",0)
    c0=phases.get("phase_passive_1_to_2",0)-phases.get("phase_passive_2_to_1",0)
    durations=[end-start for tag,start,end,left,right in episodes]
    observed=[end-start for tag,start,end,left,right in episodes if not left and not right]
    return dict(integral_B=integral["B"],integral_X=integral["X"],integral_Y=integral["Y"],integral_templates=integral["templates"],
        net_X_activity_copy_count=cx,net_Y_activity_copy_count=cy,net_background_A_to_B_count=c0,
        births={key:value for key,value in counts.items() if value},losses=losses,birth_species=births,death_species=deaths,
        copy_channel_counts=copycounts,phase_channel_counts=phases,first_new_Y_beyond_initial=first_y,
        focal_paid_Y_from_initial_X_appeared_by_deadline=focal is not None,focal_tag=focal,focal_appearance_time=focal_time,
        focal_established=established is not None,focal_establishment_time=established,maximum_live_focal_subtree=peak,
        focal_extinction_time=extinction,focal_survives_to_horizon=bool(from_focal),final_focal_subtree=len(from_focal),
        final_focal_Y=sum(live[tag]==2 for tag in from_focal),integral_focal_subtree=integral["focal"],integral_focal_Y=integral["focal_Y"],
        focal_Y_first_loss_time=first_y_loss,focal_Y_restorations_after_loss=restorations,
        focal_Y_uninterrupted_to_horizon=focal is not None and first_y_loss is None,
        initial_X_founders=len(initial_x),initial_X_founder_cohort_present=bool(initial_x),
        original_X_founders_cleared_after_presence=clearance is not None,original_X_founder_clearance_time=clearance,
        final_original_X_founders=len(set(live)&initial_x),final_initial_X_lineage=len(from_x),
        integral_original_X_founders=integral["originals"],integral_initial_X_lineage=integral["lineage"],
        integral_initial_X_lineage_Y=integral["lineage_Y"],integral_templates_after_original_X_clearance=integral["post_templates"],
        integral_initial_X_lineage_after_original_X_clearance=integral["post_lineage"],
        forward_paid_copies_after_original_X_clearance=postcopies,
        forward_paid_copies_by_initial_X_descendants_after_clearance=postdesc,
        reverse_consumed_initial_X_descendants_after_clearance=postreverse,
        first_paid_copy_by_initial_X_descendant_after_clearance=first_post,
        final_X=sum(r==1 for r in live.values()),final_Y=sum(r==2 for r in live.values()),final_B=a["material"].count(2),
        B_episode_summary=dict(episodes=len(episodes),left_censored=sum(e[3] for e in episodes),right_censored=sum(e[4] for e in episodes),
            fully_observed_episodes=len(observed),truncated_duration_sum=math.fsum(durations),
            truncated_duration_mean=math.fsum(durations)/len(durations) if durations else None,
            fully_observed_duration_mean=math.fsum(observed)/len(observed) if observed else None,
            truncated_duration_min=min(durations) if durations else None,truncated_duration_max=max(durations) if durations else None,
            entry_causes=entrycauses,exit_causes=exitcauses),
        episodes=sorted(episodes))


def check_renewal():
    from renewal_fixtures import histories,nested
    from renewal_readouts import measure
    cases=histories()
    for name,run,manual in cases:
        actual=measure(run,include_details=True)
        expected=reconstruct_renewal_independent(run)
        episodes=expected.pop("episodes")
        compare(expected,actual)
        observed=sorted((e["carrier"],e["start"],e["end"],e["left_censored"],e["right_censored"]) for e in actual["B_episode_records"])
        compare(episodes,observed)
        for key,value in manual.items():compare(value,nested(actual,key))
    invalid=load(ROOT/"renewal-fixture-results.json")["rejections"]
    assert all(row["rejected"] for row in invalid)
    return dict(manual_histories=len(cases),independent_global_reconstruction=True,
                independently_compared_fields="all numeric genealogy and occupation integrals, origins, signed currents and complete carrier episodes",
                invalid_histories_in_author_report=len(invalid),
                note="Manual adversarial histories, not new physical simulations. Later physical replay independently checks event admissibility and energy.")


def preflight():
    modelcheck=check_model(); constructed=check_constructed_instances(); spatial=check_spatial(); renewal=check_renewal()
    smoke=load(ROOT/"software-preflight.json"); rows=[]
    for row in smoke["traces"]:
        path=ROOT/row["file"]; assert digest(path)==row["sha256"]
        rows.append(dict(path=row["file"],sha256=digest(path),**verify_trajectory(load(path))))
    sources=("productive_model.py","model-spec.md","balance-and-cost.md","generalist-comparison.md",
             "study-design.md","spatial_readouts.py","dependencies/morphology.py","renewal_readouts.py",
             "renewal_fixtures.py","run_study.py","analyze_study.py","verify_productive.py",
             "../active-phase/verify_active_phase.py")
    output=dict(checks_passed=True,model_sha256=digest(ROOT/"productive_model.py"),source_sha256={name:digest(ROOT/name) for name in sources},
        scalar_model_checks=modelcheck,constructed_transition_checks=constructed,spatial_readout_checks=spatial,
        renewal_readout_checks=renewal,software_smoke_records=rows,
        totals=dict(software_records=len(rows),accepted_events=sum(r["events"] for r in rows),
                    proposals=sum(r["proposals"] for r in rows),samples=sum(r["samples"] for r in rows)),
        prospective_review=dict(records=64,seeds=list(range(730100,730108)),sides=[16,32],arms=["specialist","generalist","X_only","Y_only"],
            claim_A_required_gates=20,claim_B_required_gates=2,seed_level_t_intervals=True,
            undefined_required_values_retained=True,bootstrap_sensitivity_only=True,
            source_rates_have_no_cluster_or_ancestry_readout=True,
            no_blocking_source_defect_found=True),
        scope="Independent software/structural preflight only. All five saved traces are transport-only; constructed fixtures additionally cover paid/thermal births, losses, unrelated reverse catalyst and the complete driven cycle. No prospective physical outcomes have been generated by this audit, and no phase, population, generalist advantage, individuality or moral conclusion is established.")
    (ROOT/"preflight-verification.json").write_text(json.dumps(output,indent=2,sort_keys=True,allow_nan=False)+"\n")
    return output


if __name__=="__main__":
    result=preflight()
    print(json.dumps(dict(checks_passed=result["checks_passed"],model_sha256=result["model_sha256"],totals=result["totals"])))

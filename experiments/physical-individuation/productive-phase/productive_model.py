"""Joined local material conversion and productive copying; no group feedback.

One finite physical occupancy process and one Hamiltonian. Copying-instance IDs
and cumulative ledgers are observer histories, never coordinates in the rates.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Any

import numpy as np

S, A, B = 0, 1, 2
D, X, Y = 0, 1, 2
MATERIAL_SWAP, TEMPLATE_SWAP = "material_swap", "template_swap"
TEMPLATE_PASSIVE, PHASE_PASSIVE, COPY = "template_passive", "phase_passive", "copy"


@dataclass(frozen=True)
class Params:
    side: int = 16
    theta: float = 1.0
    epsilon_B: float = 2.0
    J_BB: float = 2.4
    E_X: float = 6.0
    E_Y: float = 6.0
    g_X: float = -2.5
    g_Y: float = 2.5
    mu: float = 12.0
    copy_rate: float = 0.1
    loss_rate: float = 0.02
    background_AB_rate: float = 0.0
    material_mobility: float = 1.0
    template_mobility: float = 0.5
    mutation_error: float = 0.05
    material_count: int = 77
    initial_X: int = 8
    initial_Y: int = 0


@dataclass(frozen=True)
class Flags:
    # Names refer to productive ACTIVITIES, also in generalist mode.
    enable_X: bool = True
    enable_Y: bool = True


@dataclass
class State:
    material: np.ndarray
    templates: np.ndarray
    material_tags: np.ndarray
    instance_tags: np.ndarray
    next_tag: int
    time: float = 0.0
    cumulative_work: float = 0.0
    cumulative_heat: float = 0.0
    net_fuel_to_waste: int = 0
    event_count: int = 0


@dataclass(frozen=True)
class Event:
    kind: str
    i: int
    j: int = -1
    before: int = -1
    after: int = -1
    activity: int = 0
    child: int = 0
    direction: int = 0


class Model:
    def __init__(self, params: Params | None = None, flags: Flags | dict[str, bool] | None = None,
                 activity_mode: str = "specialist"):
        self.params = params or Params()
        self.flags = Flags(**flags) if isinstance(flags, dict) else (flags or Flags())
        p = self.params
        if not isinstance(p.side, int) or p.side < 4:
            raise ValueError("side must be an integer >=4")
        self.n = p.side ** 2
        if activity_mode not in ("specialist", "generalist"):
            raise ValueError("activity_mode must be specialist or generalist")
        self.activity_mode = activity_mode
        for name in ("material_count", "initial_X", "initial_Y"):
            value = getattr(p, name)
            if not isinstance(value, int) or not 0 <= value <= self.n:
                raise ValueError(f"{name} must be an integer in [0,N]")
        if p.initial_X + p.initial_Y > self.n:
            raise ValueError("too many initial templates")
        for name in ("theta", "epsilon_B", "J_BB", "E_X", "E_Y", "g_X", "g_Y", "mu",
                     "copy_rate", "loss_rate", "background_AB_rate", "material_mobility",
                     "template_mobility", "mutation_error"):
            if not math.isfinite(getattr(p, name)):
                raise ValueError(f"{name} must be finite")
        if p.theta <= 0 or not 0 <= p.mutation_error <= .5:
            raise ValueError("temperature positive and mutation_error in [0,.5] required")
        if min(p.copy_rate, p.loss_rate, p.background_AB_rate, p.material_mobility, p.template_mobility) < 0:
            raise ValueError("mobilities must be nonnegative")
        if not all(isinstance(value, bool) for value in asdict(self.flags).values()):
            raise ValueError("activity flags must be booleans")
        self.E = np.array([0., p.E_X, p.E_Y])
        self.g = np.array([0., p.g_X, p.g_Y])
        neighbors, bonds = [], []
        for i in range(self.n):
            row, col = divmod(i, p.side)
            right, left = row*p.side+(col+1)%p.side, row*p.side+(col-1)%p.side
            down, up = ((row+1)%p.side)*p.side+col, ((row-1)%p.side)*p.side+col
            neighbors.append((right, down, left, up))
            bonds.extend(((i,right),(i,down)))
        self.neighbors = np.asarray(neighbors, dtype=np.int64)
        self.bonds = np.asarray(bonds, dtype=np.int64)
        self.bond_i, self.bond_j = self.bonds.T
        events, bounds = [], []
        for kind, mobility in ((MATERIAL_SWAP,p.material_mobility),(TEMPLATE_SWAP,p.template_mobility)):
            for i,j in bonds:
                events.append(Event(kind,i,j)); bounds.append(mobility)
        for i in range(self.n):
            for child in (X,Y):
                events.extend((Event(TEMPLATE_PASSIVE,i,before=D,after=child),
                               Event(TEMPLATE_PASSIVE,i,before=child,after=D)))
                bounds.extend((p.loss_rate,p.loss_rate))
        for i in range(self.n):
            events.extend((Event(PHASE_PASSIVE,i,before=A,after=B),Event(PHASE_PASSIVE,i,before=B,after=A)))
            bounds.extend((p.background_AB_rate,p.background_AB_rate))
        copy_bound = p.copy_rate*(1-p.mutation_error)
        for i in range(self.n):
            for j in self.neighbors[i]:
                for activity in (X,Y):
                    for child in (X,Y):
                        events.extend((Event(COPY,i,int(j),activity=activity,child=child,direction=1),
                                       Event(COPY,i,int(j),activity=activity,child=child,direction=-1)))
                        # Independent of mode/flags: controls share proposal clocks;
                        # disabled physical channels simply reject proposals.
                        bounds.extend((copy_bound,copy_bound))
        self.events = tuple(events)
        self._lookup = {edge:k for k,edge in enumerate(events)}
        self.reverse_index = np.r_[np.arange(4*self.n), np.arange(4*self.n,42*self.n)^1]
        self.proposal_mobilities = np.asarray(bounds,dtype=float)
        if len(self.events) != 42*self.n:
            raise RuntimeError("unexpected channel catalogue size")
        self.family_starts = np.array([0,2*self.n,4*self.n,8*self.n,10*self.n])
        self.family_counts = np.array([2*self.n,2*self.n,4*self.n,2*self.n,32*self.n])
        self.family_mobilities = np.array([p.material_mobility,p.template_mobility,p.loss_rate,
                                           p.background_AB_rate,copy_bound])
        self._proposal_family_cumulative = np.cumsum(self.family_counts*self.family_mobilities)
        self.total_proposal_rate = float(self._proposal_family_cumulative[-1])

    def initial(self, seed: int = 0) -> State:
        material_rng = np.random.default_rng(np.random.SeedSequence([seed,85601]))
        template_rng = np.random.default_rng(np.random.SeedSequence([seed,85602]))
        material = np.zeros(self.n,dtype=np.int8)
        material_tags = np.zeros(self.n,dtype=np.int64)
        sites = material_rng.choice(self.n,self.params.material_count,replace=False)
        material[sites] = A
        material_tags[sites] = np.arange(1,len(sites)+1)
        templates = np.zeros(self.n,dtype=np.int8)
        instance_tags = np.zeros(self.n,dtype=np.int64)
        sites = template_rng.choice(self.n,self.params.initial_X+self.params.initial_Y,replace=False)
        templates[sites[:self.params.initial_X]] = X
        templates[sites[self.params.initial_X:]] = Y
        instance_tags[sites] = np.arange(1,len(sites)+1)
        return State(material,templates,material_tags,instance_tags,len(sites)+1)

    @staticmethod
    def copy_state(state: State) -> State:
        return State(state.material.copy(),state.templates.copy(),state.material_tags.copy(),
                     state.instance_tags.copy(),state.next_tag,state.time,state.cumulative_work,
                     state.cumulative_heat,state.net_fuel_to_waste,state.event_count)

    @staticmethod
    def state_to_dict(state: State) -> dict[str, Any]:
        return dict(material=state.material.tolist(),templates=state.templates.tolist(),
                    material_tags=state.material_tags.tolist(),instance_tags=state.instance_tags.tolist(),
                    next_tag=int(state.next_tag),time=float(state.time),cumulative_work=float(state.cumulative_work),
                    cumulative_heat=float(state.cumulative_heat),net_fuel_to_waste=int(state.net_fuel_to_waste),
                    event_count=int(state.event_count))

    @staticmethod
    def state_from_dict(data: dict[str, Any]) -> State:
        return State(np.asarray(data["material"],dtype=np.int8),np.asarray(data["templates"],dtype=np.int8),
                     np.asarray(data["material_tags"],dtype=np.int64),np.asarray(data["instance_tags"],dtype=np.int64),
                     int(data["next_tag"]),float(data.get("time",0.)),float(data.get("cumulative_work",0.)),
                     float(data.get("cumulative_heat",0.)),int(data.get("net_fuel_to_waste",0)),
                     int(data.get("event_count",0)))

    def validate_state(self,state: State) -> None:
        arrays = (state.material,state.templates,state.material_tags,state.instance_tags)
        if any(a.shape != (self.n,) for a in arrays):
            raise ValueError("state arrays must have N entries")
        if not np.isin(state.material,(S,A,B)).all() or not np.isin(state.templates,(D,X,Y)).all():
            raise ValueError("invalid occupancy")
        carriers, templates = state.material != S, state.templates != D
        if int(carriers.sum()) != self.params.material_count:
            raise ValueError("material inventory changed")
        for present,tags in ((carriers,state.material_tags),(templates,state.instance_tags)):
            if np.any(tags[~present] != 0) or np.any(tags[present] <= 0) or len(np.unique(tags[present])) != int(present.sum()):
                raise ValueError("invalid observer tag assignment")
        if np.any(state.instance_tags >= state.next_tag) or state.next_tag <= 0:
            raise ValueError("new instance counter must exceed all live IDs")
        if not math.isfinite(state.time) or state.time < 0:
            raise ValueError("invalid observer time")

    def energy(self,state: State) -> float:
        b = state.material == B
        return float(self.params.epsilon_B*b.sum()-self.params.J_BB*np.sum(b[self.bond_i]&b[self.bond_j])
                     +self.E[state.templates].sum()-np.dot(self.g[state.templates],b))

    def event_index(self,event: Event | int) -> int:
        return int(event) if isinstance(event,(int,np.integer)) else self._lookup[event]

    def reverse(self,event: Event | int) -> Event:
        return self.events[int(self.reverse_index[self.event_index(event)])]

    def work(self,event: Event | int) -> float:
        edge = self.events[self.event_index(event)]
        return edge.direction*self.params.mu if edge.kind == COPY else 0.

    def delta_energy(self,state: State,event: Event | int) -> float:
        edge = self.events[self.event_index(event)]
        p,i,j = self.params,edge.i,edge.j
        b_i = int(state.material[i] == B)
        if edge.kind == MATERIAL_SWAP:
            change = int(state.material[j] == B)-b_i
            near_i = np.count_nonzero(state.material[self.neighbors[i]] == B)
            near_j = np.count_nonzero(state.material[self.neighbors[j]] == B)
            return float(-p.J_BB*change*(near_i-near_j-change)
                         -change*(self.g[state.templates[i]]-self.g[state.templates[j]]))
        if edge.kind == TEMPLATE_SWAP:
            return float(-(self.g[state.templates[j]]-self.g[state.templates[i]])
                         *(b_i-int(state.material[j] == B)))
        if edge.kind == TEMPLATE_PASSIVE:
            if state.templates[i] != edge.before:
                raise ValueError("template conversion source mismatch")
            return float(self.E[edge.after]-self.E[edge.before]
                         -(self.g[edge.after]-self.g[edge.before])*b_i)
        neighbors_b = np.count_nonzero(state.material[self.neighbors[i]] == B)
        if edge.kind == PHASE_PASSIVE:
            if state.material[i] != edge.before:
                raise ValueError("material conversion source mismatch")
            change = int(edge.after == B)-b_i
            return float(change*(p.epsilon_B-p.J_BB*neighbors_b-self.g[state.templates[i]]))
        if edge.activity == X:
            forward = self.E[edge.child]+p.epsilon_B-self.g[edge.child]-p.J_BB*neighbors_b
        else:
            forward = self.E[edge.child]-p.epsilon_B+p.J_BB*neighbors_b
        return float(edge.direction*forward)

    def event_rate(self,state: State,event: Event | int) -> float:
        k = self.event_index(event)
        edge,p = self.events[k],self.params
        if edge.kind == MATERIAL_SWAP:
            if state.material[edge.i] == state.material[edge.j]: return 0.
            mobility = p.material_mobility
        elif edge.kind == TEMPLATE_SWAP:
            if state.templates[edge.i] == state.templates[edge.j]: return 0.
            mobility = p.template_mobility
        elif edge.kind == TEMPLATE_PASSIVE:
            if state.templates[edge.i] != edge.before: return 0.
            mobility = p.loss_rate
        elif edge.kind == PHASE_PASSIVE:
            if state.material[edge.i] != edge.before: return 0.
            mobility = p.background_AB_rate
        else:
            if (edge.activity == X and not self.flags.enable_X) or (edge.activity == Y and not self.flags.enable_Y):
                return 0.
            parent = int(state.templates[edge.j])
            if parent == D or (self.activity_mode == "specialist" and parent != edge.activity):
                return 0.
            material_before,material_after = (A,B) if edge.activity == X else (B,A)
            target_template = D if edge.direction == 1 else edge.child
            target_material = material_before if edge.direction == 1 else material_after
            if state.templates[edge.i] != target_template or state.material[edge.i] != target_material:
                return 0.
            weight = 1-p.mutation_error if parent == edge.child else p.mutation_error
            mobility = p.copy_rate*weight
        if mobility == 0: return 0.
        exponent = min(0.,(self.work(k)-self.delta_energy(state,k))/p.theta)
        rate = mobility*math.exp(exponent)
        if rate == 0:
            raise FloatingPointError("enabled reverse rate underflow; parameter scale is not representable")
        return rate

    def rates(self,state: State):
        return self.events,np.fromiter((self.event_rate(state,k) for k in range(len(self.events))),
                                       dtype=float,count=len(self.events))

    @staticmethod
    def _snapshot(state: State,sites: list[int]) -> dict[str,list[int]]:
        return {name:[int(getattr(state,name)[i]) for i in sites]
                for name in ("material","templates","material_tags","instance_tags")}

    def apply(self,state: State,event: Event | int) -> dict[str,Any]:
        k = self.event_index(event)
        edge = self.events[k]
        if self.event_rate(state,k) <= 0:
            raise ValueError("cannot apply disabled event")
        sites = [edge.i] if edge.j == -1 else [edge.i,edge.j]
        before = self._snapshot(state,sites)
        delta,work = self.delta_energy(state,k),self.work(k)
        helper = int(state.templates[edge.j]) if edge.kind == COPY else 0
        catalyst = int(state.instance_tags[edge.j]) if edge.kind == COPY else 0
        created,removed,parent = 0,0,0
        origin = None
        if edge.kind == MATERIAL_SWAP:
            state.material[edge.i],state.material[edge.j] = state.material[edge.j],state.material[edge.i]
            state.material_tags[edge.i],state.material_tags[edge.j] = state.material_tags[edge.j],state.material_tags[edge.i]
        elif edge.kind == TEMPLATE_SWAP:
            state.templates[edge.i],state.templates[edge.j] = state.templates[edge.j],state.templates[edge.i]
            state.instance_tags[edge.i],state.instance_tags[edge.j] = state.instance_tags[edge.j],state.instance_tags[edge.i]
        elif edge.kind == PHASE_PASSIVE:
            state.material[edge.i] = edge.after
        else:
            birth = edge.after != D if edge.kind == TEMPLATE_PASSIVE else edge.direction == 1
            if birth:
                created = state.next_tag
                state.next_tag += 1
                state.instance_tags[edge.i] = created
                origin = "thermal" if edge.kind == TEMPLATE_PASSIVE else "paid_copy"
                parent = catalyst if edge.kind == COPY else 0
            else:
                removed = int(state.instance_tags[edge.i])
                state.instance_tags[edge.i] = 0
            if edge.kind == TEMPLATE_PASSIVE:
                state.templates[edge.i] = edge.after
            else:
                state.templates[edge.i] = edge.child if edge.direction == 1 else D
                substrate,product = (A,B) if edge.activity == X else (B,A)
                state.material[edge.i] = product if edge.direction == 1 else substrate
        heat = work-delta
        fuel = edge.direction if edge.kind == COPY else 0
        state.cumulative_work += work
        state.cumulative_heat += heat
        state.net_fuel_to_waste += fuel
        state.event_count += 1
        return dict(time=float(state.time),event_index=k,kind=edge.kind,i=edge.i,j=edge.j,
                    activity=edge.activity,direction=edge.direction,offspring_species=edge.child,
                    helper_species=helper,catalyst_instance=catalyst,parent_instance=parent,
                    instance_created=created,instance_removed=removed,birth_origin=origin,
                    before=before,after=self._snapshot(state,sites),delta_H=delta,chemical_work=work,
                    heat_to_bath=heat,fuel_to_waste=fuel,cumulative_work=float(state.cumulative_work),
                    cumulative_heat=float(state.cumulative_heat),net_fuel_to_waste=int(state.net_fuel_to_waste),
                    event_count=state.event_count)

    def summary(self,state: State,include_state: bool = False) -> dict[str,Any]:
        b = state.material == B
        result = dict(time=float(state.time),H=self.energy(state),S=int((state.material == S).sum()),
                      A=int((state.material == A).sum()),B=int(b.sum()),D=int((state.templates == D).sum()),
                      X=int((state.templates == X).sum()),Y=int((state.templates == Y).sum()),
                      X_on_B=int(np.sum((state.templates == X)&b)),Y_on_B=int(np.sum((state.templates == Y)&b)),
                      BB_contacts=int(np.sum(b[self.bond_i]&b[self.bond_j])),chemical_work=float(state.cumulative_work),
                      heat_to_bath=float(state.cumulative_heat),net_fuel_to_waste=int(state.net_fuel_to_waste),
                      events=state.event_count,next_tag=state.next_tag)
        if include_state:
            result.update({name:getattr(state,name).tolist()
                           for name in ("material","templates","material_tags","instance_tags")})
        return result


def simulate(params: Params | None,seed: int,horizon: float,flags: Flags | dict[str,bool] | None = None,
             activity_mode: str = "specialist",*,initial_state: State | None = None,sample_dt: float = 5.,
             max_proposals: int = 10000000,max_events: int = 5000000):
    if not math.isfinite(horizon) or horizon < 0 or not math.isfinite(sample_dt) or sample_dt <= 0:
        raise ValueError("invalid observation times")
    if max_proposals <= 0 or max_events <= 0:
        raise ValueError("safeguards must be positive")
    model = Model(params,flags,activity_mode)
    state = model.initial(seed) if initial_state is None else model.copy_state(initial_state)
    model.validate_state(state)
    if state.time != 0 or state.cumulative_work != 0 or state.cumulative_heat != 0 or state.net_fuel_to_waste != 0 or state.event_count != 0:
        raise ValueError("initial observer counters must be zero; retain actual physical occupancy and IDs")
    initial,initial_h = model.state_to_dict(state),model.energy(state)
    rng = np.random.default_rng(np.random.SeedSequence([seed,85603]))
    grid = np.r_[np.arange(0.,horizon,sample_dt),horizon]
    samples,ledger = [],[]
    index,proposals,rejected_since_last = 0,0,0
    proposed_by_family = np.zeros(5,dtype=np.int64)
    accepted_by_family = np.zeros(5,dtype=np.int64)
    terminal_censored_time = None
    while state.time < horizon:
        next_time = (state.time+float(rng.exponential(1/model.total_proposal_rate))
                     if model.total_proposal_rate > 0 else math.inf)
        while index < len(grid) and grid[index] <= min(next_time,horizon):
            sample = model.summary(state,include_state=True); sample["time"] = float(grid[index])
            samples.append(sample); index += 1
        if next_time > horizon:
            terminal_censored_time = next_time if math.isfinite(next_time) else None
            state.time = float(horizon)
            break
        state.time = next_time
        proposals += 1
        family = int(np.searchsorted(model._proposal_family_cumulative,
                                     rng.random()*model.total_proposal_rate,side="right"))
        offset = int(rng.integers(model.family_counts[family]))
        k = int(model.family_starts[family]+offset)
        proposed_by_family[family] += 1
        probability = model.event_rate(state,k)/model.proposal_mobilities[k]
        if rng.random() < probability:
            record = model.apply(state,k)
            record.update(proposal_number=proposals,rejected_since_last_accepted=rejected_since_last)
            ledger.append(record)
            accepted_by_family[family] += 1
            rejected_since_last = 0
        else:
            rejected_since_last += 1
        if proposals >= max_proposals or len(ledger) >= max_events:
            raise RuntimeError("event safeguard exceeded; no truncated result returned")
    if horizon == 0:
        samples.append(model.summary(state,include_state=True))
    model.validate_state(state)
    first_law_error = abs(model.energy(state)-initial_h-state.cumulative_work+state.cumulative_heat)
    if first_law_error > 2e-6:
        raise RuntimeError("cumulative first law failed")
    families = (MATERIAL_SWAP,TEMPLATE_SWAP,TEMPLATE_PASSIVE,PHASE_PASSIVE,COPY)
    return dict(
        params=asdict(model.params),flags=asdict(model.flags),activity_mode=activity_mode,seed=seed,horizon=float(horizon),
        initial_state=initial,final_state=model.state_to_dict(state),ledger=ledger,samples=samples,
        final_summary=model.summary(state),first_law_error=first_law_error,
        proposal_statistics=dict(total_bound=model.total_proposal_rate,proposals=proposals,accepted=len(ledger),
            rejected=proposals-len(ledger),terminal_rejected_since_last_accepted=rejected_since_last,
            proposed_by_family=dict(zip(families,map(int,proposed_by_family))),
            accepted_by_family=dict(zip(families,map(int,accepted_by_family))),terminal_censored_proposal_time=terminal_censored_time),
        random_streams=dict(initial_material=[seed,85601],initial_templates=[seed,85602],
                            proposal_clock_event_acceptance=[seed,85603]),
        algorithm="Exact constant-bound thinning, one channel uniformly within mobility-weighted family; all rejected waiting times included. Proposal bounds independent of activity mode/enable flags.",
        supplied_resources="Two-layer lattice, randomized A carriers and pure-X founders, energies/affinities, finite catalytic alphabet and paired mutation grammar, thermal bath and ideal homogeneous F/W chemostats. Initial molecular preparation and reservoir construction excluded from chemical-work ledger.",
        scope="Joined productive transformations only. No observer phase/group/lineage feedback. Generalists have both activities at full mobility; no added specificity penalty. Finite units, stable coexistence, ancestry, useful differentiation, heredity, aging, individuality and ethical status are separate unproved outcomes.",
    )

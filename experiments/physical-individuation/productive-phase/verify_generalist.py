"""Symbolic specialist/generalist comparison; no trajectories or old-file edits.

Requires SymPy. Full state alphabet, mutation and thermal appearance are kept.
The homogeneous proof is not a spatial/finite-population persistence theorem.
"""
from __future__ import annotations

import hashlib
from itertools import product
import json
import math
from pathlib import Path
import sys

sys.dont_write_bytecode = True
import sympy as s


def certificate() -> dict:
    x,y,b,k,ell,delta,eta,M,S,m,u = s.symbols(
        "x y b k ell delta eta M S m u", positive=True
    )
    c = [x,y]
    d,a = M-x-y,S-b
    mutation = [[1-m,m],[m,1-m]]

    def build(weights):
        ca = s.Matrix(2,2,lambda i,j: weights[i][0]*mutation[i][j]*c[i]*(k*d*a-ell*c[j]*b))
        cb = s.Matrix(2,2,lambda i,j: weights[i][1]*mutation[i][j]*c[i]*(k*d*b-ell*c[j]*a))
        f = s.Matrix([
            sum(ca[i,0]+cb[i,0] for i in range(2))-delta*x+eta*d,
            sum(ca[i,1]+cb[i,1] for i in range(2))-delta*y+eta*d,
            sum(ca)-sum(cb),
        ])
        return ca,cb,f

    specialist = build([[1,0],[0,1]])
    half = s.Rational(1,2)
    generalist = build([[half,half],[half,half]])
    full = build([[1,1],[1,1]])
    symmetric = {x:u,y:u,b:S/2}
    assert s.simplify((specialist[2]-generalist[2]).subs(symmetric)) == s.zeros(3,1)
    for direction in [0,1]:
        assert s.simplify((sum(specialist[direction])-sum(generalist[direction])).subs(symmetric)) == 0

    U,V = S*(2*k+ell)/2,k*M*S/2-delta-2*eta
    symmetric_rhs = -U*u*u+V*u+eta*M
    assert s.simplify(generalist[2][0].subs(symmetric)-symmetric_rhs) == 0
    P = s.Matrix([[1,1,0],[1,-1,0],[0,0,1]])
    J = s.simplify(P*generalist[2].jacobian([x,y,b]).subs(symmetric)*P.inv())
    r,c0,t,Q = 1-2*m,ell*u*S/2,eta*(M-2*u)/u,2*u*(k*(M-2*u)+ell*u)
    damping = c0+(1-r)*delta+r*t
    expected = s.diag(s.diff(symmetric_rhs,u),-damping,-Q)
    delta_fp = k*(M-2*u)*S/2-ell*u*S/2+t
    assert s.simplify((J-expected).subs(delta,delta_fp)) == s.zeros(3,3)
    assert all(s.simplify(J[i,j])==0 for i in range(3) for j in range(3) if i!=j)

    # Doubling both generalist pathway mobilities is not a free change of the
    # instantaneous lattice hazard cap: the eligible target states are exclusive.
    full_rhs = s.simplify(full[2][0].subs(symmetric))
    assert s.simplify(full_rhs-(S*u*(k*(M-2*u)-ell*u)-delta*u+eta*(M-2*u))) == 0

    # Same energies and alphabet imply the same zero-affinity thermal reference,
    # even with the extra paired kinetic pathways and positive mutation.
    q = s.symbols("q", positive=True)
    thermal_u = q*M/(1+2*q)
    equilibrium = {x:thermal_u,y:thermal_u,b:S/2,eta:delta*q,ell:k/q}
    for candidate in [specialist,generalist,full]:
        assert s.simplify(candidate[2].subs(equilibrium)) == s.zeros(3,1)
        assert s.simplify(candidate[0].subs(equilibrium)) == s.zeros(2,2)
        assert s.simplify(candidate[1].subs(equilibrium)) == s.zeros(2,2)

    # Algebraic opportunity allocation: alpha*A=(1-alpha)*B gives AB/(A+B).
    aa,bb,N = s.symbols("a0 b0 N",positive=True)
    alpha = bb/(aa+bb)
    assert s.simplify(alpha*aa-(1-alpha)*bb)==0
    assert s.simplify(alpha*aa-aa*bb/(aa+bb))==0
    pp = s.sqrt(aa)/(s.sqrt(aa)+s.sqrt(bb))
    denom = 1/(aa*(1-pp))+1/(bb*pp)
    assert s.simplify(1/denom-aa*bb/(s.sqrt(aa)+s.sqrt(bb))**2)==0
    adaptive_p = aa/(aa+bb)
    assert s.simplify(aa*(1-adaptive_p)-aa*bb/(aa+bb))==0
    assert s.simplify(bb*adaptive_p-aa*bb/(aa+bb))==0

    # Exhaustive four-neighbor microscopic eligibility screen. Every neighbor
    # has exactly one material state S/A/B and one template state D/X/Y.
    # Reactions only act on A/B. At most one family can make a forward copy at
    # a D target; acceptance <=1 and summed child mutation weights equal1.
    states = list(product(["S","A","B"],["D","X","Y"]))
    maxima = {"specialist":0.,"half_generalist":0.,"full_generalist":0.}
    fixtures=0
    for neighbors in product(states,repeat=4):
        for parent in [0,1]:
            for name,weights in [
                ("specialist",[[1,0],[0,1]]),
                ("half_generalist",[[.5,.5],[.5,.5]]),
                ("full_generalist",[[1,1],[1,1]]),
            ]:
                hazard=0.
                for material,target in neighbors:
                    if target=="D" and material in ["A","B"]:
                        hazard+=weights[parent][0 if material=="A" else 1]
                assert hazard<=4.
                maxima[name]=max(maxima[name],hazard)
            fixtures+=1
    assert maxima=={"specialist":4.,"half_generalist":2.,"full_generalist":4.}

    params={M:2.,S:2.,k:1.,ell:math.exp(-4),delta:.1,eta:.1*math.exp(-6),m:.02}
    uv=float(((V+s.sqrt(V*V+4*U*eta*M))/(2*U)).subs(params))
    sub=params|{u:uv}
    half_eigen=[float(s.diff(symmetric_rhs,u).subs(sub)), -float(damping.subs(sub)),-float(Q.subs(sub))]
    full_U = S*(2*k+ell)
    full_V = k*M*S-delta-2*eta
    full_u=float(((full_V+s.sqrt(full_V**2+4*full_U*eta*M))/(2*full_U)).subs(params))
    net=2*(.1*uv-params[eta]*(2-2*uv))
    return {
        "scope":"Exact homogeneous concentration algebra and eligibility bounds; no simulation",
        "checks":{
            "same_alphabet_sources_reverse_and_mutation":True,
            "specialist_half_generalist_symmetric_vector_fields_equal":True,
            "two_material_currents_equal_on_symmetric_submanifold":True,
            "generalist_diagonal_jacobian_at_fixed_point":True,
            "same_positive_thermal_equilibrium_at_zero_affinity":True,
            "local_opportunity_allocation_identities":True,
            "four_neighbor_eligibility_parent_fixtures":fixtures,
            "max_forward_hazard_per_unit_k":maxima,
            "full_generalist_not_excluded_by_four_k_hazard_cap":True,
        },
        "theorem_domain":{"positive":"M,S,k,ell,delta,eta","mutation_interval":[0,.5]},
        "example":{
            "theta":1.,"E":6.,"Delta":10.,"M":2.,"S":2.,"mutation":.02,
            "k":1.,"ell":params[ell],"delta":.1,"eta":params[eta],
            "specialist_and_half_generalist_each_type":uv,
            "half_generalist_eigenvalues":half_eigen,
            "common_net_copying":net,
            "common_chemical_power":10*net,
            "full_generalist_each_type":full_u,
            "full_generalist_net_copying":2*(.1*full_u-params[eta]*(2-2*full_u)),
        },
        "not_claimed":["finite_path_law_equality_for_generalist_control","spatial_unit_formation","taxonomic_necessity","universal_biodiversity_optimum","exploiter_resistance"],
    }


if __name__=="__main__":
    here=Path(__file__).resolve().parent
    result=certificate()
    result["source_sha256"]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result["note_sha256"]=hashlib.sha256((here/"generalist-comparison.md").read_bytes()).hexdigest()
    output=here/"generalist-certificate.json"
    output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"passed":True,"output":str(output),"example":result["example"]},indent=2))

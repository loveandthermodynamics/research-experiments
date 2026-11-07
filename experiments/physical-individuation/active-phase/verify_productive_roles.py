"""Exact algebra checks for productive-roles-bridge.md; no trajectory simulation.

Requires SymPy. Writes only productive-role-certificate.json alongside this file.
The ideal ODE, finite-state balance identities and spatial proposal have separate
scopes. This certificate checks the ideal algebra and stoichiometric identities;
it does not certify a spatial implementation or statistical lineage experiment.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import sys

sys.dont_write_bytecode = True
import sympy as s


def certificate() -> dict:
    x, y, b, k, ell, delta, eta, M, S, m, u = s.symbols(
        "x y b k ell delta eta M S m u", positive=True
    )
    d, a = M-x-y, S-b
    currents = s.Matrix([
        (1-m)*x*(k*d*a-ell*x*b),
        m*x*(k*d*a-ell*y*b),
        m*y*(k*d*b-ell*x*a),
        (1-m)*y*(k*d*b-ell*y*a),
    ])
    jxx, jxy, jyx, jyy = currents
    f = s.Matrix([
        jxx+jyx-delta*x+eta*d,
        jxy+jyy-delta*y+eta*d,
        jxx+jxy-jyx-jyy,
    ])
    # Transformation to x+y, x-y, b. The total and antisymmetric modes decouple
    # only at the stated symmetric fixed point.
    P = s.Matrix([[1,1,0], [1,-1,0], [0,0,1]])
    J = s.simplify(P*f.jacobian([x,y,b]).subs({x:u,y:u,b:S/2})*P.inv())
    for i, j in [(0,1), (0,2), (1,0), (2,0)]:
        assert s.simplify(J[i,j]) == 0
    U = S*(2*k+ell)/2
    V = k*M*S/2-delta-2*eta
    sym_rhs = -U*u*u+V*u+eta*M
    assert s.simplify(f[0].subs({x:u,y:u,b:S/2})-sym_rhs) == 0
    assert s.simplify(J[0,0]-s.diff(sym_rhs,u)) == 0
    d0, r = M-2*u, 1-2*m
    c0, t, Q = ell*u*S/2, eta*(M-2*u)/u, 2*u*(k*(M-2*u)+ell*u)
    Jminus = s.Matrix([
        [-(c0+(1-r)*delta+r*t), -r*Q],
        [delta-r*c0-t, -Q],
    ])
    # This substitution is exactly the symmetric fixed-point relation.
    delta_at_fixed_point = k*d0*S/2-ell*u*S/2+t
    assert s.simplify((J[1:3,1:3]-Jminus).subs(delta,delta_at_fixed_point)) == s.zeros(2)
    determinant = Q*(delta+4*m*(1-m)*c0)
    assert s.simplify(J[1:3,1:3].det()-determinant) == 0
    # At zero fuel affinity the kinetic channels are retained, but their paired
    # rates have the same equilibrium as passive template conversion.
    q = s.symbols("q", positive=True)
    thermal_u = q*M/(1+2*q)
    equilibrium_sub = {x:thermal_u,y:thermal_u,b:S/2,eta:delta*q,ell:k/q}
    assert s.simplify(currents.subs(equilibrium_sub)) == s.zeros(4,1)
    assert s.simplify(f.subs(equilibrium_sub)) == s.zeros(3,1)

    # Species order: X,Y,D,A,B,F,W. Copy reactions are XX,XY,YX,YY;
    # the parent catalytic instance cancels from each net vector.
    nu = s.Matrix([
        [1,0,1,0,-1,0],
        [0,1,0,1,0,-1],
        [-1,-1,-1,-1,1,1],
        [-1,-1,1,1,0,0],
        [1,1,-1,-1,0,0],
        [-1,-1,-1,-1,0,0],
        [1,1,1,1,0,0],
    ])
    inventories = s.Matrix([
        [1,1,1,0,0,0,0],
        [0,0,0,1,1,0,0],
        [0,0,0,0,0,1,1],
    ])
    assert inventories*nu == s.zeros(3,6)
    template_count = s.Matrix([[1,1,0,0,0,0,0]])
    phase_count = s.Matrix([[0,0,0,0,1,0,0]])
    assert template_count*nu == s.Matrix([[1,1,1,1,-1,-1]])
    assert phase_count*nu == s.Matrix([[1,1,-1,-1,0,0]])

    # Work/standard-energy ledger for explicit closed fuel: every copy consumes
    # one fuel, creates one E-energy template, and changes no ideal A/B energy.
    E, Delta = s.symbols("E Delta", positive=True)
    energy = s.Matrix([[E,E,0,0,0,Delta,0]])
    assert energy*nu == s.Matrix([[E-Delta]*4+[-E,-E]])

    pv = {M:2., S:2., k:1., ell:math.exp(-4), delta:.1,
          eta:.1*math.exp(-6), m:.02}
    uv = float(((V+s.sqrt(V*V+4*U*eta*M))/(2*U)).subs(pv))
    dv = 2-2*uv
    sub = pv | {u:uv}
    jm = Jminus.subs(sub)
    tr, det = float(s.trace(jm)), float(jm.det())
    disc = tr*tr-4*det
    assert disc < 0 and tr < 0 and det > 0
    family = .1*uv-pv[eta]*dv
    return {
        "scope":"Symbolic ideal-concentration algebra and stoichiometry only; no simulation",
        "checks":{
            "four_net_channels_reconstructed":True,
            "symmetric_fixed_point_polynomial":True,
            "symmetric_antisymmetric_decoupling":True,
            "antisymmetric_jacobian_exact_at_fixed_point":True,
            "positive_determinant_identity":str(determinant),
            "zero_affinity_copy_channels_have_thermal_fixed_point_and_zero_net_currents":True,
            "three_material_inventories_conserved":True,
            "copy_and_passive_loss_count_balances":True,
            "closed_fuel_standard_energy_ledger":True,
        },
        "theorem_domain":{
            "positive":"M,S,k,ell,delta,eta",
            "mutation_interval":[0,.5],
            "claim":"Local asymptotic stability of the symmetric physical fixed point",
            "not_claimed":["global attraction","spatial formation","finite-population lineage persistence","universal diversity optimum"],
        },
        "example":{
            "theta":1.,"template_standard_energy":6.,"fuel_chemical_potential_gap":10.,
            "M":2.,"S":2.,"k":1.,"ell":pv[ell],"delta":.1,"eta":pv[eta],"mutation":.02,
            "x_equals_y":uv,"D":dv,"A_equals_B":1.,
            "thermal_only_each_template":pv[eta]*2/(.1+2*pv[eta]),
            "net_copy_current_each_family":family,
            "total_net_copy_current":2*family,
            "chemical_power_per_normalized_volume":20*family,
            "symmetric_eigenvalue":-float(s.sqrt(V*V+4*U*eta*M).subs(pv)),
            "antisymmetric_eigenvalues":{"real":tr/2,"imaginary_magnitudes":math.sqrt(-disc)/2},
            "antisymmetric_trace":tr,"antisymmetric_determinant":det,
        },
    }


if __name__ == "__main__":
    here = Path(__file__).resolve().parent
    report = certificate()
    report["source_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    note = here/"productive-roles-bridge.md"
    report["note_sha256"] = hashlib.sha256(note.read_bytes()).hexdigest()
    output = here/"productive-role-certificate.json"
    output.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"passed":True,"output":str(output),"example":report["example"]},indent=2))

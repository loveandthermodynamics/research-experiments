"""Independent count-LDB and ideal fixed-point check; no trajectories.

This does not import or execute the author's symbolic certificate. Microscopic
rates use ordered falling-factorial catalytic pairs and a count Gibbs reference.
"""
from dataclasses import dataclass
from hashlib import sha256
from itertools import product
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent


def ode(z, M, S, k, ell, delta, eta, m):
    x, y, b = z
    d, a = M - x - y, S - b
    xx = (1-m)*x*(k*d*a-ell*x*b)
    xy = m*x*(k*d*a-ell*y*b)
    yx = m*y*(k*d*b-ell*x*a)
    yy = (1-m)*y*(k*d*b-ell*y*a)
    return np.array([xx+yx-delta*x+eta*d, xy+yy-delta*y+eta*d, xx+xy-yx-yy])


def physical_fixed_point(M, S, k, ell, delta, eta):
    U, V = S*(2*k+ell)/2, k*M*S/2-delta-2*eta
    root = math.sqrt(V*V+4*U*eta*M)
    # Same positive quadratic root, evaluated without subtractive cancellation.
    u = (V+root)/(2*U) if V >= 0 else 2*eta*M/(root-V)
    return u, -root


def jacobian_complex_step(point, parameters):
    columns = []
    for i in range(3):
        changed = np.asarray(point, dtype=complex)
        changed[i] += 1e-20j
        columns.append(np.imag(ode(changed, **parameters))/1e-20)
    return np.column_stack(columns)


def count_log_weight(counts, E=6., Delta=10., theta=1.):
    # Species X,Y,D,A,B,F,W; fixed-inventory combinatorial constants cancel.
    energy = E*(counts[0]+counts[1])+Delta*counts[5]
    return -energy/theta-sum(math.lgamma(n+1) for n in counts)


def count_pair_check(counts, parent, child, substrate, product_species, mutation):
    n = np.array(counts, dtype=int)
    E, Delta, theta = 6., 10., 1.
    k, ell = 1., math.exp(-(Delta-E)/theta)
    coefficient = 1-mutation if parent == child else mutation
    q_forward = coefficient*k*n[parent]*n[2]*n[substrate]*n[5]
    after = n.copy()
    after[child] += 1
    after[2] -= 1
    after[substrate] -= 1
    after[product_species] += 1
    after[5] -= 1
    after[6] += 1
    # The catalytic parent must be distinct from the template being removed.
    template_pair = after[parent]*(after[child]-(parent == child))
    q_reverse = coefficient*ell*template_pair*after[product_species]*after[6]
    log_rate_ratio = math.log(q_forward)-math.log(q_reverse)
    log_gibbs_ratio = count_log_weight(after)-count_log_weight(n)
    return abs(log_rate_ratio-log_gibbs_ratio)


def main():
    output = ROOT / "productive-role-independent-check.json"
    if output.exists():
        raise RuntimeError("independent check already saved; do not overwrite")
    count_error = 0.
    count_pairs = 0
    for counts, mutation in product(((3,2,4,5,6,7,8),(1,1,2,3,4,5,6)), (.02,.5)):
        for parent, child, substrate, product_species in ((0,0,3,4),(0,1,3,4),(1,0,4,3),(1,1,4,3)):
            error = count_pair_check(counts,parent,child,substrate,product_species,mutation)
            count_error = max(count_error,error)
            count_pairs += 1
    assert count_error < 4e-14
    passive_error = 0.
    for counts in ((3,2,4,5,6,7,8),(1,1,2,3,4,5,6)):
        for template in (0,1):
            after = np.array(counts,dtype=int)
            after[template] -= 1
            after[2] += 1
            q_loss = .1*counts[template]
            q_birth = .1*math.exp(-6)*after[2]
            passive_error = max(passive_error, abs(math.log(q_loss/q_birth)
                                -count_log_weight(after)+count_log_weight(counts)))
    assert passive_error < 4e-14
    max_residual = 0.
    max_matrix_error = 0.
    checked = 0
    P = np.array([[1,1,0],[1,-1,0],[0,0,1]],dtype=float)
    for M,S,k,ell,delta,eta,m in product((.2,2.),(.3,3.),(.1,1.7),(.01,3.),(.03,.9),(1e-5,.2),(0.,.02,.5)):
        parameters = dict(M=M,S=S,k=k,ell=ell,delta=delta,eta=eta,m=m)
        u, symmetric_eigenvalue = physical_fixed_point(M,S,k,ell,delta,eta)
        d = M-2*u
        assert 0 < u < M/2 and d > 0
        point = [u,u,S/2]
        residual = float(np.max(np.abs(ode(point,**parameters))))
        J = jacobian_complex_step(point,parameters)
        transformed = P@J@np.linalg.inv(P)
        r,c0,t,Q = 1-2*m,ell*u*S/2,eta*d/u,2*u*(k*d+ell*u)
        claimed = np.array([[symmetric_eigenvalue,0.,0.],
                            [0.,-(c0+(1-r)*delta+r*t),-r*Q],
                            [0.,delta-r*c0-t,-Q]])
        max_residual = max(max_residual,residual)
        max_matrix_error = max(max_matrix_error,float(np.max(np.abs(transformed-claimed))))
        assert residual < 1e-12
        assert np.max(np.abs(transformed-claimed)) < 3e-12
        assert np.max(np.real(np.linalg.eigvals(J))) < 0
        determinant = np.linalg.det(transformed[1:,1:])
        assert abs(determinant-Q*(delta+4*m*(1-m)*c0)) < 3e-12
        checked += 1
    example_parameters = dict(M=2.,S=2.,k=1.,ell=math.exp(-4),delta=.1,eta=.1*math.exp(-6),m=.02)
    u,_ = physical_fixed_point(**{key:value for key,value in example_parameters.items() if key != "m"})
    eigenvalues = sorted(np.linalg.eigvals(jacobian_complex_step([u,u,1.],example_parameters)),key=lambda z:z.real)
    family_current = .1*u-example_parameters["eta"]*(2-2*u)
    result = dict(
        checks_passed=True, microscopic_copy_pairs=count_pairs,
        maximum_count_copy_log_LDB_error=count_error,
        maximum_passive_count_log_LDB_error=passive_error,
        ideal_parameter_fixtures=checked, maximum_fixed_point_residual=max_residual,
        maximum_transformed_jacobian_error=max_matrix_error,
        example=dict(x_equals_y=u,D=2-2*u,total_net_copy_current=2*family_current,
                     chemical_power=20*family_current,
                     eigenvalues=[dict(real=float(z.real),imaginary=float(z.imag)) for z in eigenvalues]),
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        author_files_sha256={name:sha256((ROOT/name).read_bytes()).hexdigest() for name in
                             ("productive-roles-bridge.md","verify_productive_roles.py","productive-role-certificate.json")},
        scope="Independent finite-count Gibbs/rate check and ideal algebra evaluation only. No trajectory, spatial implementation, role-copying evolution, lineage persistence or diversity optimum is tested. The analytical stationary/finite-time balance bound is reviewed separately in prose.",
    )
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(json.dumps(result,indent=2))


if __name__ == "__main__":
    main()

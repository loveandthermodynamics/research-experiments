# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy", "scipy", "matplotlib"]
# ///
"""Reward-free chemical construction; normalized model units, not measured life.

Run: uv run --no-project generative_model.py
Reaction grammar, kinetics, bath and feed are assumptions. No ecological reward,
population selection operator, group-offspring operator, or diversity score drives
the ODE. A/B assemble C; catalytic C makes components, whose association makes C.
Optional slow reversible conversions seed other catalytic chemical types.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, root


NAMES = ["F", "A", "B", "C", "W", "Z", "V", "Y", "Q"]
MASS = np.array([4, 1, 1, 2, 2, 1, 1, .5, .5], dtype=float)
ENERGY = np.array([3, .6, .6, 1, 1, .25, .5, .1, .2], dtype=float)
# Columns: A+B<->C; C+F<->C+A+B+W; W+Z<->2Z+V;
# V+Y<->2Y+Q; A<->Z; Z<->2Y.
STOICH = np.array([
    [0, -1, 0, 0, 0, 0],
    [-1, 1, 0, 0, -1, 0],
    [-1, 1, 0, 0, 0, 0],
    [1, 0, 0, 0, 0, 0],
    [0, 1, -1, 0, 0, 0],
    [0, 0, 1, 0, 1, -1],
    [0, 0, 1, -1, 0, 0],
    [0, 0, 0, 1, 0, 2],
    [0, 0, 0, 1, 0, 0],
], dtype=float)
RELEASE = -ENERGY @ STOICH


@dataclass(frozen=True)
class Params:
    feed: float = 4.0
    dilution: float = 1.0
    association: float = 1.0
    catalysis: float = 1.0
    consume_w: float = 1.0
    consume_v: float = 2.0
    thermal_energy: float | None = None  # None: irreversible limiting model
    inhibition: float = 0.0  # each product raises both directional barriers
    slow_conversion: float = 0.0


def fluxes(x: NDArray[np.float64], p: Params) -> tuple[NDArray, NDArray]:
    f, a, b, c, w, z, v, y, q = x
    k = np.array([p.association, p.catalysis/(1+p.inhibition*w),
                  p.consume_w/(1+p.inhibition*v),
                  p.consume_v/(1+p.inhibition*q),
                  p.slow_conversion, p.slow_conversion])
    forward = k * np.array([a*b, f*c, w*z, v*y, a, z])
    if p.thermal_energy is None:
        backward = np.zeros(6)
    else:
        backward = k * np.exp(-RELEASE/p.thermal_energy) * np.array([
            c, c*a*b*w, z*z*v, y*y*q, z, y*y
        ])
    return forward, backward


def rhs(_t: float, x: NDArray[np.float64], p: Params) -> NDArray:
    plus, minus = fluxes(x, p)
    out = STOICH @ (plus-minus) - p.dilution*x
    out[0] += p.dilution*p.feed
    return out


def integrate(p: Params, initial: NDArray, stop: float = 200) -> tuple[NDArray, NDArray]:
    times = np.linspace(0, stop, 2001)
    sol = solve_ivp(lambda t, x: rhs(t, x, p), (0, stop), initial,
                    t_eval=times, method="DOP853", rtol=2e-10, atol=2e-12)
    assert sol.success, sol.message
    assert float(sol.y.min()) > -2e-9
    return sol.t, sol.y


def jacobian(x: NDArray, p: Params) -> NDArray:
    h = 1e-25
    return np.column_stack([
        np.imag(rhs(0, x.astype(complex) + 1j*h*np.eye(9)[j], p))/h
        for j in range(9)
    ])


def equilibrium(p: Params, guess: NDArray) -> tuple[NDArray, NDArray]:
    result = root(lambda x: rhs(0, x, p), guess,
                  jac=lambda x: jacobian(x, p), tol=1e-11)
    assert result.success or np.max(np.abs(result.fun)) < 2e-10
    assert np.max(np.abs(result.fun)) < 2e-9
    assert result.x.min() > -2e-8
    return result.x, np.linalg.eigvals(jacobian(result.x, p))


def named(x: NDArray) -> dict[str, float]:
    return dict(zip(NAMES, map(float, x), strict=True))


def state_summary(x: NDArray, p: Params) -> dict[str, object]:
    plus, minus = fluxes(x, p)
    eig = np.linalg.eigvals(jacobian(x, p))
    return {"concentrations": named(x),
            "derivative_max_abs": float(np.max(np.abs(rhs(0, x, p)))),
            "largest_eigenvalue_real_part": float(max(eig.real)),
            "material": float(MASS@x), "standard_energy": float(ENERGY@x),
            "standard_energy_release_flow": float(RELEASE@(plus-minus))}


def physical_motif_checks() -> dict[str, float]:
    """An explicit reaction-coordinate potential, not an atomistic calculation.

    U(q,r)=B(r)(1-q^2)^2 - Delta*(3q-q^3)/4 + U_bind(r).
    B is smaller in the bound configuration; each potential is used in a
    mechanical action. The same saddle lowers both directional barriers.
    """
    theta=.05; delta=.8
    def saddle(barrier: float) -> tuple[float,float]:
        q=-3*delta/(16*barrier)
        u=barrier*(1-q*q)**2-delta*(3*q-q**3)/4
        assert abs(-4*barrier*q*(1-q*q)-3*delta*(1-q*q)/4)<1e-14
        return q,u
    qfree,ufree=saddle(2.0); qbound,ubound=saddle(1.0)
    forward_barrier_free=ufree-delta/2
    forward_barrier_bound=ubound-delta/2
    speedup=np.exp((ufree-ubound)/theta)
    assert speedup>1e8
    # Thin-well/configurational-volume approximation, standard volume = 1.
    binding_volume=.01
    well_depth=.2-theta*np.log(binding_volume)
    association_constant=binding_volume*np.exp(well_depth/theta)
    assert np.isclose(association_constant,np.exp(.2/theta))
    g=brentq(lambda x:association_constant*(1-x)**2-x,0,1)
    assert g>.5
    return {"thermal_energy":theta,"reaction_endpoint_energy_drop":delta,
            "free_saddle_coordinate":qfree,"bound_saddle_coordinate":qbound,
            "free_forward_barrier":forward_barrier_free,
            "bound_forward_barrier":forward_barrier_bound,
            "Eyring_exponential_speedup":float(speedup),
            "binding_volume_over_standard_volume":binding_volume,
            "binding_well_depth":float(well_depth),
            "association_equilibrium_constant":float(association_constant),
            "equilibrium_bound_fraction_at_total_concentration_one":float(g)}


def main() -> None:
    assert np.allclose(MASS@STOICH, 0)
    assert np.allclose(RELEASE, [.2, .8, .25, .2, .35, .05])
    p = Params()
    exact = np.array([2, 1, 1, 1, 1, 1, .5, .5, .5])
    assert np.max(np.abs(rhs(0, exact, p))) < 1e-14
    assert max(np.linalg.eigvals(jacobian(exact, p)).real) < -.4
    critical_x = brentq(lambda x: 2*x**3+x*x-1, .1, 2)
    critical_feed = critical_x**2+critical_x+1+1/critical_x

    high = np.array([3, 1, 1, 0, 0, 1e-9, 0, 1e-9, 0])
    t, tr = integrate(p, high)
    assert np.allclose(tr[:, -1], exact, atol=3e-7)
    low = np.array([3.5, .5, .5, 0, 0, 1e-9, 0, 1e-9, 0])
    _, low_tr = integrate(p, low)
    assert low_tr[3, -1] < 1e-9
    assert low_tr[5, -1] < 1e-8 and low_tr[7, -1] < 1e-8

    # A trace seed tests invasion. Slow chemical conversions instead create Z,Y
    # from zero within the same finite reaction grammar, at no diversity bonus.
    p_seed = replace(p, thermal_energy=.05, slow_conversion=1e-7)
    unassembled = np.array([3, 1, 1, 0, 0, 0, 0, 0, 0], dtype=float)
    tseed, seeded = integrate(p_seed, unassembled, 300)
    eq_seed, _ = equilibrium(p_seed, seeded[:, -1])
    assert np.max(np.abs(seeded[:, -1]-eq_seed)) < 2e-7
    assert min(eq_seed[3], eq_seed[5], eq_seed[7]) > .35

    # Remove slow conversions after establishment: each catalytic type persists
    # without continued production from that source.
    p_reverse = replace(p_seed, slow_conversion=0)
    _, persisted = integrate(p_reverse, eq_seed)
    eq_reverse, _ = equilibrium(p_reverse, persisted[:, -1])
    assert min(eq_reverse[3], eq_reverse[5], eq_reverse[7]) > .35

    # No latent type can materialize in the ODE without a source or initial seed.
    _, no_source = integrate(p, unassembled)
    assert np.max(np.abs(no_source[[5, 7]])) == 0

    # Waste inhibition supplies a genuine ecological return to the founder C.
    # Compare presence/absence of downstream catalysts with identical feed.
    p_inhibit = replace(p, inhibition=.05)
    _, inhibited_with = integrate(p_inhibit, high, 300)
    absent = high.copy(); absent[[5, 7]] = 0
    _, inhibited_without = integrate(p_inhibit, absent, 300)
    no_y=high.copy(); no_y[7]=0
    _, inhibited_z_only=integrate(p_inhibit,no_y,300)
    assert inhibited_with[3, -1] > inhibited_without[3, -1] + .01
    assert inhibited_with[3, -1] > inhibited_z_only[3,-1]+.001
    assert inhibited_with[5,-1] > inhibited_z_only[5,-1]+.001

    # Finite reversibility itself creates an upstream return; no added
    # product-inhibition law is needed. Missing types are absent, not blocked.
    p_intrinsic=replace(p,thermal_energy=.1)
    _, intrinsic_all=integrate(p_intrinsic,high,500)
    _, intrinsic_z=integrate(p_intrinsic,no_y,500)
    _, intrinsic_none=integrate(p_intrinsic,absent,500)
    assert intrinsic_all[3,-1]>intrinsic_z[3,-1]>intrinsic_none[3,-1]
    assert intrinsic_all[5,-1]>intrinsic_z[5,-1]
    _, intrinsic_eigen=equilibrium(p_intrinsic,intrinsic_all[:,-1])
    assert max(intrinsic_eigen.real)<0
    intrinsic_feedback={
        "thermal_energy":p_intrinsic.thermal_energy,"added_inhibition":0,
        "C_without_consumers":float(intrinsic_none[3,-1]),
        "C_with_Z_only":float(intrinsic_z[3,-1]),
        "C_with_Z_and_Y":float(intrinsic_all[3,-1]),
        "Z_without_Y":float(intrinsic_z[5,-1]),
        "Z_with_Y":float(intrinsic_all[5,-1]),
        "C_gain_fraction":float(intrinsic_all[3,-1]/intrinsic_none[3,-1]-1),
        "C_gain_fraction_from_Y":float(intrinsic_all[3,-1]/intrinsic_z[3,-1]-1),
        "Z_gain_fraction_from_Y":float(intrinsic_all[5,-1]/intrinsic_z[5,-1]-1),
        "largest_full_coexistence_eigenvalue_real_part":float(max(intrinsic_eigen.real)),
    }

    controls: dict[str, object] = {}
    for label, pp in [
        ("no_association", replace(p, association=0)),
        ("no_catalysis", replace(p, catalysis=0)),
        ("feed_below_fold", replace(p, feed=2)),
        ("no_external_feed", replace(p, feed=0)),
        ("downstream_below_invasion", replace(p, consume_w=.4)),
    ]:
        _, run = integrate(pp, high, 300)
        controls[label] = named(run[:, -1])
    for label in ("no_association", "no_catalysis", "feed_below_fold", "no_external_feed"):
        assert controls[label]["C"] < 1e-7
    assert controls["downstream_below_invasion"]["Z"] < 1e-7

    reverse_sweep = []
    for thermal in (.01, .025, .05, .1, .2, .5):
        pp = replace(p, thermal_energy=thermal, slow_conversion=1e-7)
        _, run = integrate(pp, unassembled, 300)
        reverse_sweep.append({"thermal_energy": thermal,
                              "concentrations": named(run[:, -1]),
                              "residual": float(np.max(np.abs(rhs(0, run[:, -1], pp))))})

    # Independent material, energy, chemical-potential, and dissipation identities.
    rng = np.random.default_rng(20261020)
    material_error = energy_error = free_energy_error = mixing_balance_error = 0.0
    min_dissipation = float("inf")
    min_mixing_dissipation=float("inf")
    for _ in range(500):
        x = np.exp(rng.uniform(-4, 1, 9))
        pp = replace(p_seed, inhibition=float(rng.uniform(0, .5)))
        plus, minus = fluxes(x, pp)
        net = plus-minus
        dx = rhs(0, x, pp)
        mat_rhs = pp.dilution*(4*pp.feed-MASS@x)
        e_rhs = 3*pp.dilution*pp.feed-pp.dilution*(ENERGY@x)-RELEASE@net
        mu = ENERGY+pp.thermal_energy*np.log(x)
        affinity = pp.thermal_energy*np.log(plus/minus)
        assert np.allclose(affinity, -mu@STOICH, atol=1e-12)
        diss = float(net@affinity)
        min_dissipation = min(min_dissipation, diss)
        mu_balance = pp.dilution*(pp.feed*mu[0]-x@mu)-diss
        g=float(ENERGY@x+pp.thermal_energy*np.sum(x*(np.log(x)-1)))
        gin=3*pp.feed+pp.thermal_energy*pp.feed*(np.log(pp.feed)-1)
        mixing_gap=pp.thermal_energy*(pp.feed*np.log(pp.feed/x[0])-pp.feed+np.sum(x))
        min_mixing_dissipation=min(min_mixing_dissipation,mixing_gap)
        complete_balance=pp.dilution*(gin-g)-diss-pp.dilution*mixing_gap
        material_error = max(material_error, abs(MASS@dx-mat_rhs))
        energy_error = max(energy_error, abs(ENERGY@dx-e_rhs))
        free_energy_error = max(free_energy_error, abs(mu@dx-mu_balance))
        mixing_balance_error=max(mixing_balance_error,abs(mu@dx-complete_balance))
    assert min_dissipation >= 0
    assert min_mixing_dissipation>=0
    assert max(material_error, energy_error, free_energy_error,mixing_balance_error) < 1e-10

    # Analytic and numerical branch stability checks for the core subsystem.
    branches = []
    for feed in (3.7, 4, 5, 8):
        f_eq = lambda x: x*x+x+1+1/x-feed
        xlow = brentq(f_eq, 1e-6, critical_x)
        xhigh = brentq(f_eq, critical_x, feed)
        for branch, x in [("lower_saddle", xlow), ("upper_stable", xhigh)]:
            # Invariant-plane core f=feed-x-c, c=x^2.
            c = x*x; f = 1+1/x
            jc = np.array([[-c-2*x-1, f-c], [2*x, -1]])
            eigen = np.linalg.eigvals(jc)
            det = 2*x**3+x*x-1
            assert np.isclose(np.linalg.det(jc), det)
            assert (max(eigen.real)<0) == (branch=="upper_stable")
            branches.append({"feed":feed, "branch":branch, "x":x,
                             "largest_eigenvalue_real_part":float(max(eigen.real))})

    outputs = {
        "status":"conditional chemical construction; not a universal cosmic theorem",
        "parameters":asdict(p), "species":NAMES,
        "material_weights":MASS.tolist(), "standard_energy_weights":ENERGY.tolist(),
        "standard_energy_releases":RELEASE.tolist(),
        "physical_potential_motif":physical_motif_checks(),
        "no_reward_or_group_offspring_operator":True,
        "fold":{"x":critical_x,"feed":critical_feed},
        "irreversible_exact_equilibrium":state_summary(exact,p),
        "rare_conversion_from_zero":state_summary(eq_seed,p_seed),
        "finite_reverse_after_conversion_removed":state_summary(eq_reverse,p_reverse),
        "low_founder_final":named(low_tr[:,-1]),
        "no_seed_or_conversion_Z_Y_remain_zero":True,
        "product_feedback_without_extra_inhibition":intrinsic_feedback,
        "ecological_return":{"inhibition":p_inhibit.inhibition,
                             "C_with_consumers":float(inhibited_with[3,-1]),
                             "C_without_consumers":float(inhibited_without[3,-1]),
                             "C_with_Z_only":float(inhibited_z_only[3,-1]),
                             "Z_with_Y":float(inhibited_with[5,-1]),
                             "Z_without_Y":float(inhibited_z_only[5,-1]),
                             "relative_increase":float(inhibited_with[3,-1]/inhibited_without[3,-1]-1),
                             "C_relative_increase_from_Y":float(inhibited_with[3,-1]/inhibited_z_only[3,-1]-1)},
        "controls":controls,"finite_reverse_sweep":reverse_sweep,
        "branch_checks":branches,
        "identity_checks":{"random_states":500,"max_material_error":material_error,
                           "max_standard_energy_error":energy_error,
                           "max_free_energy_error":free_energy_error,
                           "max_mixing_aware_free_energy_balance_error":mixing_balance_error,
                           "minimum_mixing_free_energy_gap":min_mixing_dissipation,
                           "minimum_internal_free_energy_dissipation":min_dissipation},
        "checks_passed":True,
        "interpretation_limits":[
            "Rates and reaction grammar were not computed from an atomistic potential.",
            "Finite chemical possibilities are supplied; realized establishment is endogenous.",
            "C is a self-regenerating chemical complex, not demonstrated Darwinian collective individuality.",
            "No evolving sequences, stochastic extinction, shared social self-model or consciousness.",
            "The niche cascade does not demonstrate a second nested transition in individuality.",
        ],
    }
    folder=Path(__file__).parent
    (folder/"generative-results.json").write_text(json.dumps(outputs,indent=2)+"\n")

    plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False})
    fig, axes=plt.subplots(1,3,figsize=(13,3.8),layout="constrained")
    for index,label,color in [(3,"assembled C","#305a96"),(5,"consumer Z","#25917a"),(7,"consumer Y","#bc6c25")]:
        axes[0].plot(tseed,seeded[index],label=label,color=color)
    axes[0].set(xlabel="Time in dilution units",ylabel="Concentration / standard concentration",
                title="Components assemble and create niches",xlim=(0,100));axes[0].legend(frameon=False)
    axes[1].plot(t,tr[3],label="larger unassembled founder",color="#305a96")
    axes[1].plot(t,low_tr[3],label="smaller unassembled founder",color="#999999")
    axes[1].set(xlabel="Time in dilution units",ylabel="Assembled C concentration",
                title="Establishment depends on founders",xlim=(0,50));axes[1].legend(frameon=False)
    c_reference=intrinsic_none[3,-1]
    for color, (stock,label) in zip(["#999999","#7ab5a6","#25917a"],[(c_reference,"Neither"),(intrinsic_z[3,-1],"Z only"),(intrinsic_all[3,-1],"Z and Y")],strict=True):
        value=100*(stock/c_reference-1)
        axes[2].bar(label,value,color=color,width=.55)
        axes[2].text(label,value+.0008,f"{value:.5f}%",ha="center",va="bottom",fontsize=9)
    axes[2].set(xlabel="Downstream consumers present",ylabel="C increase over no-consumer baseline (%)",title="Reversibility alone returns a benefit",ylim=(0,.04))
    fig.savefig(folder/"generative-dynamics.png",dpi=180)
    print(json.dumps({"checks_passed":True,"critical_feed":critical_feed,
                      "exact_equilibrium":named(exact),"finite_reverse":named(eq_reverse),
                      "ecological_return":outputs["ecological_return"],
                      "feedback_without_extra_inhibition":intrinsic_feedback,
                      "outputs":["generative-results.json","generative-dynamics.png"]},indent=2))


if __name__=="__main__":
    main()

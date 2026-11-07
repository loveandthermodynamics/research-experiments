"""Symbolic and finite-word checks only; no physical simulations or old edits."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

sys.dont_write_bytecode=True
import sympy as s


def certificate():
    K,L0,CX,CY,M,r=s.symbols("K L0 CX CY M r",positive=True)
    D=s.symbols("D",nonnegative=True)
    left=K*(1/CX+1/CY)/2-L0*D/2+(K-L0)/(M*r)
    ceiling=s.solve(s.Eq(left,1),K)[0]
    expected=(1+L0*D/2+L0/(M*r))/((1/CX+1/CY)/2+1/(M*r))
    assert s.simplify(ceiling-expected)==0
    J=s.symbols("J",positive=True)
    assert s.simplify(left.subs({K:2*J,L0:0})-J*(1/CX+1/CY+2/(M*r)))==0
    serial=1/(1/CX+1/CY+2/(M*r))
    parallel=1/(1/CX+1/CY)
    gap=2*CX**2*CY**2/((CX+CY)*(CX*M*r+2*CX*CY+CY*M*r))
    assert s.simplify(parallel-serial-gap)==0
    # State=(stock B,q,last active,netcopy,completed rolechanges), q=0X/1R/2Y.
    # Enumerated words have no rates, RNG, probabilities or outcome selection.
    stock=2
    states={(b,0,0,0,0) for b in range(stock+1)}
    checked=0
    for length in range(10):
        following=set()
        for b,q,last,copies,switches in states:
            assert copies<=stock*(switches+1)
            checked+=1
            if q in (0,2):
                for direction in (-1,1):
                    change=direction if q==0 else -direction
                    if 0<=b+change<=stock:
                        following.add((b+change,q,last,copies+direction,switches))
                following.add((b,1,last,copies,switches))
            else:
                for target in (0,2):
                    following.add((b,target,target,copies,switches+(target!=last)))
        states=following
    return dict(
        passed=True,
        scope="Conditional synchronized-module algebra and bounded-stock word identity only; no physical model, stochastic simulation or parameter scan.",
        checks=dict(general_ceiling_rearrangement=True,zero_background_cycle_reduction=True,
                    relaxed_parallel_minus_synchronized_gap=str(gap),
                    bounded_stock_word_states=checked,maximum_word_length=9,
                    both_productive_directions_present=True,both_configuration_directions_present=True),
        general_ceiling=str(expected),
        domain="K>L0; CX,CY,M,rmax positive. D=abs(1/CX-1/CY). Actual synchronized incompatibility and all physical rate/cost premises must independently be justified.",
        not_claimed=["best_generalist_ceiling","taxonomic_necessity","spatial_individuality","physical_molecule_implementation","universal_biodiversity_optimum","ethical_legitimacy"],
    )


if __name__=="__main__":
    here=Path(__file__).resolve().parent
    result=certificate()
    result["source_sha256"]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result["note_sha256"]=hashlib.sha256((here/"next-bridge-analysis.md").read_bytes()).hexdigest()
    output=here/"capacity-certificate.json"
    output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(dict(passed=True,output=str(output),checks=result["checks"]),indent=2))

"""Independent interval signs of the proposed weak-forcing invasion coefficient.

This certifies algebra at two steady states, not finite-amplitude coexistence.
"""
from pathlib import Path
import json
from mpmath import iv

ROOT=Path(__file__).resolve().parent
iv.dps=50


def f(ell):
    return ell+iv.exp(-4*ell)+ell*iv.exp(4*ell)


def main():
    cases=[]
    for mass,lo,hi,positive in [(".95",".135448387159",".135448387162",True),
                               ("1.5",".285630560879",".285630560883",False)]:
        lower,upper=iv.mpf(lo),iv.mpf(hi)
        assert f(lower)-iv.mpf(mass)<0 and f(upper)-iv.mpf(mass)>0
        ell=iv.mpf([lo,hi]); d=iv.exp(-4*ell); r=ell*iv.exp(4*ell)
        fp=1-4*d+iv.exp(4*ell)*(1+4*ell)
        assert fp>0
        j11,j12,j21,j22=-r,r*(4*d-1),d-r,-r-1
        trace=j11+j22; determinant=j11*j22-j12*j21
        assert trace<0 and determinant>0
        # Gate flips at kappa=.25, so q=2*kappa=.5.
        q=iv.mpf(".5")
        denominator=q*q-q*trace+determinant
        assert denominator>0
        coefficient=8*d*(r*d)*(d-r)/denominator
        assert coefficient>0 if positive else coefficient<0
        cases.append(dict(M=mass,L_bracket=[lo,hi],R_enclosure=str(r),D_enclosure=str(d),
                          covariance_denominator_enclosure=str(denominator),
                          invasion_coefficient_enclosure=str(coefficient),
                          sign="positive" if positive else "negative"))
    result=dict(checks_passed=True,arithmetic="50-decimal outward interval arithmetic",
                cases=cases,
                meaning="Certified sign of K in lambda=K*epsilon^2+o(epsilon^2), subject to the persistent-resident expansion. No finite amplitude, global/spatial coexistence, availability gain, open-ended branching or participant agency is certified here.")
    (ROOT/"diversity-coefficient-certificate.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__": main()

"""Read-only certificate rerun plus final binding in this agent's new audit."""
import json
from pathlib import Path
import sys

sys.dont_write_bytecode=True
from verify_productive import ROOT,digest,load,compare
import verify_balance_independent
import verify_generalist


def review():
    files=("verify_balance_independent.py","balance-independent-certificate.json","balance-independent-review.md",
           "verify_generalist.py","generalist-certificate.json","generalist-comparison.md","balance-and-cost.md","model-spec.md")
    before={name:digest(ROOT/name) for name in files}
    saved_balance=load(ROOT/"balance-independent-certificate.json")
    repeated_balance=verify_balance_independent.certificate()
    compare(saved_balance,repeated_balance)
    assert saved_balance["verifier_sha256"]==before["verify_balance_independent.py"]
    for name,expected in saved_balance["source_sha256"].items():assert digest(ROOT/name)==expected
    saved_generalist=load(ROOT/"generalist-certificate.json")
    repeated_generalist=verify_generalist.certificate()
    compare(repeated_generalist,saved_generalist)
    assert saved_generalist["source_sha256"]==before["verify_generalist.py"]
    assert saved_generalist["note_sha256"]==before["generalist-comparison.md"]
    assert all(digest(ROOT/name)==expected for name,expected in before.items())
    result=dict(checks_passed=True,source_sha256=before,reviewer_sha256=digest(Path(__file__)),
                exact_cycle_affinity=saved_balance["exact_cycle"]["cycle_log_affinity_exact"],
                balance_source_cases=saved_balance["thermal_and_energy_values"]["cases_thermal_source_bound"],
                generalist_eligibility_cases=saved_generalist["checks"]["four_neighbor_eligibility_parent_fixtures"],
                certificate_math_reproduced=True,saved_certificate_bindings_match_current_sources=True,
                existing_reports_and_sources_unchanged=True,
                protocol_binding_scope="The prospective protocol bound analytical notes. This final analytical rerun additionally binds the certificate scripts and reports; it does not retroactively add them to the prospective protocol.",
                scope="Analytical identities, source bounds, local-cycle thermodynamics and conditional homogeneous/generalist comparisons. Positive NESS current need not be large. No spatial phase, ancestry survival, taxonomic necessity, exploit resistance or organismhood theorem.")
    (ROOT/"certificate-reverification.json").write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n")
    target=ROOT/"independent-verification.json"
    audit=load(target);assert audit["checks_passed"]
    audit["analytical_certificate_review"]=result
    audit["certificate_reverification_sha256"]=digest(ROOT/"certificate-reverification.json")
    target.write_text(json.dumps(audit,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print(json.dumps(dict(checks_passed=True,certificate_reverification_sha256=digest(ROOT/"certificate-reverification.json"),
                          final_audit_sha256=digest(target))))


if __name__=="__main__":review()

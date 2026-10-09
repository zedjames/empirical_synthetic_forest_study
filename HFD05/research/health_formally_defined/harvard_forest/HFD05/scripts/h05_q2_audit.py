"""Independent output auditor: exact fractions and a separately written interval."""
from fractions import Fraction
from statistics import NormalDist
import math
import h05_common as h

def interval(n,K):
    z=NormalDist().inv_cdf(.975);z2=z*z
    root=z*math.sqrt(n*(K-n)/K+z2/4)
    return (0 if n==0 else (n+z2/2-root)/(K+z2),1 if n==K else (n+z2/2+root)/(K+z2))
def run():
    h.guard("q2")
    rows=h.read_csv(h.ROOT/"q2_correction/finite_bank_identity.csv")
    transitions={};evaluations=0
    for r in rows:
        v,n,K=map(int,[r["nv"],r["nr"],r["K"]]);t=Fraction(r["theta"]);present=r["present"]=="True"
        assert 0<=n<=v<=K
        joint=present and Fraction(v,K)>=t and Fraction(n,K)>=t
        assert joint==(present and Fraction(n,K)>=t)==(r["simplified_point"]=="True")
        assert r["original_point"]==r["simplified_point"]
        lo,hi=interval(n,K)
        assert abs(lo-float(r["corrected_lower"]))<2e-15 and abs(hi-float(r["corrected_upper"]))<2e-15
        status="FALSE" if not present else ("TRUE" if t==0 else ("FALSE" if hi<float(t) else ("TRUE" if lo>float(t) else "MC_UNRESOLVED")))
        assert status==r["corrected_status"]
        assert float(r["corrected_width"])<=float(r["old_width"])+1e-14
        key=r["old_status"]+" -> "+status
        transitions[key]=transitions.get(key,0)+int(r["multiplicity"]);evaluations+=int(r["multiplicity"])
    for r in h.read_csv(h.ROOT/"q2_correction/primary_surface_identity.csv"):
        assert r["original_health_fraction"]==r["corrected_health_fraction"] and r["finite_identical"]=="True"
    h.write_json("mathematical_audit/q2_independent_audit.json",dict(status="PASS",count_cells=len(rows),weighted_evaluations=evaluations,
        source=h.sha(__file__),transitions=transitions,interval_formula="Independently evaluated count-form Wilson95",point_labels="EXACTLY_PRESERVED"))
    print("PASS independent Q2 count/formula audit",len(rows),evaluations,flush=True)
if __name__=="__main__":run()

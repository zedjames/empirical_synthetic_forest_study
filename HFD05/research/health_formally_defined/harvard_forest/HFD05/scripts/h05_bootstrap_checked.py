"""Direct-contraction replay of unchanged bootstrap laws and scalar reference audit."""
import itertools,json,math,warnings
from collections import defaultdict
import numpy as np
import h05_common as h

def product(weights,totals):
    values=np.einsum("ij,jk->ik",weights,totals,dtype=np.float64,optimize=False)
    if not np.isfinite(values).all():raise ValueError("Nonfinite checked bootstrap aggregation")
    return values
def units():
    rows=[]
    for n in [15,64,1024]:
        g=h.rng("fixture",31,n);weights=g.multinomial(n,np.full(n,1/n),100);totals=g.uniform(0,100,size=(n,5));totals[:,0]+=1
        reference=np.array([[math.fsum(float(w)*float(t) for w,t in zip(row,totals[:,j])) for j in range(5)] for row in weights]);direct=product(weights,totals)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always");blas=weights@totals
        error=float(abs(direct-reference).max());blas_error=float(abs(blas-reference).max())
        if error>1e-8 or blas_error>1e-8 or not np.isfinite(blas).all():raise ValueError("Numerical bootstrap product disagrees with scalar reference")
        rows.append(dict(columns=n,replicates=100,scalar_reference="math.fsum",direct_max_error=error,BLAS_max_error=blas_error,BLAS_warnings=[str(w.message) for w in caught],finite_all=True))
    h.write_json("mathematical_audit/bootstrap_numeric_units.json",dict(status="PASS",cases=rows,source_sha256=h.sha(__file__),warning_suppression_used_for_scientific_outputs=False))

def adult():
    raw=h.read_csv(h.ROOT/"external_sources/HF453/mortality_predictions.csv");rows=[]
    for r in raw:
        r=dict(r)
        for key in ["start_year","end_year","taxon","sector_proxy"]:r[key]=int(r[key])
        for key in ["observed_survival","predicted_survival","brier","negative_log_score"]:r[key]=float(r[key])
        rows.append(r)
    windows=[("CUMULATIVE",2021,end) for end in [2022,2023,2024]]+[("ANNUAL",start,start+1) for start in [2021,2022,2023]];axes=[("ALL",["ALL"]),("taxon",list(range(4))),("hemlock",["HEMLOCK","OTHER"]),("sector_proxy",list(range(4))),("size_group",["10_TO_30","30_PLUS"]),("damage",["RECORDED_DAMAGE","NO_RECORDED_DAMAGE","UNKNOWN"])];intervals=[]
    for wi,(kind,start,end) in enumerate(windows):
        for si,scheme in enumerate(["STRICT","DN_DEATH","DN_ALIVE"]):
            selected=[r for r in rows if (r["target"],r["start_year"],r["end_year"],r["scheme"])==(kind,start,end,scheme)]
            for ai,(axis,labels) in enumerate(axes):
                for li,label in enumerate(labels):
                    chosen=selected if axis=="ALL" else [r for r in selected if r[axis]==label]
                    for ui,unit in enumerate(["spatial_cluster","external_tag"]):
                        keys=sorted({r[unit] for r in chosen});by={k:i for i,k in enumerate(keys)};totals=np.zeros((len(keys),5))
                        for r in chosen:totals[by[r[unit]]]+=[1,1-r["observed_survival"],1-r["predicted_survival"],r["brier"],r["negative_log_score"]]
                        if keys:
                            g=h.rng("external_bootstrap",0,wi,si,ai,li,ui);weights=g.multinomial(len(keys),np.full(len(keys),1/len(keys)),size=1000);sample=product(weights,totals)
                            if np.any(sample[:,0]<=0):raise ValueError("Empty sampled cluster support")
                            rates=(sample[:,1]-sample[:,2])/sample[:,0];brier=sample[:,3]/sample[:,0];lo,hi=np.quantile(rates,[.025,.975]);bl,bh=np.quantile(brier,[.025,.975]);lo,hi,bl,bh=map(float,[lo,hi,bl,bh])
                        else:lo=hi=bl=bh=None
                        intervals.append(dict(target=kind,start_year=start,end_year=end,scheme=scheme,axis=axis,label=label,bootstrap_unit=unit,replicates=1000,cluster_units=len(keys),mortality_gap_lower=lo,mortality_gap_upper=hi,brier_lower=bl,brier_upper=bh,
                            scope="Conditional source-support cluster bootstrap; descriptive, not unbiased ecological population interval"))
        print("Checked adult bootstrap",kind,start,end,flush=True)
    h.write_csv("uncertainty/HF453_cluster_intervals_checked.csv",intervals)
    return compare("uncertainty/HF453_cluster_intervals.csv",intervals)

def seedlings():
    raw=h.read_csv(h.ROOT/"external_sources/HF355/transition_records.csv");rows=[]
    for r in raw:
        r=dict(r);r["observed_survival"]=int(r["observed_survival"]);r["flagged"]=r["flagged"]=="True";rows.append(r)
    # The full original taxon grid includes zero-transition taxa; keep them.
    taxa=[r["label"] for r in h.read_csv(h.ROOT/"external_sources/HF355/transition_summary.csv") if r["axis"]=="taxon"]
    panels=[("ALL","ALL",rows)]+[("taxon",code,[r for r in rows if r["taxonCode"]==code]) for code in taxa]+[("end_year",year,[r for r in rows if r["end_year"]==year]) for year in ["2018","2019","2020","2021"]]+[("UNFLAGGED_ONLY","0",[r for r in rows if not r["flagged"]])];intervals=[]
    for pi,(axis,label,selected) in enumerate(panels):
        keys=sorted({r["coordinate"] for r in selected});by={k:i for i,k in enumerate(keys)};totals=np.zeros((len(keys),2))
        for r in selected:totals[by[r["coordinate"]]]+=[1,r["observed_survival"]]
        if keys:
            g=h.rng("external_bootstrap",1,pi);weights=g.multinomial(len(keys),np.full(len(keys),1/len(keys)),1000);sample=product(weights,totals);lo,hi=np.quantile(sample[:,1]/sample[:,0],[.025,.975]);lo,hi=float(lo),float(hi)
        else:lo=hi=None
        intervals.append(dict(axis=axis,label=label,subplot_units=len(keys),replicates=1000,survival_lower=lo,survival_upper=hi,scope="Whole-subplot repeat-history bootstrap, descriptive sampled support only"))
    h.write_csv("uncertainty/HF355_subplot_intervals_checked.csv",intervals)
    return compare("uncertainty/HF355_subplot_intervals.csv",intervals)

def compare(original,checked):
    old=h.read_csv(h.ROOT/original)
    if len(old)!=len(checked):raise ValueError("Bootstrap rows changed")
    maximum=0
    for a,b in zip(old,checked):
        for key in [k for k in b if k.endswith("lower") or k.endswith("upper")]:
            if b[key] is None:
                if a[key]!="":raise ValueError("Undefined interval changed")
            else:
                if not math.isfinite(float(a[key])):raise ValueError("Initial bootstrap nonfinite")
                maximum=max(maximum,abs(float(a[key])-b[key]))
    h.write_json("mathematical_audit/"+("HF453" if "HF453" in original else "HF355")+"_bootstrap_comparison.json",dict(status="PASS" if maximum<=1e-12 else "NUMERICAL_DISCREPANCY_REQUIRES_ATTRIBUTION",original_sha256=h.sha(h.ROOT/original),checked_rows=len(checked),max_interval_endpoint_discrepancy=maximum,
        law_seeds_and_records_unchanged=True,original_preserved=True,checked_aggregation="Unoptimized einsum, independently checked scalar fsum"))
    if maximum>1e-12:raise ValueError("Initial/checked bootstrap discrepancy exceeds contract")
    print("PASS checked bootstrap",original,maximum,flush=True)

if __name__=="__main__":
    import sys
    h.guard("bootstrap_numeric");{"units":units,"adult":adult,"seedlings":seedlings}[sys.argv[1]]()

"""Exact unconditional response census over a declared finite semantic grid."""
import numpy as np
import h05_common as h
from h05_inference import wilson_array
def count(bank,lawful,model,cells):
    if len(bank)!=len(lawful) or len(bank)==0:raise ValueError("Invalid candidate bank")
    K=len(bank);cache={};horizons={c[0] for c in cells};tau_values={c[3] for c in cells};rhos={c[5] for c in cells}
    initial=bank[:,0,2];present={}
    for alpha in sorted({c[1] for c in cells}):
        for beta in sorted({c[2] for c in cells}):
            realizes=(bank[:,:,1]/model["baseline_ba"]>=alpha)&(bank[:,:,2]/model["baseline_juveniles"]>=beta)
            current=np.zeros(K,int);longest=current.copy();present[(alpha,beta)]=bool(realizes[0,0])
            if np.any(realizes[:,0]!=realizes[0,0]):raise ValueError("Conditional state varies across candidate histories")
            for year in range(bank.shape[1]):
                current=np.where(realizes[:,year],0,current+1);longest=np.maximum(longest,current)
                if year not in horizons:continue
                reserve=np.divide(bank[:,year,2],initial,out=np.zeros(K),where=initial>0)
                for tau in tau_values:
                    viable=lawful&realizes[:,year]&(longest<=tau)
                    for rho in rhos:
                        rv=viable&(reserve>=rho)
                        if np.any(rv&~viable) or np.any(viable&~lawful):raise ValueError("Reserve/lawful subset violation")
                        cache[(year,alpha,beta,tau,rho)]=(int(viable.sum()),int(rv.sum()))
    nv=np.array([cache[(c[0],c[1],c[2],c[3],c[5])][0] for c in cells],np.int64);nr=np.array([cache[(c[0],c[1],c[2],c[3],c[5])][1] for c in cells],np.int64)
    p=np.array([present[(c[1],c[2])] for c in cells]);return nv,nr,p
def infer(nv,nr,p,K,cells):
    if np.any(nr>nv) or np.any(nv>K):raise ValueError("Invalid unconditional census")
    q2=np.array([c[4]=="Q2" for c in cells]);theta=np.array([c[6] for c in cells]);n=np.where(q2,nr,nv);lo,hi=wilson_array(n,K)
    finite=p&(n/K>=theta);status=np.where(~p,0,np.where(theta==0,1,np.where(hi<theta,0,np.where(lo>theta,1,2)))).astype(np.uint8)
    return finite,status

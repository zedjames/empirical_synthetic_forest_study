"""Exact finite census and single-component pointwise Wilson inference."""
import math
from fractions import Fraction
from statistics import NormalDist

def census(nv,nr,K):
    if any(type(v)is not int for v in [nv,nr,K]) or K<=0 or not 0<=nr<=nv<=K:
        raise ValueError("Unconditional reserve subset census required")
def finite(present,nv,nr,K,theta,query):
    census(nv,nr,K);t=Fraction(str(theta))
    if query not in {"Q1","Q2"} or not 0<=t<=1:raise ValueError("Unknown licensed query/threshold")
    n=nv if query=="Q1" else nr
    return bool(present and n*t.denominator>=K*t.numerator)
def wilson(n,K,confidence=.95):
    if type(n)is not int or type(K)is not int or K<=0 or not 0<=n<=K or confidence!=.95:
        raise ValueError("Single registered Wilson95% census required")
    z=NormalDist().inv_cdf((1+confidence)/2);p=n/K;den=1+z*z/K
    center=(p+z*z/(2*K))/den
    radius=z*math.sqrt(p*(1-p)/K+z*z/(4*K*K))/den
    return (0.0 if n==0 else max(0.,center-radius),1.0 if n==K else min(1.,center+radius))
def classify(present,nv,nr,K,theta,query):
    point=finite(present,nv,nr,K,theta,query)
    lo,hi=wilson(nv if query=="Q1" else nr,K)
    status="FALSE" if not present else ("TRUE" if theta==0 else ("FALSE" if hi<theta else ("TRUE" if lo>theta else "MC_UNRESOLVED")))
    return dict(FINITE_BANK_HEALTH=point,KERNEL_MC_STATUS=status,lower=lo,upper=hi,confidence=.95,query=query)
def wilson_array(n,K):
    import numpy as np
    n=np.asarray(n)
    if K<=0 or np.any(n<0) or np.any(n>K) or np.any(n!=np.floor(n)):raise ValueError("Invalid census")
    z=NormalDist().inv_cdf(.975);p=n/K;den=1+z*z/K
    center=(p+z*z/(2*K))/den;radius=z*np.sqrt(p*(1-p)/K+z*z/(4*K*K))/den
    return np.where(n==0,0,np.maximum(0,center-radius)),np.where(n==K,1,np.minimum(1,center+radius))

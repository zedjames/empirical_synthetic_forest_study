"""Separately addressable corrected nonnegative growth, without historical mutation."""
import numpy as np
import h05_common as h

def laws(mu,sd,corrected):
    mu,sd=np.broadcast_arrays(np.asarray(mu,float),np.asarray(sd,float))
    if np.any(~np.isfinite(mu)) or np.any(~np.isfinite(sd)) or np.any(mu<0) or np.any(sd<0):raise ValueError("Invalid growth moments")
    if corrected:
        positive=(mu>0)&(sd>0)
        a=np.maximum(np.divide(mu,sd,out=np.zeros_like(mu),where=positive)**2,1e-6)
        b=np.divide(mu,a,out=np.zeros_like(mu),where=positive)
        mean=np.where(sd==0,mu,np.where(mu==0,0,a*b))
        var=np.where(positive,a*b*b,0)
    else:
        a=np.maximum((mu/np.maximum(sd,1e-9))**2,1e-6)
        b=np.where(mu>0,sd**2/np.maximum(mu,1e-9),0)
        mean=a*b;var=a*b*b
    return dict(shape=a,scale=b,mean=mean,variance=var)

def sample(mu,sd,size,generator,corrected):
    law=laws(mu,sd,corrected)
    x=generator.gamma(law["shape"],law["scale"],size=size)
    # Consume the declared Gamma draw before deterministic replacement; this
    # makes the coupling explicit, not a claim of universal shared path values.
    if corrected:x=np.where(np.asarray(sd)==0,np.asarray(mu),np.where(np.asarray(mu)==0,0,x))
    return x

class GrowthRNG:
    """Local adapter owns the changed law; original simulator globals untouched."""
    def __init__(self,generator,params,scenario,corrected=True,recorder=None):
        self.generator=generator;self.params=params;self.scenario=scenario;self.corrected=corrected;self.recorder=recorder
    def __getattr__(self,name):return getattr(self.generator,name)
    def gamma(self,shape,scale,size=None):
        length=size[1];cells=np.concatenate([np.arange(64),np.tile(np.arange(16)*4,(length-64)//16)])
        mu=self.params["growth"][cells]*self.scenario["growth_factor"];sd=self.params["sd"][cells]*self.scenario["growth_factor"]
        old=laws(mu,sd,False)
        if not np.array_equal(shape,old["shape"]) or not np.array_equal(scale,old["scale"]):raise ValueError("Actual historical call differs from audited formula")
        if self.recorder is not None:self.recorder(mu,sd,shape,scale,size)
        return sample(mu,sd,size,self.generator,self.corrected)

def simulate(state,params,scenario,years,count,generator,variant="gamma_growth",corrected=True,recorder=None):
    adapter=GrowthRNG(generator,params,scenario,corrected,recorder) if variant=="gamma_growth" else generator
    return h.original.demography.simulate(state,params,scenario,years,count,adapter,variant)

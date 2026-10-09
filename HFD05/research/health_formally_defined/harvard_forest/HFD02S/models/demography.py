"""Cohort-resolution state reconstruction and unconditioned demographic kernel.

No individual missing tree positions are generated. Growth is a cohort-level
random effect; all size/stratum detail is therefore explicitly approximate.
"""
import math
import numpy as np
from contracts import load


def parameters(model, rng):
    hazard = np.array([rng.choice(model["hazard_grid"], p=p) for p in model["hazard_posterior"]])
    growth = np.clip(rng.normal(model["growth_mean"], model["growth_se"]), 0, 0.7)
    cfg = load("config/design.json")["fit"]
    shapes = np.minimum(model["recruitment_mean"]*cfg["recruitment_exposure_years"]
                        + cfg["recruitment_pseudo_count"], cfg["parameter_effective_sample_cap"])
    # Zero observed entry is not a structural zero; positive pseudo-exposure.
    means = model["recruitment_mean"] + cfg["recruitment_pseudo_count"]/cfg["recruitment_exposure_years"]
    recruit = rng.gamma(shapes, means/shapes)
    return dict(hazard=hazard, growth=growth, sd=model["growth_sd"], recruit=recruit)


def reconstruct(core, model, rng):
    cfg = load("config/design.json")
    rc = cfg["reconstruction"]
    stats = core["stats"]
    p = parameters(model, rng)
    wet = (np.arange(64)//4)%4 == 3
    mnar_shift = rng.normal(0, rc["mnar_log_hazard_shift_sd"])
    wet_shift = rng.normal(0, rc["wet_log_hazard_shift_sd"])
    tail = np.divide(stats["tail_years"], stats["n1"], out=np.zeros(64), where=stats["n1"]>0)
    measured_n = rng.binomial(stats["n1"], np.exp(-p["hazard"]*tail))
    measured_d = np.clip(model["diameter1"]+p["growth"]*tail
                         + rng.normal(0, cfg["fit"]["measurement_dbh_sd_cm"], 64), 1, 300)
    years = np.divide(stats["unknown_dt"], stats["unknown"], out=np.zeros(64), where=stats["unknown"]>0)
    unknown_h = p["hazard"]*np.exp(mnar_shift+wet*wet_shift)
    imputed_n = rng.binomial(stats["unknown"], np.exp(-unknown_h*years))
    imputed_d = np.clip(model["diameter0"]+p["growth"]*years
                       + rng.normal(0, rc["synthetic_dbh_prior_sd_cm"], 64), 1, 300)
    unseen_n = rng.poisson(p["recruit"]*cfg["fit"]["recruitment_exposure_years"]
                          *wet*np.exp(rng.normal(0, rc["unseen_recruit_factor_log_sd"])))
    # Dated alive but unmeasured diameter: retain alive evidence and infer size.
    missing_tail = np.divide(stats["alive_missing_size_tail"],stats["alive_missing_size"],
                             out=np.zeros(64),where=stats["alive_missing_size"]>0)
    missing_n = rng.binomial(stats["alive_missing_size"],np.exp(-p["hazard"]*missing_tail))
    missing_d = np.clip(model["diameter0"]+p["growth"]*missing_tail
                         +rng.normal(0, rc["synthetic_dbh_prior_sd_cm"], 64), 1, 300)
    n = measured_n+imputed_n+unseen_n+missing_n
    # Exact aggregate basal area within the declared four-component RMS model.
    d2 = measured_n*measured_d**2+imputed_n*imputed_d**2+unseen_n*1.7**2+missing_n*missing_d**2
    branches = [0,0]
    branch_ba = 0.
    for conflict in core["conflicts"]:
        branch = int(rng.integers(2))
        branches[branch] += 1
        c = conflict["old_cell"] if branch==0 else conflict["new_cell"]
        if branch==1 and conflict["new_fate"] is not None:
            alive = bool(conflict["new_fate"]) and bool(rng.random()<math.exp(-p["hazard"][c]*conflict["tail"]))
        else:
            alive = bool(rng.random()<math.exp(-unknown_h[c]*conflict["years"]))
        diameter = conflict["new_diameter"] if branch==1 else None
        diameter = diameter if diameter is not None else conflict["old_diameter"]+p["growth"][c]*conflict["years"]
        if alive:
            n[c] += 1
            square = np.clip(diameter,1,300)**2
            d2[c] += square
            branch_ba += square*math.pi/40000
    d = np.sqrt(np.divide(d2, n, out=model["diameter0"]**2, where=n>0))
    components = dict(measured_projected_ba=float(np.sum(measured_n*measured_d**2)*math.pi/40000),
                      imputed_ba=float(np.sum(imputed_n*imputed_d**2+unseen_n*1.7**2
                                               +missing_n*missing_d**2)*math.pi/40000+branch_ba),
                      wet_imputed_stems=int((imputed_n+unseen_n+missing_n)[wet].sum()),
                      association_branches=branches,
                      mnar_shift=float(mnar_shift), wet_shift=float(wet_shift))
    return dict(n=n, d=d, components=components,ancestry=core["identity_ledger_sha256"])


def check_cohorts(n,d,identities):
    unique = len(set(identities))==len(identities)==n.shape[1]
    return unique & np.all(np.isfinite(n)&(n>=0)&(n==np.floor(n)),axis=1) & np.all(
        np.isfinite(d)&(d>=1)&(d<=300),axis=1)


def metrics(n, d, cells, deaths, recruits):
    ba = n*d*d*math.pi/40000
    return np.stack([n.sum(axis=1), ba.sum(axis=1),
                     (n*(d<10)).sum(axis=1), (ba*(cells//16==0)).sum(axis=1),
                     deaths, recruits], axis=1)


def check_history(bank):
    """Independent archived-history checker, not a viability filter."""
    finite = np.isfinite(bank).all(axis=(1,2))
    nonnegative = (bank >= 0).all(axis=(1,2))
    n, juvenile, deaths, recruits = (bank[:,:,i] for i in (0,2,4,5))
    balance = np.all(np.diff(n,axis=1) == -np.diff(deaths,axis=1)+np.diff(recruits,axis=1), axis=1)
    monotone = np.all(np.diff(deaths,axis=1)>=0,axis=1)&np.all(np.diff(recruits,axis=1)>=0,axis=1)
    return finite&nonnegative&balance&monotone&np.all(juvenile<=n,axis=1)


def simulate(state, params, scenario, years, count, rng, variant="gamma_growth"):
    width = 64+16*years
    n = np.zeros((count,width), dtype=np.int64)
    d = np.ones((count,width))
    n[:,:64] = state["n"]
    d[:,:64] = state["d"]
    cells = np.concatenate([np.arange(64), np.tile(np.arange(16)*4,years)])
    bank = np.zeros((count,years+1,6))
    deaths, recruits = np.zeros(count), np.zeros(count)
    bank[:,0] = metrics(n[:,:64], d[:,:64], cells[:64], deaths, recruits)
    kernel_lawful = np.ones(count,dtype=bool)
    for year in range(1, years+1):
        length = 64+16*(year-1)
        idx = cells[:length]
        h = params["hazard"][idx]*scenario["hazard_factor"]+(idx//16==0)*scenario["hemlock_extra_hazard"]
        surviving = rng.binomial(n[:,:length], np.exp(-h))
        deaths += (n[:,:length]-surviving).sum(axis=1)
        n[:,:length] = surviving
        mean = params["growth"][idx]*scenario["growth_factor"]
        sd = params["sd"][idx]*scenario["growth_factor"]
        if variant == "gamma_growth":
            shape = np.maximum((mean/np.maximum(sd,1e-9))**2,1e-6)
            increment = rng.gamma(shape, np.where(mean>0,sd**2/np.maximum(mean,1e-9),0),
                                  size=(count,length))
        elif variant == "normal_growth":
            increment = np.maximum(0,rng.normal(mean,sd,size=(count,length)))
        else:
            raise ValueError("Unknown generator variant")
        d[:,:length] = np.clip(d[:,:length]+increment,1,300)
        birth = rng.poisson(params["recruit"].reshape(16,4).sum(axis=1)
                            *scenario["recruitment_factor"],size=(count,16))
        n[:,length:length+16] = birth
        d[:,length:length+16] = 1.7
        recruits += birth.sum(axis=1)
        identities = [("initial",c) for c in range(64)]+[
            ("birth",t,g) for t in range(1,year+1) for g in range(16)]
        kernel_lawful &= check_cohorts(n[:,:length+16],d[:,:length+16],identities)
        bank[:,year] = metrics(n[:,:length+16],d[:,:length+16],cells[:length+16],deaths,recruits)
    return bank, kernel_lawful & check_history(bank)


def response(bank, lawful, model, horizon, alpha, beta, tau, reserve_ratio):
    history = bank[:,:horizon+1]
    structure = history[:,:,1]/model["baseline_ba"]
    regeneration = history[:,:,2]/model["baseline_juveniles"]
    realizes = (structure>=alpha)&(regeneration>=beta)
    longest, current = np.zeros(len(bank),int), np.zeros(len(bank),int)
    for year in range(horizon+1):
        current = np.where(realizes[:,year],0,current+1)
        longest = np.maximum(longest,current)
    continues = realizes[:,-1] & (longest<=tau)
    viable = lawful&continues
    reserve = np.divide(history[:,-1,2],history[:,0,2],out=np.zeros(len(bank)),where=history[:,0,2]>0)
    adequate_reserve = viable&(reserve>=reserve_ratio)
    return dict(present=bool(realizes[0,0]), viable=viable,
                reserve_viable=adequate_reserve, longest_excursion=longest,
                structure=structure[:,-1], regeneration=regeneration[:,-1], reserve=reserve)


def mass(viable, lawful, weights):
    weights = np.asarray(weights)
    if len(weights)!=len(viable) or np.any(weights<0) or not np.isclose(weights.sum(),1,atol=1e-12):
        raise ValueError("Unconditional candidate measure required")
    if np.any(viable & ~lawful):
        raise ValueError("Nonlawful continuation counted as capacity")
    return dict(candidate_mass=float(weights.sum()), unlawful_mass=float(weights[~lawful].sum()),
                lawful_failure_mass=float(weights[lawful&~viable].sum()),
                viable_mass=float(weights[viable].sum()))

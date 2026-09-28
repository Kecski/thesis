"""
Simulate wastewater measurements from a phylogeny.

Generative chain (following the general structure of wastewater
renewal/shedding models such as Lison et al.):

  phylogeny --> infection times --> daily incidence I_t
            --> individual shedding (kernel x per-person load multiplier)
            --> fecal load entering sewer, with decay + transport delay
            --> concentration = load / flow
            --> noisy measurement (lognormal, replicates, limit of detection)
            --> values at discrete sampling days

All default parameters are PLACEHOLDERS: replace them with values from
the paper / your pathogen.

Requirements: numpy, pandas, scipy, biopython (only for reading trees)
"""
from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd
from scipy import stats


# --------------------------------------------------------------------------
# 1. Phylogeny -> infection times -> daily incidence
# --------------------------------------------------------------------------
def infection_times_from_tree(newick_path: str,
                              last_tip_day: float,
                              mode: str = "nodes",
                              fmt: str = "newick") -> np.ndarray:
    """
    Extract infection times (in days) from a time-scaled phylogeny
    (branch lengths in days).

    mode = "nodes": internal nodes = transmission events (typical for
                    transmission trees / birth-death trees where a
                    branching event is an infection).
    mode = "tips" : tip times are treated as infection times.
    mode = "all"  : both.

    Time axis: the most recent tip is placed at `last_tip_day`.
    """
    from Bio import Phylo

    tree = Phylo.read(newick_path, fmt)
    depths = tree.depths()  # distance from root, requires branch lengths
    max_depth = max(depths[t] for t in tree.get_terminals())

    def abs_time(clade):
        return last_tip_day - (max_depth - depths[clade])

    times = []
    if mode in ("nodes", "all"):
        times += [abs_time(c) for c in tree.get_nonterminals()]
    if mode in ("tips", "all"):
        times += [abs_time(c) for c in tree.get_terminals()]
    return np.array(times)


def daily_incidence(infection_times: np.ndarray,
                    n_days: int,
                    t0: float = 0.0,
                    sampling_fraction: float = 1.0,
                    rng: Optional[np.random.Generator] = None) -> np.ndarray:
    """
    Bin infection times into daily counts. The phylogeny usually contains
    only a fraction of all infections, so counts are scaled up by
    1/sampling_fraction (Poisson draw to keep them integer).
    """
    rng = rng or np.random.default_rng(seed=42)
    days = np.floor(infection_times - t0).astype(int)
    days = days[(days >= 0) & (days < n_days)]
    counts = np.bincount(days, minlength=n_days).astype(float)
    if sampling_fraction < 1.0:
        counts = rng.poisson(counts / sampling_fraction).astype(float)
    return counts


# --------------------------------------------------------------------------
# 2. Shedding / transport model
# --------------------------------------------------------------------------
def discretised_gamma(mean: float, sd: float, max_days: int,
                      offset: float = 0.0) -> np.ndarray:
    """Discretised gamma pmf on days 0..max_days-1 (mass normalised to 1)."""
    shape = (mean / sd) ** 2
    scale = sd ** 2 / mean
    edges = np.arange(0, max_days + 1, dtype=float)
    cdf = stats.gamma.cdf(edges - offset, a=shape, scale=scale)
    pmf = np.diff(cdf)
    return pmf / pmf.sum()


@dataclass
class WastewaterParams:
    # --- shedding ---
    total_load_per_infection: float = 1e9    # gene copies shed in total per person
    shed_mean: float = 8.0                   # days since infection, mean of shedding profile
    shed_sd: float = 5.0
    shed_max_days: int = 40
    indiv_shape: float = 1.0                 # gamma shape of per-person load multiplier
                                             # (mean 1; small = heavy heterogeneity; np.inf = none)
    # --- sewer transport ---
    decay_fraction: float = 0.5              # fraction of load surviving to the sampler
    transport_mean: float = 0.5              # mean delay in sewer/sampling (days)
    transport_sd: float = 0.3
    # --- flow ---
    flow_mean: float = 1e8                   # litres / day at the plant
    flow_cv: float = 0.10                    # day-to-day variability of flow
    # --- measurement ---
    meas_cv: float = 0.3                     # CV of multiplicative lab noise (lognormal)
    n_replicates: int = 2
    lod: float = 1.0                         # limit of detection (copies / L); below -> censored
    sample_days_every: int = 1               # sampling interval (days)


def expected_load(incidence: np.ndarray, p: WastewaterParams,
                  rng: np.random.Generator) -> np.ndarray:
    """
    Load (gene copies/day) arriving at the sampling point.

    For each infection day s, the cohort of n_s infected people has a summed
    load multiplier  M_s = sum of n_s iid Gamma(k, 1/k)  ~ Gamma(n_s*k, 1/k)
    (mean n_s), which captures between-person variability exactly since all
    people in a cohort share the same shedding-time kernel.
    """
    n = len(incidence)
    if np.isinf(p.indiv_shape):
        mult = incidence.copy()
    else:
        mult = np.zeros(n)
        pos = incidence > 0
        mult[pos] = rng.gamma(shape=incidence[pos] * p.indiv_shape,
                              scale=1.0 / p.indiv_shape)

    shed = discretised_gamma(p.shed_mean, p.shed_sd, p.shed_max_days)
    transport = discretised_gamma(p.transport_mean, p.transport_sd, 10)
    kernel = np.convolve(shed, transport) * p.decay_fraction * p.total_load_per_infection

    load = np.convolve(mult, kernel)[:n]
    return load


def simulate_wastewater(incidence: np.ndarray,
                        p: WastewaterParams = WastewaterParams(),
                        seed: Optional[int] = None) -> pd.DataFrame:
    """Return a DataFrame with true and measured concentrations at sampling days."""
    rng = np.random.default_rng(seed)
    n = len(incidence)

    load = expected_load(incidence, p, rng)

    flow_sd = np.sqrt(np.log(1 + p.flow_cv ** 2))
    flow = p.flow_mean * rng.lognormal(-0.5 * flow_sd ** 2, flow_sd, size=n)
    conc = load / flow  # copies / L

    # multiplicative lognormal noise with mean 1
    s = np.sqrt(np.log(1 + p.meas_cv ** 2))
    noise = rng.lognormal(-0.5 * s ** 2, s, size=(n, p.n_replicates))
    meas = conc[:, None] * noise
    censored = meas < p.lod
    meas = np.where(censored, np.nan, meas)  # NaN = below LOD

    df = pd.DataFrame({
        "day": np.arange(n),
        "incidence": incidence,
        "load": load,
        "flow": flow,
        "true_concentration": conc,
    })
    for r in range(p.n_replicates):
        df[f"measurement_{r+1}"] = meas[:, r]
        df[f"below_lod_{r+1}"] = censored[:, r]

    # keep only discrete sampling days
    return df[df["day"] % p.sample_days_every == 0].reset_index(drop=True)


# --------------------------------------------------------------------------
# Demo
# --------------------------------------------------------------------------
if __name__ == "__main__":
    import sys

    n_days = 120
    rng = np.random.default_rng(1)

    if len(sys.argv) > 1:
        # python simulate_wastewater.py tree.nwk
        times = infection_times_from_tree(sys.argv[1], last_tip_day=n_days - 1,
                                          mode="nodes")
        inc = daily_incidence(times, n_days, sampling_fraction=0.05, rng=rng)
    else:
        # fallback: fake infection times from an epidemic with exponential growth
        r = 0.05
        t = np.arange(n_days)
        inc = rng.poisson(5 * np.exp(r * t) / 3).astype(float)

    df = simulate_wastewater(inc, WastewaterParams(sample_days_every=2), seed=42)
    print(df.head(10).to_string())
    df.to_csv("simulated_wastewater.csv", index=False)

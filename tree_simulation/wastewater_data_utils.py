"""
Utilities to simulate wastewater measurements from a FULL phylogeny
(every infection is a lineage in the tree).

Generative chain:

  full phylogeny --> infection times (days) --> daily incidence I_t
      --> shedding (kernel x per-person load multiplier)
      --> sewer decay + transport delay
      --> concentration = load / flow
      --> noisy lab measurement (lognormal, replicates, limit of detection)
      --> values on discrete sampling dates

All default parameter values are PLACEHOLDERS: replace them with the values
from the paper / your pathogen.

Requirements: numpy, pandas, scipy, biopython (for reading trees)
"""
from dataclasses import dataclass
from typing import Dict, Optional, Union

import numpy as np
import pandas as pd
from scipy import stats
from Bio import Phylo


# --------------------------------------------------------------------------
# 1. Phylogeny -> infection times -> daily incidence
# --------------------------------------------------------------------------
def infection_times_from_tree(tree: Union[str, "Bio.Phylo.BaseTree.Tree"],
                              root_day: float = 0.0,
                              mode: str = "transmission_nodes",
                              include_index_case: bool = True,
                              fmt: str = "newick",
                              branch_length_unit_days: float = 365.25
                              ) -> np.ndarray:
    """
    Extract one infection time per infected individual from a full,
    time-scaled phylogeny.

    Parameters
    ----------
    tree : str or Bio.Phylo tree
        Path to a tree file, or an already-loaded Bio.Phylo tree. Branch
        lengths must be time-scaled, in the unit given by
        `branch_length_unit_days`.
    root_day : float
        Day on the simulation axis assigned to the root of the tree
        (default 0, i.e. the tree start is day 0 and corresponds to
        `start_date` in `simulate_wastewater`). All other node times are
        root_day + their distance from the root, in days.
    mode : {"transmission_nodes", "tip_parent"}
        How infections are read from the tree.

        - "transmission_nodes": every internal node is a transmission event
          creating a new infection (a node with k children creates k-1 new
          infections, so a binary node creates 1). Use this for trees from a
          birth-death / transmission process where branching = transmission.
        - "tip_parent": every tip is one infected host, and its infection
          time is the time of its parent node (where its branch splits off).
    include_index_case : bool
        Only for mode="transmission_nodes". If True, one extra infection at
        the root time is added for the index case.
    fmt : str
        Tree file format understood by Bio.Phylo (default "newick";
        also e.g. "nexus", "phyloxml").
    branch_length_unit_days : float
        Number of days per branch-length unit: 365.25 if the tree is in
        years (default), 30.44 for months, 7 for weeks, 1 for days.

    Returns
    -------
    np.ndarray
        Infection times in days, one entry per infection (may be
        fractional).
    """
    if isinstance(tree, str):
        tree = Phylo.read(tree, fmt)

    depths = tree.depths()  # distance from root; needs branch lengths
    depths = {c: d * branch_length_unit_days for c, d in depths.items()}

    def abs_time(clade):
        return root_day + depths[clade]

    if mode == "transmission_nodes":
        times = []
        for node in tree.get_nonterminals():
            times += [abs_time(node)] * (len(node.clades) - 1)
        if include_index_case:
            times.append(abs_time(tree.root))
        return np.array(times)

    if mode == "tip_parent":
        parent: Dict = {}
        for node in tree.find_clades(order="level"):
            for child in node.clades:
                parent[child] = node
        return np.array([abs_time(parent[t]) for t in tree.get_terminals()
                         if t in parent])

    raise ValueError("mode must be 'transmission_nodes' or 'tip_parent'")


def daily_incidence(infection_times: np.ndarray,
                    n_days: int,
                    start_day: Optional[int] = None) -> np.ndarray:
    """
    Count infections per day (exact; no thinning, since the tree is full).

    Parameters
    ----------
    infection_times : np.ndarray
        Infection times in days, e.g. from `infection_times_from_tree`.
    n_days : int
        Length of the simulated time axis (days).
    start_day : int, optional
        Day corresponding to index 0 of the output. Infections before it
        are dropped (not shifted). Default: the floor of the earliest
        infection time. Use 0 to align with `root_day=0`.

    Returns
    -------
    np.ndarray of shape (n_days,)
        Number of new infections on each day.
    """
    if start_day is None:
        start_day = int(np.floor(infection_times.min()))
    days = np.floor(infection_times - start_day).astype(int)
    days = days[(days >= 0) & (days < n_days)]
    return np.bincount(days, minlength=n_days).astype(float)


# --------------------------------------------------------------------------
# 2. Shedding / transport model
# --------------------------------------------------------------------------
def discretised_gamma(mean: float, sd: float, max_days: int) -> np.ndarray:
    """
    Gamma probability mass function discretised to whole days.

    Parameters
    ----------
    mean : float
        Mean of the continuous gamma distribution (days).
    sd : float
        Standard deviation of the gamma distribution (days).
    max_days : int
        Number of daily bins; the pmf covers days 0 .. max_days-1
        (bin d = probability mass in [d, d+1)).

    Returns
    -------
    np.ndarray of shape (max_days,)
        Probabilities summing to 1 (mass beyond max_days is renormalised).
    """
    shape = (mean / sd) ** 2
    scale = sd ** 2 / mean
    edges = np.arange(0, max_days + 1, dtype=float)
    pmf = np.diff(stats.gamma.cdf(edges, a=shape, scale=scale))
    return pmf / pmf.sum()


@dataclass
class WastewaterParams:
    """
    Parameters of the wastewater model.

    Shedding
    --------
    total_load_per_infection : float
        Total gene copies shed by one infected person over the whole
        infection (mean over people).
    shed_mean, shed_sd : float
        Mean and sd (days) of the gamma-shaped distribution of shedding
        over days since infection.
    shed_max_days : int
        Maximum number of days after infection that shedding is modelled.
    indiv_shape : float
        Shape k of the per-person load multiplier ~ Gamma(k, scale=1/k)
        (mean 1). Small k = strong between-person heterogeneity, large k =
        nearly identical people, inf = no heterogeneity.

    Sewer transport
    ---------------
    decay_fraction : float
        Fraction of the shed load that survives (decay, losses) until it
        reaches the sampling point, between 0 and 1.
    transport_mean, transport_sd : float
        Mean and sd (days) of the gamma-distributed delay between shedding
        and sampling.

    Flow
    ----
    flow_mean : float
        Mean wastewater flow at the plant (litres per day).
    flow_cv : float
        Coefficient of variation of day-to-day flow (lognormal, mean
        `flow_mean`). Set 0 for constant flow.

    Measurement
    -----------
    meas_cv : float
        Coefficient of variation of multiplicative lognormal laboratory
        noise (mean 1). Set 0 for no noise.
    n_replicates : int
        Number of technical replicates measured per sampling day.
    lod : float
        Limit of detection (copies per litre). Replicate values below it
        are reported as NaN and flagged in `below_lod_<r>`.
    sample_days_every : int
        Sampling interval in days (1 = daily, 2 = every second day, ...),
        counted from `start_date`.
    """
    total_load_per_infection: float = 1e9
    shed_mean: float = 8.0
    shed_sd: float = 5.0
    shed_max_days: int = 21
    indiv_shape: float = 1.0
    decay_fraction: float = 0.5
    transport_mean: float = 7
    transport_sd: float = 2
    flow_mean: float = 1e8
    flow_cv: float = 0.10
    meas_cv: float = 0.3
    n_replicates: int = 1
    lod: float = 1.0
    sample_days_every: int = 1


def expected_load(incidence: np.ndarray, p: WastewaterParams,
                  rng: np.random.Generator) -> np.ndarray:
    """
    Gene-copy load (copies/day) arriving at the sampling point.

    For each infection day s, the cohort of n_s infected people has a summed
    load multiplier M_s = sum of n_s iid Gamma(k, 1/k), which is exactly
    Gamma(n_s * k, 1/k) with mean n_s. All people in a cohort share the same
    shedding-time kernel, so this reproduces between-person variability
    without simulating individuals. Load = M convolved with
    (shedding kernel * transport kernel * decay * total load).

    Parameters
    ----------
    incidence : np.ndarray of shape (n_days,)
        New infections per day.
    p : WastewaterParams
        Model parameters (see class docstring).
    rng : np.random.Generator
        Random number generator (used for individual heterogeneity).

    Returns
    -------
    np.ndarray of shape (n_days,)
        Load arriving at the sampler per day (gene copies/day).
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
    kernel = (np.convolve(shed, transport)
              * p.decay_fraction * p.total_load_per_infection)
    return np.convolve(mult, kernel)[:n]


def simulate_wastewater(incidence: np.ndarray,
                        start_date: Union[str, pd.Timestamp],
                        sequencing_delay: int=0,
                        p: Optional[WastewaterParams] = None,
                        seed: Optional[int] = None) -> pd.DataFrame:
    """
    Simulate wastewater measurements from a daily incidence series.

    Parameters
    ----------
    incidence : np.ndarray of shape (n_days,)
        Number of new infections per day (index 0 = `start_date`).
    start_date : str or pd.Timestamp
        Calendar date of the first entry of `incidence`, e.g. "2020-03-01".
        With `root_day=0`, this is the date of the tree root.
    p : WastewaterParams, optional
        Model parameters; defaults to `WastewaterParams()`.
    seed : int, optional
        Seed for reproducibility.

    Returns
    -------
    pd.DataFrame
        One row per sampling day with columns: `date`, `incidence`,
        `load` (copies/day arriving at the sampler), `flow` (L/day),
        `true_concentration` (noise-free, copies/L), and for each replicate
        r: `measurement_<r>` (copies/L, NaN if below LOD) and
        `below_lod_<r>` (bool).
    """
    p = p or WastewaterParams()
    rng = np.random.default_rng(seed)
    n = len(incidence)

    load = expected_load(incidence, p, rng)

    flow_sd = np.sqrt(np.log(1 + p.flow_cv ** 2))
    flow = p.flow_mean * rng.lognormal(-0.5 * flow_sd ** 2, flow_sd, size=n)
    conc = load / flow

    s = np.sqrt(np.log(1 + p.meas_cv ** 2))
    noise = rng.lognormal(-0.5 * s ** 2, s, size=(n, p.n_replicates))
    meas = conc[:, None] * noise
    censored = meas < p.lod
    meas = np.where(censored, np.nan, meas)

    dates = pd.date_range(start=start_date, periods=n, freq="D")
    df = pd.DataFrame({"date": dates, "incidence": incidence,
                       "load": load, "flow": flow,
                       "true_concentration": conc})
    for r in range(p.n_replicates):
        df[f"measurement_{r+1}"] = meas[:, r]
        df[f"below_lod_{r+1}"] = censored[:, r]

    keep = np.arange(n) % p.sample_days_every == 0
    return df[keep].reset_index(drop=True)


# --------------------------------------------------------------------------
# 3. Export for EpiSewer (R)
# --------------------------------------------------------------------------
def export_for_episewer(df: pd.DataFrame, output_dir: str="./") -> None:
    """
    Write EpiSewer-ready CSV files from the output of `simulate_wastewater`.

    Parameters
    ----------
    df : pd.DataFrame
        Output of `simulate_wastewater`.
    prefix : str
        Path prefix; writes `<prefix>_measurements.csv` (columns `date`,
        `concentration` in gc/mL, one row per replicate, below-LOD
        replicates dropped) and `<prefix>_flows.csv` (columns `date`,
        `flow` in mL/day).
    """
    rep_cols = [c for c in df.columns if c.startswith("measurement_")]
    meas = (df.melt(id_vars="date", value_vars=rep_cols,
                    value_name="concentration")[["date", "concentration"]]
              .dropna().sort_values("date"))
    meas["concentration"] /= 1000      # copies/L -> gc/mL

    flows = df[["date", "flow"]].copy()
    flows["flow"] *= 1000              # L/day -> mL/day

    meas.to_csv(f"{output_dir}/measurements.csv", index=False)
    flows.to_csv(f"{output_dir}/flows.csv", index=False)
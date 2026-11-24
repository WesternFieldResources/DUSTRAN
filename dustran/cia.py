"""Collision-induced absorption (CIA), e.g. O2-O2 and O2-N2, from HITRAN .cia files.

Written for DUSTRAN. HITRAN CIA files (Karman et al. 2019 and updates) are downloaded by the
user from hitran.org/cia and never bundled.

File format: blocks, each with a header line
    pair  numin  numax  npts  T[K]  max  res  comments...
followed by npts lines "wavenumber [cm-1]  k [cm^5 molecule^-2]".
Optical depth of one layer:  tau = k(nu, T) * n_a * n_b * dz   (number densities in cm^-3).

Blocks sharing the same (numin, numax) are one data set at several temperatures; we interpolate
linearly in T (clamped). Where data sets overlap, the first in the file wins.
"""
from pathlib import Path
import numpy as np


def read_cia(path):
    """List of data sets: dict(range=(lo, hi), T=array, nu=[arrays], k=[arrays])."""
    lines = Path(path).read_text().splitlines()
    sets, i = {}, 0
    while i < len(lines):
        h = lines[i].split()
        if len(h) < 5:
            i += 1
            continue
        lo, hi, n, T = float(h[1]), float(h[2]), int(h[3]), float(h[4])
        a = np.array([l.split()[:2] for l in lines[i + 1:i + 1 + n]], float)
        key = (round(lo, 1), round(hi, 1))
        s = sets.setdefault(key, dict(range=(lo, hi), T=[], nu=[], k=[]))
        s["T"].append(T); s["nu"].append(a[:, 0]); s["k"].append(np.maximum(a[:, 1], 0.0))
        i += 1 + n
    out = []
    for s in sets.values():
        o = np.argsort(s["T"])
        out.append(dict(range=s["range"], T=np.array(s["T"])[o],
                        nu=[s["nu"][j] for j in o], k=[s["k"][j] for j in o]))
    return out


def coefficient(sets, wn, T):
    """k(wn) [cm^5 molecule^-2] at temperature T, zero where no data set covers wn."""
    k = np.zeros_like(wn, dtype=float)
    done = np.zeros(len(wn), bool)
    for s in sets:
        lo, hi = s["range"]
        m = (wn >= lo) & (wn <= hi) & ~done
        if not m.any():
            continue
        ks = [np.interp(wn[m], nu, kk, left=0, right=0) for nu, kk in zip(s["nu"], s["k"])]
        Ts = s["T"]
        if len(Ts) == 1:
            k[m] = ks[0]
        else:
            t = np.clip(T, Ts[0], Ts[-1]); j = min(np.searchsorted(Ts, t), len(Ts) - 1); j = max(j, 1)
            f = (t - Ts[j - 1]) / (Ts[j] - Ts[j - 1])
            k[m] = (1 - f) * ks[j - 1] + f * ks[j]
        done |= m
    return k


def optical_depth(L, files, wn, groups, coarse_step=1.0):
    """CIA optical depth on grid `wn`, summed over layer groups (dict name -> layer mask).

    files  dict: path -> (gas_a, gas_b) using layer column names, with "N2" taken as 0.781 x air.
    Computed on a coarse grid (CIA is smooth) and interpolated to `wn`.
    """
    grid = np.arange(wn.min(), wn.max() + coarse_step, coarse_step)
    out = {g: np.zeros(len(grid)) for g in groups}
    dens = lambda g, i: (0.781 * L["air"][i] if g == "N2" else L["col"][g][i]) / L["dz"][i]
    for path, (a, b) in files.items():
        sets = read_cia(path)
        for i, T in enumerate(L["T"]):
            tau = coefficient(sets, grid, T) * dens(a, i) * dens(b, i) * L["dz"][i]
            for g, m in groups.items():
                if m[i]:
                    out[g] += tau
    return {g: np.interp(wn, grid, v) for g, v in out.items()}

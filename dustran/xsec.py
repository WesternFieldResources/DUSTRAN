"""Temperature-dependent absorption cross sections (e.g. ozone), as tables of sigma(wl, T).

Written for DUSTRAN. Default reader is for the IUP Bremen ozone data set (Serdyuchenko,
Gorshelev, Weber et al. 2014, AMT 7, 625-636), file `serdyuchenkogorshelev5digits.dat`:
45 header lines, then vacuum wavelength [nm] and sigma [cm^2 molecule^-1] at
293, 283, ..., 193 K. Downloaded by the user from the IUP Bremen site; never bundled.

Layer optical depth: tau = sigma(wl, T_layer) * column. Interpolated linearly in T (clamped).
Ozone bands are broad, so this is evaluated on a coarse grid and interpolated to the fine one.
"""
from pathlib import Path
import numpy as np

BREMEN_T = np.array([293, 283, 273, 263, 253, 243, 233, 223, 213, 203, 193], float)


def read_bremen(path, wl_min=380.0, wl_max=1100.0, step_nm=0.1):
    """(wl [nm], T [K] ascending, sigma (n_T, n_wl)), thinned to `step_nm`."""
    a = np.loadtxt(Path(path), skiprows=45)
    a = a[(a[:, 0] >= wl_min) & (a[:, 0] <= wl_max)]
    keep = np.r_[0, np.where(np.diff(np.floor(a[:, 0] / step_nm)) > 0)[0] + 1]
    a = a[keep]
    order = np.argsort(BREMEN_T)
    return a[:, 0], BREMEN_T[order], np.maximum(a[:, 1:][:, order].T, 0.0)


def sigma_at(table, T):
    wl, Ts, sig = table
    t = np.clip(T, Ts[0], Ts[-1]); j = int(np.clip(np.searchsorted(Ts, t), 1, len(Ts) - 1))
    f = (t - Ts[j - 1]) / (Ts[j] - Ts[j - 1])
    return (1 - f) * sig[j - 1] + f * sig[j]


def optical_depth(L, table, gas, wn, groups, column_scale=1.0, coarse_step=2.0):
    """Optical depth of `gas` (a layer column name, e.g. "O3") on grid `wn` [cm^-1], summed
    over layer groups. `column_scale` rescales the profile (e.g. to a measured total column)."""
    grid = np.arange(wn.min(), wn.max() + coarse_step, coarse_step)
    wl_tab = table[0]
    out = {g: np.zeros(len(grid)) for g in groups}
    for i, T in enumerate(L["T"]):
        s = np.interp(1e7 / grid, wl_tab, sigma_at(table, T), left=0, right=0)
        tau = s * L["col"][gas][i] * column_scale
        for g, m in groups.items():
            if m[i]:
                out[g] += tau
    return {g: np.interp(wn, grid, v) for g, v in out.items()}


def dobson(L, gas="O3"):
    """Total column above the surface in Dobson units (1 DU = 2.6867e16 molecules cm^-2)."""
    return L["col"][gas].sum() / 2.6867e16

"""Layered atmosphere: per-layer pressure, temperature and gas columns from a level profile.

Written for DUSTRAN from first principles.

Input is a level profile, top of atmosphere first, in the AFGL column layout
    z [km], p [hPa], T [K], air, o3, o2, h2o, co2, no2 [molecules cm^-3]
(Anderson et al. 1986, AFGL-TR-86-0110; US government work). DUSTRAN does not bundle these
files: pass the path of your own copy.

Layer columns assume each gas's number density varies exponentially with height between two
levels (the hydrostatic shape), so the column is integrated exactly for that shape:
    N = (n1 - n2) * dz / ln(n1 / n2)      (falls back to the trapezoid when n1 ~ n2 or a level is 0)
"""
from pathlib import Path
import numpy as np

AFGL_COLS = ["z", "p", "T", "air", "o3", "o2", "h2o", "co2", "no2"]

# Well-mixed gases not in the AFGL main table, as volume mixing ratios (override per scene).
DEFAULT_VMR = {"CH4": 1.9e-6, "N2O": 3.3e-7, "CO": 1.0e-7}


def load_afgl(path):
    """Dict of level arrays from an AFGL-layout text file (comment lines start with #)."""
    a = np.loadtxt(Path(path), comments="#")
    prof = {k: a[:, i] for i, k in enumerate(AFGL_COLS)}
    if prof["z"][0] < prof["z"][-1]:                     # enforce top of atmosphere first
        prof = {k: v[::-1] for k, v in prof.items()}
    return prof


def _exp_column(n1, n2, dz_cm):
    n1, n2 = np.asarray(n1, float), np.asarray(n2, float)
    trap = 0.5 * (n1 + n2) * dz_cm
    with np.errstate(divide="ignore", invalid="ignore"):
        r = n1 / n2
        exp = (n1 - n2) * dz_cm / np.log(r)
    ok = (n1 > 0) & (n2 > 0) & (np.abs(r - 1) > 1e-6)
    return np.where(ok, exp, trap)


def insert_level(prof, z_km):
    """Profile with an extra level at z_km (log-linear in pressure and densities, linear in T)."""
    z = np.asarray(prof["z"], float)
    if np.any(np.isclose(z, z_km)):
        return prof
    zz = z[::-1]
    new = {}
    for k, v in prof.items():
        v = np.asarray(v, float)
        if k == "z":
            val = z_km
        elif k == "T":
            val = np.interp(z_km, zz, v[::-1])
        else:
            val = np.exp(np.interp(z_km, zz, np.log(np.maximum(v[::-1], 1e-30))))
        new[k] = np.append(v, val)
    order = np.argsort(-new["z"])
    return {k: v[order] for k, v in new.items()}


AFGL_CO2_PPM = 330.0      # the AFGL (1986) tables assume ~330 ppm CO2; real air was ~413 ppm in 2020


def co2_ppm(year):
    """Approximate global-mean CO2 [ppm] for a year (NOAA global trend: ~339 ppm in 1980,
    ~370 in 2000, ~390 in 2010, ~412 in 2020, ~420 in 2023), linear between those anchors."""
    yrs = [1980, 1990, 2000, 2010, 2020, 2023, 2030]
    ppm = [339, 354, 369.5, 389.9, 412.2, 419.3, 436]
    return float(np.interp(year, yrs, ppm))


def layers(prof, surface_km=0.0, vmr=None, h2o_scale=1.0, sensor_km=None, co2_ppm=None):
    """Layer table from a level profile, truncated at the surface elevation.

    If `sensor_km` is given, a level is inserted there and the table gets a boolean
    `below_sensor` per layer (the layers the upward path crosses).
    `co2_ppm` rescales the profile's CO2 to that mixing ratio (AFGL tables carry ~330 ppm).

    Returns dict with, per layer (top first): p [bar], T [K], dz [cm], and per gas
    col [molecules cm^-2] and vmr (layer column / air column).
    Gases: H2O, CO2, O2, O3 from the profile; CH4, N2O, CO well mixed from `vmr`.
    `h2o_scale` multiplies the water profile (total column scales by the same factor).
    """
    vmr = {**DEFAULT_VMR, **(vmr or {})}
    prof = insert_level(prof, surface_km)
    if sensor_km is not None:
        prof = insert_level(prof, sensor_km)
    keep = np.asarray(prof["z"], float) >= surface_km - 1e-9
    P = {k: np.asarray(v, float)[keep] for k, v in prof.items()}

    dz = (P["z"][:-1] - P["z"][1:]) * 1e5
    p_bar = _exp_column(P["p"][:-1], P["p"][1:], dz) / dz / 1000.0   # pressure averaged like density
    T = 0.5 * (P["T"][:-1] + P["T"][1:])
    air = _exp_column(P["air"][:-1], P["air"][1:], dz)

    col = {}
    for g, key in (("H2O", "h2o"), ("CO2", "co2"), ("O2", "o2"), ("O3", "o3")):
        scale = h2o_scale if g == "H2O" else 1.0
        col[g] = scale * _exp_column(P[key][:-1], P[key][1:], dz)
    if co2_ppm is not None:
        col["CO2"] = col["CO2"] * co2_ppm / (col["CO2"].sum() / air.sum() * 1e6)
    for g, x in vmr.items():
        col[g] = x * air
    out = dict(p=p_bar, T=T, dz=dz, air=air, col=col, vmr={g: c / air for g, c in col.items()},
               z_top=P["z"][:-1], z_bot=P["z"][1:])
    if sensor_km is not None:
        out["below_sensor"] = P["z"][:-1] <= sensor_km + 1e-9
    return out

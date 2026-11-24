"""Line-by-line gas optical depth from HITRAN, with RADIS computing the Voigt cross sections.

RADIS is LGPL-3 and is used unmodified as an imported library. HITRAN line data are downloaded
by RADIS to the user's machine on first use (~/.radisdb); DUSTRAN does not redistribute them.

Optical depths are returned per gas, so water can be rescaled per pixel without recomputing:
    tau_total = sum_g tau_g,   tau_H2O(w) = w * tau_H2O(1)   (fixed profile shape and layer p, T)
The water continuum (absorption between lines) is NOT included yet.
"""
from pathlib import Path
import numpy as np

K_B = 1.380649e-16          # erg/K
HITRAN_GASES = ("H2O", "CO2", "CH4", "O2", "N2O", "CO")


def _factory(gas, wn_min, wn_max, wstep):
    from radis import SpectrumFactory
    sf = SpectrumFactory(wavenum_min=wn_min, wavenum_max=wn_max, molecule=gas, isotope="all",
                         wstep=wstep, truncation=25, verbose=0,
                         warnings={"MissingSelfBroadeningWarning": "ignore",
                                   "NegativeEnergiesWarning": "ignore"})
    sf.fetch_databank("hitran")
    return sf


def optical_depth(L, wl_min_nm, wl_max_nm, gases=("H2O", "CO2", "CH4", "O2"), wstep=0.002, cache=None):
    """Per-gas layer optical depths on a uniform wavenumber grid.

    L      layer table from dustran.atmosphere.layers
    wstep  grid step [cm^-1]; 0.002 resolves upper-layer Doppler widths (~0.007 cm^-1 near 2.3 um)
    Returns (wn [cm^-1], {gas: (n_layers, n_wn) optical depth}).
    """
    if cache and Path(cache).exists():
        z = np.load(cache)
        return z["wn"], {g: z[g] for g in gases}
    wn_min, wn_max = 1e7 / wl_max_nm, 1e7 / wl_min_nm
    wn = np.arange(wn_min, wn_max, wstep)
    out = {}
    for g in gases:
        sf = _factory(g, wn_min, wn_max, wstep)
        tau = np.zeros((len(L["p"]), len(wn)))
        for i, (p, T) in enumerate(zip(L["p"], L["T"])):
            if L["col"][g][i] <= 0:
                continue
            x = max(float(L["vmr"][g][i]), 1e-12)
            s = sf.eq_spectrum(Tgas=float(T), pressure=float(p), mole_fraction=x, path_length=1)
            w, k = s.get("abscoeff", wunit="cm-1")              # cm^-1 at mole fraction x
            n_gas = x * p * 1e6 / (K_B * T)                     # molecules cm^-3 (bar -> dyn cm^-2)
            tau[i] = np.interp(wn, w, k / n_gas, left=0, right=0) * L["col"][g][i]
        out[g] = tau
    if cache:
        np.savez_compressed(cache, wn=wn, **out)
    return wn, out


def column_tau(L, gas, wn_min, wn_max, groups, wstep=0.002, chunk=200.0, log=None):
    """One gas's optical depth summed over groups of layers, over a wide range, in chunks.

    groups  dict name -> boolean layer mask (e.g. all layers, layers below the sensor)
    Returns (wn, {name: float32 array}). Memory stays at one chunk per layer at a time;
    RADIS includes lines within its truncation distance of each chunk, so chunks join cleanly.
    """
    n = int(round((wn_max - wn_min) / wstep))
    wn = wn_min + wstep * np.arange(n)
    out = {k: np.zeros(n, np.float32) for k in groups}
    for c0 in np.arange(wn_min, wn_max, chunk):
        c1 = min(c0 + chunk, wn_max)
        i0, i1 = np.searchsorted(wn, [c0, c1])
        if i1 <= i0:
            continue
        try:
            sf = _factory(gas, c0 - 1.0, c1 + 1.0, wstep)
        except Exception as e:                # RADIS EmptyDatabaseError: no lines in this chunk
            if type(e).__name__ != "EmptyDatabaseError":
                raise
            if log:
                log(f"{gas} {c0:.0f}-{c1:.0f} cm-1 no lines")
            continue
        for i, (p, T) in enumerate(zip(L["p"], L["T"])):
            if L["col"][gas][i] <= 0 or not any(m[i] for m in groups.values()):
                continue
            x = max(float(L["vmr"][gas][i]), 1e-12)
            s = sf.eq_spectrum(Tgas=float(T), pressure=float(p), mole_fraction=x, path_length=1)
            w, k = s.get("abscoeff", wunit="cm-1")
            n_gas = x * p * 1e6 / (K_B * T)
            tau = np.interp(wn[i0:i1], w, k / n_gas, left=0, right=0) * L["col"][gas][i]
            for name, m in groups.items():
                if m[i]:
                    out[name][i0:i1] += tau.astype(np.float32)
        if log:
            log(f"{gas} {c0:.0f}-{c1:.0f} cm-1 done")
    return wn, out


def two_way_airmass(sza_deg, vza_deg):
    """Geometric airmass sun->ground->sensor (plane parallel; sensor above the atmosphere)."""
    return 1 / np.cos(np.radians(sza_deg)) + 1 / np.cos(np.radians(vza_deg))

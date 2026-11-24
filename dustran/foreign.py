"""Water vapour foreign continuum (water colliding with N2/O2) as a collision band under each line.

Physical picture (the same one used for quasibound water dimers in dustran/dimer.py): while a water
molecule is colliding with an air molecule, the pair absorbs like a short-lived complex. Its
spectrum follows the monomer lines, smeared by the collision lifetime into broad profiles:

    Cf(nu, T) = A * sum_j S_j(T) * f(nu - nu_j; gamma)          [cm2 molec-1 atm-1]
    k(nu)     = Cf * n_H2O * P_air,foreign [atm]                 [cm-1]
    tau       = Cf * P_air,foreign * column_H2O

    S_j(T)   HITRAN line intensities at the layer temperature (dimer.line_strength)
    f        Lorentzian, half width gamma [cm-1], cut at +-250 cm-1
    A        effective pairing constant [atm-1]: the fraction of water molecules in a collision
             complex per atm of foreign gas

A and gamma are not taken from any model. They are measured from playa scenes with independent
photometer water (experiments/foreign_continuum_fit.py) and frozen. No laboratory or field
measurement of the foreign continuum exists inside the 0.94 um band, so this measurement is the
constraint. It scales linearly with water (unlike the self-continuum, which scales with water^2).
"""
import numpy as np
from dustran import dimer

CUT = 250.0


def cross_section_unit(wn, T, lines, gamma):
    """Sum_j S_j(T) f(nu - nu_j; gamma) on grid wn [cm2 molec-1], i.e. Cf for A = 1 atm-1."""
    nu, S296, El = lines[:3]
    m = S296 > 1e-27
    nu, S = nu[m], dimer.line_strength(S296[m], El[m], nu[m], T)
    out = np.zeros(len(wn))
    for k in range(len(wn)):
        d = wn[k] - nu
        sel = np.abs(d) < CUT
        out[k] = np.sum(S[sel] * gamma / np.pi / (d[sel] ** 2 + gamma ** 2))
    return out


def optical_depth_unit(L, wn, groups, gamma, lines=None, coarse_step=1.0):
    """Foreign-continuum optical depth for A = 1 atm-1, water at the layer table's amounts."""
    lines = lines if lines is not None else dimer.water_lines()
    grid = np.arange(wn.min() - 1, wn.max() + coarse_step + 1, coarse_step)
    out = {g: np.zeros(len(grid)) for g in groups}
    kB = 1.380649e-16
    for i, T in enumerate(L["T"]):
        n_w = L["col"]["H2O"][i] / L["dz"][i]
        n_air = L["air"][i] / L["dz"][i]
        p_for_atm = (n_air - n_w) * kB * T / 1.01325e6
        tau = cross_section_unit(grid, T, lines, gamma) * p_for_atm * L["col"]["H2O"][i]
        for g, m in groups.items():
            if m[i]:
                out[g] += tau
    return {g: np.interp(wn, grid, v) for g, v in out.items()}

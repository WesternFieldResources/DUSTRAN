"""Rayleigh scattering by air: cross section, depolarisation and layer optical depth.

Written for DUSTRAN from the paper:
    Bodhaine, Wood, Dutton & Slusser (1999), "On Rayleigh optical depth calculations",
    J. Atmos. Oceanic Technol. 16, 1854-1861.

    sigma = 24 pi^3 (n^2 - 1)^2 / (lambda^4 Ns^2 (n^2 + 2)^2) * F_air
with the refractive index of dry air from Peck & Reeder (1972), scaled for CO2, and the King
factor F_air as the volume-weighted mean of N2, O2, Ar and CO2.
"""
import numpy as np

NS = 2.546899e19            # molecules cm^-3 at 288.15 K, 1013.25 hPa


def refractive_index(wl_nm, co2_ppm=400.0):
    lam = np.asarray(wl_nm, float) * 1e-3                       # um
    s2 = lam ** -2
    n300 = 1 + 1e-8 * (8060.51 + 2480990 / (132.274 - s2) + 17455.7 / (39.32957 - s2))
    c = co2_ppm * 1e-6
    return 1 + (n300 - 1) * (1 + 0.54 * (c - 0.0003))


def king_factor(wl_nm, co2_ppm=400.0):
    lam = np.asarray(wl_nm, float) * 1e-3
    s2 = lam ** -2
    f_n2 = 1.034 + 3.17e-4 * s2
    f_o2 = 1.096 + 1.385e-3 * s2 + 1.448e-4 * s2 ** 2
    c = co2_ppm * 1e-4                                          # percent by volume
    return (78.084 * f_n2 + 20.946 * f_o2 + 0.934 * 1.00 + c * 1.15) / (78.084 + 20.946 + 0.934 + c)


def cross_section(wl_nm, co2_ppm=400.0):
    """Rayleigh scattering cross section per air molecule [cm^2]."""
    lam_cm = np.asarray(wl_nm, float) * 1e-7
    n = refractive_index(wl_nm, co2_ppm)
    return 24 * np.pi ** 3 * (n ** 2 - 1) ** 2 / (lam_cm ** 4 * NS ** 2 * (n ** 2 + 2) ** 2) * king_factor(wl_nm, co2_ppm)


def depolarisation(wl_nm, co2_ppm=400.0):
    """Depolarisation ratio rho from the King factor F = (6 + 3 rho) / (6 - 7 rho)."""
    F = king_factor(wl_nm, co2_ppm)
    return 6 * (F - 1) / (3 + 7 * F)


def phase_moments(wl_nm, co2_ppm=400.0):
    """Legendre moments [1, 0, beta2] of the Rayleigh phase function P = sum beta_l P_l(cos),
    normalised so that (1/4pi) * integral P dOmega = 1 (beta2 = 0.5 without depolarisation)."""
    rho = depolarisation(wl_nm, co2_ppm)
    g = rho / (2 - rho)
    return np.array([1.0, 0.0, (1 - g) / (2 * (1 + 2 * g))])


def layer_tau(L, wl_nm, co2_ppm=400.0):
    """(n_layers,) Rayleigh optical depth of each layer = cross section x air column."""
    return cross_section(wl_nm, co2_ppm) * L["air"]

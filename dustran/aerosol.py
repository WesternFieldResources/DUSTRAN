"""Dust aerosol: optical properties from Mie theory, and a vertical profile.

Written for DUSTRAN.

Refractive index (provisional): real part 1.53; imaginary part falling from 0.0055 at 400 nm
to 0.0006 beyond 1000 nm. These are representative of mineral dust with ~1-2 % iron oxides,
within the range of published laboratory and field values (e.g. Balkanski et al. 2007, ACP 7, 81;
Di Biagio et al. 2019, ACP 19, 15503). TO REPLACE with a tabulated, cited data set
(Di Biagio 2019 is CC-BY) once downloaded.

Size: one lognormal mode. Its median radius is chosen to reproduce a measured Angstrom exponent
(440-870 nm), so the spectral shape of the optical depth follows the observation.

Vertical profile: optical depth density falls exponentially above the ground with a scale height
(default 2 km, a boundary-layer value).
"""
import numpy as np
from scipy.optimize import brentq
from dustran import mie

K_WL = np.array([350, 400, 450, 500, 550, 600, 700, 1000, 2500.0])
K_VAL = np.array([0.0075, 0.0055, 0.0040, 0.0025, 0.0018, 0.0012, 0.0008, 0.0006, 0.0006])
N_REAL = 1.53


def refractive_index(wl_nm):
    return N_REAL + 1j * np.exp(np.interp(wl_nm, K_WL, np.log(K_VAL)))


class Dust:
    def __init__(self, r_g=0.5, sigma_g=2.0, n_moments=64):
        self.r_g, self.sigma_g, self.n_moments = r_g, sigma_g, n_moments
        self._cache = {}

    def optics(self, wl_nm):
        key = round(float(wl_nm), 3)
        if key not in self._cache:
            self._cache[key] = mie.lognormal_optics(wl_nm / 1000, refractive_index(wl_nm), self.r_g,
                                                     self.sigma_g, n_moments=self.n_moments)
        return self._cache[key]

    def angstrom(self, wl1=440.0, wl2=870.0):
        e1, e2 = self.optics(wl1)["ext"], self.optics(wl2)["ext"]
        return -np.log(e1 / e2) / np.log(wl1 / wl2)

    @classmethod
    def fit_angstrom(cls, alpha, sigma_g=2.0, n_moments=64):
        """Dust whose median radius reproduces Angstrom exponent `alpha` (440-870 nm)."""
        f = lambda lr: cls(np.exp(lr), sigma_g, 4).angstrom() - alpha
        r = np.exp(brentq(f, np.log(0.02), np.log(3.0), xtol=1e-3))
        return cls(r, sigma_g, n_moments)


class Measured(Dust):
    """Aerosol from a measured volume size distribution and refractive index, e.g. an AERONET
    inversion: radii [um], dV/dlnr [um^3 um^-2] (any scale), refractive index real/imag parts at
    a few wavelengths [nm] (interpolated linearly, held constant outside)."""

    def __init__(self, radii, dvdlnr, m_wl, m_real, m_imag, n_moments=64, n_fine=150):
        self.n_moments = n_moments
        self._cache = {}
        lr = np.log(radii)
        self.lnr = np.linspace(lr[0], lr[-1], n_fine)
        dv = np.maximum(np.interp(self.lnr, lr, dvdlnr), 0)
        r = np.exp(self.lnr)
        self.weights = dv / (4 / 3 * np.pi * r ** 3)          # number per dlnr (relative)
        self.weights /= self.weights.sum()
        self.m = (np.asarray(m_wl, float), np.asarray(m_real, float), np.asarray(m_imag, float))

    def refractive_index(self, wl_nm):
        w, n, k = self.m
        return np.interp(wl_nm, w, n) + 1j * np.interp(wl_nm, w, k)

    def optics(self, wl_nm):
        key = round(float(wl_nm), 3)
        if key not in self._cache:
            self._cache[key] = mie.distribution_optics(wl_nm / 1000, self.refractive_index(wl_nm),
                                                        np.exp(self.lnr), self.weights, self.n_moments)
        return self._cache[key]


def layer_tau(L, aod_ref, ext_ratio, scale_km=2.0):
    """(n_layers,) aerosol optical depth: total `aod_ref * ext_ratio` (ext_ratio = ext(wl)/ext(ref)),
    distributed exponentially in height above the ground."""
    zg = L["z_bot"][-1]
    f = np.exp(-(L["z_bot"] - zg) / scale_km) - np.exp(-(L["z_top"] - zg) / scale_km)
    return aod_ref * ext_ratio * f / f.sum()

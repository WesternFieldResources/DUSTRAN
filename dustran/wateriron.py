"""Water or iron? Separate a water-correction error from iron absorption by shape.

A wrong water amount leaves a residue in ln(reflectance) with a fixed shape: the derivative of
the band gas transmittance with respect to water, g(band) = d ln<E T>/d w. That shape has the
whole water family in it (720, 820, 940, 1140 nm), each band at its own depth. Iron does not:
ferric iron (hematite, goethite, jarosite) is one broad dip near 860-920 nm, and ferrous iron
(olivine, pyroxene) is one broad dip near 1000-1100 nm.

Per pixel, over FIT_NM (one spectrometer, so no seam), fit

    ln rho = c0 + c1 x + c2 x^2  +  a g  +  sum_k b_k ln(rho_k)

with rho_k the iron references. A pure-water error shows up as `a`; an iron dip shows up in the
b_k and leaves `a` alone. Because ln rho_est = ln rho_true - g (w_est - w_true), the corrected
water is w_est + a. The water fit before this step (dustran.batch.fit_water) draws a straight
line from 850-880 to 1000-1040 nm: both ends sit inside iron bands, so iron biases it. This
step removes that bias.

References: USGS Spectral Library 7 (Kokaly et al. 2017, USGS Data Series 1035), convolved to
AVIRIS Classic 2014 bands by the library's authors and resampled here to the scene's band
centres. Stored under DATA/refs_iron (not in the repository).
"""
from pathlib import Path
import numpy as np

from dustran import instrument
from dustran.paths import DATA

FIT_NM = (700.0, 1260.0)
REFS = {
    "ferric: hematite GDS27": "s07_AV14_Hematite_GDS27_BECKa_AREF.txt",
    "ferric: goethite GDS134": "s07_AV14_Goethite_GDS134_ASDFRb_AREF.txt",
    "ferric: jarosite GDS99": "s07_AV14_Jarosite_GDS99_K_200C_Syn_BECKa_AREF.txt",
    "ferrous: olivine GDS70": "s07_AV14_Olivine_GDS70.c_Fo89_70um_BECKb_AREF.txt",
    "ferrous: bronzite HS9.3B": "s07_AV14_Bronzite_HS9.3B_ASDFRc_AREF.txt",
}
WAVES = "s07_AV14_Wavelengths_in_microns_224_ch_AVIRIS14a.txt"
RIDGE = 1e-6
ITER = 2


def load_refs(centers, folder=None):
    """ln(reflectance) of each iron reference at the scene's band centres: (K, B)."""
    d = Path(folder) if folder else DATA / "refs_iron"
    wl = np.loadtxt(d / WAVES, skiprows=1) * 1000.0
    out = []
    for f in REFS.values():
        r = np.loadtxt(d / f, skiprows=1)
        ok = (r > 0) & (r < 2) & np.isfinite(wl)
        o = np.argsort(wl[ok])
        out.append(np.log(np.interp(centers, wl[ok][o], r[ok][o])))
    return np.array(out)


class WaterIron:
    def __init__(self, centers, good, refs=None):
        c = np.asarray(centers, float)
        # bands of one spectrometer only: B on Classic (655-1262 nm) holds the whole fit range
        b = [0] + instrument.seams(c) + [len(c)]
        one = np.zeros(len(c), bool)                                # the spectrometer that holds 940 nm
        i940 = int(np.argmin(np.abs(c - 940)))
        for lo, hi in zip(b[:-1], b[1:]):
            if lo <= i940 < hi:
                one[lo:hi] = True
        self.sel = good & one & (c >= FIT_NM[0]) & (c <= FIT_NM[1])
        self.c = c
        x = (c[self.sel] - c[self.sel].mean()) / 200.0
        self.base = np.c_[np.ones_like(x), x, x ** 2]
        R = load_refs(c) if refs is None else refs
        Rs = R[:, self.sel]
        self.iron = (Rs - Rs.mean(1, keepdims=True)).T              # (b, K)
        self.names = list(REFS)
        self.iterations = ITER

    def fit(self, ln_rho, g):
        """ln_rho (P, B) and g (P, B) = d ln<ET>/dw at each pixel's water. Returns a (P,), b (P, K)
        and the rms log residual (P,) after the fit."""
        y = ln_rho[:, self.sel]
        gs = g[:, self.sel]
        P = len(y)
        X = np.concatenate([np.broadcast_to(self.base, (P,) + self.base.shape), gs[:, :, None],
                            np.broadcast_to(self.iron, (P,) + self.iron.shape)], 2)     # (P, b, M)
        XtX = np.einsum("pbm,pbn->pmn", X, X)
        XtX += RIDGE * np.trace(XtX, axis1=1, axis2=2)[:, None, None] * np.eye(X.shape[2])
        Xty = np.einsum("pbm,pb->pm", X, y)
        beta = np.linalg.solve(XtX, Xty[:, :, None])[:, :, 0]
        res = y - np.einsum("pbm,pm->pb", X, beta)
        return beta[:, 3], beta[:, 4:], np.sqrt(np.mean(res ** 2, 1))

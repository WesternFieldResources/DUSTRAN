"""Step 1 (go/no-go), part B: correct the playa radiance for gases only, then score the saw-tooth.

    python experiments/playa_gas_score.py TAUDIR [--plot out.png]

What it does, for the median radiance of the Railroad Valley playa square (AVIRIS-NG
ang20200712t214704, NASA L1B radiance, 100 x 100 px around the RadCalNet RVUS site):

    rho = pi * L * d^2 / (cos(sza) * band<E_sun * T_gas>)

T_gas is DUSTRAN's line-by-line two-way transmittance (sun path through the whole atmosphere,
view path from the ground up to the sensor at 5.31 km), multiplied by the TSIS solar spectrum
on the fine grid and only then averaged over each band. Scattering is NOT modelled. It is
smooth with wavelength, so it changes the overall level and slope but not the band-to-band
saw-tooth that this test measures. Free parameters: water scale w (x AFGL US standard) and
wavelength shift.

Scoring: departure from a smoothing spline through clean bands (rms %) and roughness (rms of
the 2nd difference of the spectrum, % of median) per gas window.
"""
import argparse, sys
from pathlib import Path
import numpy as np
import sys as _sys; _sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1]))
from dustran.paths import DATA
from scipy.interpolate import UnivariateSpline

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from dustran import instrument

D = DATA
SCENE = D / "ng_20200712_playa/scene.npz"        # frozen from NASA L1B radiance + obs; band table from its header
TSIS = D / "tsis/tsis_f0_0p1.txt"                  # TSIS-1 HSRS (Coddington et al. 2021), 0.1 nm
R, C, HW = 149, 541, 50
DOY = 194
GAS = {"O2 760": (750, 775), "H2O 940": (880, 1000), "H2O 1140": (1080, 1180),
       "CO2 2010/2060": (1990, 2090), "CH4 2280-2380": (2280, 2380)}


def load_scene():
    z = np.load(SCENE)
    return {k: (z[k].item() if z[k].ndim == 0 else z[k]) for k in z.files}


def score(wl, rho, bbl):
    good = (bbl > 0) & (wl > 420) & (wl < 2460)
    clean = good.copy()
    for lo, hi in GAS.values():
        clean &= ~((wl >= lo) & (wl <= hi))
    clean &= ~((wl > 1300) & (wl < 1500)) & ~((wl > 1750) & (wl < 2000))
    sm = UnivariateSpline(wl[clean], rho[clean], s=clean.sum() * (0.002 * np.median(rho[clean])) ** 2)(wl)
    res = 100 * (rho / sm - 1)
    out = {"clean": np.sqrt(np.mean(res[clean] ** 2))}
    for k, (lo, hi) in GAS.items():
        m = good & (wl >= lo) & (wl <= hi); idx = np.where(m)[0]
        out[k] = (np.sqrt(np.mean(res[m] ** 2)),
                  100 * np.sqrt(np.mean(np.diff(rho[idx], n=2) ** 2)) / np.median(rho[idx]))
    return out, res, good


class Model:
    """Gas-only forward model on one or more line-by-line segments (each a TAUDIR with
    tau_<gas>.npz on its own grid; gases missing from a segment are treated as absent there)."""

    def __init__(self, taudirs, S, cia_files=None, o3_file=None, o3_du=None, gas_scale=None):
        """gas_scale: multiply a gas's optical depth, e.g. {"CO2": 412.2 / 330} to bring the
        AFGL CO2 (~330 ppm) to the flight year (optical depth is linear in column at fixed p, T)."""
        if gas_scale is None:              # tau files hold AFGL CO2 (~330 ppm); this scene is July 2020
            from dustran.atmosphere import co2_ppm
            gas_scale = {"CO2": co2_ppm(2020) / 330.2}
        if isinstance(taudirs, (str, Path)):
            taudirs = [taudirs]
        mu_s, mu_v = np.cos(np.radians(S["sza"])), np.cos(np.radians(S["vza"]))
        t = np.loadtxt(TSIS)
        self.segs = []
        for d in taudirs:
            z = {g: np.load(f) for g in ("H2O", "CO2", "CH4", "O2")
                 if (f := Path(d) / f"tau_{g}.npz").exists()}
            wn = z["H2O"]["wn"]
            dry = sum(gas_scale.get(g, 1.0) * (z[g]["all"] / mu_s + z[g]["below"] / mu_v).astype(np.float64)
                      for g in z if g != "H2O")
            wet = (z["H2O"]["all"] / mu_s + z["H2O"]["below"] / mu_v).astype(np.float64)
            E = np.interp(1e7 / wn, t[:, 0], t[:, 1])
            self.segs.append(dict(wn=wn, E=E, Edry=E * np.exp(-dry), wet=wet))
        if len(self.segs) > 1:            # join segments into one grid so bands across a seam are covered
            self.segs.sort(key=lambda s: s["wn"][0])
            keep = [self.segs[0]]
            for s in self.segs[1:]:
                m = s["wn"] > keep[-1]["wn"][-1]
                keep.append({k: v[m] for k, v in s.items()})
            self.segs = [{k: np.concatenate([s[k] for s in keep]) for k in keep[0]}]
        if cia_files:                      # collision-induced absorption, added to the dry gases
            from dustran import atmosphere, cia
            import playa_gas_tau as T
            L = atmosphere.layers(atmosphere.load_afgl(T.AFGL), surface_km=T.GROUND_KM, sensor_km=T.SENSOR_KM)
            groups = {"all": np.ones(len(L["p"]), bool), "below": L["below_sensor"]}
            for s in self.segs:
                c = cia.optical_depth(L, cia_files, s["wn"], groups)
                s["Edry"] = s["Edry"] * np.exp(-(c["all"] / mu_s + c["below"] / mu_v))
                s["cia"] = c
        if o3_file:                        # ozone cross sections, scaled to a measured column if given
            from dustran import atmosphere, xsec
            import playa_gas_tau as T
            L = atmosphere.layers(atmosphere.load_afgl(T.AFGL), surface_km=T.GROUND_KM, sensor_km=T.SENSOR_KM)
            groups = {"all": np.ones(len(L["p"]), bool), "below": L["below_sensor"]}
            scale = (o3_du / xsec.dobson(L)) if o3_du else 1.0
            self.o3_du = xsec.dobson(L) * scale
            tab = xsec.read_bremen(o3_file)
            for s in self.segs:
                c = xsec.optical_depth(L, tab, "O3", s["wn"], groups, column_scale=scale)
                s["Edry"] = s["Edry"] * np.exp(-(c["all"] / mu_s + c["below"] / mu_v))
        self.wn, self.E, self.wet = self.segs[0]["wn"], self.segs[0]["E"], self.segs[0]["wet"]
        self.S = S
        d = 1 - 0.01672 * np.cos(np.radians(0.9856 * (S.get("doy", DOY) - 4)))   # Earth-Sun distance [AU]
        self.k = np.pi * d ** 2 / np.cos(np.radians(S["sza"]))

    def rho(self, w=1.0, shift=0.0, gases=True):
        c, f = self.S["centers"], self.S["fwhm"]
        M = np.full(len(c), np.nan)
        for s in self.segs:
            y = s["Edry"] * np.exp(-w * s["wet"]) if gases else s["E"]
            lo, hi = 1e7 / s["wn"].max(), 1e7 / s["wn"].min()
            inside = (c + shift - 2.5 * f >= lo) & (c + shift + 2.5 * f <= hi) & np.isnan(M)
            if inside.any():
                m = instrument.convolve(s["wn"], y, c, f, shift)
                M[inside] = m[inside]
        out = ~np.isfinite(M)     # bands outside every line-by-line segment: sun only, no gas correction
        if out.any():
            if not hasattr(self, "_tsis"):
                t = np.loadtxt(TSIS); self._tsis = (1e7 / t[::-1, 0], t[::-1, 1])
            M[out] = instrument.convolve(self._tsis[0], self._tsis[1], c, f, shift)[out]
        return self.k * self.S["L"] * 10.0 / M      # uW cm-2 -> mW m-2 (TSIS file units)


def fit(model, wl, bbl, windows, ws, shifts):
    best = None
    for sh in shifts:
        for w in ws:
            sc, _, _ = score(wl, model.rho(w, sh), bbl)
            cost = np.sqrt(np.mean([sc[k][0] ** 2 for k in windows]))
            if best is None or cost < best[0]:
                best = (cost, w, sh)
    return best


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("taudir")
    ap.add_argument("--plot")
    a = ap.parse_args()
    S = load_scene()
    wl = S["centers"]
    M = Model(a.taudir, S)
    ws = np.round(np.arange(0.5, 3.51, 0.05), 3)
    shifts = np.round(np.arange(-1.0, 0.51, 0.1), 2)

    runs = {"no gas correction": M.rho(gases=False)}
    c, w940, _ = fit(M, wl, S["bbl"], ["H2O 940"], ws, [0.0])
    _, w1140, _ = fit(M, wl, S["bbl"], ["H2O 1140"], ws, [0.0])
    runs[f"no shift, w fit on 940 (w={w940})"] = M.rho(w940, 0.0)
    _, wj, shj = fit(M, wl, S["bbl"], list(GAS), ws, shifts)
    runs[f"w + shift fit on all windows (w={wj}, shift={shj} nm)"] = M.rho(wj, shj)

    print(f"playa square: {S['npx']} px, sza {S['sza']:.1f}, vza {S['vza']:.1f}")
    print(f"water fitted separately: 940 -> w={w940}, 1140 -> w={w1140}  "
          f"(x AFGL US std above the playa; agreement = consistent water)")
    print("\nrms departure from smooth % / roughness %   (lower = better)")
    hdr = f"{'window':16s}"
    keys = list(runs)
    print(hdr + "".join(f" | {k[:34]:>34s}" for k in keys))
    print(f"joint fit: w={wj}, shift={shj} nm")
    scored = {k: score(wl, r, S["bbl"]) for k, r in runs.items()}
    print(f"{'clean bands':16s}" + "".join(f" | {scored[k][0]['clean']:34.2f}" for k in keys))
    for g in GAS:
        print(f"{g:16s}" +
              "".join(f" | {scored[k][0][g][0]:24.2f} / {scored[k][0][g][1]:6.2f}" for k in keys))
    print(f"\nrho at 1600 nm (no scattering correction): {runs[keys[-1]][np.argmin(abs(wl - 1600))]:.3f}")

    if a.plot:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(13, 4.5))
        for k in keys[1:]:
            _, res, good = scored[k]
            ax.plot(wl[good], res[good], ".-", ms=2, lw=.7, label="DUSTRAN gas: " + k)
        for lo, hi in GAS.values():
            ax.axvspan(lo, hi, color="orange", alpha=.15)
        ax.axhline(0, c="k", lw=.5); ax.set_ylim(-10, 10); ax.legend(fontsize=8)
        ax.set_xlabel("nm"); ax.set_ylabel("departure from smooth, %")
        fig.savefig(a.plot, dpi=110, bbox_inches="tight")

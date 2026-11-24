"""AVIRIS Classic over Ivanpah Playa vs Sentinel-2 L2A, several years.

    python experiments/ivanpah_three_way.py

Per flight: playa radiance + geometry from classic_site.py (site.npz), DUSTRAN gases (line-by-line,
one tau set for ground 0.79 km / sensor 20.15 km), O2 CIA, ozone (AFGL column), Rayleigh +
lognormal dust scaled to Sentinel-2's own AOD (Sen2Cor) with our SOS solver. Water fitted on 940
(coarse grid, no shift fit: speed). No RadCalNet here; no time-of-day transfer (Sentinel 18:2x UTC,
flights 17.0-19.4 UTC, noted per row). The ratio DUSTRAN / Sentinel per band is the brightness check.
"""
import sys
from pathlib import Path
import numpy as np
import sys as _sys; _sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1]))
from dustran.paths import DATA
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import playa_gas_tau as T
T.GROUND_KM, T.SENSOR_KM = 0.79, 20.15          # CIA/ozone layers inside Model use these
import playa_gas_score as P
from dustran import atmosphere, rayleigh, sos, aerosol

D = (DATA / "ivanpah")
XS = (DATA / "hitran_xsec")
FLIGHTS = {"f171207t01p00r07": ("s2_20171207", 2017, 341), "f180906t01p00r07": ("s2_20180908", 2018, 249),
           "f201005t01p00r09": ("s2_20201007", 2020, 279), "f210224t01p00r07": ("s2_20210224", 2021, 55)}
BANDS = {"B02": (492.4, 66), "B03": (559.8, 36), "B04": (664.6, 31), "B05": (704.1, 15), "B06": (740.5, 15),
         "B07": (782.8, 20), "B8A": (864.7, 21), "B11": (1613.7, 91), "B12": (2202.4, 175)}


def s2(dirname, b):
    a = np.fromfile(D / dirname / f"{b}.img", "<f4").reshape(202, 202)[91:111, 91:111]
    return np.median(a)


def band_avg(wl, y, c, w):
    g = np.exp(-0.5 * ((wl - c) / (w / 2.355)) ** 2); m = np.isfinite(y)
    return np.sum(g[m] * y[m]) / np.sum(g[m])


def run(f, s2dir, year, doy):
    z = np.load(D / f / "site.npz")
    wl, o = z["centers"], z["obs"]
    L = z["L"]
    bbl = ((L > 0.01) & ~((wl > 1340) & (wl < 1450)) & ~((wl > 1790) & (wl < 1960)) & (wl > 395)).astype(float)
    S = dict(L=L, centers=wl, fwhm=z["fwhm"], sza=o[4], vza=o[2], saz=o[3], vaz=o[1], bbl=bbl, doy=doy)
    M = P.Model([D / "tau_swir", D / "tau_vnir"], S,
                cia_files={XS / "O2-O2_2024.cia": ("O2", "O2"), XS / "O2-N2_2024.cia": ("O2", "N2")},
                o3_file=str(DATA / "bremen_o3/serdyuchenkogorshelev5digits.dat"),
                gas_scale={"CO2": atmosphere.co2_ppm(year) / 330.2})
    prof = atmosphere.load_afgl(T.AFGL)
    cm1 = atmosphere.layers(prof, 0.79)["col"]["H2O"].sum() / 3.3428e22
    good = bbl > 0

    def cost(w):
        q = M.rho(w)
        inw = good & (wl >= 880) & (wl <= 1000)
        e = good & (((wl >= 850) & (wl <= 880)) | ((wl >= 1000) & (wl <= 1040)))
        b = np.polyfit(wl[e], q[e], 1)
        return np.sqrt(np.mean((q[inw] / np.polyval(b, wl[inw]) - 1) ** 2))
    ws = np.arange(0.05, 2.0, 0.05)
    c = [cost(w) for w in ws]; i = int(np.argmin(c))
    if 0 < i < len(ws) - 1:                     # parabola through the three best points
        a, b_, _ = np.polyfit(ws[i - 1:i + 2], c[i - 1:i + 2], 2); w = -b_ / (2 * a)
    else:
        w = ws[i]
    r_app = M.rho(w)

    Ly = atmosphere.layers(prof, surface_km=0.79, sensor_km=float(z["sensor_km"]))
    ks = int(np.argmax(Ly["below_sensor"]))
    mu0, muv = np.cos(np.radians(S["sza"])), np.cos(np.radians(S["vza"]))
    dphi = S["vaz"] - (S["saz"] + 180.0)
    aod = s2(s2dir, "AOT") / 1000.0
    dust = aerosol.Dust(r_g=0.08, sigma_g=1.8)
    grid = np.r_[np.arange(400, 1000, 50.0), np.arange(1000, 2501, 250.0)]
    terms = []
    for wv in grid:
        op = dust.optics(wv)
        tR = rayleigh.layer_tau(Ly, wv); bR = np.zeros(len(op["beta"])); bR[:3] = rayleigh.phase_moments(wv)
        tA = aerosol.layer_tau(Ly, aod, op["ext"] / dust.optics(550)["ext"])
        tau = tR + tA; sca = tR + op["ssa"] * tA
        beta = (tR[:, None] * bR + (op["ssa"] * tA)[:, None] * op["beta"]) / sca[:, None]
        terms.append(sos.surface_terms(tau, sca / tau, beta, mu0, muv, dphi, ks, n_streams=8))
    it = {k: np.interp(wl, grid, np.array([t[k] for t in terms])) for k in terms[0]}
    x = r_app - it["R_path"]
    rho = x / (it["T_dn"] * it["T_up"] + it["S"] * x)
    np.save(D / f / "dustran_rho.npy", np.vstack([wl, rho]))
    o_ = np.argsort(wl)
    rat = {b: band_avg(wl[o_], np.where(good, rho, np.nan)[o_], cw, fw) / ((s2(s2dir, b) - 1000) / 1e4)
           for b, (cw, fw) in BANDS.items()}
    return w * cm1, s2(s2dir, "WVP") / 1000, aod, o[9], S["sza"], rat


if __name__ == "__main__":
    print(f"{'flight':18s} {'UTC':>5s} {'sza':>4s} {'H2O fit/S2':>11s} {'AOD':>5s}   DUSTRAN / Sentinel-2: "
          + " ".join(f"{b:>5s}" for b in BANDS))
    out = {}
    for f, (s2dir, year, doy) in FLIGHTS.items():
        h2o, h2o_s2, aod, utc, sza, rat = run(f, s2dir, year, doy)
        out[f] = [rat[b] for b in BANDS]
        print(f"{f:18s} {utc:5.2f} {sza:4.0f} {h2o:5.2f}/{h2o_s2:4.2f} {aod:5.3f}   "
              + " ".join(f"{rat[b]:5.3f}" for b in BANDS), flush=True)
    np.savez(D / "three_way_ivanpah.npz", bands=list(BANDS), **out)

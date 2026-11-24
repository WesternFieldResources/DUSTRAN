"""Playa with Rayleigh (+ optional aerosol) from DUSTRAN's own solver, vs RadCalNet in absolute terms.

    python experiments/playa_rayleigh.py TAUDIRS RCN_INPUT_FILE [--plot out.png] [--aod 0.053 0.234]

Gases: line by line + O2 CIA + ozone (as playa_radcalnet.py --cia --o3).
Scattering: dustran/sos.py on a coarse wavelength grid, interpolated to the bands.
--aod: aerosol optical depths at 550 nm to try (e.g. RadCalNet's measured value and a high alternative). Aerosol size is
fitted to RadCalNet's Angstrom exponent (dustran/aerosol.py); 16 streams with delta-M.
Surface reflectance from
    r = R_path + T_dn T_up rho / (1 - S rho)   ->   rho = (r - R_path) / (T_dn T_up + S (r - R_path))
where r is the gas-corrected apparent reflectance (gas transmittance divided out of everything; Chandrasekhar's surface-coupling form).
"""
import argparse, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import playa_gas_score as P
import playa_gas_tau as T
from dustran import atmosphere, rayleigh, sos, aerosol

ap = argparse.ArgumentParser(); ap.add_argument("taudirs"); ap.add_argument("rcn"); ap.add_argument("--plot")
ap.add_argument("--aod", nargs="*", type=float, default=[])
ap.add_argument("--aeronet", action="store_true", help="add a run with AERONET-measured aerosol (sizes, index, AOD)")
a = ap.parse_args()
rows = [l.split("\t") for l in open(a.rcn).read().splitlines()]
head = {}
for r in rows:
    if r and r[0].strip() and not r[0].strip().isdigit():
        head.setdefault(r[0].rstrip(":"), r[1:])
col = head["UTC"].index("21:30")
spec = np.array([(float(r[0]), float(r[1 + col])) for r in rows if r and r[0].strip().isdigit()][:191])
spec = spec[spec[:, 1] < 9]; wr, rr = spec.T

S = P.load_scene(); wl = S["centers"]
X = Path("D:/workspace/dustran_data/hitran_xsec")
M = P.Model(a.taudirs.split(","), S, cia_files={X / "O2-O2_2024.cia": ("O2", "O2"), X / "O2-N2_2024.cia": ("O2", "N2")},
            o3_file="D:/workspace/dustran_data/bremen_o3/serdyuchenkogorshelev5digits.dat",
            o3_du=float(head["O3"][col]))
r_app = M.rho(2.15, 0.389 + 0.063e-3 * (wl - 1500))

prof = atmosphere.load_afgl(T.AFGL)
L = atmosphere.layers(prof, surface_km=T.GROUND_KM, sensor_km=T.SENSOR_KM)
p_afgl = np.exp(np.interp(T.GROUND_KM, prof["z"][::-1], np.log(prof["p"][::-1])))
L["air"] = L["air"] * float(head["P"][col]) / p_afgl       # scale air column to the measured surface pressure
ks = int(np.argmax(L["below_sensor"]))
mu0, muv = np.cos(np.radians(S["sza"])), np.cos(np.radians(S["vza"]))
dphi = S["vaz"] - (S["saz"] + 180.0)
grid = np.r_[np.arange(400, 1000, 25.0), np.arange(1000, 2501, 150.0)]
dust = aerosol.Dust.fit_angstrom(float(head["Ang"][col])) if a.aod else None
if dust:
    print(f"aerosol: lognormal r_g {dust.r_g:.3f} um, sigma_g {dust.sigma_g}, fitted to Angstrom {head['Ang'][col]}")


def solve(aod, aer=None, ref=550):
    aer = aer or dust
    t0 = time.time(); terms = []
    for w in grid:
        tR = rayleigh.layer_tau(L, w)
        if aod:
            o = aer.optics(w)
            bR = np.zeros(len(o["beta"])); bR[:3] = rayleigh.phase_moments(w)
            tA = aerosol.layer_tau(L, aod, o["ext"] / aer.optics(ref)["ext"])
            tau = tR + tA; sca = tR + o["ssa"] * tA
            beta = (tR[:, None] * bR + (o["ssa"] * tA)[:, None] * o["beta"]) / sca[:, None]
            ssa = sca / tau; n = 16
        else:
            tau, ssa, beta, n = tR, np.ones_like(tR), np.tile(rayleigh.phase_moments(w), (len(tR), 1)), 8
        terms.append(sos.surface_terms(tau, ssa, beta, mu0, muv, dphi, ks, n_streams=n))
    tab = {k: np.array([t[k] for t in terms]) for k in terms[0]}
    it = {k: np.interp(wl, grid, v) for k, v in tab.items()}
    x = r_app - it["R_path"]
    print(f"AOD {aod}: {len(grid)} wavelengths in {time.time() - t0:.0f}s; at 450 nm R_path "
          f"{np.interp(450, grid, tab['R_path']):.4f} T_dn {np.interp(450, grid, tab['T_dn']):.4f} "
          f"T_up {np.interp(450, grid, tab['T_up']):.4f} S {np.interp(450, grid, tab['S']):.4f}")
    return x / (it["T_dn"] * it["T_up"] + it["S"] * x)


runs = {"gas only": r_app, "Rayleigh": solve(0.0)}
for aod in a.aod:
    runs[f"Ray+aer {aod}"] = solve(aod)
if a.aeronet:
    import aeronet_aerosol
    aer, meas = aeronet_aerosol.load()
    runs[f"AERONET {meas[500]:.3f}"] = solve(meas[500], aer, 500)

np.savez("D:/workspace/dustran_data/ng_2020_playa_runs.npz", wl=wl, **{k.replace(" ", "_"): v for k, v in runs.items()})
good = S["bbl"] > 0
to10 = lambda y: np.array([np.nanmean(np.where(good & (abs(wl - c) <= 5), y, np.nan)) for c in wr])
ok = np.isfinite(to10(r_app))
print("\nabsolute ratio to RadCalNet ground (mean / spread %):")
print(f"{'range':12s} " + " ".join(f"{k:>15s}" for k in runs))
for lo, hi in ((400, 450), (450, 500), (500, 700), (700, 1000), (1000, 1300), (1500, 1750), (2000, 2300)):
    m = ok & (wr >= lo) & (wr <= hi)
    cells = []
    for y in runs.values():
        q = to10(y)[m] / rr[m]; cells.append(f"{q.mean():.3f} / {100 * q.std():.1f}")
    print(f"{lo}-{hi:<7d} " + " ".join(f"{c:>15s}" for c in cells))

if a.plot:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(12, 4.5))
    for lab, y in runs.items():
        ax.plot(wr[ok], to10(y)[ok] / rr[ok], label=lab)
    ax.axhline(1, c="gray"); ax.set_ylim(0.85, 1.2); ax.legend(); ax.set_xlabel("nm"); ax.set_ylabel("ratio to RadCalNet ground")
    fig.savefig(a.plot, dpi=110, bbox_inches="tight")

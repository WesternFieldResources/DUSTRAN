"""AVIRIS Classic over the Railroad Valley playa vs RadCalNet (second validation day).

    python experiments/classic_playa.py

Scene: f190813t01p00r07 (2019-08-13, ~22:13 UTC, ER-2, sensor ~20.5 km, 15.2 m pixels).
Playa pixels: ortho pixels whose raw pixel lies within +-12 raw pixels of the RadCalNet RVUS
point (found from the IGM/GLT; ~370 m square). Ground truth: RadCalNet RVUS 2019 day 225 at
22:00 UTC; atmosphere (water, AOD, Angstrom, pressure, ozone) from the same RadCalNet record.

Same chain as the 2020 AVIRIS-NG run: line-by-line gases (CO2 scaled to 2019), O2 CIA, ozone
(Bremen), Rayleigh + fine aerosol (lognormal fitted to the Angstrom exponent) with our SOS
solver. Water and a wavelength shift are fitted; everything else is from measurements.
"""
import sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import playa_gas_score as P
import playa_gas_tau as T
from dustran import atmosphere, rayleigh, sos, aerosol

D = Path("D:/workspace/dustran_data")
SC = D / "classic_f190813t01p00r07"
RDN = SC / "f190813t01p00r07rdn_e"
RCN = D / "radcalnet/RVUS00_2019_225_v00.06.input"
GROUND, SENSOR, DOY, YEAR, UTC = 1.435, 20.53, 225, 2019, 22.2208   # flight time, decimal UTC hours


def radcalnet(path, utc):
    """RadCalNet record at `utc`: "HH:MM" picks that 30-min column; a float (decimal hours)
    interpolates linearly between the two neighbouring columns (the playa's reflectance changes
    ~1 % per 15 min in the afternoon as the sun drops)."""
    rows = [l.split("\t") for l in open(path).read().splitlines()]
    head = {}
    for r in rows:
        if r and r[0].strip() and not r[0].strip().isdigit():
            head.setdefault(r[0].rstrip(":"), r[1:])
    hours = np.array([int(t[:2]) + int(t[3:5]) / 60 for t in head["UTC"]])
    if isinstance(utc, str):
        cols, wts = [head["UTC"].index(utc)], [1.0]
    else:
        j = int(np.clip(np.searchsorted(hours, utc), 1, len(hours) - 1))
        f = (utc - hours[j - 1]) / (hours[j] - hours[j - 1])
        cols, wts = [j - 1, j], [1 - f, f]
    first, second = {}, {}                  # second block repeats wavelengths as 1-sigma uncertainties
    for r in rows:
        if r and r[0].strip().isdigit():
            vals = [float(r[1 + c]) for c in cols]
            v = 9999.0 if any(x > 9 for x in vals) else sum(w * x for w, x in zip(wts, vals))
            (second if int(r[0]) in first else first)[int(r[0])] = v
    w = np.array(sorted(first)); v = np.array([first[k] for k in w]); u = np.array([second.get(k, 9999) for k in w])
    ok = (v < 9) & (u < 9)
    atm = {k: sum(wt * float(head[k][c]) for wt, c in zip(wts, cols)) for k in ("WV", "AOD", "Ang", "P", "O3")}
    return atm, w[ok].astype(float), v[ok], u[ok]


def scene():
    spc = np.loadtxt(RDN / "f190813t01p00r07rdn_e.spc")
    L = np.load(SC / "playa_rdn_median.npy")
    obs = np.memmap(RDN / "f190813t01p00r07rdn_e_obs_ort", "<f4", "r", shape=(12251, 10, 1126))
    oy, ox = np.load(SC / "playa_ortho_idx.npy")
    o = np.median(np.stack([np.asarray(obs[y, :, x]) for y, x in zip(oy, ox)]), 0)
    bbl = ((L > 0.01) & ~((spc[:, 0] > 1340) & (spc[:, 0] < 1450)) & ~((spc[:, 0] > 1790) & (spc[:, 0] < 1960))
           & (spc[:, 0] > 395)).astype(float)
    if "--gain" in sys.argv:                 # smooth vicarious gain on radiance (classic_gain.py)
        gw, gg = np.load(SC / "classic_2019_gain.npy")
        L = L * np.interp(spc[:, 0], gw, gg)
    return dict(L=L, centers=spc[:, 0], fwhm=spc[:, 1], sza=o[4], vza=o[2], saz=o[3], vaz=o[1],
                bbl=bbl, npx=len(oy), doy=DOY)


if __name__ == "__main__":
    atm, wr, rr, ur = radcalnet(RCN, UTC)
    print("RadCalNet", UTC, atm)
    S = scene(); wl = S["centers"]
    print(f"playa: {S['npx']} px, sza {S['sza']:.1f}, vza {S['vza']:.1f}")
    XS = Path("D:/workspace/dustran_data/hitran_xsec")
    M = P.Model([SC / "tau_swir", SC / "tau_vnir"], S,
                cia_files={XS / "O2-O2_2024.cia": ("O2", "O2"), XS / "O2-N2_2024.cia": ("O2", "N2")},
                o3_file="D:/workspace/dustran_data/bremen_o3/serdyuchenkogorshelev5digits.dat",
                o3_du=atm["O3"], gas_scale={"CO2": atmosphere.co2_ppm(YEAR) / 330.2})
    prof = atmosphere.load_afgl(T.AFGL)
    cm1 = atmosphere.layers(prof, GROUND)["col"]["H2O"].sum() / 3.3428e22

    # water + shift: fit water on 940 with the RadCalNet surface as baseline (as playa_water.py)
    ground = np.interp(wl, wr, rr)
    good = S["bbl"] > 0
    def water_cost(w, shift, lo=880, hi=1000, e=((850, 880), (1000, 1040))):
        q = M.rho(w, shift) / ground
        inw = good & (wl >= lo) & (wl <= hi)
        edge = good & (((wl >= e[0][0]) & (wl <= e[0][1])) | ((wl >= e[1][0]) & (wl <= e[1][1])))
        b = np.polyfit(wl[edge], q[edge], 1)
        return np.sqrt(np.mean((q[inw] / np.polyval(b, wl[inw]) - 1) ** 2))
    best = min(((water_cost(w, s), w, s) for s in np.arange(-1.5, 1.51, 0.25) for w in np.arange(0.6, 2.6, 0.05)))
    _, w, sh = best
    w1140 = min((water_cost(x, sh, 1080, 1180, ((1040, 1075), (1230, 1260))), x) for x in np.arange(0.6, 2.6, 0.05))[1]
    print(f"fitted: water {w * cm1:.3f} cm (940), {w1140 * cm1:.3f} cm (1140); RadCalNet {atm['WV']} cm; shift {sh:+.2f} nm")
    r_app = M.rho(w, sh)

    L = atmosphere.layers(prof, surface_km=GROUND, sensor_km=SENSOR)
    p_afgl = np.exp(np.interp(GROUND, prof["z"][::-1], np.log(prof["p"][::-1])))
    L["air"] = L["air"] * atm["P"] / p_afgl
    ks = int(np.argmax(L["below_sensor"]))
    mu0, muv = np.cos(np.radians(S["sza"])), np.cos(np.radians(S["vza"]))
    dphi = S["vaz"] - (S["saz"] + 180.0)
    import aeronet_aerosol                 # AERONET Railroad_Valley2: inversion 21:50, direct sun 22:12:30
    dust, meas = aeronet_aerosol.load("21:50:02", "22:12:30", "13:08:2019", D / "aeronet_2019_225",
                                      "aod_Railroad_Valley2.txt")
    print(f"AERONET at 22:12:30: AOD500 {meas[500]:.3f}, water {meas['PW']:.3f} cm")
    grid = np.r_[np.arange(400, 1000, 25.0), np.arange(1000, 2501, 150.0)]
    t0 = time.time(); terms = []
    for wv in grid:
        o = dust.optics(wv)
        tR = rayleigh.layer_tau(L, wv); bR = np.zeros(len(o["beta"])); bR[:3] = rayleigh.phase_moments(wv)
        tA = aerosol.layer_tau(L, meas[500], o["ext"] / dust.optics(500)["ext"])
        tau = tR + tA; sca = tR + o["ssa"] * tA
        beta = (tR[:, None] * bR + (o["ssa"] * tA)[:, None] * o["beta"]) / sca[:, None]
        terms.append(sos.surface_terms(tau, sca / tau, beta, mu0, muv, dphi, ks, n_streams=16))
    tab = {k: np.array([t[k] for t in terms]) for k in terms[0]}
    it = {k: np.interp(wl, grid, v) for k, v in tab.items()}
    x = r_app - it["R_path"]
    rho = x / (it["T_dn"] * it["T_up"] + it["S"] * x)
    # scene-wide terms, so any pixel's radiance can be turned into reflectance in closed form
    np.savez(SC / ("scene_terms_gain.npz" if "--gain" in sys.argv else "scene_terms.npz"), wl=wl,
             r_per_L=r_app / S["L"], R_path=it["R_path"], T=it["T_dn"] * it["T_up"], S=it["S"])
    print(f"scattering solved in {time.time() - t0:.0f}s")

    order = np.argsort(wl)
    to10 = lambda y: np.array([np.nanmean(np.where(good & (abs(wl - c) <= 5), y, np.nan)) for c in wr])
    q = to10(rho) / rr; qa = to10(r_app) / rr
    print("\nratio to RadCalNet ground (mean / spread %), RadCalNet 1-sigma:")
    for lo, hi in ((400, 450), (450, 500), (500, 700), (700, 1000), (1000, 1300), (1500, 1750), (2000, 2300)):
        m = np.isfinite(q) & (wr >= lo) & (wr <= hi)
        print(f"  {lo}-{hi:<5d} gas only {np.mean(qa[m]):.3f}   full {np.mean(q[m]):.3f} / {100 * np.std(q[m]):.1f}%"
              f"   (RadCalNet +-{100 * np.mean(ur[m] / rr[m]):.1f}%)")
    np.save(SC / ("dustran_playa_rho_gain.npy" if "--gain" in sys.argv else "dustran_playa_rho.npy"), np.vstack([wl, rho]))

    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(12, 4.5))
    ok = np.isfinite(q)
    ax.fill_between(wr, 1 - ur / rr, 1 + ur / rr, color="gray", alpha=.2, label="RadCalNet 1-sigma")
    ax.plot(wr[ok], q[ok], label="DUSTRAN / RadCalNet (AVIRIS Classic 2019-08-13)")
    ax.axhline(1, c="k", lw=.5); ax.set_ylim(0.85, 1.15); ax.legend(); ax.set_xlabel("nm"); ax.set_ylabel("ratio")
    fig.savefig(Path(__file__).resolve().parents[1] / f"docs/classic_2019_playa{'_gain' if '--gain' in sys.argv else ''}.png", dpi=110, bbox_inches="tight")

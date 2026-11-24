"""Three independent looks at the Railroad Valley playa, 2019: Sentinel-2, RadCalNet, AVIRIS Classic.

    python experiments/three_way.py

No source is treated as truth. Each has its own calibration chain:
  Sentinel-2A L2A (ESA Sen2Cor, N0500 reprocessing), 2019-08-12 18:29 UTC, 20 m bands
  RadCalNet RVUS (ground instruments + model), 2019-08-13, interpolated to the time needed
  AVIRIS Classic f190813t01p00r07, 2019-08-13 22:13 UTC: DUSTRAN and NASA reflectance

Time of day: the playa's reflectance changes ~10 % between late morning and mid-afternoon (sun
angle). AVIRIS is moved to Sentinel's time of day with RadCalNet's own diurnal curve,
k(wl) = RCN(18:29) / RCN(22:13). That uses RadCalNet's *relative* shape over the day, not its
absolute level. Day offset: Sentinel-2 is from the day before (no Aug 13 pass); RadCalNet day 224
was not available, so day 225 at 18:29 stands in for it (the playa's change overnight is assumed small).

Band responses: Gaussians at the S2A centre wavelengths and widths (the official spectral
response functions aren't downloaded; fine for a smooth playa spectrum).
"""
import sys
from pathlib import Path
import numpy as np
import sys as _sys; _sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1]))
from dustran.paths import DATA
sys.path.insert(0, str(Path(__file__).resolve().parent))
import classic_playa as C

DAY = "2020" if "--2020" in sys.argv else "2019"
S2 = Path(f"{DATA}/{'s2_20200712' if DAY == '2020' else 's2_20190812'}/crop")
BANDS = {"B02": (492.4, 66), "B03": (559.8, 36), "B04": (664.6, 31), "B05": (704.1, 15), "B06": (740.5, 15),
         "B07": (782.8, 20), "B8A": (864.7, 21), "B11": (1613.7, 91), "B12": (2202.4, 175)}
T_S2, T_AV = (18 + 29.3 / 60, "21:30") if DAY == "2020" else (18 + 29 / 60, 22.2208)   # 2020: RadCalNet 22:00 missing; flight 21:50
HW = 9                                   # +-9 px x 20 m = ~370 m box, like the AVIRIS box


def s2_band(b, hw=HW):
    a = np.fromfile(S2 / f"{b}.img", "<f4").reshape(201, 201)
    box = a[100 - hw:100 + hw + 1, 100 - hw:100 + hw + 1]
    return (np.median(box) - 1000) / 10000, np.std(box) / 10000


def band_avg(wl, y, c, w):
    g = np.exp(-0.5 * ((wl - c) / (w / 2.355)) ** 2)
    m = np.isfinite(y)
    return np.sum(g[m] * y[m]) / np.sum(g[m])


if __name__ == "__main__":
    scl = np.fromfile(S2 / "SCL.img", "<f4").reshape(201, 201)[100 - HW:100 + HW + 1, 100 - HW:100 + HW + 1]
    print("Sentinel-2 scene classes in the box:", dict(zip(*np.unique(scl, return_counts=True))), "(5 = bare soil, 8-10 = cloud)")
    rcn = (DATA / "radcalnet") / ("RVUS00_2020_194_v00.06.input" if DAY == "2020" else C.RCN.name)
    _, wr, r_s2t, u_s2t = C.radcalnet(rcn, T_S2)
    _, wr2, r_avt, _ = C.radcalnet(rcn, T_AV)
    k = r_s2t / np.interp(wr, wr2, r_avt)          # RadCalNet diurnal factor, Sentinel time / AVIRIS time
    if DAY == "2020":                              # AVIRIS-NG: DUSTRAN with AERONET aerosol (no NASA product here)
        import playa_gas_score as P
        z = np.load(str(DATA / "ng_2020_playa_runs.npz"))
        wl, dus = z["wl"], z["AERONET_0.062"]; nasa = np.full(len(wl), np.nan)
        good = P.load_scene()["bbl"] > 0
        other = "-"
    else:
        S = C.scene(); wl = S["centers"]; good = S["bbl"] > 0
        dus = np.load(C.SC / "dustran_playa_rho.npy")[1]
        nasa = np.load(C.SC / "nasa_playa_rfl.npy")
        other = "NASA"
    o = np.argsort(wl)
    kk = np.interp(wl, wr, k)
    print(f"day {DAY}: diurnal factor (Sentinel time / AVIRIS time) {kk[np.argmin(abs(wl - 865))]:.3f} at 865 nm")
    rows = []
    print(f"\n{'band':5s} {'nm':>6s} {'S2':>7s} {'RadCalNet':>10s} {'DUSTRAN':>8s} {other:>7s}   ratios to S2: RCN / DUSTRAN / {other}")
    for b, (c, w) in BANDS.items():
        s2, sd = s2_band(b)
        rc = band_avg(wr, r_s2t, c, w)
        du = band_avg(wl[o], np.where(good, dus * kk, np.nan)[o], c, w)
        na = band_avg(wl[o], np.where(good, nasa * kk, np.nan)[o], c, w)
        rows.append((b, c, s2, rc, du, na))
        print(f"{b:5s} {c:6.0f} {s2:7.4f} {rc:10.4f} {du:8.4f} {na:7.4f}   {rc / s2:.3f} / {du / s2:.3f} / {na / s2:.3f}")
    np.save(C.SC / "three_way.npy" if DAY == "2019" else (DATA / "three_way_2020.npy"), np.array([r[1:] for r in rows]))

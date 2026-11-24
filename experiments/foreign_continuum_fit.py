"""Measure the water foreign continuum from playa days with independent photometer water.

    python experiments/foreign_continuum_fit.py

Two days so far (both Railroad Valley, both with AERONET water):
  AVIRIS Classic 2019-08-13 (sensor 20.5 km, photometer 1.00 cm)
  AVIRIS-NG     2020-07-12 (sensor 5.3 km,  photometer 1.19 cm)
For each collision-band width gamma and strength A (dustran/foreign.py), water is fitted on 940
and on 1140 nm against RadCalNet's ground shape, as before. A is chosen to bring the fitted water
to the photometer on one day, then applied blind to the other day (and the reverse). A continuum
that is real physics should transfer between instruments and altitudes.
"""
import sys
from pathlib import Path
import numpy as np
R = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(R)); sys.path.insert(0, str(R / "experiments"))
from dustran import atmosphere, dimer, foreign
import playa_gas_score as P
import playa_gas_tau as T
import classic_playa as C

XS = Path("D:/workspace/dustran_data/hitran_xsec")
CIA = {XS / "O2-O2_2024.cia": ("O2", "O2"), XS / "O2-N2_2024.cia": ("O2", "N2")}
O3 = "D:/workspace/dustran_data/bremen_o3/serdyuchenkogorshelev5digits.dat"
B940 = (880, 1000, ((850, 880), (1000, 1040)))
B1140 = (1080, 1180, ((1040, 1075), (1230, 1260)))
GAMMAS = (3.0, 10.0, 30.0)


class ModelF(P.Model):
    A = 0.0

    def rho(self, w=1.0, shift=0.0, gases=True):
        saved = [s["Edry"] for s in self.segs]
        for s in self.segs:
            s["Edry"] = s["Edry"] * np.exp(-w * self.A * s["funit"])
        try:
            return super().rho(w, shift, gases)
        finally:
            for s, e in zip(self.segs, saved):
                s["Edry"] = e


def scenes():
    out = {}
    atm, wr, rr, _ = C.radcalnet(C.RCN, C.UTC)
    S = C.scene()
    M = ModelF([C.SC / "tau_swir", C.SC / "tau_vnir"], S, cia_files=CIA, o3_file=O3, o3_du=atm["O3"],
               gas_scale={"CO2": atmosphere.co2_ppm(C.YEAR) / 330.2})
    out["Classic 2019"] = dict(M=M, S=S, ground=np.interp(S["centers"], wr, rr), shift=0.0, photo=1.00,
                               surface=C.GROUND, sensor=C.SENSOR)
    atm, wr, rr, _ = C.radcalnet("D:/workspace/dustran_data/radcalnet/RVUS00_2020_194_v00.06.input", "21:30")
    S = P.load_scene()
    M = ModelF("D:/workspace/dustran_data/playa_rrv", S)
    out["NG 2020"] = dict(M=M, S=S, ground=np.interp(S["centers"], wr, rr),
                          shift=0.389 + 0.063e-3 * (S["centers"] - 1500), photo=1.19,
                          surface=T.GROUND_KM, sensor=T.SENSOR_KM)
    return out


def water_fit(d, band):
    M, S, wl, good = d["M"], d["S"], d["S"]["centers"], d["S"]["bbl"] > 0
    lo, hi, e = band

    def cost(w):
        q = M.rho(w, d["shift"]) / d["ground"]
        inw = good & (wl >= lo) & (wl <= hi)
        ed = good & (((wl >= e[0][0]) & (wl <= e[0][1])) | ((wl >= e[1][0]) & (wl <= e[1][1])))
        b = np.polyfit(wl[ed], q[ed], 1)
        return np.sqrt(np.mean((q[inw] / np.polyval(b, wl[inw]) - 1) ** 2))
    ws = np.arange(0.3, 3.0, 0.1); c = [cost(w) for w in ws]; i = int(np.argmin(c))
    w2 = np.linspace(ws[max(i - 1, 0)], ws[min(i + 1, len(ws) - 1)], 11); c2 = [cost(w) for w in w2]
    j = int(np.argmin(c2))
    return w2[j] * d["cm1"], c2[j]


if __name__ == "__main__":
    lines = dimer.water_lines()
    D = scenes()
    prof = atmosphere.load_afgl(T.AFGL)
    for name, d in D.items():
        d["cm1"] = atmosphere.layers(prof, d["surface"])["col"]["H2O"].sum() / 3.3428e22
        L = atmosphere.layers(prof, surface_km=d["surface"], sensor_km=d["sensor"])
        d["L"] = L; d["groups"] = {"all": np.ones(len(L["p"]), bool), "below": L["below_sensor"]}
    res = {}
    for g in GAMMAS:
        for name, d in D.items():
            mu_s, mu_v = np.cos(np.radians(d["S"]["sza"])), np.cos(np.radians(d["S"]["vza"]))
            for s in d["M"].segs:
                t = foreign.optical_depth_unit(d["L"], s["wn"], d["groups"], g, lines)
                s["funit"] = t["all"] / mu_s + t["below"] / mu_v
        # strength grid: A0 gives two-way optical depth ~0.05 at the 940 band peak (Classic)
        s0 = D["Classic 2019"]["M"].segs[0]
        m = (1e7 / s0["wn"] > 930) & (1e7 / s0["wn"] < 950)
        A0 = 0.05 / s0["funit"][m].max()
        for f in (0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0):
            for name, d in D.items():
                d["M"].A = f * A0
                w9, r9 = water_fit(d, B940); w11, r11 = water_fit(d, B1140)
                res[(g, f, name)] = (w9, r9, w11, r11)
                print(f"gamma {g:4.0f}  A {f * A0:.3e} ({f:4.2f})  {name:12s} 940 {w9:.3f} cm ({100 * r9:.2f}%)  "
                      f"1140 {w11:.3f} cm ({100 * r11:.2f}%)   photometer {d['photo']:.2f}", flush=True)
        res[(g, "A0")] = A0
    np.save(Path("D:/workspace/dustran_data/foreign_fit.npy"), np.array([res], dtype=object))

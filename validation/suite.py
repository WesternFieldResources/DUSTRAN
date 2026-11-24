"""Real-scene regression suite: every validated day, every reference, one table.

    python validation/suite.py              evaluate saved outputs against references and baseline
    python validation/suite.py --rerun      rerun the scene pipelines first (slow: minutes per scene)
    python validation/suite.py --update     accept the current numbers as the new baseline

Cases (data on D:, never in the repo; a case whose data is missing is skipped, not failed):
  ng2020        AVIRIS-NG, Railroad Valley 2020-07-12: DUSTRAN / RadCalNet ground by range
  classic2019   AVIRIS Classic, Railroad Valley 2019-08-13: DUSTRAN / RadCalNet; seams
  ivanpah_*     AVIRIS Classic, Ivanpah 2017/2020/2021: DUSTRAN / Sentinel-2 per band; seams
  water         fitted water at 940 and 1140 nm vs photometer, both Railroad Valley days
Each number is compared with validation/baseline.json; moves larger than TOL are flagged.
Seam ratios are reported against 1 (flag >1 %), never corrected.
"""
import json, subprocess, sys
from pathlib import Path
import numpy as np
import sys as _sys; _sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1]))
from dustran.paths import DATA
R = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(R)); sys.path.insert(0, str(R / "experiments"))

D = DATA
PY = sys.executable
BASE = R / "validation/baseline.json"
TOL = {"ratio": 0.005, "water": 0.02, "seam": 0.002}
RANGES = ((400, 450), (450, 500), (500, 700), (700, 1000), (1000, 1300), (1500, 1750), (2000, 2300))
IVANPAH = {"ivanpah_2017": "f171207t01p00r07", "ivanpah_2020": "f201005t01p00r09", "ivanpah_2021": "f210224t01p00r07"}


def rcn_ratio(wl, rho, good, wr, rr):
    to10 = np.array([np.nanmean(np.where(good & (abs(wl - c) <= 5), rho, np.nan)) for c in wr])
    q = to10 / rr
    return {f"{lo}-{hi}": float(np.nanmean(q[(wr >= lo) & (wr <= hi)])) for lo, hi in RANGES}


def case_ng2020():
    import playa_gas_score as P, classic_playa as C
    z = np.load(D / "ng_2020_playa_runs.npz")
    _, wr, rr, _ = C.radcalnet(D / "radcalnet/RVUS00_2020_194_v00.06.input", "21:30")
    return {"ratio": rcn_ratio(z["wl"], z["AERONET_0.062"], P.load_scene()["bbl"] > 0, wr, rr)}


def case_classic2019():
    import classic_playa as C
    from dustran import instrument
    import playa_gas_score as P
    wl, rho = np.load(C.SC / "dustran_playa_rho.npy")
    _, wr, rr, _ = C.radcalnet(C.RCN, C.UTC)
    S = C.scene()
    t = np.loadtxt(P.TSIS)
    sun = instrument.convolve(1e7 / t[::-1, 0], t[::-1, 1], S["centers"], S["fwhm"])
    seams = instrument.seam_check(S["centers"], S["fwhm"], S["L"], sun)
    return {"ratio": rcn_ratio(wl, rho, S["bbl"] > 0, wr, rr),
            "seam": {f"{s['ends_nm']:.0f}|{s['starts_nm']:.0f}": s["ratio"] for s in seams if not s["low_signal"]}}


def case_ivanpah(key):
    f = IVANPAH[key]
    z = np.load(D / "ivanpah/three_way_ivanpah.npz", allow_pickle=True)
    s = np.load(D / "ivanpah" / f / "site.npz")
    seams = {f"seam{i}": float(r) for i, (r, low) in enumerate(zip(s["seam_ratio"], s["seam_low"])) if not low}
    return {"ratio": dict(zip([str(b) for b in z["bands"]], map(float, z[f]))), "seam": seams}


def case_water():
    import foreign_continuum_fit as F
    from dustran import atmosphere
    import playa_gas_tau as T
    out = {}
    Ds = F.scenes(); prof = atmosphere.load_afgl(T.AFGL)
    for name, d in Ds.items():
        d["cm1"] = atmosphere.layers(prof, d["surface"])["col"]["H2O"].sum() / 3.3428e22
        for s in d["M"].segs:
            s["funit"] = 0 * s["wn"]          # continuum off until a measured one is adopted
        for band, b in (("940", F.B940), ("1140", F.B1140)):
            out[f"{name} {band} - photometer"] = F.water_fit(d, b)[0] - d["photo"]
    return {"water": out}


CASES = {"ng2020": case_ng2020, "classic2019": case_classic2019,
         **{k: (lambda k=k: case_ivanpah(k)) for k in IVANPAH}, "water": case_water}
RERUN = [["experiments/playa_rayleigh.py", f"{DATA}/playa_rrv_vnir,{DATA}/playa_rrv",
          str(DATA / "radcalnet/RVUS00_2020_194_v00.06.input"), "--aeronet"],
         ["experiments/classic_playa.py"], ["experiments/ivanpah_three_way.py"]]


def flatten(res):
    return {f"{case}/{kind}/{k}": v for case, kinds in res.items() for kind, d in kinds.items() for k, v in d.items()}


if __name__ == "__main__":
    if "--rerun" in sys.argv:
        for cmd in RERUN:
            print("running", cmd[0], flush=True)
            subprocess.run([PY, *cmd], cwd=R, check=True, stdout=subprocess.DEVNULL)
    res = {}
    for name, fn in CASES.items():
        try:
            res[name] = fn()
        except FileNotFoundError as e:
            print(f"skip {name}: missing {Path(e.filename).name if e.filename else e}")
    now = flatten(res)
    base = json.loads(BASE.read_text()) if BASE.exists() else {}
    print(f"\n{'check':52s} {'baseline':>9s} {'now':>9s}  ")
    nflag = 0
    for k, v in now.items():
        kind = k.split("/")[1]
        b = base.get(k)
        moved = b is not None and abs(v - b) > TOL[kind]
        seam_bad = kind == "seam" and abs(v - 1) > 0.01
        nflag += moved
        print(f"{k:52s} {b if b is not None else float('nan'):9.4f} {v:9.4f}  "
              f"{'MOVED ' if moved else ''}{'seam>1% ' if seam_bad else ''}")
    print(f"\n{len(now)} checks, {nflag} moved beyond tolerance")
    if "--update" in sys.argv:
        BASE.write_text(json.dumps(now, indent=1))
        print("baseline updated")

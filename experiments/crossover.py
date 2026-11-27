"""Brightness check between two passes over the same ground: is the difference a level per spectrometer?

    python experiments/crossover.py CUBE_A_HDR CUBE_B_HDR [--json OUT] [--step 3]

Both cubes are DUSTRAN (or any ENVI) reflectance in the same UTM zone, with ENVI map info (rotation
allowed). Every `step`-th pixel of A is mapped to B by map coordinates (nearest pixel). The rotation
sign convention is not taken on trust: both are tried and the one whose matched 860 nm images
correlate better is used (both correlations are reported).

Matched pixels kept: valid in both, bare (NDVI < 0.15 in both), not the darkest 20 % at 1650 nm.
Per band: median of B/A. Per spectrometer (split at the seams): the mean ratio and the slope of the
ratio across the spectrometer (per 100 nm), outside strong gas bands. The read-out:

  * flat inside each spectrometer, a different level in each: an instrument (calibration) level
    difference, and the per-spectrometer levels are the gains that would match B to A;
  * a smooth slope that runs across the seams: atmosphere, illumination or the surface itself;
    a per-spectrometer gain would hide it, so none is suggested.
Nothing is corrected here.
"""
import argparse, json
import numpy as np
from dustran.batch import read_header
from dustran import instrument

GAS = [(750, 775), (920, 980), (1100, 1170), (1320, 1500), (1780, 2000), (2000, 2080), (2350, 2500)]
FLAT_PER_100NM = 0.005          # |slope| below this per 100 nm counts as flat


def open_cube(p):
    _, h, arr = read_header(p)
    nl, nb, ns = int(h["lines"]), int(h["bands"]), int(h["samples"])
    mi = [x.strip() for x in h["map info"].strip("{}").split(",")]
    rot = 0.0
    for x in mi:
        if x.startswith("rotation="):
            rot = float(x.split("=")[1])
    g = dict(px=float(mi[1]) - 1, py=float(mi[2]) - 1, x0=float(mi[3]), y0=float(mi[4]),
             dx=float(mi[5]), dy=float(mi[6]), zone=mi[7], rot=np.radians(rot))
    return np.memmap(p[:-4], "<f4", "r", shape=(nl, nb, ns)), arr("wavelength"), g


def to_map(g, col, row, sign):
    i, j = col - g["px"], row - g["py"]
    c, s = np.cos(g["rot"]), sign * np.sin(g["rot"])
    return g["x0"] + g["dx"] * (i * c + j * s), g["y0"] + g["dx"] * (i * s) - g["dy"] * (j * c)


def to_pix(g, x, y, sign):
    c, s = np.cos(g["rot"]), sign * np.sin(g["rot"])
    u, v = (x - g["x0"]) / g["dx"], (g["y0"] - y) / g["dy"]
    # invert [u, -v'] : u = i c + j s ; v = -(i s) + j c
    i = u * c - v * s
    j = u * s + v * c
    return np.rint(i + g["px"]).astype(int), np.rint(j + g["py"]).astype(int)


def sample(cube, rows, cols, bands=None):
    out = []
    order = np.argsort(rows)
    r, cl = rows[order], cols[order]
    for rr in np.unique(r):
        m = r == rr
        x = np.asarray(cube[rr][:, cl[m]] if bands is None else cube[rr][bands][:, cl[m]], np.float32).T
        out.append(x)
    X = np.concatenate(out)
    back = np.empty_like(order); back[order] = np.arange(len(order))
    return X[back]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("a"); ap.add_argument("b"); ap.add_argument("--json"); ap.add_argument("--step", type=int, default=3)
    a = ap.parse_args()
    A, c, ga = open_cube(a.a)
    B, cb, gb = open_cube(a.b)
    assert ga["zone"] == gb["zone"] and len(c) == len(cb)
    nl, nb, ns = A.shape
    rr, cc = np.meshgrid(np.arange(0, nl, a.step), np.arange(0, ns, a.step), indexing="ij")
    rr, cc = rr.ravel(), cc.ravel()
    i860 = int(np.argmin(np.abs(c - 860)))
    best = None
    tried = {}
    for sign in (+1, -1):
        x, y = to_map(ga, cc, rr, sign)
        bi, bj = to_pix(gb, x, y, sign)
        ok = (bi >= 0) & (bi < B.shape[2]) & (bj >= 0) & (bj < B.shape[0])
        if ok.sum() < 1000:
            tried[sign] = dict(overlap=int(ok.sum()), r=None)
            continue
        k = np.flatnonzero(ok)
        k = k[np.random.default_rng(0).choice(len(k), min(len(k), 20000), replace=False)]
        va = sample(A, rr[k], cc[k], [i860])[:, 0]
        vb = sample(B, bj[k], bi[k], [i860])[:, 0]
        v = (va > 0) & (vb > 0)
        r = float(np.corrcoef(va[v], vb[v])[0, 1]) if v.sum() > 100 else None
        tried[sign] = dict(overlap=int(ok.sum()), r=r)
        if r is not None and (best is None or r > best[0]):
            best = (r, sign, ok, bi, bj)
    print("rotation sign tried:", tried)
    if best is None:
        print("no overlap"); return
    r, sign, ok, bi, bj = best
    k = np.flatnonzero(ok)
    XA = sample(A, rr[k], cc[k]); XB = sample(B, bj[k], bi[k])
    i = lambda w: int(np.argmin(np.abs(c - w)))
    nd = lambda X: (X[:, i(860)] - X[:, i(660)]) / (X[:, i(860)] + X[:, i(660)] + 1e-9)
    v = (XA[:, i(1650)] > 0) & (XB[:, i(1650)] > 0) & (XA[:, i(860)] > 0) & (XB[:, i(860)] > 0)
    v &= (nd(XA) < 0.15) & (nd(XB) < 0.15)
    v &= XA[:, i(1650)] > np.percentile(XA[v, i(1650)], 20)
    print(f"matched bare pixels: {v.sum()} (860 nm correlation {r:.3f}, rotation sign {sign:+d})")
    ratio = np.nanmedian(XB[v] / np.where(XA[v] > 0.005, XA[v], np.nan), 0)
    gas = np.zeros(len(c), bool)
    for lo, hi in GAS:
        gas |= (c >= lo) & (c <= hi)
    bounds = [0] + instrument.seams(c) + [len(c)]
    specs = {}
    for name, (lo, hi) in zip("ABCD", zip(bounds[:-1], bounds[1:])):
        m = np.zeros(len(c), bool); m[lo:hi] = True; m &= ~gas & np.isfinite(ratio) & (c > 420)
        if m.sum() < 4:
            continue
        p = np.polyfit((c[m] - c[m].mean()) / 100, ratio[m], 1)
        specs[name] = dict(range_nm=[float(c[lo]), float(c[hi - 1])], level=float(np.median(ratio[m])),
                           slope_per_100nm=float(p[0]), flat=bool(abs(p[0]) < FLAT_PER_100NM),
                           rms_about_line=float(np.std(ratio[m] - np.polyval(p, (c[m] - c[m].mean()) / 100))))
    seams = []
    for s in bounds[1:-1]:
        seams.append(dict(at_nm=float(c[s]), ratio_before=float(ratio[s - 1]), ratio_after=float(ratio[s]),
                          step=float(ratio[s] / ratio[s - 1])))
    allflat = all(v["flat"] for v in specs.values())
    lv = [v["level"] for v in specs.values()]
    verdict = ("level per spectrometer (flat inside each): instrument-type difference; levels = B/A gains"
               if allflat else
               "ratio slopes inside at least one spectrometer: atmosphere, illumination or surface; no gain suggested")
    for k_, v_ in specs.items():
        print(f"  {k_} {v_['range_nm'][0]:.0f}-{v_['range_nm'][1]:.0f} nm: B/A level {v_['level']:.3f}, "
              f"slope {v_['slope_per_100nm']:+.4f}/100 nm {'flat' if v_['flat'] else 'SLOPED'}, rms {v_['rms_about_line']:.4f}")
    for s in seams:
        print(f"  seam {s['at_nm']:.0f} nm: B/A {s['ratio_before']:.3f} -> {s['ratio_after']:.3f} (step {s['step']:.3f})")
    print("verdict:", verdict, f"(levels spread {max(lv) - min(lv):.3f})")
    if a.json:
        open(a.json, "w").write(json.dumps(dict(a=a.a, b=a.b, matched=int(v.sum()), corr860=r, rotation_sign=sign,
                                                 tried={str(k): t for k, t in tried.items()},
                                                 spectrometers=specs, seams=seams, verdict=verdict,
                                                 ratio=[None if not np.isfinite(x) else float(x) for x in ratio],
                                                 wavelength=c.tolist()), indent=1))


if __name__ == "__main__":
    main()

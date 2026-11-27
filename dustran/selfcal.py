"""Spectral self-calibration from the scene: band-centre offset, slope and width per spectrometer.

    python -m dustran.selfcal RADIANCE OUT_JSON --taudirs DIR[,DIR] --ground KM --sensor KM
                              [--obs OBS_ORT_HDR] [--year 2018]

The scene's median radiance over bare ground contains features at wavelengths that are known
exactly: solar lines (Fraunhofer), O2, water, CO2 and CH4. Inside each window the measured
radiance is modelled as

    L(band) = P(wavelength) * < E_sun * T_gas(water) >_band(centre + d, width * s)

where P is a quadratic continuum (surface and scattering, smooth over ~50 nm), < >_band is the
Gaussian band response, d the centre offset and s the width scale. For each window, d, s and the
water scale are chosen by grid search to minimise the rms log residual. Then, per spectrometer:

  * offset + slope: a straight line d(wavelength) = a + b (wavelength - mid), weighted by each
    window's sharpness (inverse curvature of its cost), through the window offsets;
  * width scale: weighted median of the window values.

A shift and a width error leave different marks: a shift leaves residuals with opposite signs on
the two sides of a line (antisymmetric), and a width error leaves residuals with the same sign on
both sides (symmetric; Green 1998, Applied Optics 37:683). Fitting both on several lines per
spectrometer is the in-flight method of Gao, Montes and Davis (2004, RSE 90:424) and Guanter,
Richter and Moreno (2006, Applied Optics 45:2360).

Nothing is corrected here. The JSON lists every window, the per-spectrometer fit, and flags:
  * a window whose minimum sits on the grid edge (no fit);
  * windows in one spectrometer that disagree by more than 1 nm after the line is fitted;
  * a spectrometer with fewer than 2 usable windows (slope not fitted, offset only).
dustran.batch --selfcal JSON applies the result (centres + d, widths * s) to the band table it models.
"""
import argparse, json
from pathlib import Path
import numpy as np

from dustran import atmosphere, instrument
from dustran.batch import Scene, AFGL, open_radiance, read_header

# (name, lo, hi, has_water). Window edges include continuum on both sides of each feature.
WINDOWS = [
    ("Sun 430-486 G band, H-gamma, H-beta", 415, 500, False),
    ("Sun 517 Mg b", 495, 545, False),
    ("Sun 589 Na D", 560, 620, False),
    ("Sun 656 H-alpha", 625, 668, False),
    ("O2 687 + H2O 720", 665, 750, True),
    ("O2 760", 740, 800, True),
    ("H2O 820", 790, 865, True),
    ("H2O 1140", 1075, 1205, True),
    # H2O 940 is left out on purpose: the ferric iron band (860-920 nm) sits under it, so on
    # iron-bearing ground its best width runs to the grid edge. It is used for water, not calibration.
    ("O2 1270", 1230, 1310, True),
    ("CO2 1575/1605", 1540, 1640, True),
    ("CH4 1666", 1630, 1710, True),
    ("CO2 2010/2060", 1980, 2090, True),
    ("CH4/H2O 2320", 2280, 2370, True),
    ("CH4/H2O 2380", 2350, 2430, True),
]
MIN_BANDS = 6
# Fitted and reported, but kept out of the per-spectrometer line: on all six AVIRIS Classic lines
# tested (2014-2018) this window sits +1.3 to +1.7 nm off the line through the others, the same in
# every year, so the offset is in our model of the window (O2 B band at the spectrometer's blue
# end), not in the instrument.
REPORT_ONLY = {"O2 687 + H2O 720"}
D_GRID = np.round(np.arange(-5.0, 5.001, 0.2), 2)
S_GRID = np.round(np.arange(0.80, 1.301, 0.05), 3)
W_GRID = np.linspace(0.1, 5.0, 50)
DISAGREE_NM = 1.0


class Fine:
    """The fine-grid sun x gas spectrum of a scene (no gas table, no scattering)."""

    def __init__(self, taudirs, ground_km, sensor_km, year, o3_du=None):
        s = Scene.__new__(Scene)
        s.L = atmosphere.layers(atmosphere.load_afgl(str(AFGL)), surface_km=ground_km, sensor_km=sensor_km)
        groups = {"all": np.ones(len(s.L["p"]), bool), "below": s.L["below_sensor"]}
        self.wn, self.E, self.da, self.db, self.wa, self.wb = s._fine(
            taudirs, groups, atmosphere.co2_ppm(year) / 330.2, o3_du)
        self.wl = 1e7 / self.wn

    def window(self, lo, hi, mu_s, mu_v, pad=40.0, bin_nm=0.05, ws=None):
        """Fine spectrum over the window, averaged into bin_nm bins (exact for band averages to well
        under 0.1 %: a 10 nm Gaussian is flat across 0.05 nm). Returns bin centres, bin widths and,
        per water scale in ws, the bin-mean of E_sun * T_gas."""
        m = (self.wl >= lo - pad) & (self.wl <= hi + pad)
        wl = self.wl[m]
        base = (self.E * np.exp(-(self.da / mu_s + self.db / mu_v)))[m]
        wet = (self.wa / mu_s + self.wb / mu_v)[m]
        edges = np.arange(lo - pad, hi + pad + bin_nm, bin_nm)
        k = np.clip(np.digitize(wl, edges) - 1, 0, len(edges) - 2)
        cnt = np.bincount(k, minlength=len(edges) - 1)
        ok = cnt > 0
        ws = np.array([1.0]) if ws is None else ws
        Y = np.stack([np.bincount(k, base * np.exp(-w * wet), len(edges) - 1)[ok] / cnt[ok] for w in ws], 1)
        ctr = 0.5 * (edges[:-1] + edges[1:])[ok]
        return ctr, np.full(ok.sum(), bin_nm), Y


def _gauss(wl, dw, c, f):
    s = f / (2 * np.sqrt(2 * np.log(2)))
    g = np.exp(-0.5 * ((wl[None, :] - c[:, None]) / s[:, None]) ** 2) * dw[None, :]
    return g / g.sum(1, keepdims=True)


def fit_window(fine, c, f, L, lo, hi, water, mu_s, mu_v):
    """Grid search over offset d, width scale s, water w for bands with centres in [lo, hi]."""
    sel = np.where((c >= lo) & (c <= hi) & (L > 0))[0]
    if len(sel) < MIN_BANDS:
        return None
    ws = W_GRID if water else np.array([1.0])
    wl, dw, Y = fine.window(lo, hi, mu_s, mu_v, ws=ws)                  # bins, widths, bins x water
    if wl.size < 50 or wl[0] > c[sel].min() - 3 * f[sel].max() or wl[-1] < c[sel].max() + 3 * f[sel].max():
        return None
    x = (c[sel] - c[sel].mean()) / 50.0
    A = np.c_[np.ones_like(x), x, x ** 2]
    Hres = np.eye(len(sel)) - A @ np.linalg.pinv(A)                    # removes the quadratic continuum
    lnL = np.log(L[sel])
    cost = np.empty((len(D_GRID), len(S_GRID), len(ws)))
    for i, d in enumerate(D_GRID):
        for j, s in enumerate(S_GRID):
            M = _gauss(wl, dw, c[sel] + d, f[sel] * s) @ Y                # bands x water
            r = Hres @ (lnL[:, None] - np.log(np.maximum(M, 1e-30)))
            cost[i, j] = np.sqrt(np.mean(r ** 2, 0))
    i, j, k = np.unravel_index(np.argmin(cost), cost.shape)
    edge = i in (0, len(D_GRID) - 1)                      # offset not bracketed: no fit
    width_edge = j in (0, len(S_GRID) - 1)                 # reported; the shared width decides
    water_edge = bool(water and k in (0, len(ws) - 1))
    # parabolic refinement along d and s; sharpness = curvature of cost along d
    def refine(v, grid, idx):
        if 0 < idx < len(grid) - 1:
            den = v[idx - 1] - 2 * v[idx] + v[idx + 1]
            if den > 0:
                return grid[idx] + np.clip(0.5 * (v[idx - 1] - v[idx + 1]) / den, -1, 1) * (grid[1] - grid[0]), den
        return grid[idx], 0.0
    d, curv = refine(cost[:, j, k], D_GRID, i)
    s, _ = refine(cost[i, :, k], S_GRID, j)
    # the residual pattern at the best fit vs at d=0, s=1 (what the delivered table leaves)
    i0, j0 = int(np.argmin(np.abs(D_GRID))), int(np.argmin(np.abs(S_GRID - 1)))
    k0 = int(np.argmin(cost[i0, j0]))
    return dict(bands=[int(b) for b in sel], d_nm=float(d), width_scale=float(s), water_scale=float(ws[k]),
                rms_pct=float(100 * cost[i, j, k]), rms_pct_delivered_table=float(100 * cost[i0, j0, k0]),
                water_scale_delivered_table=float(ws[k0]),
                sharpness=float(curv / (0.2 ** 2) / max(cost[i, j, k], 1e-9)), edge=bool(edge), width_edge=bool(width_edge),
                water_edge=water_edge, _cost=cost)


def at_width(r, js):
    """Offset (refined), water and rms of a window at a fixed width-scale index js."""
    cost = r["_cost"][:, js, :]
    v = cost.min(1)
    i = int(np.argmin(v)); k = int(np.argmin(cost[i]))
    d, curv = D_GRID[i], 0.0
    if 0 < i < len(D_GRID) - 1:
        den = v[i - 1] - 2 * v[i] + v[i + 1]
        if den > 0:
            d = D_GRID[i] + np.clip(0.5 * (v[i - 1] - v[i + 1]) / den, -1, 1) * (D_GRID[1] - D_GRID[0])
            curv = den
    edge = i in (0, len(D_GRID) - 1)
    return float(d), int(k), float(100 * v[i]), float(curv / (0.2 ** 2) / max(v[i], 1e-9)), bool(edge)


def bare_median(rd, n=6000, seed=0):
    """Median radiance of bright, non-vegetated pixels (NDVI < 0.2, brightness 30-95th pct)."""
    rng = np.random.default_rng(seed)
    rows = np.sort(rng.choice(rd.nl, min(rd.nl, 300), replace=False))
    c = rd.centers
    ir, ired = int(np.argmin(np.abs(c - 860))), int(np.argmin(np.abs(c - 660)))
    isw = int(np.argmin(np.abs(c - 1650)))
    picks, geo = [], []
    for r in rows:
        L, g, ok = rd.block(r, r + 1)
        if ok.any():
            L, g = L[ok], g[ok]
            k = rng.choice(len(L), min(len(L), max(1, n // len(rows))), replace=False)
            picks.append(L[k]); geo.append(g[k])
    L, g = np.concatenate(picks), np.concatenate(geo)
    ndvi = (L[:, ir] - L[:, ired]) / (L[:, ir] + L[:, ired] + 1e-9)
    b = L[:, isw]
    lo, hi = np.percentile(b, [30, 95])
    m = (ndvi < 0.2) & (b >= lo) & (b <= hi)
    return np.median(L[m], 0), g[m], int(m.sum())


def run(radiance, taudirs, ground_km, sensor_km, year, obs=None, o3_du=None, log=print):
    rd = open_radiance(radiance, obs)
    c, f = rd.centers.astype(float), rd.fwhm.astype(float)
    L, g, npx = bare_median(rd)
    sza, vza = float(np.median(g[:, 4])), float(np.median(g[:, 2]))
    mu_s, mu_v = np.cos(np.radians(sza)), np.cos(np.radians(vza))
    log(f"{rd.name}: {npx} bare bright pixels; sun zenith {sza:.1f}, view zenith {vza:.1f}")
    fine = Fine(taudirs, ground_km, sensor_km, year, o3_du)
    bounds = [0] + instrument.seams(c) + [len(c)]
    out = dict(radiance=str(radiance), pixels=npx, sza=sza, vza=vza, spectrometers={})
    for si, (a, b) in enumerate(zip(bounds[:-1], bounds[1:])):
        name = "ABCD"[si] if len(bounds) == 5 else str(si)
        cc, ff, LL = c[a:b], f[a:b], L[a:b]
        wins = {}
        for wname, lo, hi, water in WINDOWS:
            if lo < cc.min() or hi > cc.max():
                continue
            r = fit_window(fine, cc, ff, LL, lo, hi, water, mu_s, mu_v)
            if r is None:
                continue
            r["bands"] = [x + a for x in r["bands"]]
            r["centre_nm"] = float(np.mean(c[r["bands"]]))
            wins[wname] = r
            log(f"  {name} {wname[:22]:22s} d {r['d_nm']:+5.2f} nm  width x{r['width_scale']:.3f}  "
                f"water x{r['water_scale']:.2f}  rms {r['rms_pct']:.2f}% (table {r['rms_pct_delivered_table']:.2f}%)"
                f"{'  OFFSET-EDGE' if r['edge'] else ''}{'  water-edge' if r['water_edge'] else ''}")
        # one width scale per spectrometer, shared by its windows: minimise the sum over windows of
        # (cost / own best cost)^2, each window's cost minimised over offset and water at that width
        cand = {k: v for k, v in wins.items() if not v["edge"]}
        flags = []
        if cand:
            tot = sum((v["_cost"].min((0, 2)) / v["_cost"].min()) ** 2 for v in cand.values())
            js = int(np.argmin(tot))
            width = float(S_GRID[js])
            if js in (0, len(S_GRID) - 1):
                flags.append(f"shared width on the grid edge (x{width:.2f})")
            for k, v in wins.items():
                d, kw, rms, sh, e = at_width(v, js)
                v.update(d_at_shared_width_nm=d, rms_pct_at_shared_width=rms, sharpness_at_shared_width=sh,
                         edge_at_shared_width=e)
        else:
            js, width = int(np.argmin(np.abs(S_GRID - 1))), 1.0
        use = {k: v for k, v in cand.items() if not v["edge_at_shared_width"] and v["sharpness_at_shared_width"] > 0
               and k not in REPORT_ONLY}
        mid = float(np.mean(cc))
        if len(use) == 2:                                  # two points cannot test a line: offset only
            y = np.array([v["d_at_shared_width_nm"] for v in use.values()])
            w = np.array([v["sharpness_at_shared_width"] for v in use.values()])
            off, slope = float(np.sum(w * y) / w.sum()), 0.0
            flags.append(f"two windows: offset only (they differ by {y[1] - y[0]:+.2f} nm)")
            for k, yy in zip(use, y):
                use[k]["resid_nm"] = float(yy - off)
        elif len(use) >= 3:
            x = np.array([v["centre_nm"] - mid for v in use.values()])
            y = np.array([v["d_at_shared_width_nm"] for v in use.values()])
            w = np.sqrt(np.array([v["sharpness_at_shared_width"] for v in use.values()]))
            A = np.c_[np.ones_like(x), x] * w[:, None]
            off, slope = np.linalg.lstsq(A, y * w, rcond=None)[0]
            for k, rr in zip(use, y - (off + slope * x)):
                use[k]["resid_nm"] = float(rr)
                if abs(rr) > DISAGREE_NM:
                    flags.append(f"{k} disagrees with the line by {rr:+.2f} nm")
        elif len(use) == 1:
            v = next(iter(use.values()))
            off, slope = v["d_at_shared_width_nm"], 0.0
            flags.append("one usable window: offset only, no slope")
        else:
            off, slope, width = 0.0, 0.0, 1.0
            flags.append("no usable window: delivered table kept")
        for v in wins.values():
            v.pop("_cost")
        out["spectrometers"][name] = dict(
            bands=[a, b], range_nm=[float(cc.min()), float(cc.max())], mid_nm=mid,
            offset_nm=float(off), slope_nm_per_100nm=float(100 * slope), width_scale=width,
            n_windows=len(use), flags=flags, windows=wins)
        log(f"{name} {cc.min():.0f}-{cc.max():.0f} nm: offset {off:+.2f} nm, slope {100 * slope:+.2f} nm/100 nm, "
            f"width x{width:.3f} ({len(use)} windows){'; ' + '; '.join(flags) if flags else ''}")
    return out


def apply(selfcal, centers, fwhm):
    """Per-band shift [nm] and corrected FWHM from a selfcal JSON (dict)."""
    c = np.asarray(centers, float)
    shift, f = np.zeros_like(c), np.asarray(fwhm, float).copy()
    for v in selfcal["spectrometers"].values():
        a, b = v["bands"]
        shift[a:b] = v["offset_nm"] + v["slope_nm_per_100nm"] / 100 * (c[a:b] - v["mid_nm"])
        f[a:b] *= v["width_scale"]
    return shift, f


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("radiance"); ap.add_argument("out_json")
    ap.add_argument("--obs"); ap.add_argument("--taudirs", required=True)
    ap.add_argument("--ground", type=float, required=True); ap.add_argument("--sensor", type=float, required=True)
    ap.add_argument("--year", type=int, default=2020); ap.add_argument("--o3-du", type=float)
    a = ap.parse_args(argv)
    out = run(a.radiance, a.taudirs.split(","), a.ground, a.sensor, a.year, a.obs, a.o3_du,
              log=lambda s: print(s, flush=True))
    Path(a.out_json).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out_json).write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()

"""Whole-scene atmospheric correction: every pixel of a flight line at once, in numpy.

    python -m dustran.batch RADIANCE_HDR OBS_ORT OUT_PREFIX --taudirs DIR[,DIR] --ground KM --sensor KM
                            [--o3-du 282] [--aeronet DIR --inv-time T --sun-time T --date DD:MM:YYYY]
                            [--aod 0.05] [--pressure-ratio 1.0] [--shift A B] [--rows N]

Same physics as the validated playa chain (experiments/playa_rayleigh.py, classic_playa.py),
reorganised so that the per-pixel work is array arithmetic:

Once per scene
  * Line-by-line gas optical depths (from experiments/playa_gas_tau.py at the scene's ground
    elevation and sensor altitude), O2 collision-induced absorption, ozone (Bremen cross sections
    scaled to a measured column), TSIS-1 solar spectrum.
  * A small gas table built on the fly (no stored lookup tables): ln band<E_sun T_gas> for
    ~40 water scales x a few view zenith angles, with the fitted Gaussian band shift.
  * Scattering terms R_path, T_dn T_up, S from dustran/sos.py on a coarse wavelength grid, for a
    few view geometries (view zenith x relative azimuth); Rayleigh + aerosol (AERONET-measured
    sizes and index, or a lognormal fitted to an Angstrom exponent).

Per pixel (vectorised, in row blocks)
  * Water: the scale that makes the gas-corrected spectrum flattest across the 940 nm band,
    relative to a straight line through clean bands on either side (as experiments/playa_water.py),
    picked from the table and refined with a parabola. 1140 nm is fitted the same way as a
    cross-check; both go into the water map, and their disagreement is the scene's water flag.
  * r = gas-corrected apparent reflectance at that water; x = r - R_path;
    rho = x / (T_dn T_up + S x).

Output: an ENVI BIL float32 reflectance cube with the input's band table and map info, a
3-band water map (940 cm, 1140 cm, 940 fit rms %), and a JSON summary.
Ground elevation, solar zenith for the gas paths and the aerosol are single values per scene.
"""
import argparse, json, re, sys, time
from pathlib import Path
import numpy as np

from dustran import atmosphere, aerosol, cia, instrument, rayleigh, sos, xsec
from dustran.paths import DATA

AFGL = DATA / "afgl/afglus.dat"
TSIS = DATA / "tsis/tsis_f0_0p1.txt"
CIA = {DATA / "hitran_xsec/O2-O2_2024.cia": ("O2", "O2"), DATA / "hitran_xsec/O2-N2_2024.cia": ("O2", "N2")}
O3 = DATA / "bremen_o3/serdyuchenkogorshelev5digits.dat"
H2O_CM = 1 / 3.3428e22                    # molecules cm^-2 -> g cm^-2 (precipitable cm)
WATER = {"940": ((880, 1000), ((850, 880), (1000, 1040))),
         "1140": ((1080, 1180), ((1040, 1075), (1230, 1260)))}
SCAT_GRID = np.r_[np.arange(400, 1000, 25.0), np.arange(1000, 2501, 150.0)]


# ---------------------------------------------------------------- ENVI helpers
def read_header(path):
    text = Path(path).read_text()
    hdr = {m.group(1).strip().lower(): m.group(2).strip()
           for m in re.finditer(r"^\s*([^=\n]+?)\s*=\s*(\{[^}]*\}|[^\n]*)", text, re.M)}
    arr = lambda k: np.array([float(v) for v in hdr[k].strip("{}").split(",")]) if k in hdr else None
    return text, hdr, arr


def open_bil(hdr_path, bands=None):
    _, hdr, _ = read_header(hdr_path)
    shape = (int(hdr["lines"]), bands or int(hdr["bands"]), int(hdr["samples"]))
    dt = {"4": "<f4", "5": "<f8"}[hdr["data type"]]
    img = str(hdr_path)[:-4]
    il = hdr.get("interleave", "bil").lower()
    if il == "bil":
        return np.memmap(img, dt, "r", shape=shape), "bil"
    if il == "bip":
        return np.memmap(img, dt, "r", shape=(shape[0], shape[2], shape[1])), "bip"
    raise ValueError(f"interleave {il} not supported")


def write_header(path, template_text, description, bands=None, band_names=None, extra=None):
    """Copy a header, replacing description (and band count/names for non-cube outputs)."""
    t = re.sub(r"description\s*=\s*\{[^}]*\}", "description = {\n  " + description + "}", template_text, count=1)
    if bands is not None:
        keep = []
        for line in t.splitlines():
            k = line.split("=")[0].strip().lower()
            if k in ("wavelength", "fwhm", "bbl", "band names", "smoothing factors", "wavelength file",
                     "wavelength units", "data gain values", "default bands"):
                continue
            keep.append(line)
        t = "\n".join(keep)
        t = re.sub(r"^bands\s*=.*$", f"bands = {bands}", t, flags=re.M)
        t += "\nband names = {" + ", ".join(band_names) + "}"
    t = re.sub(r"^data type\s*=.*$", "data type = 4", t, flags=re.M)
    t = re.sub(r"^interleave\s*=.*$", "interleave = bil", t, flags=re.M)
    if "data ignore value" not in t:
        t += "\ndata ignore value = -9999"
    for k, v in (extra or {}).items():
        t += f"\n{k} = {v}"
    Path(path).write_text(t.rstrip() + "\n")


# ---------------------------------------------------------------- once per scene
class Scene:
    """Everything that does not depend on the pixel."""

    def __init__(self, centers, fwhm, taudirs, ground_km, sensor_km, sza, doy_esd, vza_nodes, dphi_nodes,
                 shift=(0.389, 0.063e-3), o3_du=None, aer=None, aod=None, aod_ref=550.0,
                 pressure_ratio=1.0, ws=np.linspace(0.1, 2.5, 40), co2_ppm=None, log=print):
        self.c, self.f = np.asarray(centers, float), np.asarray(fwhm, float)
        self.shift = shift[0] + shift[1] * (self.c - 1500)
        self.ws, self.vza = np.asarray(ws, float), np.asarray(vza_nodes, float)
        self.dphi = np.asarray(dphi_nodes, float)
        self.mu_s = np.cos(np.radians(sza))
        self.esd = doy_esd
        prof = atmosphere.load_afgl(str(AFGL))
        self.L = atmosphere.layers(prof, surface_km=ground_km, sensor_km=sensor_km)
        self.L["air"] = self.L["air"] * pressure_ratio
        self.cm1 = atmosphere.layers(prof, ground_km)["col"]["H2O"].sum() * H2O_CM
        groups = {"all": np.ones(len(self.L["p"]), bool), "below": self.L["below_sensor"]}
        if co2_ppm is None:
            co2_ppm = atmosphere.co2_ppm(2020)
        t0 = time.time()
        wn, E, dry_all, dry_below, wet_all, wet_below = self._fine(taudirs, groups, co2_ppm / 330.2, o3_du)
        log(f"fine grid: {len(wn)} points, {time.time() - t0:.0f}s")
        # gas table: ln band<E T> over (view node, water scale, band)
        t0 = time.time()
        lnET = np.empty((len(self.vza), len(self.ws), len(self.c)))
        lo, hi = 1e7 / wn.max(), 1e7 / wn.min()
        inside = (self.c + self.shift - 2.5 * self.f >= lo) & (self.c + self.shift + 2.5 * self.f <= hi)
        t = np.loadtxt(TSIS)
        sun_only = instrument.convolve(1e7 / t[::-1, 0], t[::-1, 1], self.c, self.f, self.shift)
        self.E_band = sun_only
        for i, v in enumerate(self.vza):
            mu_v = np.cos(np.radians(v))
            dry = np.exp(-(dry_all / self.mu_s + dry_below / mu_v)) * E
            wet = wet_all / self.mu_s + wet_below / mu_v
            for j, w in enumerate(self.ws):
                m = instrument.convolve(wn, dry * np.exp(-w * wet), self.c, self.f, self.shift)
                m = np.where(inside & np.isfinite(m), m, sun_only)
                lnET[i, j] = np.log(m)
        self.lnET = lnET
        log(f"gas table {lnET.shape}: {time.time() - t0:.0f}s")
        # scattering terms per view node
        t0 = time.time()
        self.scat = self._scattering(aer, aod, aod_ref)
        log(f"scattering {len(SCAT_GRID)} wavelengths x {len(self.vza)} x {len(self.dphi)} geometries: "
            f"{time.time() - t0:.0f}s")

    def _fine(self, taudirs, groups, co2_scale, o3_du):
        segs = []
        for d in taudirs:
            z = {g: np.load(f) for g in ("H2O", "CO2", "CH4", "O2") if (f := Path(d) / f"tau_{g}.npz").exists()}
            wn = z["H2O"]["wn"]
            sc = lambda g: co2_scale if g == "CO2" else 1.0
            segs.append(dict(wn=wn,
                             da=sum(sc(g) * z[g]["all"].astype(np.float64) for g in z if g != "H2O"),
                             db=sum(sc(g) * z[g]["below"].astype(np.float64) for g in z if g != "H2O"),
                             wa=z["H2O"]["all"].astype(np.float64), wb=z["H2O"]["below"].astype(np.float64)))
        segs.sort(key=lambda s: s["wn"][0])
        keep = [segs[0]]
        for s in segs[1:]:
            m = s["wn"] > keep[-1]["wn"][-1]
            keep.append({k: v[m] for k, v in s.items()})
        s = {k: np.concatenate([x[k] for x in keep]) for k in keep[0]}
        wn = s["wn"]
        c = cia.optical_depth(self.L, CIA, wn, groups)
        s["da"] += c["all"]; s["db"] += c["below"]
        scale = (o3_du / xsec.dobson(self.L)) if o3_du else 1.0
        o = xsec.optical_depth(self.L, xsec.read_bremen(str(O3)), "O3", wn, groups, column_scale=scale)
        s["da"] += o["all"]; s["db"] += o["below"]
        t = np.loadtxt(TSIS)
        E = np.interp(1e7 / wn, t[:, 0], t[:, 1])
        return wn, E, s["da"], s["db"], s["wa"], s["wb"]

    def _scattering(self, aer, aod, ref):
        L = self.L
        ks = int(np.argmax(L["below_sensor"]))
        out = np.empty((len(self.vza), len(self.dphi), 4, len(SCAT_GRID)))
        for k, wv in enumerate(SCAT_GRID):
            tR = rayleigh.layer_tau(L, wv)
            if aod:
                o = aer.optics(wv)
                bR = np.zeros(len(o["beta"])); bR[:3] = rayleigh.phase_moments(wv)
                tA = aerosol.layer_tau(L, aod, o["ext"] / aer.optics(ref)["ext"])
                tau = tR + tA; sca = tR + o["ssa"] * tA
                beta = (tR[:, None] * bR + (o["ssa"] * tA)[:, None] * o["beta"]) / sca[:, None]
                ssa, n = sca / tau, 16
            else:
                tau, ssa, beta, n = tR, np.ones_like(tR), np.tile(rayleigh.phase_moments(wv), (len(tR), 1)), 8
            for i, v in enumerate(self.vza):
                for j, p in enumerate(self.dphi):
                    r = sos.surface_terms(tau, ssa, beta, self.mu_s, np.cos(np.radians(v)), p, ks, n_streams=n)
                    out[i, j, :, k] = r["R_path"], r["T_dn"] * r["T_up"], r["S"], r["T_dn"]
        # to bands: (vza, dphi, term, band)
        return np.stack([[[np.interp(self.c, SCAT_GRID, out[i, j, t]) for t in range(4)]
                          for j in range(len(self.dphi))] for i in range(len(self.vza))])


# ---------------------------------------------------------------- per pixel
def _bracket(nodes, x):
    """Index i and fraction f with x ~ nodes[i] + f (nodes[i+1]-nodes[i]); clipped to the ends."""
    x = np.clip(x, nodes[0], nodes[-1])
    i = np.clip(np.searchsorted(nodes, x) - 1, 0, len(nodes) - 2)
    return i, (x - nodes[i]) / (nodes[i + 1] - nodes[i])


def _line_operator(wl, inw, edge):
    """Matrix H so that H @ y[edge] is the straight line through y[edge] evaluated at wl[inw]."""
    A = np.c_[np.ones(edge.sum()), wl[edge]]
    return np.c_[np.ones(inw.sum()), wl[inw]] @ np.linalg.pinv(A)


def fit_water(sc, r0, iv, fv, window, good):
    """Per-pixel water scale from one band. r0: (P, B) radiance-based numerator k*L*10 (so that
    r = r0 / band<E T>). Returns scale w (P,) and the rms % at the minimum."""
    (lo, hi), ((a0, a1), (b0, b1)) = window
    wl = sc.c
    inw = good & (wl >= lo) & (wl <= hi)
    edge = good & (((wl >= a0) & (wl <= a1)) | ((wl >= b0) & (wl <= b1)))
    sel = inw | edge
    H = _line_operator(wl[sel], inw[sel], edge[sel])
    ii, ie = inw[sel], edge[sel]
    ln_r0 = np.log(np.maximum(r0[:, sel], 1e-12))
    lnET = sc.lnET[:, :, sel]
    cost = np.empty((len(r0), len(sc.ws)))
    for j in range(len(sc.ws)):
        lnm = lnET[iv, j] * (1 - fv[:, None]) + lnET[iv + 1, j] * fv[:, None]
        r = np.exp(ln_r0 - lnm)
        line = r[:, ie] @ H.T
        cost[:, j] = np.sqrt(np.mean((r[:, ii] / line - 1) ** 2, axis=1))
    j = np.clip(np.argmin(cost, 1), 1, len(sc.ws) - 2)
    P = np.arange(len(r0))
    c0, c1, c2 = cost[P, j - 1], cost[P, j], cost[P, j + 1]
    den = c0 - 2 * c1 + c2
    d = np.where(den > 0, 0.5 * (c0 - c2) / np.where(den > 0, den, 1), 0.0)
    d = np.clip(d, -1, 1)
    step = sc.ws[1] - sc.ws[0]
    return sc.ws[j] + d * step, 100 * np.min(cost, 1)


def correct_block(sc, Lrad, vza, dphi, sza_px, good):
    """Lrad (P, B) radiance [uW cm-2 sr-1 nm-1]; vza, dphi, sza_px (P,). Returns rho (P,B),
    water scale at 940 and 1140 (P,), 940 fit rms %."""
    k = np.pi * sc.esd ** 2 / np.cos(np.radians(sza_px))
    r0 = k[:, None] * Lrad * 10.0                      # uW cm-2 -> mW m-2 (TSIS units)
    iv, fv = _bracket(sc.vza, vza)
    w940, rms940 = fit_water(sc, r0, iv, fv, WATER["940"], good)
    w1140, _ = fit_water(sc, r0, iv, fv, WATER["1140"], good)
    jw, fw = _bracket(sc.ws, w940)
    T = sc.lnET
    lnm = ((T[iv, jw] * (1 - fw[:, None]) + T[iv, jw + 1] * fw[:, None]) * (1 - fv[:, None])
           + (T[iv + 1, jw] * (1 - fw[:, None]) + T[iv + 1, jw + 1] * fw[:, None]) * fv[:, None])
    r = r0 * np.exp(-lnm)
    ip, fp = _bracket(sc.dphi, dphi)
    S4 = sc.scat
    wv, wp = fv[:, None, None], fp[:, None, None]
    terms = ((S4[iv, ip] * (1 - wp) + S4[iv, ip + 1] * wp) * (1 - wv)
             + (S4[iv + 1, ip] * (1 - wp) + S4[iv + 1, ip + 1] * wp) * wv)     # (P, 4, B)
    x = r - terms[:, 0]
    rho = x / (terms[:, 1] + terms[:, 2] * x)
    return rho, w940, w1140, rms940


def relative_azimuth(vaz, saz):
    """|view azimuth - (sun azimuth + 180)| folded into [0, 180] degrees."""
    d = np.abs(((vaz - (saz + 180.0)) + 180.0) % 360.0 - 180.0)
    return d


def run(rdn_hdr, obs_hdr, out_prefix, scene_kw, rows=64, log=print, max_lines=None):
    text, hdr, arr = read_header(rdn_hdr)
    centers, fwhm = arr("wavelength"), arr("fwhm")
    nl, nb, ns = int(hdr["lines"]), int(hdr["bands"]), int(hdr["samples"])
    rdn, _ = open_bil(rdn_hdr)
    obs, il = open_bil(obs_hdr)
    if max_lines:                                  # quick test on the first lines
        nl = min(nl, max_lines); rdn, obs = rdn[:nl], obs[:nl]
        text = re.sub(r"^lines\s*=.*$", f"lines = {nl}", text, flags=re.M)
    assert il == "bip"
    # scene geometry from the ortho obs (valid pixels)
    o = np.asarray(obs[::7, ::7, :]).reshape(-1, obs.shape[2])
    o = o[o[:, 0] > 0]
    sza, esd = float(np.median(o[:, 4])), float(np.median(o[:, 10]))
    vmax = float(np.max(o[:, 2]))
    vza_nodes = np.linspace(0, max(vmax, 1.0) + 0.5, 3)
    dphi_nodes = np.array([0.0, 90.0, 180.0])
    log(f"scene: {nl} x {ns} px, {nb} bands; sza {sza:.2f}, view zenith up to {vmax:.1f}, esd {esd:.5f}")
    sc = Scene(centers, fwhm, sza=sza, doy_esd=esd, vza_nodes=vza_nodes, dphi_nodes=dphi_nodes, log=log, **scene_kw)
    good = np.ones(nb, bool)
    if "bbl" in hdr:
        good = arr("bbl") > 0
    good &= (centers > 395)
    # seams (AVIRIS Classic has three; reported, never corrected)
    med = np.nanmedian(np.where(np.asarray(rdn[nl // 2, :, :]).T > 0, np.asarray(rdn[nl // 2, :, :]).T, np.nan), 0)
    seams = instrument.seam_check(centers, fwhm, med, sc.E_band)
    log("seams: " + ("none (single spectrometer)" if not seams else
                     "; ".join(f"{s['ends_nm']:.0f}|{s['starts_nm']:.0f} ratio {s['ratio']:.4f}"
                               f"{' FLAG' if s['flagged'] else ''}{' low-signal' if s['low_signal'] else ''}"
                               for s in seams)))

    Path(out_prefix).parent.mkdir(parents=True, exist_ok=True)
    out_img = Path(str(out_prefix) + "_rfl_img")
    wat_img = Path(str(out_prefix) + "_h2o_img")
    rfl = np.memmap(out_img, "<f4", "w+", shape=(nl, nb, ns))
    wat = np.memmap(wat_img, "<f4", "w+", shape=(nl, 3, ns))
    t0 = time.time()
    stats = []
    for r0 in range(0, nl, rows):
        r1 = min(nl, r0 + rows)
        L = np.asarray(rdn[r0:r1], np.float64).transpose(0, 2, 1).reshape(-1, nb)   # (P, B)
        g = np.asarray(obs[r0:r1], np.float64).reshape(-1, obs.shape[2])
        ok = (L[:, 0] > -9000) & (g[:, 0] > 0) & (np.nanmax(L, 1) > 0)
        R = np.full((len(L), nb), -9999.0, np.float32)
        W = np.full((len(L), 3), -9999.0, np.float32)
        if ok.any():
            rho, w940, w1140, rms = correct_block(sc, L[ok], g[ok, 2], relative_azimuth(g[ok, 1], g[ok, 3]),
                                                  g[ok, 4], good)
            R[ok] = rho
            W[ok] = np.c_[w940 * sc.cm1, w1140 * sc.cm1, rms]
            stats.append(np.c_[w940, w1140])
        rfl[r0:r1] = R.reshape(r1 - r0, ns, nb).transpose(0, 2, 1)
        wat[r0:r1] = W.reshape(r1 - r0, ns, 3).transpose(0, 2, 1)
        if (r0 // rows) % 8 == 0:
            el = time.time() - t0
            log(f"rows {r1}/{nl}  {el / 60:.1f} min  eta {el / r1 * (nl - r1) / 60:.1f} min")
    rfl.flush(); wat.flush()
    extra = None
    if "bbl" not in hdr:       # radiance headers often carry none: mark the opaque water bands and < 395 nm
        bbl = (centers > 395) & ~((centers > 1350) & (centers < 1440)) & ~((centers > 1805) & (centers < 1985))
        extra = {"bbl": "{ " + " , ".join(f"{float(b):.1f}" for b in bbl) + " }"}
    write_header(str(out_img) + ".hdr", text, "DUSTRAN surface reflectance (Lambertian, plane-parallel)",
                 extra=extra)
    write_header(str(wat_img) + ".hdr", text, "DUSTRAN water column", bands=3,
                 band_names=["water 940 (cm)", "water 1140 (cm)", "940 fit rms (%)"])
    s = np.concatenate(stats)
    rel = (s[:, 1] - s[:, 0]) / s[:, 0]
    summary = dict(
        radiance=str(rdn_hdr), reflectance=str(out_img) + ".hdr", water=str(wat_img) + ".hdr",
        sza=sza, esd=esd, vza_nodes=vza_nodes.tolist(), afgl_water_cm=sc.cm1,
        water_940_cm_pct=dict(zip(["p5", "p50", "p95"], (np.percentile(s[:, 0], [5, 50, 95]) * sc.cm1).tolist())),
        water_1140_cm_pct=dict(zip(["p5", "p50", "p95"], (np.percentile(s[:, 1], [5, 50, 95]) * sc.cm1).tolist())),
        water_1140_vs_940_median_rel=float(np.median(rel)),
        water_flag_frac_gt10pct=float(np.mean(np.abs(rel) > 0.10)),
        water_flag=bool(abs(np.median(rel)) > 0.10),
        seams=[{k: v for k, v in x.items() if k != "pairs"} for x in seams],
        minutes=(time.time() - t0) / 60,
        **{k: (v if isinstance(v, (int, float, str, list)) else str(v)) for k, v in scene_kw.items()
           if k not in ("aer",)})
    Path(str(out_prefix) + "_summary.json").write_text(json.dumps(summary, indent=1))
    log(json.dumps({k: summary[k] for k in ("water_940_cm_pct", "water_1140_cm_pct",
                                             "water_1140_vs_940_median_rel", "water_flag", "minutes")}))
    return summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rdn_hdr"); ap.add_argument("obs_hdr"); ap.add_argument("out_prefix")
    ap.add_argument("--taudirs", required=True, help="comma-separated line-by-line tau folders")
    ap.add_argument("--ground", type=float, required=True, help="ground elevation [km] (gas tau run)")
    ap.add_argument("--sensor", type=float, required=True, help="sensor altitude [km] (gas tau run)")
    ap.add_argument("--o3-du", type=float, default=None, help="measured ozone column [DU]")
    ap.add_argument("--pressure-ratio", type=float, default=1.0, help="scale the air column (Rayleigh, CIA)")
    ap.add_argument("--shift", nargs=2, type=float, default=[0.389, 0.063e-3],
                    help="band shift A + B (nm - 1500) [nm]")
    ap.add_argument("--aeronet", help="AERONET folder (SIZ_hyb.txt, RIN_hyb.txt, aod file)")
    ap.add_argument("--aeronet-aod-file", default="aod.txt")
    ap.add_argument("--inv-time"); ap.add_argument("--sun-time"); ap.add_argument("--date")
    ap.add_argument("--aod", type=float, default=None, help="AOD at 550 nm (lognormal dust, Angstrom below)")
    ap.add_argument("--angstrom", type=float, default=1.2)
    ap.add_argument("--rows", type=int, default=64)
    ap.add_argument("--max-lines", type=int, default=None, help="process only the first N lines (testing)")
    a = ap.parse_args()
    kw = dict(taudirs=a.taudirs.split(","), ground_km=a.ground, sensor_km=a.sensor, o3_du=a.o3_du,
              pressure_ratio=a.pressure_ratio, shift=tuple(a.shift))
    if a.aeronet:
        sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments"))
        import aeronet_aerosol
        aer, meas = aeronet_aerosol.load(a.inv_time, a.sun_time, a.date, a.aeronet, a.aeronet_aod_file)
        kw.update(aer=aer, aod=meas[500], aod_ref=500.0)
        print(f"AERONET: AOD500 {meas[500]:.3f}, photometer water {meas['PW']:.3f} cm")
    elif a.aod:
        kw.update(aer=aerosol.Dust.fit_angstrom(a.angstrom), aod=a.aod, aod_ref=550.0)
    run(a.rdn_hdr, a.obs_hdr, a.out_prefix, kw, rows=a.rows, log=lambda s: print(s, flush=True), max_lines=a.max_lines)

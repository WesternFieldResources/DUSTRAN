"""Instrument spectral response: convolve a fine spectrum to sensor bands.

Gaussian band responses from centre wavelengths and FWHM, with a wavelength shift as a free
parameter (AVIRIS-NG on Railroad Valley showed -0.4 to -0.7 nm). Written for DUSTRAN.
"""
import numpy as np


def seams(centers_nm):
    """Indices where band centres step backwards: the start of each new spectrometer/detector
    in instruments built from several (AVIRIS Classic has 3 seams; AVIRIS-NG has none)."""
    c = np.asarray(centers_nm, float)
    return [i for i in range(1, len(c)) if c[i] < c[i - 1]]


def seam_check(centers_nm, fwhm_nm, radiance, solar, max_dnm=6.0, n=4, flag=0.01):
    """Free relative calibration from the scene itself.

    Where two spectrometers overlap, their bands see the same light, so radiance per unit
    sunlight must agree. For each seam, pair each band of the earlier spectrometer with the
    nearest band of the next (within `max_dnm`) and report the ratio next/previous.
    `solar` is the band-averaged solar irradiance (removes the sun's own spectral shape).
    Returns a list of dicts: seam index, wavelength range, pair ratios, median ratio, flagged.
    A median ratio more than `flag` from 1 means the spectrometers disagree for this scene.
    No correction is applied here; the caller decides, and the numbers are reported.
    """
    c = np.asarray(centers_nm, float)
    r = np.asarray(radiance, float) / np.asarray(solar, float)
    out = []
    for s in seams(c):
        prev = np.arange(max(0, s - n), s)
        nxt = np.arange(s, min(len(c), s + n))
        pairs = []
        for i in prev:
            j = nxt[np.argmin(np.abs(c[nxt] - c[i]))]
            if abs(c[j] - c[i]) <= max_dnm and r[i] > 0 and r[j] > 0:
                pairs.append((float(c[i]), float(c[j]), float(r[j] / r[i])))
        med = float(np.median([p[2] for p in pairs])) if pairs else np.nan
        # low signal: the seam sits in an opaque absorption band (e.g. AVIRIS Classic C|D at
        # ~1870 nm, inside the 1.9 um water band), so the ratio is noise-dominated
        level = np.median(r[prev]) / np.median(r[r > 0])
        low = bool(level < 0.25)
        out.append(dict(seam=s, ends_nm=float(c[s - 1]), starts_nm=float(c[s]), pairs=pairs,
                        ratio=med, signal_level=float(level), low_signal=low,
                        flagged=bool(np.isfinite(med) and abs(med - 1) > flag and not low)))
    return out


def convolve(wn, y, centers_nm, fwhm_nm, shift_nm=0.0, weight=None):
    """Band averages of `y` (on wavenumber grid `wn` [cm^-1]) for Gaussian bands.

    weight  optional fine-grid weighting (e.g. solar irradiance) applied inside each band.
    """
    wn = np.asarray(wn, float)
    y = np.asarray(y)
    c = np.asarray(centers_nm, float) + shift_nm
    s = np.asarray(fwhm_nm, float) / (2 * np.sqrt(2 * np.log(2)))
    asc = wn[0] < wn[-1]
    out = np.full(len(c), np.nan)
    for k in range(len(c)):
        lo, hi = 1e7 / (c[k] + 4 * s[k]), 1e7 / (c[k] - 4 * s[k])     # band window in wavenumber
        if asc:
            i0, i1 = np.searchsorted(wn, [lo, hi])
        else:
            m = np.where((wn >= lo) & (wn <= hi))[0]
            i0, i1 = (m[0], m[-1] + 1) if len(m) else (0, 0)
        if i1 - i0 < 3:
            continue
        wl = 1e7 / wn[i0:i1]
        g = np.exp(-0.5 * ((wl - c[k]) / s[k]) ** 2)
        if weight is not None:
            g = g * weight[i0:i1]
        out[k] = np.sum(g * y[i0:i1]) / np.sum(g)
    return out

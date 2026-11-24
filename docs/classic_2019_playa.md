# Second day: AVIRIS Classic over the playa, 2019-08-13 (f190813t01p00r07)

Run 2026-11-23. Script: `experiments/classic_playa.py`. Data on `D:/workspace/dustran_data`.
**No corrections applied to the radiance:** it is used exactly as delivered (int16 / gain).

## Setup
- **Scene:** ER-2, sensor at ~20.5 km, 15.2 m pixels, 22:13 UTC, sun zenith 40°, view nadir.
- **Playa pixels:** 779 pixels within ±180 m of the RadCalNet point. The pixel-to-pixel spread
  is 0.6 %.
- **Atmosphere (measured):**
  - AERONET Railroad_Valley2 at 22:12:30: AOD(500) = 0.050, water = 1.00 cm
  - aerosol sizes and refractive index from the 21:50 inversion
  - RadCalNet at 22:00: pressure 859 hPa, ozone 295 DU
- **CO2:** set to 2019.
- **Ground truth:** RadCalNet RVUS 2019 day 225, 22:00 UTC. The playa was darker than in 2020
  (0.28 at 550 nm vs 0.37).

## Result: DUSTRAN / RadCalNet

| Range | gases only | full DUSTRAN | spread | RadCalNet 1σ |
|---|---|---|---|---|
| 400–450 | 1.228 | 0.902 | 7.0 % | 6.4 % |
| 450–500 | 1.087 | 0.958 | 0.5 % | 5.6 % |
| 500–700 | 0.958 | 0.935 | 0.9 % | 4.6 % |
| 700–1000 | 0.901 | 0.903 | 2.7 % | 4.2 % |
| 1000–1300 | 0.825 | **0.829** | 5.0 % | 4.2 % |
| 1500–1750 | 0.869 | 0.873 | 1.0 % | 4.0 % |
| 2000–2300 | 0.950 | 0.953 | 4.9 % | 3.7 % |

**Water:** 1.12 cm fitted at 940 nm, 1.42 cm at 1140 nm. The photometer measured 1.00 cm.

## The darkness is in the radiance, not in DUSTRAN
- **The check:** RadCalNet also publishes a *predicted top-of-atmosphere* reflectance. The ER-2
  flies at 20 km, so the raw AVIRIS radiance, divided by the sun with no atmosphere model at
  all, can be compared with it directly.

| | 450 | 550 | 650 | 870 | 1050 | 1250 | 1650 | 2200 nm |
|---|---|---|---|---|---|---|---|---|
| AVIRIS / RadCalNet TOA | 0.94 | 0.95 | 0.95 | 0.90 | 0.85 | 0.84 | 0.88 | 0.95 |

- **The same deficit, with the same shape, appears before any of our physics runs.** DUSTRAN
  carries it through faithfully.
- **Location doesn't explain it.** Moving the sample box up to 1.8 km around the point changes
  brightness by at most ±5 %.
- ~~**The shape steps at Classic's spectrometer boundaries**~~ (WITHDRAWN, see the seam test below) (A ≈ 365–700, B ≈ 660–1270,
  C ≈ 1250–1870, D ≈ 1840–2500 nm):
  - A: −5 %
  - B: −10 to −16 %
  - C: −12 %
  - D: −5 %

  Per-spectrometer structure like this is the fingerprint of radiometric calibration.
- **The 940 vs 1140 water disagreement** (1.12 vs 1.42 cm) fits the same picture: 1140 nm sits
  in spectrometer B, near its upper edge, where the deficit changes fastest.

## Open questions (nothing changed yet)
- Is Classic's 2019 calibration really 5–16 % low, or is RadCalNet's 2019 prediction high?
  NASA's own reflectance product for this flight (downloading) will show which side NASA lands
  on. An independent satellite measurement the same day (Landsat 8 / Sentinel-2 surface
  reflectance) would be a third opinion.
- How water behaves once the radiometry question is settled. At 940 nm (spectrometer B, but in
  its middle) the fit is 1.12 vs 1.00 measured: +12 %, against +28 % on the 2020 NG day.

## NASA's own reflectance for the same flight (f190813t01p00r07_rfl_v1l1)
- **NASA's method** (from its README): a gas and scattering correction, then an **empirical
  multiplier per band** fitted to in-situ invariant targets. The factors are in the header:
  mostly 0.97–1.03 in clean bands, far from 1 in gas bands.

| Range | NASA / RadCalNet | NASA without its gains | DUSTRAN / RadCalNet | NASA / DUSTRAN |
|---|---|---|---|---|
| 400–450 | 0.830 | 0.897 | 0.902 | 0.919 |
| 450–500 | 0.932 | 0.942 | 0.958 | 0.973 |
| 500–700 | 0.894 | 0.897 | 0.935 | 0.956 |
| 700–1000 | 0.861 | 0.877 | 0.903 | 0.954 |
| 1000–1300 | 0.807 | 0.795 | 0.829 | 0.978 |
| 1500–1750 | 0.829 | 0.825 | 0.873 | 0.950 |
| 2000–2300 | 0.922 | 0.909 | 0.953 | 0.970 |

- **NASA's playa is also 7–19 % darker than RadCalNet, with the same spectrometer-step
  pattern,** and on average 2–8 % darker than DUSTRAN.
- **Two independent processors agree with each other and disagree with RadCalNet the same
  way.** That confirms the deficit is upstream of atmospheric correction: in the delivered
  radiance, or in RadCalNet's 2019 values.
- **NASA water at the playa: 1.45 cm.** DUSTRAN gets 1.12 cm (940 nm); AERONET measured
  1.00 cm. DUSTRAN is closer.

## Can it be corrected programmatically?
That depends on whether the deficit is **consistent**: the same per-spectrometer factors across
Classic flights from the same season, or explainable by something like detector temperature or
date. One flight can't show that. Next evidence:
- Sentinel-2 surface reflectance of the playa near 2019-08-13: an independent calibration chain,
  to test RadCalNet.
- Other Classic flights over a calibration site with ground truth, to see whether the factors
  repeat.

## Correction: it is NOT per-spectrometer steps (seam test)
- **The test:** Classic's spectrometers overlap at their seams. Compare signal per unit sunlight
  in the overlapping bands.
- **Result:**
  - A→B (655–668 nm): bands agree to 0.2–0.9 %
  - B→C (1253–1263 nm): bands agree to 0.1–0.4 %
  - C→D (1866–1876 nm): bands agree to 2–3 %, inside an opaque water band, so noisy
- **Meaning:** a per-spectrometer calibration step would show as a jump at the seam. There is
  none. The earlier "spectrometer fingerprint" reading came from too few sample wavelengths and
  is withdrawn.
- **What the deficit actually is:** a smooth, broad dip relative to RadCalNet. About −5 % in the
  visible, deepening to −15/−16 % around 1050–1250 nm, then −12 % at 1650 nm and −5 % at
  2200 nm. It runs continuously across the seams.
- **Candidates, in order:**
  1. RadCalNet's spectral shape. Not checked: the ground measurement may be sparse in wavelength
     and filled in by a model, which would be least constrained in exactly this region.
  2. A smooth spectral error in Classic's calibration.
  3. The playa's own state (moisture or crust) differing between the ground measurement and the
     flight.
- **Decider:** Sentinel-2 L2A on 2019-08-12/14. B8A (865 nm) and B11 (1610 nm) sit in the dip.

## 940 vs 1140 water on Classic: band position and width ruled out
Fit water at each band with a free wavelength shift (−2 to +2 nm) and band width (FWHM ×0.9 to
×1.3), against the RadCalNet surface, on the gain-adjusted radiance:

| FWHM × | 940 nm: water / shift / rms | 1140 nm: water / shift / rms |
|---|---|---|
| 0.9 | 1.11 cm / +0.00 / 4.3 % | 1.38 cm / +0.75 / 2.4 % |
| 1.0 | 1.14 cm / −0.25 / 3.9 % | 1.41 cm / +0.75 / 1.5 % |
| 1.1 | 1.17 cm / −0.25 / 3.7 % | 1.44 cm / +1.00 / 1.8 % |
| 1.2 | 1.20 cm / −0.25 / 3.7 % | 1.44 cm / +1.00 / 2.7 % |

- **No combination brings the two bands together.** 1140 stays 0.25–0.3 cm above 940. So the
  disagreement isn't a simple band shift or band width.
- **The 940 fit is the poorer one here** (rms ~3.7–4.3 % vs 1.5 % at 1140). That's the opposite
  of what I assumed; 940 may be the suspect band on Classic.
- **The photometer says 1.00 cm,** below both.
- **What differs from the NG 2020 day, where the two bands agreed:** the sensor at 20 km (vs 5 km),
  sun zenith 40° (vs 30°), 10 nm Classic bands (vs 5.5 nm NG). Open; the Ivanpah flights
  (different altitudes and seasons) should separate these.

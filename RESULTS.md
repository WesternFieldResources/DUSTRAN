# Results

What has been validated so far, and against what. Ratios are DUSTRAN / reference.

## 1. AVIRIS-NG, Railroad Valley, 2020-07-12 (vs RadCalNet RVUS ground reflectance)

Atmosphere measured on the day: AERONET Railroad_Valley2 (sizes, refractive index, AOD 0.053,
water 1.19 cm), RadCalNet (pressure, ozone). Free parameters: water and a wavelength shift.

| Range (nm) | DUSTRAN / RadCalNet (mean / spread %) |
|---|---|
| 400–450 | 1.036 / 2.3 |
| 450–500 | 1.016 / 0.8 |
| 500–700 | 1.002 / 1.0 |
| 700–1000 | 1.027 / 1.7 |
| 1000–1300 | 1.028 / 1.1 |
| 1500–1750 | 1.017 / 0.6 |

- **Gases:** the wavelength shift rises from +0.30 nm in the near-infrared to +0.50 nm at
  2330 nm. Scaling CO2 to 2020 (412 ppm) cut the spread against the ground at 2000–2300 nm
  from 4.6 % to 1.3 %.
- **The 2300 nm dip is real surface:** RadCalNet's ground reflectance drops 8.3 % from 2270 to
  2300 nm.
- **Aerosol barely matters over a bright playa:** AOD 0.053 vs 0.234 changes reflectance by
  under 1 %.
- **Sentinel-2B the same day:** RadCalNet is 0–5 % brighter than Sentinel-2, and DUSTRAN tracks
  RadCalNet.

## 2. AVIRIS Classic, Railroad Valley, 2019-08-13

Details: `docs/classic_2019_playa.md`.
- Raw radiance is 5–16 % below what RadCalNet predicts, as a smooth hump deepest at
  1050–1600 nm. NASA's reflectance product is also 7–19 % low.
- Sentinel-2A (2019-08-12) agrees with RadCalNet within ~2 %, so the deficit is in the Classic
  radiance.
- Seams agree (A|B 0.994, B|C 1.001): the deficit is not a spectrometer step.

## 3. AVIRIS Classic, Ivanpah Playa, four flights vs Sentinel-2

| Flight | A\|B seam | B\|C seam | DUSTRAN / Sentinel-2 (B02 … B12) |
|---|---|---|---|
| 2017-12-07 | 0.993 | 1.005 | 0.86–0.90 |
| 2018-09-06 | 0.992 | 1.003 | cloud shadows and a ~1 km location error; excluded |
| 2020-10-05 | **0.954** | 0.997 | 0.72–0.92 |
| 2021-02-24 | **0.960** | **1.012** | 0.92–0.97 |

- **Classic is consistently dark,** by 5–12 % in most bands on every usable flight, but the size
  and spectral shape change from flight to flight. One coefficient does not fit them all.
- **Pixel-by-pixel regression** against Sentinel-2 across a 4 km box shows the darkness is mostly
  multiplicative (gain). Exception: 2020-10-05 has an additive deficit of about 8 % at 1610 nm.
- **The A|B seam moved by ~4 % between the 2019 and 2020 seasons** and stayed moved in 2021.
- No calibration correction has been adopted.

## 4. Water

| Scene | 940 nm | 1140 nm | Photometer |
|---|---|---|---|
| NG 2020, lines only | 1.53 cm | 1.56 cm | 1.19 cm |
| NG 2020, + foreign continuum (yardstick) | 1.36 | 1.32 | 1.19 |
| Classic 2019, lines only | 1.14 | 1.44 | 1.00 |
| Classic 2019, + foreign continuum (yardstick) | 1.03 | 1.21–1.24 | 1.00 |

- No band shift (±2 nm) or width change (90–130 %) reconciles 940 with 1140 on Classic.
- The water self-continuum (water dimers, `dimer.py`) changes nothing over dry deserts: it
  scales with water squared.
- A foreign continuum of the expected size moves both scenes 10–15 % toward the photometer, but
  the yardstick version worsens the in-band fit shape in three of four fits.
- Lorentzian line wings beyond 25 cm⁻¹ supply only a seventh of that absorption, so far wings
  are not the mechanism.
- **Measured foreign continuum (`dustran/foreign.py`, `experiments/foreign_continuum_fit.py`):**
  a collision band under every HITRAN water line, widths 3, 10 and 30 cm⁻¹, strength A fitted to
  the photometer on one day and applied blind to the other.

  | Strength that puts 940 water on the photometer | Classic 2019 | NG 2020 |
  |---|---|---|
  | width 30 cm⁻¹ | A ≈ 0.024 atm⁻¹ | A ≈ 0.049 atm⁻¹ |

  - **It does not transfer.** NG needs twice the strength Classic needs. Fitted on Classic, NG's
    940 water stays 12 % high; fitted on NG, Classic's drops 9 % low.
  - **The two bands disagree differently on each day:** Classic needs more continuum at 1140
    than at 940; NG needs more at 940 than at 1140.
  - **The in-band fit shape gets worse** as the continuum grows, at every width.
  - **Conclusion:** a single universal absorber can't explain the water bias. What is left
    differs by day or instrument, which points to the water profile shape (the standard
    atmosphere is wrong differently on each day) or the band response, not to missing
    spectroscopy. Not adopted.
- Remaining suspects: the water profile shape (no per-scene profile yet), the line shape inside
  25 cm⁻¹, and the photometer's own ~10 % uncertainty.

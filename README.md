# DUSTRAN

Line-by-line atmospheric correction for imaging spectrometers over desert scenes, in plain
Python/numpy. No containers, no stored lookup tables, and every step can be inspected.

- **Gases:** HITRAN line by line (Voigt profiles, 0.002 cm⁻¹), per gas, so water can be
  rescaled per pixel without recomputing.
- **Collision-induced and cross-section absorbers:** O2–O2 and O2–N2 (HITRAN CIA), ozone
  (Serdyuchenko et al. 2014).
- **Scattering:** our own successive-orders-of-scattering solver (Fourier azimuth modes, Gauss
  streams, delta-M), Rayleigh from Bodhaine et al. (1999), dust as the only aerosol with our own
  Mie code.
- **Instrument:** Gaussian band responses with a fitted wavelength shift, and a seam check that
  measures the relative calibration of multi-spectrometer sensors from the scene itself.
- **Validation:** real observations only. AVIRIS-NG and AVIRIS Classic over Railroad Valley and
  Ivanpah playas, against RadCalNet ground reflectance, sun photometer measurements and
  Sentinel-2 surface reflectance.

License: MIT (`LICENSE`). Data sources keep their own terms (`PROVENANCE.md`).

## Getting started

### 1. Install

Python 3.11 or newer.

```
git clone https://github.com/aseever/DUSTRAN.git
cd DUSTRAN
python -m venv .venv
.venv/Scripts/activate            # Windows;  source .venv/bin/activate on macOS/Linux
pip install -r requirements.txt
python -m pytest tests            # 7 physics tests: energy conservation, Mie cases, delta-M ...
```

### 2. Choose a data folder

DUSTRAN never ships data. Pick a folder for downloads and point `DUSTRAN_DATA` at it (default:
`~/dustran_data`):

```
setx DUSTRAN_DATA "D:\dustran_data"          # Windows (new shells)
export DUSTRAN_DATA=~/dustran_data           # macOS/Linux
```

### 3. Get the reference data

| What | Where to get it | Put it in `DUSTRAN_DATA/` |
|---|---|---|
| HITRAN line lists | automatic: downloaded on first use to `~/.radisdb` | — |
| O2–O2 and O2–N2 collision-induced absorption | hitran.org/cia (`O2-O2_2024.cia`, `O2-N2_2024.cia`) | `hitran_xsec/` |
| Ozone cross sections | IUP Bremen spectroscopy page (Serdyuchenko et al. 2014), `serdyuchenkogorshelev5digits.dat` | `bremen_o3/` |
| Solar spectrum | TSIS-1 HSRS, LASP (Coddington et al. 2021), 0.1 nm, two columns: nm, irradiance | `tsis/tsis_f0_0p1.txt` |
| Standard atmosphere | AFGL US standard (Anderson et al. 1986), columns z p T air o3 o2 h2o co2 no2 | `afgl/afglus.dat` |

### 4. First run: check an AVIRIS Classic scene's spectrometer seams

Download any AVIRIS Classic radiance flight (NASA AVIRIS data portal, `fYYMMDDt01pXXrYY.tar.gz`),
unpack it, and pick a ground point in the scene's UTM zone:

```
python experiments/classic_site.py <unpacked>/fYYMMDDt01pXXrYY_rdn_X <utm_easting> <utm_northing>
```

It prints the sun and view geometry and the ratio at each spectrometer seam. Seams within 1 % of
1.0 mean the spectrometers agree; larger ratios are flagged.

### 5. Full validation

The experiments in `experiments/` reproduce every result in `RESULTS.md`. They need the
validation data listed in `PROVENANCE.md` (RadCalNet, AERONET, Sentinel-2, the AVIRIS scenes in
`docs/validation_downloads.md`). Line-by-line gas optical depths are computed once per site and
sensor altitude (`experiments/playa_gas_tau.py`, 15–30 min). Then:

```
python validation/suite.py          # every validated day vs its references, compared with the baseline
```

Cases whose data are missing are skipped.

### 6. Correct a whole flight line

`dustran/batch.py` corrects every pixel of an orthorectified radiance cube at once (numpy, in
row blocks; no per-pixel Python loops). First compute the line-by-line gas optical depths for the
scene's median ground elevation and the sensor altitude (step 5), then:

```
python -m dustran.batch RDN_img.hdr RDN_obs_ort.hdr OUT/prefix --taudirs TAU_VNIR,TAU_SWIR \
    --ground 1.52 --sensor 4.64 --o3-du 282 \
    --aeronet DUSTRAN_DATA/aeronet_... --inv-time HH:MM:SS --sun-time HH:MM:SS --date DD:MM:YYYY
```

Per scene it builds a small gas table (band-averaged sunlight x gas transmittance for ~40 water
amounts and 3 view angles) and the scattering terms for 9 view geometries; per pixel it fits water
at 940 nm (1140 nm as a cross-check), then inverts the radiance equation. Writes an ENVI
reflectance cube with the input's band table and map info, a water map and a JSON summary
(water percentiles, 940/1140 agreement flag, seams). An AVIRIS-NG line (3681 x 737 px) takes
~4 min of set-up and ~7 min of correction.

Status: research code. See `STACK.md` for the design and each layer's status, `RESULTS.md` for
what has been validated, `PROVENANCE.md` for where every piece came from and its license, and
`BIBLIOGRAPHY.md` for the published work it is built on.

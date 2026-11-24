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

Status: research code. See `STACK.md` for the design and each layer's status, `RESULTS.md` for
what has been validated, `PROVENANCE.md` for where every piece came from and its license, and
`BIBLIOGRAPHY.md` for the published work it is built on.

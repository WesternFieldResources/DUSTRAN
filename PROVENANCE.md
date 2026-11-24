# Provenance and licenses

Every piece of DUSTRAN, where it came from, and its license status. All code in this repository
was written for DUSTRAN from published papers and first principles (see `BIBLIOGRAPHY.md`). No
third-party code is vendored, and no third-party data are bundled.

## Code

| Module | Origin |
|---|---|
| `dustran/atmosphere.py` | Layering of standard-atmosphere profiles; exponential-profile column integral derived from first principles. CO2 scaled to the flight year from published global-mean records. |
| `dustran/gas.py` | Line-by-line optical depth per gas and layer. Voigt cross sections are computed by RADIS (imported, unmodified) from HITRAN line lists. |
| `dustran/cia.py` | HITRAN `.cia` reader and collision-induced optical depth. |
| `dustran/xsec.py` | Temperature-dependent cross-section tables (ozone). |
| `dustran/rayleigh.py` | Bodhaine et al. (1999). |
| `dustran/sos.py` | Successive orders of scattering with Fourier azimuth decomposition and delta-M (Hansen & Travis 1974; Wiscombe 1977). |
| `dustran/mie.py` | Bohren & Huffman (1983), ch. 4; checked against published cases. |
| `dustran/aerosol.py` | Lognormal and measured size distributions. Refractive index values are provisional representative numbers, to be replaced with a cited table (Di Biagio et al. 2019). |
| `dustran/instrument.py` | Gaussian band response, wavelength shift, spectrometer seam check. |
| `dustran/dimer.py` | Water-dimer self-continuum model of Simonova et al. (eqs. 3 and 5, Tables 2 and 3). The bound-dimer temperature extrapolation is our own assumption and is flagged in the module. |

## Dependencies

| Package | License | Use |
|---|---|---|
| numpy | BSD-3 | everywhere |
| scipy | BSD-3 | Legendre functions, splines |
| RADIS | LGPL-3 | imported unmodified for Voigt cross sections. LGPL places no conditions on DUSTRAN's own license as long as RADIS itself is not modified or vendored. |
| netCDF4 | MIT | reading coefficient files |
| matplotlib | PSF-based | plots in experiments only |

## Data (never bundled; each user downloads from the source)

| Data | Source | Terms |
|---|---|---|
| HITRAN line lists | hitran.org (Gordon et al. 2022), downloaded by RADIS to `~/.radisdb` | free to use; cite HITRAN |
| HITRAN CIA (O2–O2, O2–N2) | hitran.org/cia (Karman et al. 2019) | as HITRAN |
| Ozone cross sections | IUP Bremen (Serdyuchenko et al. 2014) | cite the paper; confirm terms before release |
| TSIS-1 solar spectrum | LASP (Coddington et al. 2021) | public; cite the paper |
| AFGL standard atmospheres | Anderson et al. (1986), US government work | public domain |
| AERONET v3 | NASA GSFC (Holben et al. 1998) | public; acknowledge site PIs |
| RadCalNet RVUS | radcalnet.org (Bouvet et al. 2019) | user account; RadCalNet terms |
| AVIRIS-NG / AVIRIS Classic | NASA | public |
| Sentinel-2 L2A | Copernicus | free and open |

## Before public release
- Confirm the Bremen ozone and HITRAN terms in writing.
- Fetch the TSIS-1 HSRS directly from LASP and the AFGL tables from the original report, and
  record checksums.


"""Crop a Sentinel-2 L2A zip to a 201 x 201 box of 20 m pixels around a point (ENVI float32).

    python experiments/s2_crop.py <SAFE.zip> <out_dir> <utm_e> <utm_n>

Reads straight from the zip with GDAL (/vsizip/); bands B02-B12 at 20 m plus SCL, AOT, WVP.
Values are left as stored: reflectance = (DN - 1000) / 10000 for N0400+ processing.
"""
import os, subprocess, sys, zipfile
from pathlib import Path

GDAL = r"C:\Program Files\QGIS 3.44.12\bin"
BANDS = ["B02", "B03", "B04", "B05", "B06", "B07", "B8A", "B11", "B12", "SCL", "AOT", "WVP"]

if __name__ == "__main__":
    z, out, e, n = sys.argv[1], Path(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
    out.mkdir(parents=True, exist_ok=True)
    names = zipfile.ZipFile(z).namelist()
    env = dict(os.environ, PATH=GDAL + os.pathsep + os.environ["PATH"])
    h = 100.5 * 20
    for b in BANDS:
        src = next(x for x in names if "/R20m/" in x and x.endswith(f"_{b}_20m.jp2"))
        subprocess.run([GDAL + r"\gdal_translate.exe", "-q", "-of", "ENVI", "-ot", "Float32",
                        "-projwin", str(e - h), str(n + h), str(e + h), str(n - h),
                        f"/vsizip/{z}/{src}", str(out / f"{b}.img")], check=True, env=env)
    print("cropped", len(BANDS), "bands to", out)

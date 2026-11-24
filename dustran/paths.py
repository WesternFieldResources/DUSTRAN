"""Where DUSTRAN keeps downloaded data (never in the repository).

Set the environment variable DUSTRAN_DATA to choose the folder. Without it, DUSTRAN uses
D:/workspace/dustran_data if that exists, else ~/dustran_data.
"""
import os
from pathlib import Path


def _data_dir():
    if os.environ.get("DUSTRAN_DATA"):
        return Path(os.environ["DUSTRAN_DATA"])
    legacy = Path("D:/workspace/dustran_data")
    return legacy if legacy.exists() else Path.home() / "dustran_data"


DATA = _data_dir()

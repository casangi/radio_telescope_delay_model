"""FITS-IDI extraction and uvw comparison against real EHT archive data.

The dataset is a compact extract (3000 UV rows) from the DiFX-correlated
EHT 2021 track e21d15 (job E21D15.1, 3C273, band 1); its uvw were written by
difx2fits from the difxcalc11 ``.im`` model in float32 light-seconds, so
this pins the embedded difxcalc11 pipeline against a real archive at the
file's own precision floor (measured: median 0.056 m, max 0.168 m).
"""

import os

import numpy as np
import pytest

from radio_telescope_delay_model.delay_model_cpp import calc_available
from radio_telescope_delay_model.fits_idi import compare_uvw, read_fits_idi

DATA = os.path.join(os.path.dirname(__file__), "data", "eht_e21d15_3c273_uvw.npz")


def test_read_fits_idi_extract():
    data = read_fits_idi(DATA)
    assert data["station_name"] == ["AA", "AX", "GL", "KT", "MG", "MM", "SW"]
    assert data["antenna_position"].shape == (7, 3)
    # Mixed alt-az and Nasmyth (4/5) mounts all map to AZEL.
    assert set(data["mount_type"]) == {"AZEL"}
    assert data["source_name"] == ["3C273"]
    assert len(data["time"]) == len(data["uvw"]) == len(data["antenna1"])
    assert np.all(data["antenna1"] != data["antenna2"])  # autos dropped
    assert 1.0e6 < np.abs(data["uvw"]).max() < 1.3e7  # Earth-scale baselines


@pytest.mark.skipif(not calc_available(), reason="CALC11 extension not built")
def test_difxcalc11_reproduces_eht_archive_uvw():
    data = read_fits_idi(DATA)
    result = compare_uvw(data, method="difxcalc11")
    residual = np.abs(result["residual"])
    assert np.median(residual) < 0.1  # the file's float32 quantization floor
    assert residual.max() < 0.5
    # The wrong orientation misses by the full uvw scale.
    flipped = np.abs(-result["model_uvw"] - result["file_uvw"]).max()
    assert flipped > 1.0e6


@pytest.mark.skipif(not calc_available(), reason="CALC11 extension not built")
def test_other_methods_show_their_conventions():
    """The geometric CALC recipe and the pure projection differ from the
    DiFX archive by their documented conventions (atmosphere-in-uvw and
    annual aberration respectively), not by orientation."""
    data = read_fits_idi(DATA)
    calc = np.abs(compare_uvw(data, method="calc")["residual"])
    astropy_res = np.abs(compare_uvw(data, method="astropy")["residual"])
    assert 1.0 < np.median(calc) < 50.0
    assert 5.0 < np.median(astropy_res) < 300.0
    with pytest.raises(ValueError, match="method"):
        compare_uvw(data, method="bogus")

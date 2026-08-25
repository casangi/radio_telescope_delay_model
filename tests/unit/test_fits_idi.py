"""FITS-IDI extraction and uvw comparison against real EHT archive data.

The dataset is a compact extract (3000 UV rows) from the DiFX-correlated
EHT 2021 track e21d15 (job E21D15.1, 3C273, band 1); its uvw were written by
difx2fits from the difxcalc11 ``.im`` model in float32 light-seconds, so
this pins the embedded difxcalc11 pipeline against a real archive at the
file's own precision floor (float32 light-seconds: ~1 m ulp on the longest
baselines, plus the ~1.8 m ALMA-position inconsistency of this delivery;
measured per-component median ~0.09 m, max ~2.1 m across all 21 baselines --
EOP choice moves these by under a centimeter).
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
    # float32 ulp reaches ~1 m at the longest-baseline uvw components, and
    # the EOP values used at correlation time differ from IERS-B finals.
    assert residual.max() < 2.5
    # The wrong orientation misses by the full uvw scale.
    flipped = np.abs(-result["model_uvw"] - result["file_uvw"]).max()
    assert flipped > 1.0e6


@pytest.mark.skipif(not calc_available(), reason="CALC11 extension not built")
def test_atmosphere_free_alma_reaches_storage_floor():
    """With the production convention -- the file's own EOPs and an
    atmosphere-free phased ALMA (its real-time phasing removes the
    troposphere, so the correlator model must not add it) -- every baseline
    reaches the float32 storage floor with no fitted parameters at all,
    including the 2.6 km AA-AX baseline that the plain model misses by the
    full ALMA line-of-sight troposphere (~1.8 m at this scan's elevation)."""
    data = read_fits_idi(DATA)
    result = compare_uvw(
        data, method="difxcalc11", eop=data["eop"], atmosphere_free_station=["AA"]
    )
    residual = np.linalg.norm(result["residual"], axis=1)
    aa_ax = ((data["antenna1"] == 0) & (data["antenna2"] == 1)) | (
        (data["antenna1"] == 1) & (data["antenna2"] == 0)
    )
    # Vector-norm metrics; the medians sit at the float32 storage floor
    # (measured 0.136 median / 0.373 max / 0.09 mm AA-AX).
    assert np.median(residual) < 0.2
    assert residual.max() < 0.45
    assert np.median(residual[aa_ax]) < 0.005  # sub-5 mm on the 2.6 km baseline


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

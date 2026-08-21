"""Tests for the difxcalc11 station catalogs and named-telescope support.

"All the telescopes" that difxcalc11 / pycalc11 support are the stations of
the catalogs both packages ship (this package packages the same files);
every one of them is exercised here -- at the parsing level, and through the
CALC delay model with its ocean loading applied.
"""

import numpy as np
import pytest
from astropy.time import Time

from radio_telescope_delay_model import (
    MOUNT_TYPE_TO_AXIS_CODE,
    calculate_delays_calc,
    telescope_names,
)
from radio_telescope_delay_model.delay_model_cpp import calc_available
from radio_telescope_delay_model.station_data import (
    axis_tilt_catalog,
    ocean_loading_catalog,
    ocean_pole_tide_catalog,
    station_geophysics,
)

TIME = Time(["2019-10-03T19:13:20.000", "2019-10-03T19:40:00.000"], scale="utc")
PHASE_CENTER = np.array([[5.233697011, -0.710938054]])  # ICRS radians

ALMA_POSITION = np.array(
    [
        [2225061.164, -5440061.789, -2481681.151],
        [2224993.209, -5440087.586, -2481675.855],
        [2225003.052, -5440024.184, -2481806.733],
    ]
)


def test_all_telescopes_catalog():
    """Every telescope in the ocean loading catalog parses with sane values."""
    names = telescope_names()
    assert len(names) > 300  # the difxcalc11 catalog carries ~364 stations
    catalog = ocean_loading_catalog()
    for name in names:
        entry = catalog[name]
        assert entry["vertical_amplitude"].shape == (11,)
        assert entry["horizontal_amplitude"].shape == (2, 11)
        assert entry["vertical_phase"].shape == (11,)
        assert entry["horizontal_phase"].shape == (2, 11)
        for key in entry:
            if key != "name":
                assert np.all(np.isfinite(entry[key])), name
        # BLQ amplitudes are metres; even extreme coastal sites are < 20 cm.
        assert entry["vertical_amplitude"].max() < 0.2, name
        assert np.abs(entry["horizontal_amplitude"]).max() < 0.2, name
        assert np.all(np.abs(entry["vertical_phase"]) <= np.pi + 1e-9), name
        assert np.all(np.abs(entry["horizontal_phase"]) <= np.pi + 1e-9), name


def test_all_telescopes_pole_tide():
    entries = ocean_pole_tide_catalog()
    assert len(entries) > 300
    for entry in entries:
        assert entry["coefficients"].shape == (6,)
        assert np.all(np.isfinite(entry["coefficients"])), entry["name"]
        assert np.abs(entry["coefficients"]).max() < 1.0, entry["name"]
        assert -90.0 <= entry["latitude_deg"] <= 90.0, entry["name"]
        assert 0.0 <= entry["longitude_deg"] <= 360.0, entry["name"]


def test_station_geophysics_known_and_unknown():
    with pytest.warns(UserWarning, match="NOSUCHSTATION"):
        arrays = station_geophysics(["ONSALA60", "NOSUCHSTATION"])
    # ONSALA60's first vertical amplitude straight from the BLQ block.
    assert arrays["ocean_vertical_amplitude"][0, 0] == pytest.approx(0.00351)
    assert arrays["ocean_pole_tide_coefficients"][0].any()
    np.testing.assert_array_equal(arrays["ocean_vertical_amplitude"][1], 0.0)
    np.testing.assert_array_equal(arrays["ocean_pole_tide_coefficients"][1], 0.0)


def test_axis_tilt_interpolation():
    """dANTILT semantics: clamped before the first epoch, linear between
    bracketing epochs, constant for single-epoch stations."""
    early = axis_tilt_catalog(1980.0)["PIETOWN"]
    assert early == pytest.approx((0.14, 0.19))  # clamped to 1988.666 entry
    fraction = (1988.733 - 1988.666) / (1988.800 - 1988.666)
    mid = axis_tilt_catalog(1988.733)["PIETOWN"]
    assert mid[0] == pytest.approx(0.14 + fraction * (0.26 - 0.14))
    assert mid[1] == pytest.approx(0.19 + fraction * (0.07 - 0.19))
    # Single-epoch stations are constant at any epoch, both name columns.
    assert axis_tilt_catalog(2050.0)["SC-VLBA"] == pytest.approx((1.04, 1.58))
    assert axis_tilt_catalog(1980.0)["SC"] == pytest.approx((1.04, 1.58))


def _positions_from_pole_tide(names):
    """Approximate geocentric positions (spherical Earth) for every name
    with latitude/longitude in the ocean pole tide catalog."""
    by_name = {}
    for entry in ocean_pole_tide_catalog():
        by_name.setdefault(entry["name"], entry)
    kept, positions = [], []
    for name in names:
        entry = by_name.get(name)
        if entry is None:
            continue
        lat = np.deg2rad(entry["latitude_deg"])
        lon = np.deg2rad(entry["longitude_deg"])
        positions.append(
            6371.0e3
            * np.array(
                [np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)]
            )
        )
        kept.append(name)
    return kept, np.asarray(positions)


@pytest.mark.skipif(not calc_available(), reason="CALC11 extension not built")
def test_all_telescopes_through_delay_model():
    """Every catalog telescope with coordinates runs through CALC with its
    ocean loading applied, producing finite delays, and the loading moves
    the geometric delay by a nonzero sub-decimetre amount."""
    names, positions = _positions_from_pole_tide(telescope_names())
    assert len(names) > 300
    effects = []
    chunk = 50
    for start in range(0, len(names), chunk):
        chunk_names = names[start : start + chunk]
        chunk_positions = positions[start : start + chunk]
        reference = chunk_positions[0]
        with_loading = calculate_delays_calc(
            chunk_positions,
            TIME,
            PHASE_CENTER,
            reference_position=reference,
            station_name=chunk_names,
        )
        without_loading = calculate_delays_calc(
            chunk_positions, TIME, PHASE_CENTER, reference_position=reference
        )
        for key in ("geometric_delay", "dry_delay", "wet_delay"):
            assert np.all(np.isfinite(with_loading[key]))
        delta = 299792458.0 * np.abs(
            with_loading["geometric_delay"] - without_loading["geometric_delay"]
        )
        effects.extend(delta.max(axis=0))  # per antenna, max over epochs
    effects = np.asarray(effects)
    assert effects.max() < 0.5  # ocean loading is a centimetre-level effect
    assert effects.max() > 1.0e-4
    assert np.mean(effects > 1.0e-8) > 0.8  # applied for nearly every station


@pytest.mark.skipif(not calc_available(), reason="CALC11 extension not built")
def test_mount_types():
    """All difxcalc11 mount types run; with a zero axis offset every mount is
    exactly the alt-az result (the axis-offset module is inert), and with a
    nonzero offset each mount type moves the delay its own way."""
    azel = calculate_delays_calc(ALMA_POSITION, TIME, PHASE_CENTER)
    for mount in MOUNT_TYPE_TO_AXIS_CODE:
        results = calculate_delays_calc(
            ALMA_POSITION, TIME, PHASE_CENTER, mount_type=mount
        )
        np.testing.assert_array_equal(
            results["geometric_delay"], azel["geometric_delay"]
        )

    offset = np.full(len(ALMA_POSITION), 2.0)
    azel_offset = calculate_delays_calc(
        ALMA_POSITION, TIME, PHASE_CENTER, axis_offset_metres=offset
    )
    assert (
        299792458.0
        * np.abs(azel_offset["geometric_delay"] - azel["geometric_delay"]).max()
        > 1.0e-4
    )
    for mount in ("EQUA", "XYNS", "XYEW", "RICH"):
        results = calculate_delays_calc(
            ALMA_POSITION,
            TIME,
            PHASE_CENTER,
            axis_offset_metres=offset,
            mount_type=mount,
        )
        difference = 299792458.0 * np.abs(
            results["geometric_delay"] - azel_offset["geometric_delay"]
        )
        assert difference.max() > 1.0e-4, mount

    with pytest.raises(ValueError, match="mount_type"):
        calculate_delays_calc(ALMA_POSITION, TIME, PHASE_CENTER, mount_type="BOGUS")

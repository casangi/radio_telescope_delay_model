"""The package is array-general: connected arrays through global VLBI.

The delay model takes every telescope property per antenna (ITRF position,
station name for the geophysical catalogs, mount type, axis offset,
weather), so nothing is ALMA-specific. These tests exercise the array
classes the package must serve -- a compact connected array (ALMA/VLA
style), a continental VLBI network (VLBA/EVN style, real catalog stations),
and an Earth-spanning network (EHT scale) -- through the astropy method,
the CALC method (geocenter reference, the VLBI convention) and the embedded
difxcalc11 pipeline, checking mutual consistency at scale-appropriate
levels.

Station coordinates for the VLBI cases are derived from the packaged ocean
pole tide catalog's latitude/longitude (spherical Earth): metre-exact
coordinates are irrelevant here -- the point is real cataloged stations at
real geographic scales.
"""

import numpy as np
import pytest
from astropy.time import Time

from radio_telescope_delay_model import (
    calculate_uvw_astropy,
    calculate_uvw_calc,
)
from radio_telescope_delay_model.delay_model_cpp import calc_available
from radio_telescope_delay_model.station_data import ocean_pole_tide_catalog

pytestmark = pytest.mark.skipif(
    not calc_available(), reason="CALC11 extension not built"
)

TIME = Time(["2019-10-03T19:13:20.000", "2019-10-03T19:40:00.000"], scale="utc")
PHASE_CENTER = np.array([[5.233697011, -0.710938054]])  # ICRS radians
NORTHERN_PHASE_CENTER = np.array([[1.5, 0.7]])  # for the northern networks

VLBA_STATIONS = ["MK-VLBA", "PIETOWN", "LA-VLBA", "FD-VLBA", "HN-VLBA"]
EVN_STATIONS = ["EFLSBERG", "MEDICINA", "ONSALA60", "WETTZELL", "YEBES40M"]
GLOBAL_STATIONS = ["MK-VLBA", "WESTFORD", "WETTZELL", "HARTRAO"]


def _catalog_positions(names):
    by_name = {entry["name"]: entry for entry in ocean_pole_tide_catalog()}
    positions = []
    for name in names:
        entry = by_name[name]
        lat = np.deg2rad(entry["latitude_deg"])
        lon = np.deg2rad(entry["longitude_deg"])
        positions.append(
            6371.0e3
            * np.array(
                [np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)]
            )
        )
    return np.asarray(positions)


def _check_array(names, phase_center, min_baseline):
    """Run one array through all three methods and check consistency.

    calc vs astropy is asserted RELATIVELY at 2e-4: the two carry different
    uvw conventions at the annual-aberration level (v/c ~ 1e-4; the CALC
    methods derive uvw from aberrated delays, the astropy method is a pure
    geometric projection), measured at 0.4-1.0e-4 relative across these
    arrays. The difxcalc11 mode is checked for finiteness and convention
    only: its exact-mode uvw carry the atmosphere, which is unbounded for
    stations that see the source below the horizon -- real difxcalc
    behaviour for global networks (correlators drop those rows).
    """
    positions = _catalog_positions(names)
    geocenter = np.zeros(3)

    uvw_calc, antenna1, antenna2 = calculate_uvw_calc(
        positions,
        TIME,
        phase_center,
        reference_position=geocenter,
        station_name=names,
    )
    uvw_astropy, a1, a2 = calculate_uvw_astropy(positions, TIME, phase_center)
    uvw_difx, d1, d2 = calculate_uvw_calc(
        positions, TIME, phase_center, mode="difxcalc11", station_name=names
    )
    np.testing.assert_array_equal(antenna1, a1)
    np.testing.assert_array_equal(antenna2, d2)

    scale = np.abs(uvw_calc).max()
    assert scale > min_baseline  # the scenario really is this large
    assert np.all(np.isfinite(uvw_calc))
    assert np.all(np.isfinite(uvw_difx))
    np.testing.assert_allclose(uvw_calc, uvw_astropy, atol=2.0e-4 * scale)


def test_continental_vlba_network():
    _check_array(VLBA_STATIONS, PHASE_CENTER, min_baseline=1.0e6)


def test_european_evn_network():
    _check_array(EVN_STATIONS, NORTHERN_PHASE_CENTER, min_baseline=5.0e5)


def test_earth_spanning_network():
    """EHT-scale baselines (Hawaii - Europe - South Africa)."""
    _check_array(GLOBAL_STATIONS, PHASE_CENTER, min_baseline=8.0e6)


def test_difxcalc11_mode_high_elevation_network():
    """With the source high above a regional network, the difxcalc11 mode
    and the geocenter-referenced geometric CALC method agree to the
    atmosphere-in-uvw / recipe level."""
    from astropy.coordinates import ITRS, SkyCoord

    positions = _catalog_positions(EVN_STATIONS)
    centroid = positions.mean(axis=0)
    up = centroid / np.linalg.norm(centroid)
    zenith = SkyCoord(
        ITRS(
            x=up[0], y=up[1], z=up[2], representation_type="cartesian", obstime=TIME[0]
        )
    ).transform_to("icrs")
    phase_center = np.array([[zenith.ra.rad, zenith.dec.rad]])

    uvw_calc, _, _ = calculate_uvw_calc(
        positions,
        TIME,
        phase_center,
        reference_position=np.zeros(3),
        station_name=EVN_STATIONS,
    )
    uvw_difx, _, _ = calculate_uvw_calc(
        positions,
        TIME,
        phase_center,
        mode="difxcalc11",
        station_name=EVN_STATIONS,
    )
    difference = np.abs(uvw_difx - uvw_calc).max()
    assert difference < 20.0  # atmosphere-in-uvw scale at ~2000 km baselines
    assert difference > 1e-4  # the recipes genuinely differ


def test_buried_reference_warning():
    """Defaulting the reference for an Earth-spanning array warns and points
    at the geocenter convention; an explicit geocenter does not warn."""
    positions = _catalog_positions(GLOBAL_STATIONS)
    with pytest.warns(UserWarning, match="geocenter"):
        calculate_uvw_calc(positions, TIME, PHASE_CENTER)
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        calculate_uvw_calc(
            positions, TIME, PHASE_CENTER, reference_position=np.zeros(3)
        )


def test_mixed_known_and_unknown_stations():
    """An EHT-style mix -- geodetic stations plus mm sites absent from the
    geodetic catalogs -- runs, with zeros (and a warning) for the unknown
    names, exactly as difxcalc11 behaves."""
    names = ["WESTFORD", "WETTZELL", "ALMA-EHT"]
    positions = np.vstack(
        [
            _catalog_positions(names[:2]),
            np.array([[2225061.164, -5440061.789, -2481681.151]]),
        ]
    )
    with pytest.warns(UserWarning, match="ALMA-EHT"):
        uvw, antenna1, antenna2 = calculate_uvw_calc(
            positions,
            TIME,
            PHASE_CENTER,
            reference_position=np.zeros(3),
            station_name=names,
        )
    assert np.all(np.isfinite(uvw))

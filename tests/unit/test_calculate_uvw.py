"""Unit tests for the uvw calculations, based on the package's example notebook.

The scenarios mirror ``example_notebooks/Radio_Telescope_Delay_Model_Example.ipynb``:
per-antenna uvw from the astropy construction, baseline uvw from CALC11 delays
(DiFX finite-difference recipe), their mutual agreement, and -- when pycalc11
is installed -- a cross-check against CALC's native UBASE/VBASE/WBASE from an
independent wrapper of the same model.

All baselines are in the archival / VLBI convention adopted by MSv4:
``uvw = P(antenna1) - P(antenna2)`` with W towards the source.
"""

import importlib.util

import numpy as np
import pytest
from astropy.time import Time

from radio_telescope_delay_model import (
    calculate_antenna_uvw_astropy,
    calculate_delays_calc,
    calculate_uvw_astropy,
    calculate_uvw_calc,
    earth_orientation_parameters,
)
from radio_telescope_delay_model.delay_model_cpp import calc_available

pycalc11_missing = importlib.util.find_spec("pycalc11") is None

# Six ALMA 12 m antennas (ITRF, metres) and a southern phase centre -- the
# example notebook's setting (TW Hya field, ALMA site).
ANTENNA_POSITION = np.array(
    [
        [2225061.164, -5440061.789, -2481681.151],
        [2224993.209, -5440087.586, -2481675.855],
        [2225003.052, -5440024.184, -2481806.733],
        [2225081.906, -5439969.101, -2481789.559],
        [2225030.837, -5440068.099, -2481654.284],
        [2224981.091, -5440131.383, -2481555.220],
    ]
)
TIME = Time(
    [
        "2019-10-03T19:13:20.000",
        "2019-10-03T19:40:00.000",
        "2019-10-03T20:06:40.000",
    ],
    scale="utc",
)
PHASE_CENTER = np.array([[5.233697011, -0.710938054]])  # ICRS radians


def test_astropy_convention_pin():
    """Baselines are the archival ``P(antenna1) - P(antenna2)`` of the
    per-antenna projections (exact by construction)."""
    antenna_uvw = calculate_antenna_uvw_astropy(ANTENNA_POSITION, TIME, PHASE_CENTER)
    uvw, antenna1, antenna2 = calculate_uvw_astropy(
        ANTENNA_POSITION, TIME, PHASE_CENTER
    )
    assert uvw.shape == (len(TIME), 15, 3)
    np.testing.assert_array_equal(
        uvw, antenna_uvw[:, antenna1, :] - antenna_uvw[:, antenna2, :]
    )
    # w of a per-antenna projection is the position component towards the
    # source: for these southern-sky times it is bounded by the antenna radius.
    assert np.all(np.abs(antenna_uvw[..., 2]) < 6.5e6)


def test_astropy_accepts_strings_and_unix():
    uvw_time_object, _, _ = calculate_uvw_astropy(ANTENNA_POSITION, TIME, PHASE_CENTER)
    uvw_strings, _, _ = calculate_uvw_astropy(
        ANTENNA_POSITION, [str(t) for t in TIME.isot], PHASE_CENTER
    )
    uvw_unix, _, _ = calculate_uvw_astropy(ANTENNA_POSITION, TIME.unix, PHASE_CENTER)
    np.testing.assert_allclose(uvw_strings, uvw_time_object, atol=1e-8)
    np.testing.assert_allclose(uvw_unix, uvw_time_object, atol=1e-6)


@pytest.mark.skipif(not calc_available(), reason="CALC11 extension not built")
def test_calc_and_astropy_methods_agree():
    """The two independent uvw constructions agree to geodetic-model level.

    CALC11 (full geodetic VLBI model) and the astropy projection differ by
    their earth models (nutation, tides, aberration handling); on these
    ~230 m baselines the measured difference is ~2 cm (1e-4 relative), far
    below the hundreds-of-metres uvw signal and consistent with the package's
    example-notebook comparison.
    """
    uvw_astropy, a1_astropy, a2_astropy = calculate_uvw_astropy(
        ANTENNA_POSITION, TIME, PHASE_CENTER
    )
    uvw_calc, a1_calc, a2_calc = calculate_uvw_calc(
        ANTENNA_POSITION, TIME, PHASE_CENTER
    )
    np.testing.assert_array_equal(a1_astropy, a1_calc)
    np.testing.assert_array_equal(a2_astropy, a2_calc)
    scale = np.abs(uvw_astropy).max()
    assert scale > 100  # metres; the comparison is not vacuous
    np.testing.assert_allclose(uvw_calc, uvw_astropy, atol=0.05)


@pytest.mark.skipif(not calc_available(), reason="CALC11 extension not built")
def test_calc_delays_shapes_and_magnitudes():
    results = calculate_delays_calc(ANTENNA_POSITION, TIME, PHASE_CENTER)
    for key in ("geometric_delay", "dry_delay", "wet_delay"):
        assert results[key].shape == (len(TIME), len(ANTENNA_POSITION))
    # Geometric delays are bounded by baseline / c.
    assert np.abs(results["geometric_delay"]).max() < 300 / 299792458.0
    # Atmospheric delays are positive and of ns order at ALMA altitude scale.
    assert np.all(results["dry_delay"] > 0)
    assert np.all(results["wet_delay"] > 0)


@pytest.mark.skipif(not calc_available(), reason="CALC11 extension not built")
@pytest.mark.skipif(pycalc11_missing, reason="pycalc11 not installed")
def test_calc_uvw_matches_pycalc11():
    """Cross-check against pycalc11's native CALC UBASE/VBASE/WBASE.

    Both packages wrap the same CALC11 Fortran; pycalc11 reads CALC's native
    per-antenna uvw (geocenter mode) while this package reconstructs uvw from
    delays and finite-difference direction derivatives (and pycalc11 samples
    on its own 24 s epoch grid), so agreement is at the recipe level
    (~1 cm on these baselines), not bit level.
    """
    import astropy.units as un
    from astropy.coordinates import EarthLocation, SkyCoord
    from pycalc11 import Calc

    stations = EarthLocation(
        x=ANTENNA_POSITION[:, 0] * un.m,
        y=ANTENNA_POSITION[:, 1] * un.m,
        z=ANTENNA_POSITION[:, 2] * un.m,
    )
    calc = Calc(
        station_names=[f"S{k:02d}" for k in range(len(ANTENNA_POSITION))],
        station_coords=stations,
        source_coords=SkyCoord(
            [PHASE_CENTER[0, 0]] * un.rad,
            [PHASE_CENTER[0, 1]] * un.rad,
            frame="icrs",
        ),  # fmt: skip
        start_time=TIME[0] - 60 * un.s,
        duration_min=int((TIME[-1] - TIME[0]).sec / 60) + 3,
        base_mode="geocenter",
        dry_atm=False,
        wet_atm=False,
    )
    calc.run_driver()
    pycalc_uvw = calc.uvw.to_value(un.m)  # (time, station, source, 3)
    pycalc_times = Time(calc.times).unix

    uvw_calc, antenna1, antenna2 = calculate_uvw_calc(
        ANTENNA_POSITION, TIME, PHASE_CENTER
    )
    n_matched = 0
    for our_index, t in enumerate(TIME.unix):
        k = int(np.argmin(np.abs(pycalc_times - t)))
        if abs(pycalc_times[k] - t) > 0.5:
            continue  # pycalc's 24 s epoch grid does not sample this time
        per_antenna = pycalc_uvw[k].reshape(len(ANTENNA_POSITION), 3)
        # pycalc's per-antenna values are delay-like (-P(r)); the archival
        # baseline is antenna2 - antenna1 of them.
        expected = per_antenna[antenna2] - per_antenna[antenna1]
        np.testing.assert_allclose(uvw_calc[our_index], expected, atol=0.03)
        n_matched += 1
    assert n_matched >= 1, "no pycalc11 epoch coincided with the test times"


def test_earth_orientation_parameters():
    eop = earth_orientation_parameters(TIME)
    for key in (
        "polar_motion_x_arcsec",
        "polar_motion_y_arcsec",
        "ut1_minus_utc_seconds",
        "leap_seconds",
    ):
        assert eop[key].shape == (len(TIME),)
        assert np.all(np.isfinite(eop[key]))
    np.testing.assert_array_equal(eop["leap_seconds"], 37.0)  # TAI - UTC in 2019
    assert np.abs(eop["polar_motion_x_arcsec"]).max() < 1.0
    assert np.abs(eop["ut1_minus_utc_seconds"]).max() < 1.0


# Epochs around the 2016-12-31 23:59:60 leap-second insertion (TAI - UTC
# steps 36 -> 37): the calc method must split such a request internally.
STRADDLE_TIME = Time(
    [
        "2016-12-31T21:00:00.000",
        "2016-12-31T23:00:00.000",
        "2017-01-01T01:00:00.000",
    ],
    scale="utc",
)


def test_earth_orientation_leap_step():
    eop = earth_orientation_parameters(STRADDLE_TIME)
    np.testing.assert_array_equal(eop["leap_seconds"], [36.0, 36.0, 37.0])


@pytest.mark.skipif(not calc_available(), reason="CALC11 extension not built")
def test_calc_straddles_leap_second_insertion():
    """A request across a leap-second insertion is split into constant-leap
    runs and stitched back in epoch order.

    Each returned epoch must be bit-identical to the same epoch computed by a
    request confined to its own run (identical run arrays by construction),
    and the stitched whole must still agree with the independent astropy
    method (a mis-stitched or mislabeled epoch would err at the 0.1 m level).
    """
    results = calculate_delays_calc(ANTENNA_POSITION, STRADDLE_TIME, PHASE_CENTER)
    np.testing.assert_array_equal(results["leap_seconds"], [36.0, 36.0, 37.0])

    uvw, antenna1, antenna2 = calculate_uvw_calc(
        ANTENNA_POSITION, STRADDLE_TIME, PHASE_CENTER
    )
    assert uvw.shape == (len(STRADDLE_TIME), 15, 3)

    # The leap-36 pair is one run of its own: identical to requesting it alone.
    uvw_before, _, _ = calculate_uvw_calc(
        ANTENNA_POSITION, STRADDLE_TIME[:2], PHASE_CENTER
    )
    np.testing.assert_array_equal(uvw[:2], uvw_before)

    # The leap-37 singleton is padded with the nearest epoch (23:00) for the
    # EOP-rate estimate: identical to requesting exactly that pair, which
    # splits into the same two runs.
    uvw_pair, _, _ = calculate_uvw_calc(
        ANTENNA_POSITION, STRADDLE_TIME[1:], PHASE_CENTER
    )
    np.testing.assert_array_equal(uvw[2], uvw_pair[1])

    # And the stitched result still matches the independent astropy method.
    uvw_astropy, _, _ = calculate_uvw_astropy(
        ANTENNA_POSITION, STRADDLE_TIME, PHASE_CENTER
    )
    np.testing.assert_allclose(uvw, uvw_astropy, atol=0.05)

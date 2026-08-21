"""Unit tests for the CALC11 pybind11 extension (thread safety, validation)."""

import threading

import numpy as np
import pytest
from astropy.time import Time

from radio_telescope_delay_model.calculate_uvw import _packaged_ephemeris_path
from radio_telescope_delay_model.delay_model_cpp import (
    calc_available,
    calc_delay_model,
)
from radio_telescope_delay_model.earth_orientation import (
    earth_orientation_parameters,
)

pytestmark = pytest.mark.skipif(
    not calc_available(), reason="CALC11 extension not built"
)

ANTENNA_POSITION = np.array(
    [
        [2225061.164, -5440061.789, -2481681.151],
        [2224993.209, -5440087.586, -2481675.855],
        [2225003.052, -5440024.184, -2481806.733],
    ]
)
TIME = Time(["2019-10-03T19:13:20.000", "2019-10-03T19:40:00.000"], scale="utc")


def _run():
    n_antenna = len(ANTENNA_POSITION)
    n_time = len(TIME)
    eop = earth_orientation_parameters(TIME)
    geometric = np.empty((n_time, n_antenna))
    dry = np.empty((n_time, n_antenna))
    wet = np.empty((n_time, n_antenna))
    calc_delay_model(
        ANTENNA_POSITION.mean(axis=0),
        ANTENNA_POSITION,
        np.zeros(n_antenna),
        np.full(n_antenna, 1013.25),
        np.full(n_antenna, 0.5),
        np.ascontiguousarray(TIME.mjd),
        np.full(n_time, 5.233697011),
        np.full(n_time, -0.710938054),
        eop["polar_motion_x_arcsec"],
        eop["polar_motion_y_arcsec"],
        eop["ut1_minus_utc_seconds"],
        float(eop["leap_seconds"][0]),
        np.zeros(n_antenna),
        _packaged_ephemeris_path(),
        geometric,
        dry,
        wet,
    )
    return geometric, dry, wet


def test_deterministic():
    first = _run()
    second = _run()
    for a, b in zip(first, second, strict=True):
        np.testing.assert_array_equal(a, b)


def test_thread_safety():
    """Concurrent calls are serialized internally (COMMON-block core) and
    produce exactly the sequential results."""
    reference = _run()
    results = [None] * 4

    def worker(k):
        results[k] = _run()

    threads = [threading.Thread(target=worker, args=(k,)) for k in range(4)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    for result in results:
        for a, b in zip(result, reference, strict=True):
            np.testing.assert_array_equal(a, b)


def test_input_validation():
    n_antenna = len(ANTENNA_POSITION)
    eop = earth_orientation_parameters(TIME)
    out = np.empty((len(TIME), n_antenna))
    with pytest.raises(Exception, match="antenna_position"):
        calc_delay_model(
            ANTENNA_POSITION.mean(axis=0),
            ANTENNA_POSITION[:, :2].copy(),
            np.zeros(n_antenna),
            np.full(n_antenna, 1013.25),
            np.full(n_antenna, 0.5),
            np.ascontiguousarray(TIME.mjd),
            np.full(len(TIME), 5.2),
            np.full(len(TIME), -0.7),
            eop["polar_motion_x_arcsec"],
            eop["polar_motion_y_arcsec"],
            eop["ut1_minus_utc_seconds"],
            float(eop["leap_seconds"][0]),
            np.zeros(n_antenna),
            _packaged_ephemeris_path(),
            out,
            out.copy(),
            out.copy(),
        )

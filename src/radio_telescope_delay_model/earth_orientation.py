"""Earth-orientation parameters for the CALC delay model, via astropy IERS."""

from __future__ import annotations

import numpy as np

__all__ = ["earth_orientation_parameters"]


def earth_orientation_parameters(time) -> dict:
    """Polar motion, UT1-UTC and leap seconds at the given times.

    Parameters
    ----------
    time : astropy.time.Time (array)

    Returns
    -------
    dict
        ``polar_motion_x_arcsec``, ``polar_motion_y_arcsec``,
        ``ut1_minus_utc_seconds`` (IERS-B, the table astropy ships) and
        ``leap_seconds`` (TAI - UTC, a step function that increments at
        leap-second insertions) -- each an ``[n_time]`` float64 array.
    """
    import erfa
    from astropy.utils import iers

    iers_b = iers.IERS_B.open()
    pm_x, pm_y = iers_b.pm_xy(time)
    dut1 = iers_b.ut1_utc(time)
    # TAI - UTC per epoch from the IAU leap-second table (steps exactly at
    # each insertion; e.g. 36 -> 37 at 2017-01-01). Not astropy's
    # ``unix_tai - unix``: astropy smears the inserted second across the whole
    # leap day, so that difference is off by up to 1 near an insertion.
    ymdhms = time.utc.ymdhms
    day_fraction = (
        ymdhms["hour"] / 24.0 + ymdhms["minute"] / 1440.0 + ymdhms["second"] / 86400.0
    )
    leap_seconds = np.ascontiguousarray(
        erfa.dat(ymdhms["year"], ymdhms["month"], ymdhms["day"], day_fraction)
    )
    return {
        "polar_motion_x_arcsec": np.ascontiguousarray(pm_x.to_value("arcsec")),
        "polar_motion_y_arcsec": np.ascontiguousarray(pm_y.to_value("arcsec")),
        "ut1_minus_utc_seconds": np.ascontiguousarray(dut1.to_value("s")),
        "leap_seconds": leap_seconds,
    }

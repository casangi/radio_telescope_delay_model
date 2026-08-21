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
        ``polar_motion_x_arcsec``, ``polar_motion_y_arcsec`` and
        ``ut1_minus_utc_seconds`` as ``[n_time]`` float64 arrays (IERS-B, the
        table astropy ships), and the scalar ``leap_seconds`` (TAI - UTC).
    """
    from astropy.utils import iers

    iers_b = iers.IERS_B.open()
    pm_x, pm_y = iers_b.pm_xy(time)
    dut1 = iers_b.ut1_utc(time)
    # TAI - UTC in integer seconds at the first time (constant per run; CALC
    # takes a single value).
    leap_seconds = float(np.round(time[0].unix_tai - time[0].unix))
    return {
        "polar_motion_x_arcsec": np.ascontiguousarray(pm_x.to_value("arcsec")),
        "polar_motion_y_arcsec": np.ascontiguousarray(pm_y.to_value("arcsec")),
        "ut1_minus_utc_seconds": np.ascontiguousarray(dut1.to_value("s")),
        "leap_seconds": leap_seconds,
    }

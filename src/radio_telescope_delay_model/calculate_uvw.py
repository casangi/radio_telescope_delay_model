"""Baseline ``(u, v, w)`` coordinates by two methods: CALC11 and astropy.

Both methods return uvw in the **archival / VLBI convention** adopted by MSv4
(and AstroVIPER): right-handed axes with ``w`` towards the source, and the
baseline ``(antenna1, antenna2)`` carrying the difference of the projected
antenna positions ``P(antenna1) - P(antenna2)`` -- the numbers that
observatory-written measurement sets, AIPS and VLBI correlators contain.

* :func:`calculate_uvw_calc` -- delays from the CALC11 geodetic VLBI model
  (the model VLBI correlators use), converted to uvw with the DiFX
  finite-difference recipe (``difxcalc.c::callCalc``): the per-antenna
  quantities are delay-like (``-P(r)``), so the archival baseline is
  ``antenna2 - antenna1`` of them.
* :func:`calculate_uvw_astropy` -- geometric projections built with astropy
  (ITRF -> GCRS, sky-offset frame of the phase centre); the same algorithm
  AstroVIPER inherited from SIRIUS. The archival baseline is
  ``antenna1 - antenna2`` of the per-antenna projections.

The two methods agree to the level of their geodetic modelling differences
(centimetres on kilometre baselines); CALC additionally provides dry and wet
tropospheric delays.
"""

from __future__ import annotations

import numpy as np

SPEED_OF_LIGHT = 299792458.0

__all__ = [
    "calculate_antenna_uvw_astropy",
    "calculate_delays_calc",
    "calculate_uvw_astropy",
    "calculate_uvw_calc",
    "baseline_antenna_pairs",
]


def baseline_antenna_pairs(n_antenna: int) -> tuple[np.ndarray, np.ndarray]:
    """Cross-correlation antenna index pairs (``antenna1 < antenna2``)."""
    antenna1, antenna2 = np.triu_indices(n_antenna, k=1)
    return antenna1.astype(np.int64), antenna2.astype(np.int64)


def calculate_antenna_uvw_astropy(
    antenna_position: np.ndarray,
    time,
    phase_center_ra_dec: np.ndarray,
    direction_frame: str = "icrs",
) -> np.ndarray:
    """Per-antenna ``(u, v, w)`` projections (metres) for each time (astropy).

    The antenna ITRF positions are transformed to GCRS and rotated into the
    sky-offset frame of the phase centre, so ``w`` points towards the phase
    centre, ``u`` east and ``v`` towards the pole.

    Parameters
    ----------
    antenna_position : np.ndarray, [n_antenna, 3], metres
        ITRF geocentric antenna positions.
    time : array-like, [n_time]
        UTC times: ISO strings (``YYYY-MM-DDTHH:MM:SS.SSS``), unix seconds, or
        an ``astropy.time.Time``.
    phase_center_ra_dec : np.ndarray, [n_time | 1, 2], radians
        Phase centre per time (or a single fixed one).
    direction_frame : str
        Astropy frame of ``phase_center_ra_dec`` (``"icrs"`` or ``"fk5"``).

    Returns
    -------
    np.ndarray, [n_time, n_antenna, 3]
    """
    import astropy.coordinates as coord
    import astropy.units as u

    astropy_time = _as_time(time)
    antenna_position = np.asarray(antenna_position, dtype=np.float64)
    n_time = len(astropy_time)
    phase_center = np.broadcast_to(
        np.asarray(phase_center_ra_dec, dtype=np.float64).reshape(-1, 2), (n_time, 2)
    )

    location = coord.EarthLocation(
        x=antenna_position[:, 0] * u.m,
        y=antenna_position[:, 1] * u.m,
        z=antenna_position[:, 2] * u.m,
    )
    site = coord.EarthLocation(*antenna_position.mean(axis=0) * u.m)

    uvw = np.zeros((n_time, len(antenna_position), 3))
    for k in range(n_time):
        t = astropy_time[k]
        site_position, site_velocity = site.get_gcrs_posvel(t)
        antenna_gcrs = coord.GCRS(
            location.get_gcrs_posvel(t)[0],
            obstime=t,
            obsgeoloc=site_position,
            obsgeovel=site_velocity,
        )
        pointing = coord.SkyCoord(
            phase_center[k, 0] * u.rad, phase_center[k, 1] * u.rad,
            frame=direction_frame,
        )  # fmt: skip
        frame_uvw = pointing.transform_to(antenna_gcrs).skyoffset_frame()
        offset = antenna_gcrs.transform_to(frame_uvw).cartesian
        # sky-offset axes: x towards the source, y east, z towards the pole.
        uvw[k] = np.stack(
            [offset.y.to_value(u.m), offset.z.to_value(u.m), offset.x.to_value(u.m)],
            axis=1,
        )
    return uvw


def calculate_uvw_astropy(
    antenna_position: np.ndarray,
    time,
    phase_center_ra_dec: np.ndarray,
    direction_frame: str = "icrs",
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Baseline ``uvw`` (metres) for every time, astropy method.

    Returns
    -------
    uvw : np.ndarray, [n_time, n_baseline, 3]
        Archival convention: ``P(antenna1) - P(antenna2)``.
    antenna1, antenna2 : np.ndarray, [n_baseline] int
    """
    antenna1, antenna2 = baseline_antenna_pairs(np.shape(antenna_position)[0])
    antenna_uvw = calculate_antenna_uvw_astropy(
        antenna_position, time, phase_center_ra_dec, direction_frame
    )
    uvw = antenna_uvw[:, antenna1, :] - antenna_uvw[:, antenna2, :]
    return np.ascontiguousarray(uvw), antenna1, antenna2


def calculate_delays_calc(
    antenna_position: np.ndarray,
    time,
    phase_center_ra_dec: np.ndarray,
    reference_position: np.ndarray | None = None,
    temperature_celsius=None,
    pressure_hpa=None,
    humidity_fraction=None,
    axis_offset_metres=None,
    ephemeris_path: str | None = None,
) -> dict:
    """CALC11 delays of every antenna relative to the array reference.

    Parameters
    ----------
    antenna_position : np.ndarray, [n_antenna, 3], metres
        ITRF geocentric antenna positions.
    time : array-like, [n_time]
        UTC times (>= 2; CALC forms EOP rates from the first and last).
    phase_center_ra_dec : np.ndarray, [n_time | 1, 2], radians (ICRS)
    reference_position : np.ndarray, [3], metres, optional
        Array reference; default: mean antenna position.
    temperature_celsius, pressure_hpa, humidity_fraction : [n_antenna], optional
        Surface weather per antenna (defaults 0 C, 1013.25 hPa, 0.5); affects
        only the dry / wet delays, not the geometric delay.
    axis_offset_metres : [n_antenna], optional
        Antenna axis offsets (default 0; alt-az mounts assumed).
    ephemeris_path : str, optional
        JPL DE421 ephemeris file; defaults to the packaged copy.

    Returns
    -------
    dict
        ``geometric_delay``, ``dry_delay``, ``wet_delay`` -- each
        ``[n_time, n_antenna]`` seconds -- plus the earth-orientation
        parameters used (``polar_motion_x_arcsec``, ``polar_motion_y_arcsec``,
        ``ut1_minus_utc_seconds``, ``leap_seconds``).
    """
    from radio_telescope_delay_model.delay_model_cpp import almacalc
    from radio_telescope_delay_model.earth_orientation import (
        earth_orientation_parameters,
    )

    astropy_time = _as_time(time)
    antenna_position = np.asarray(antenna_position, dtype=np.float64)
    n_antenna = antenna_position.shape[0]
    n_time = len(astropy_time)
    if n_time < 2:
        raise ValueError("CALC needs at least two times (EOP rate estimation).")
    phase_center = np.ascontiguousarray(
        np.broadcast_to(
            np.asarray(phase_center_ra_dec, dtype=np.float64).reshape(-1, 2),
            (n_time, 2),
        )
    )
    if reference_position is None:
        reference_position = antenna_position.mean(axis=0)
    if ephemeris_path is None:
        ephemeris_path = _packaged_ephemeris_path()

    def per_antenna(value, default):
        if value is None:
            value = default
        return np.ascontiguousarray(
            np.broadcast_to(np.asarray(value, dtype=np.float64), (n_antenna,))
        )

    temperature = per_antenna(temperature_celsius, 0.0)
    pressure = per_antenna(pressure_hpa, 1013.25)
    humidity = per_antenna(humidity_fraction, 0.5)
    axis_offset = per_antenna(axis_offset_metres, 0.0)

    eop = earth_orientation_parameters(astropy_time)

    def run(ra, dec):
        geometric = np.empty((n_time, n_antenna))
        dry = np.empty((n_time, n_antenna))
        wet = np.empty((n_time, n_antenna))
        almacalc(
            np.ascontiguousarray(reference_position, dtype=np.float64),
            antenna_position,
            temperature,
            pressure,
            humidity,
            np.ascontiguousarray(astropy_time.mjd, dtype=np.float64),
            np.ascontiguousarray(ra),
            np.ascontiguousarray(dec),
            eop["polar_motion_x_arcsec"],
            eop["polar_motion_y_arcsec"],
            eop["ut1_minus_utc_seconds"],
            eop["leap_seconds"],
            axis_offset,
            ephemeris_path,
            geometric,
            dry,
            wet,
        )
        return geometric, dry, wet

    geometric, dry, wet = run(phase_center[:, 0], phase_center[:, 1])
    return {
        "geometric_delay": geometric,
        "dry_delay": dry,
        "wet_delay": wet,
        "_run": run,
        "_phase_center": phase_center,
        **eop,
    }


def calculate_uvw_calc(
    antenna_position: np.ndarray,
    time,
    phase_center_ra_dec: np.ndarray,
    reference_position: np.ndarray | None = None,
    ephemeris_path: str | None = None,
    delta: float = 1e-5,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Baseline ``uvw`` (metres) for every time, CALC11 method.

    Follows the DiFX recipe (``difxcalc.c::callCalc``): the geometric delay is
    evaluated at the phase centre and at two offset directions, and the
    per-antenna uvw follows from the delay and its finite-difference
    direction derivatives::

        u = (c / delta) * (tau - tau_x)   with ra_x  = ra - delta / cos(dec)
        v = (c / delta) * (tau_y - tau)   with dec_y = dec + delta
        w = c * tau

    These per-antenna values are delay-like (``-P(r)``), so the archival
    baseline is assembled as ``antenna2 - antenna1``.

    Returns
    -------
    uvw : np.ndarray, [n_time, n_baseline, 3]
        Archival convention: ``P(antenna1) - P(antenna2)``.
    antenna1, antenna2 : np.ndarray, [n_baseline] int
    """
    results = calculate_delays_calc(
        antenna_position,
        time,
        phase_center_ra_dec,
        reference_position=reference_position,
        ephemeris_path=ephemeris_path,
    )
    run = results["_run"]
    phase_center = results["_phase_center"]
    tau = results["geometric_delay"]

    ra = phase_center[:, 0]
    dec = phase_center[:, 1]
    tau_x, _, _ = run(ra - delta / np.cos(dec), dec)
    tau_y, _, _ = run(ra, dec + delta)

    factor = SPEED_OF_LIGHT / delta
    antenna_uvw = np.stack(
        [factor * (tau - tau_x), factor * (tau_y - tau), SPEED_OF_LIGHT * tau],
        axis=2,
    )  # [n_time, n_antenna, 3], delay-like sign (-P(r))

    antenna1, antenna2 = baseline_antenna_pairs(np.shape(antenna_position)[0])
    uvw = antenna_uvw[:, antenna2, :] - antenna_uvw[:, antenna1, :]
    return np.ascontiguousarray(uvw), antenna1, antenna2


def _as_time(time):
    from astropy.time import Time

    if isinstance(time, Time):
        return np.atleast_1d(time)
    time = np.asarray(time)
    if np.issubdtype(time.dtype, np.number):
        return np.atleast_1d(Time(time.astype(np.float64), format="unix", scale="utc"))
    return np.atleast_1d(Time(time.astype(str), scale="utc"))


def _packaged_ephemeris_path() -> str:
    import os

    path = os.path.join(os.path.dirname(__file__), "data", "DE421_little_Endian")
    if not os.path.exists(path):
        raise FileNotFoundError(
            "Packaged JPL DE421 ephemeris not found; pass ephemeris_path."
        )
    return path

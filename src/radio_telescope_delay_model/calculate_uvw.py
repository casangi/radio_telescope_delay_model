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

# difxcalc11 mount strings (dinit.f::dSITI) -> CALC axis type codes.
MOUNT_TYPE_TO_AXIS_CODE = {
    "AZEL": 3,  # altitude-azimuth
    "EQUA": 1,  # equatorial
    "XYNS": 2,  # X/Y, X-axis north-south
    "XYEW": 4,  # X/Y, X-axis east-west
    "RICH": 5,  # Richmond special
}

__all__ = [
    "calculate_antenna_uvw_astropy",
    "calculate_delays_calc",
    "calculate_uvw_astropy",
    "calculate_uvw_calc",
    "baseline_antenna_pairs",
    "MOUNT_TYPE_TO_AXIS_CODE",
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
    station_name=None,
    mount_type=None,
) -> dict:
    """CALC11 delays of every antenna relative to the array reference.

    Parameters
    ----------
    antenna_position : np.ndarray, [n_antenna, 3], metres
        ITRF geocentric antenna positions.
    time : array-like, [n_time]
        UTC times (>= 2; CALC forms EOP rates from a run's first and last
        epoch).
    phase_center_ra_dec : np.ndarray, [n_time | 1, 2], radians (ICRS)
    reference_position : np.ndarray, [3], metres, optional
        Array reference; default: mean antenna position. Passing exactly
        ``(0, 0, 0)`` selects difxcalc11's geocenter mode: per-antenna
        delays are referenced to the geocenter and CALC's geocenter
        special-casing applies (dummy site geometry, no station-1
        atmosphere/tides) -- the convention of difxcalc11's ``.im``
        per-antenna model.
    temperature_celsius, pressure_hpa, humidity_fraction : [n_antenna], optional
        Measured surface weather per antenna; affects only the dry / wet
        delays, not the geometric delay. When omitted, CALC's internal
        height-based standard-atmosphere defaults apply (temperature
        20 C - 6.5e-3 K/m x height, pressure 1013.25 x (1 - 6.5e-3 x
        height / 293.15)^5.26 hPa, relative humidity 0.5) -- the same
        defaults difxcalc11 uses when a .calc file carries no weather.
    axis_offset_metres : [n_antenna], optional
        Antenna axis offsets (default 0).
    ephemeris_path : str, optional
        JPL DE421 ephemeris file; defaults to the packaged copy.
    station_name : sequence of str, [n_antenna], optional
        Telescope names looked up in the packaged difxcalc11 station
        catalogs (see :mod:`radio_telescope_delay_model.station_data`):
        ocean loading and ocean pole tide loading coefficients are applied
        per antenna; unknown names get zeros with a warning, exactly as
        difxcalc11 behaves.
    mount_type : str or sequence of str, [n_antenna], optional
        difxcalc11 mount strings, one of ``AZEL`` (default), ``EQUA``,
        ``XYNS``, ``XYEW``, ``RICH``; sets the CALC antenna axis type used
        by the axis-offset module.

    Returns
    -------
    dict
        ``geometric_delay``, ``dry_delay``, ``wet_delay`` -- each
        ``[n_time, n_antenna]`` seconds -- plus the earth-orientation
        parameters used (``polar_motion_x_arcsec``, ``polar_motion_y_arcsec``,
        ``ut1_minus_utc_seconds``, ``leap_seconds`` -- each ``[n_time]``).

    Notes
    -----
    TAI - UTC is evaluated per epoch; a request that straddles a leap-second
    insertion is split internally into runs of constant leap count (CALC's
    ATMUTC treats it as a run constant) and the results are stitched back in
    the caller's epoch order.
    """
    from radio_telescope_delay_model.delay_model_cpp import calc_delay_model
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
        # For continent- or Earth-spanning arrays (VLBA, EVN, EHT, ...) the
        # mean position lies far below the surface and is not a meaningful
        # site; the consensus formula's baseline nonlinearity then makes the
        # reference choice matter (~cm on 1e4 km baselines). VLBI convention
        # is the geocenter.
        if np.linalg.norm(reference_position) < 6.3e6:
            import warnings

            warnings.warn(
                "The default array reference (mean antenna position) is deep "
                "inside the Earth for this array; for VLBI-scale arrays pass "
                "reference_position=np.zeros(3) to use the geocenter (the "
                "VLBI/difxcalc11 convention), or a specific site.",
                stacklevel=2,
            )
    if ephemeris_path is None:
        ephemeris_path = _packaged_ephemeris_path()

    def per_antenna(value, default):
        if value is None:
            value = default
        return np.ascontiguousarray(
            np.broadcast_to(np.asarray(value, dtype=np.float64), (n_antenna,))
        )

    # -999 sentinels select CALC's height-based standard-atmosphere defaults
    # (catmm.f treats temperature/humidity <= -90 and pressure <= 0 as
    # "no measurement"), matching difxcalc11 when no weather is supplied.
    temperature = per_antenna(temperature_celsius, -999.0)
    pressure = per_antenna(pressure_hpa, -999.0)
    humidity = per_antenna(humidity_fraction, -999.0)
    axis_offset = per_antenna(axis_offset_metres, 0.0)

    station_kwargs = {}
    if mount_type is not None:
        if isinstance(mount_type, str):
            mount_type = [mount_type] * n_antenna
        mounts = [str(m).strip().upper() for m in mount_type]
        if len(mounts) != n_antenna:
            raise ValueError("mount_type must have one entry per antenna.")
        unknown = sorted({m for m in mounts if m not in MOUNT_TYPE_TO_AXIS_CODE})
        if unknown:
            raise ValueError(
                f"Unknown mount_type {unknown}; supported difxcalc11 mounts: "
                f"{sorted(MOUNT_TYPE_TO_AXIS_CODE)}."
            )
        station_kwargs["axis_type"] = np.ascontiguousarray(
            [MOUNT_TYPE_TO_AXIS_CODE[m] for m in mounts], dtype=np.int32
        )
    if station_name is not None:
        from radio_telescope_delay_model.station_data import station_geophysics

        if len(station_name) != n_antenna:
            raise ValueError("station_name must have one entry per antenna.")
        station_kwargs.update(station_geophysics(station_name))

    eop = earth_orientation_parameters(astropy_time)
    mjd = np.ascontiguousarray(astropy_time.mjd, dtype=np.float64)

    # CALC treats the leap count (ATMUTC) as a constant of one run, so a
    # request that straddles a leap-second insertion is split into runs of
    # constant TAI - UTC and stitched back in epoch order. A single-epoch run
    # is padded with the nearest other epoch -- discarded from the output --
    # only so the per-run EOP rates (formed from a run's first and last
    # epoch) stay defined; the pad's UT1 - UTC is re-expressed with the run's
    # leap count, keeping the rate free of the 1 s UT1 - UTC step at the
    # insertion.
    leap_seconds = eop["leap_seconds"]
    constant_leap_runs = []
    for leap_value in np.unique(leap_seconds):
        selected = np.flatnonzero(leap_seconds == leap_value)
        call_indices = selected
        if selected.size == 1:
            other = np.flatnonzero(leap_seconds != leap_value)
            pad = other[np.argmin(np.abs(mjd[other] - mjd[selected[0]]))]
            call_indices = np.sort(np.append(selected, pad))
        dut1_call = np.ascontiguousarray(
            eop["ut1_minus_utc_seconds"][call_indices]
            + (leap_value - leap_seconds[call_indices])
        )
        keep = np.isin(call_indices, selected)
        constant_leap_runs.append((call_indices, keep, float(leap_value), dut1_call))

    def run(ra, dec):
        ra = np.asarray(ra, dtype=np.float64)
        dec = np.asarray(dec, dtype=np.float64)
        geometric = np.empty((n_time, n_antenna))
        dry = np.empty((n_time, n_antenna))
        wet = np.empty((n_time, n_antenna))
        for call_indices, keep, leap_value, dut1_call in constant_leap_runs:
            call_geometric = np.empty((call_indices.size, n_antenna))
            call_dry = np.empty((call_indices.size, n_antenna))
            call_wet = np.empty((call_indices.size, n_antenna))
            calc_delay_model(
                np.ascontiguousarray(reference_position, dtype=np.float64),
                antenna_position,
                temperature,
                pressure,
                humidity,
                np.ascontiguousarray(mjd[call_indices]),
                np.ascontiguousarray(ra[call_indices]),
                np.ascontiguousarray(dec[call_indices]),
                np.ascontiguousarray(eop["polar_motion_x_arcsec"][call_indices]),
                np.ascontiguousarray(eop["polar_motion_y_arcsec"][call_indices]),
                dut1_call,
                leap_value,
                axis_offset,
                ephemeris_path,
                call_geometric,
                call_dry,
                call_wet,
                **station_kwargs,
            )
            geometric[call_indices[keep]] = call_geometric[keep]
            dry[call_indices[keep]] = call_dry[keep]
            wet[call_indices[keep]] = call_wet[keep]
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
    mode: str = "geometric",
    temperature_celsius=None,
    pressure_hpa=None,
    humidity_fraction=None,
    axis_offset_metres=None,
    station_name=None,
    mount_type=None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Baseline ``uvw`` (metres) for every time, CALC11 method.

    Two recipes are provided, selected with ``mode``:

    ``"geometric"`` (default)
        The DiFX calcserver recipe (``difxcalc.c::callCalc``): the geometric
        (vacuum) delay is evaluated at the phase centre and at two offset
        directions, and the per-antenna uvw follows from one-sided
        finite differences::

            u = (c / delta) * (tau - tau_x)   with ra_x  = ra - delta / cos(dec)
            v = (c / delta) * (tau_y - tau)   with dec_y = dec + delta
            w = c * tau

    ``"difxcalc11"``
        Runs the **embedded difxcalc11 pipeline itself** (the vendored
        dSTART/dINITL/dSCAN/dDRIVR code driven through a generated
        correlator ``.calc`` job; see
        :mod:`radio_telescope_delay_model.difxcalc11_core`). The raw
        per-antenna samples are bit-identical to what the difxcalc11
        binary computes internally ("ABERRATION CORR: EXACT", geocentric
        per-antenna model, atmosphere in the total delay and uvw,
        difxcalc11's own EOP tables and height-based weather defaults);
        arbitrary epochs are evaluated by interpolating each 2-minute
        interval's unique degree-5 sample polynomial -- the polynomial its
        ``.im`` files carry. ``reference_position``, weather and axis
        offsets must be omitted in this mode (a ``.calc`` job carries
        none). The geocentric reference matters beyond bookkeeping: the
        consensus formula is nonlinear in the baseline, so geocentric
        per-antenna differencing differs from array-referenced
        differencing by ``(K.b_geo/c)(K.(w2-w1)/c)`` terms, ~1e-4 m for a
        200 m baseline -- and DiFX products carry the geocentric
        convention.

    Either way the per-antenna values are delay-like (``-P(r)``), so the
    archival baseline is assembled as ``antenna2 - antenna1``.

    ``station_name`` / ``mount_type`` (and the weather and axis-offset
    arguments) are passed through to :func:`calculate_delays_calc`.

    Returns
    -------
    uvw : np.ndarray, [n_time, n_baseline, 3]
        Archival convention: ``P(antenna1) - P(antenna2)``.
    antenna1, antenna2 : np.ndarray, [n_baseline] int
    """
    if mode not in ("geometric", "difxcalc11"):
        raise ValueError(
            f"Unknown mode {mode!r}; expected 'geometric' or 'difxcalc11'."
        )
    if mode == "difxcalc11":
        if reference_position is not None and np.any(
            np.asarray(reference_position, dtype=np.float64) != 0.0
        ):
            raise ValueError(
                "mode='difxcalc11' uses the geocenter reference (the .im "
                "per-antenna convention); do not pass reference_position."
            )
        if any(
            value is not None
            for value in (
                temperature_celsius,
                pressure_hpa,
                humidity_fraction,
                axis_offset_metres,
            )
        ):
            raise ValueError(
                "mode='difxcalc11' runs the embedded difxcalc11 pipeline, "
                "which (like a correlator .calc job) carries no measured "
                "weather or axis offsets; its height-based defaults apply."
            )
        from radio_telescope_delay_model.difxcalc11_core import (
            difxcalc11_samples,
            evaluate_samples,
        )

        astropy_time = _as_time(time)
        phase_center = np.asarray(phase_center_ra_dec, dtype=np.float64)
        samples = difxcalc11_samples(
            antenna_position,
            astropy_time,
            phase_center,
            station_name=station_name,
            mount_type=mount_type,
        )
        epochs = np.atleast_1d(astropy_time.unix)
        antenna_uvw = np.stack(
            [evaluate_samples(samples, q, epochs) for q in ("u", "v", "w")],
            axis=2,
        )  # [n_time, n_antenna, 3], geocentric delay-like (-P(r))
        antenna1, antenna2 = baseline_antenna_pairs(np.shape(antenna_position)[0])
        uvw = antenna_uvw[:, antenna2, :] - antenna_uvw[:, antenna1, :]
        return np.ascontiguousarray(uvw), antenna1, antenna2

    results = calculate_delays_calc(
        antenna_position,
        time,
        phase_center_ra_dec,
        reference_position=reference_position,
        temperature_celsius=temperature_celsius,
        pressure_hpa=pressure_hpa,
        humidity_fraction=humidity_fraction,
        axis_offset_metres=axis_offset_metres,
        ephemeris_path=ephemeris_path,
        station_name=station_name,
        mount_type=mount_type,
    )
    run = results["_run"]
    phase_center = results["_phase_center"]
    ra = phase_center[:, 0]
    dec = phase_center[:, 1]

    tau = results["geometric_delay"]
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
    # The packaged DE421 is exactly 1715 direct-access records of 8144 bytes.
    # Any other size means the file was corrupted at rest (Dropbox sync has
    # truncated it twice), which would otherwise surface only as silent NaN
    # delays from misaligned ephemeris records.
    size = os.path.getsize(path)
    if size != 13966960:
        raise OSError(
            f"Packaged DE421 ephemeris is corrupted ({size} bytes, expected "
            "13966960). Restore it, e.g. 'git checkout -- "
            "src/radio_telescope_delay_model/data/DE421_little_Endian', and "
            "check what rewrote it (Dropbox sync has done this)."
        )
    return path

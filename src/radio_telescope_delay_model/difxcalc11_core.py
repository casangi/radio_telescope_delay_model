"""The embedded difxcalc11 pipeline: bit parity with the difxcalc binary.

``difxcalc11_samples`` writes a difxcalc ``.calc`` job, runs the vendored
difxcalc11 core (its own dSTART -> dINITL -> dSCAN -> dDRIVR sequence,
compiled unmodified except for five runtime data-path OPEN sites) and
returns the raw per-antenna model samples -- **bit-identical** to what the
difxcalc11 binary computes internally before its ``.im`` polynomial fit,
verified against an instrumented reference build compiled with the same
toolchain and reproducibility flags (see experiments/difx_uvw_comparison).

The sample grid is difxcalc11's: 6 epochs per 2-minute interval, 24 s
apart, interval boundaries repeated. ``evaluate_samples`` interpolates the
unique degree-5 polynomial through each interval's 6 samples (the same
polynomial the ``.im`` file carries, without the least-squares rounding of
its coefficient representation) at arbitrary epochs.

EOP handling is difxcalc11's own: the ``.calc`` carries 5 daily IERS-B
entries (TAI-UTC from the IAU table) and the core builds and interpolates
its internal tables from them. Surface weather is difxcalc11's height-based
default (a ``.calc`` carries no weather).
"""

from __future__ import annotations

import os
import tempfile

import numpy as np

__all__ = ["difxcalc11_samples", "evaluate_samples"]

_SAMPLES_PER_INTERVAL = 6
_SAMPLE_SPACING_SECONDS = 24.0


def _data_file(name: str) -> str:
    path = os.path.join(os.path.dirname(__file__), "data", name)
    if not os.path.exists(path):
        raise FileNotFoundError(f"Packaged data file missing: {path}")
    return path


def _calc_lines(
    antenna_position,
    station_name,
    mount_type,
    start,
    duration_s,
    ra,
    dec,
    eop_mjds,
    eop_tai_utc,
    eop_ut1_utc,
    eop_x,
    eop_y,
):
    lines = [
        ("JOB ID:", "1"),
        ("JOB START TIME:", f"{start.mjd:.7f}"),
        ("JOB STOP TIME:", f"{start.mjd + duration_s / 86400.0:.7f}"),
        ("DUTY CYCLE:", "1.000"),
        ("OBSCODE:", "RTDM"),
        ("DIFX VERSION:", "rtdm-difxcalc11"),
        ("SUBJOB ID:", "0"),
        ("SUBARRAY ID:", "0"),
        ("VEX FILE:", "none.vex"),
        ("START MJD:", f"{start.mjd:.7f}"),
        ("START YEAR:", str(start.datetime.year)),
        ("START MONTH:", str(start.datetime.month)),
        ("START DAY:", str(start.datetime.day)),
        ("START HOUR:", str(start.datetime.hour)),
        ("START MINUTE:", str(start.datetime.minute)),
        ("START SECOND:", str(start.datetime.second)),
        ("SPECTRAL AVG:", "1"),
        ("TAPER FUNCTION:", "UNIFORM"),
        ("NUM TELESCOPES:", str(len(antenna_position))),
    ]
    for k, (x, y, z) in enumerate(antenna_position):
        lines += [
            (f"TELESCOPE {k} NAME:", station_name[k]),
            (f"TELESCOPE {k} MOUNT:", mount_type[k]),
            (f"TELESCOPE {k} OFFSET (m):", "0.000000"),
            (f"TELESCOPE {k} X (m):", f"{x:.6f}"),
            (f"TELESCOPE {k} Y (m):", f"{y:.6f}"),
            (f"TELESCOPE {k} Z (m):", f"{z:.6f}"),
            (f"TELESCOPE {k} SHELF:", "NONE"),
        ]
    lines += [
        ("NUM SOURCES:", "1"),
        ("SOURCE 0 NAME:", "TARGET"),
        ("SOURCE 0 RA:", f"{ra:.16f}"),
        ("SOURCE 0 DEC:", f"{dec:.16f}"),
        ("SOURCE 0 CALCODE:", ""),
        ("SOURCE 0 QUAL:", "0"),
        ("NUM SCANS:", "1"),
        ("SCAN 0 IDENTIFIER:", "No0001"),
        ("SCAN 0 START (S):", "0"),
        ("SCAN 0 DUR (S):", str(int(duration_s))),
        ("SCAN 0 OBS MODE NAME:", "rtdm"),
        ("SCAN 0 UVSHIFT INTERVAL (NS):", "2000000000"),
        ("SCAN 0 AC AVG INTERVAL (NS):", "2000000"),
        ("SCAN 0 POINTING SRC:", "0"),
        ("SCAN 0 NUM PHS CTRS:", "1"),
        ("SCAN 0 PHS CTR 0:", "0"),
        ("NUM EOPS:", str(len(eop_mjds))),
    ]
    for k in range(len(eop_mjds)):
        lines += [
            (f"EOP {k} TIME (mjd):", str(int(eop_mjds[k]))),
            (f"EOP {k} TAI_UTC (sec):", str(int(eop_tai_utc[k]))),
            (f"EOP {k} UT1_UTC (sec):", f"{eop_ut1_utc[k]:.13f}"),
            (f"EOP {k} XPOLE (arcsec):", f"{eop_x[k]:.13f}"),
            (f"EOP {k} YPOLE (arcsec):", f"{eop_y[k]:.13f}"),
        ]
    lines += [
        ("NUM SPACECRAFT:", "0"),
        ("IM FILENAME:", "unused.im"),
        ("FLAG FILENAME:", "unused.flag"),
    ]
    return lines


def difxcalc11_samples(
    antenna_position,
    time,
    phase_center_ra_dec,
    station_name=None,
    mount_type=None,
    eop=None,
) -> dict:
    """Raw difxcalc11 model samples covering the given times.

    Parameters mirror :func:`calculate_delays_calc`; ``station_name`` feeds
    difxcalc11's own catalog lookup (unknown names get its usual warning on
    stdout and zero coefficients), ``mount_type`` its mount strings.

    ``eop`` optionally supplies the daily EOP entries for the generated
    ``.calc`` -- a dict with ``mjd``, ``tai_utc``, ``ut1_utc``,
    ``x_pole_arcsec``, ``y_pole_arcsec`` arrays (e.g. the FITS-IDI CALC
    table: the values the correlator actually used, reproducing its frame
    exactly). Default: 5 daily IERS-B entries with erfa TAI-UTC.

    Returns a dict with ``sample_unix`` [n_samples] (UTC epochs, 24 s grid,
    interval boundaries deduplicated), ``interval_start_unix``
    [n_intervals], and per-antenna arrays [n_samples, n_antenna]:
    ``delay`` (total, CALC sign, seconds; the ``.im`` DELAY is -1e6 x this),
    ``dry``, ``wet`` (seconds), ``u``, ``v``, ``w`` (geocentric per-antenna
    model, metres).
    """
    import erfa
    from astropy.time import Time
    from astropy.utils import iers

    from radio_telescope_delay_model.delay_model_cpp import _difxcalc11_ext

    antenna_position = np.asarray(antenna_position, dtype=np.float64)
    n_antenna = antenna_position.shape[0]
    if not isinstance(time, Time):
        raise TypeError("time must be an astropy Time (use _as_time upstream).")
    if station_name is None:
        station_name = [f"A{k:02d}" for k in range(n_antenna)]
    if len(station_name) != n_antenna:
        raise ValueError("station_name must have one entry per antenna.")
    station_name = [str(s).strip().upper()[:8] for s in station_name]
    if mount_type is None:
        mount_type = ["AZEL"] * n_antenna
    elif isinstance(mount_type, str):
        mount_type = [mount_type] * n_antenna
    mount_type = [str(m).strip().upper() for m in mount_type]

    phase_center = np.asarray(phase_center_ra_dec, dtype=np.float64).reshape(-1, 2)
    if phase_center.shape[0] != 1:
        raise ValueError(
            "the difxcalc11 backend supports a single phase centre per run."
        )

    # Job window: start on a whole UTC minute at least a minute before the
    # first epoch; difxcalc processes ProcMin/2 + 1 two-minute intervals, so
    # the samples always cover past the last epoch.
    t_unix = np.atleast_1d(time.unix)
    start = Time(np.floor(t_unix.min() / 60.0 - 1.0) * 60.0, format="unix", scale="utc")
    start = Time(start.isot.split(".")[0], scale="utc")  # exact whole second
    duration = int(np.ceil((t_unix.max() - start.unix) / 120.0) + 1) * 120

    if eop is not None:
        eop_mjds = np.asarray(eop["mjd"], dtype=int)
        tai_utc = np.asarray(eop["tai_utc"], dtype=np.float64)
        dut1 = np.asarray(eop["ut1_utc"], dtype=np.float64)
        pole_x = np.asarray(eop["x_pole_arcsec"], dtype=np.float64)
        pole_y = np.asarray(eop["y_pole_arcsec"], dtype=np.float64)
    else:
        # 5 daily IERS-B EOP entries around the window, exactly as a
        # correlator .calc carries; TAI-UTC per day from the IAU table.
        mjd0 = int(np.floor(start.mjd))
        eop_mjds = np.arange(mjd0 - 2, mjd0 + 3)
        eop_times = Time(eop_mjds.astype(float), format="mjd", scale="utc")
        iers_b = iers.IERS_B.open()
        pm_x, pm_y = iers_b.pm_xy(eop_times)
        pole_x = pm_x.to_value("arcsec")
        pole_y = pm_y.to_value("arcsec")
        dut1 = iers_b.ut1_utc(eop_times).to_value("s")
        ymdf = eop_times.ymdhms
        tai_utc = erfa.dat(ymdf["year"], ymdf["month"], ymdf["day"], 0.0)

    lines = _calc_lines(
        antenna_position,
        station_name,
        mount_type,
        start,
        duration,
        phase_center[0, 0],
        phase_center[0, 1],
        eop_mjds,
        tai_utc,
        dut1,
        pole_x,
        pole_y,
    )
    handle, calc_path = tempfile.mkstemp(suffix=".calc", prefix="rtdm_")
    if len(calc_path) > 127:  # CALC filename buffers are CHARACTER*128
        os.close(handle)
        os.remove(calc_path)
        handle, calc_path = tempfile.mkstemp(suffix=".calc", prefix="rtdm_", dir="/tmp")
    try:
        with os.fdopen(handle, "w") as f:
            for key, value in lines:
                f.write(f"{key:<20}{value}\n")

        max_samples = (duration // 120 + 2) * _SAMPLES_PER_INTERVAL
        ymdhms = np.zeros((max_samples, 6), dtype=np.int32)
        arrays = [np.zeros((max_samples, n_antenna)) for _ in range(6)]
        n = _difxcalc11_ext.difxcalc11_samples(
            calc_path,
            _data_file("ocean_load.coef"),
            _data_file("tilt.dat"),
            _data_file("ocean_pole_tide.coef"),
            _data_file("DE421_little_Endian"),
            _data_file("ut1ls.dat"),
            ymdhms,
            *arrays,
        )
    finally:
        os.remove(calc_path)

    ymdhms = ymdhms[:n]
    delay, dry, wet, u, v, w = (a[:n] for a in arrays)
    sample_time = Time(
        {
            "year": ymdhms[:, 0],
            "month": ymdhms[:, 1],
            "day": ymdhms[:, 2],
            "hour": ymdhms[:, 3],
            "minute": ymdhms[:, 4],
            "second": ymdhms[:, 5].astype(np.float64),
        },
        scale="utc",
    )
    return {
        "sample_unix": sample_time.unix,
        "interval_start_unix": sample_time.unix[::_SAMPLES_PER_INTERVAL],
        "n_intervals": n // _SAMPLES_PER_INTERVAL,
        "delay": delay,
        "dry": dry,
        "wet": wet,
        "u": u,
        "v": v,
        "w": w,
    }


def evaluate_samples(samples: dict, quantity: str, epochs_unix) -> np.ndarray:
    """Evaluate a sampled quantity at arbitrary epochs.

    Interpolates, per 2-minute interval, the unique degree-5 polynomial
    through the interval's 6 samples (barycentric Lagrange on the 24 s
    nodes) -- mathematically the polynomial difxcalc11 fits into its ``.im``
    file for that interval.
    """
    values = samples[quantity]  # [n_samples, n_antenna]
    n_intervals = samples["n_intervals"]
    starts = samples["interval_start_unix"]
    epochs_unix = np.atleast_1d(np.asarray(epochs_unix, dtype=np.float64))

    nodes = np.arange(_SAMPLES_PER_INTERVAL) * _SAMPLE_SPACING_SECONDS
    # Barycentric weights for the 6 nodes.
    weights = np.empty(_SAMPLES_PER_INTERVAL)
    for j in range(_SAMPLES_PER_INTERVAL):
        d = 1.0
        for m in range(_SAMPLES_PER_INTERVAL):
            if m != j:
                d *= nodes[j] - nodes[m]
        weights[j] = 1.0 / d

    out = np.empty((len(epochs_unix),) + values.shape[1:])
    for i, t in enumerate(epochs_unix):
        k = int(np.clip((t - starts[0]) // 120.0, 0, n_intervals - 1))
        block = values[k * _SAMPLES_PER_INTERVAL : (k + 1) * _SAMPLES_PER_INTERVAL]
        dt = t - starts[k]
        diff = dt - nodes
        on_node = np.isclose(diff, 0.0, atol=1e-9)
        if on_node.any():
            out[i] = block[int(np.argmax(on_node))]
        else:
            factors = weights / diff
            out[i] = (factors[:, None] * block).sum(axis=0) / factors.sum()
    return out

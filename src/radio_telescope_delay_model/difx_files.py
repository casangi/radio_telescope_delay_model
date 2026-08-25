"""Readers for DiFX correlator model files (``.calc`` and ``.im``).

A ``.calc`` file is the job description difxcalc11 consumes (telescopes,
source, EOPs); the ``.im`` file is what it produces: per 2-minute interval
and per antenna, degree-5 polynomials in seconds-from-interval-start for
the delay (microseconds, DiFX sign), dry/wet atmosphere, and the
geocentric per-antenna U/V/W (metres). These readers exist so production
correlator models can be compared against this package at full precision
(see ``tests/unit/test_difx_im.py`` and the EHT example notebook).
"""

from __future__ import annotations

import re

import numpy as np

__all__ = ["read_calc", "read_im", "evaluate_im"]

_POLY_RE = re.compile(r"SCAN\s+(\d+)\s+POLY\s+(\d+)\s+(MJD|SEC):\s*(\S+)")
_ANT_RE = re.compile(
    r"SRC\s+(\d+)\s+ANT\s+(\d+)\s+(DELAY|DRY|WET|U|V|W|AZ|EL GEOM|EL CORR|"
    r"PAR ANGLE)\s+\([^)]*\):\s*(.*)"
)


def read_calc(path) -> dict:
    """Parse a DiFX ``.calc`` file.

    Returns ``station_name`` [n], ``mount_type`` [n], ``antenna_position``
    [n, 3] ITRF metres, ``axis_offset`` [n] metres, ``source_name``,
    ``phase_center_ra_dec`` [1, 2] radians, and ``eop`` in the dict layout
    :func:`radio_telescope_delay_model.difxcalc11_core.difxcalc11_samples`
    accepts.
    """
    values = {}
    for line in open(path):
        if ":" in line:
            key, _, value = line.partition(":")
            values[key.strip()] = value.strip()
    n_ant = int(values["NUM TELESCOPES"])
    n_eop = int(values["NUM EOPS"])
    return {
        "station_name": [values[f"TELESCOPE {i} NAME"] for i in range(n_ant)],
        "mount_type": [values[f"TELESCOPE {i} MOUNT"].upper() for i in range(n_ant)],
        "antenna_position": np.array(
            [
                [float(values[f"TELESCOPE {i} {axis} (m)"]) for axis in "XYZ"]
                for i in range(n_ant)
            ]
        ),
        "axis_offset": np.array(
            [float(values.get(f"TELESCOPE {i} OFFSET (m)", 0.0)) for i in range(n_ant)]
        ),
        "source_name": values["SOURCE 0 NAME"],
        "phase_center_ra_dec": np.array(
            [[float(values["SOURCE 0 RA"]), float(values["SOURCE 0 DEC"])]]
        ),
        "eop": {
            "mjd": np.array([int(values[f"EOP {i} TIME (mjd)"]) for i in range(n_eop)]),
            "tai_utc": np.array(
                [float(values[f"EOP {i} TAI_UTC (sec)"]) for i in range(n_eop)]
            ),
            "ut1_utc": np.array(
                [float(values[f"EOP {i} UT1_UTC (sec)"]) for i in range(n_eop)]
            ),
            "x_pole_arcsec": np.array(
                [float(values[f"EOP {i} XPOLE (arcsec)"]) for i in range(n_eop)]
            ),
            "y_pole_arcsec": np.array(
                [float(values[f"EOP {i} YPOLE (arcsec)"]) for i in range(n_eop)]
            ),
        },
    }


def read_im(path) -> list:
    """Parse a DiFX ``.im`` file into a list of polynomial intervals.

    Each entry is a dict with ``mjd`` (int), ``sec`` (int, interval start,
    UTC seconds of day), ``scan`` (int) and, per quantity name (``DELAY``,
    ``DRY``, ``WET`` in microseconds; ``U``, ``V``, ``W`` in metres; ...),
    a dict mapping antenna index -> coefficient array (degree ascending, in
    seconds from the interval start).
    """
    intervals = []
    current = None
    for line in open(path):
        m = _POLY_RE.match(line)
        if m:
            if m.group(3) == "MJD":
                current = {"scan": int(m.group(1)), "mjd": int(m.group(4)), "sec": None}
                intervals.append(current)
            else:
                current["sec"] = int(m.group(4))
            continue
        m = _ANT_RE.match(line)
        if m and current is not None:
            current.setdefault(m.group(3), {})[int(m.group(2))] = np.array(
                [float(x) for x in m.group(4).split()]
            )
    return intervals


def evaluate_im(intervals, quantity, antenna, mjd, seconds):
    """Evaluate ``quantity`` for one antenna at UTC epochs (``mjd`` int
    array, ``seconds`` of day), using each epoch's covering interval.
    Returns NaN where no interval covers an epoch."""
    mjd = np.atleast_1d(np.asarray(mjd, dtype=np.int64))
    seconds = np.atleast_1d(np.asarray(seconds, dtype=np.float64))
    out = np.full(seconds.shape, np.nan)
    for interval in intervals:
        coeffs = interval.get(quantity, {}).get(antenna)
        if coeffs is None:
            continue
        t = (mjd - interval["mjd"]) * 86400.0 + seconds - interval["sec"]
        inside = (t >= -1e-6) & (t < 120.0) & np.isnan(out)
        if np.any(inside):
            out[inside] = np.polyval(coeffs[::-1], t[inside])
    return out

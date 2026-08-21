"""Named-telescope geophysics from the packaged difxcalc11 station catalogs.

The package ships the same station catalogs difxcalc11 (and pycalc11, which
vendors difxcalc11) reads at run time:

* ``data/ocean_load.coef`` -- BLQ ocean loading (11 tides: vertical and
  west/south horizontal amplitudes in metres, phases in degrees), parsed as
  difxcalc11's ``dinit.f::dOCNIN`` does;
* ``data/ocean_pole_tide.coef`` -- Desai (2002) ocean pole tide loading
  coefficients, parsed as ``dinit.f::dOPTLIN`` does;
* ``data/tilt.dat`` -- antenna fixed-axis tilts in arc-minutes, parsed and
  epoch-interpolated as ``dinit.f::dANTILT`` does. (This CALC 11 version does
  not apply tilts in its axis-offset module -- neither does difxcalc11 -- so
  the catalog is informational.)

A telescope is "supported" when it appears in these catalogs; stations
missing from a catalog get zero coefficients with a warning, exactly like
difxcalc11.
"""

from __future__ import annotations

import functools
import os
import warnings

import numpy as np

__all__ = [
    "axis_tilt_catalog",
    "ocean_loading_catalog",
    "ocean_pole_tide_catalog",
    "station_geophysics",
    "telescope_names",
]

_DATA_DIR = os.path.join(os.path.dirname(__file__), "data")


def _normalize(name) -> str:
    return str(name).strip().upper()


@functools.lru_cache(maxsize=1)
def ocean_loading_catalog() -> dict:
    """Ocean loading coefficients by station name (and aliases).

    Returns a dict mapping every name that appears in one of the three BLQ
    name columns (``dOCNIN`` matches all three) to a shared entry dict with
    ``name`` (the primary name), ``vertical_amplitude`` [11] (m),
    ``horizontal_amplitude`` [2, 11] (west, south; m), ``vertical_phase``
    [11] and ``horizontal_phase`` [2, 11] (radians).
    """
    path = os.path.join(_DATA_DIR, "ocean_load.coef")
    with open(path) as f:
        lines = f.readlines()

    catalog: dict = {}
    i = 0
    while i < len(lines):
        line = lines[i].rstrip("\n")
        i += 1
        if line.startswith("$$") or not line.startswith("  "):
            continue
        tokens = line.split()
        if not tokens:
            continue
        try:
            float(tokens[0])
            continue  # a coefficient row outside a recognised block
        except ValueError:
            pass  # a station header line
        names = [line[2:10], line[12:20], line[22:30]]
        rows = []
        while i < len(lines) and len(rows) < 6:
            data_line = lines[i]
            i += 1
            if data_line.startswith("$$"):
                continue
            values = data_line.split()
            if len(values) != 11:
                break
            rows.append([float(v) for v in values])
        if len(rows) != 6:
            continue
        table = np.asarray(rows)
        entry = {
            "name": _normalize(names[0]),
            "vertical_amplitude": np.ascontiguousarray(table[0]),
            "horizontal_amplitude": np.ascontiguousarray(table[1:3]),
            "vertical_phase": np.ascontiguousarray(np.deg2rad(table[3])),
            "horizontal_phase": np.ascontiguousarray(np.deg2rad(table[4:6])),
        }
        for name in names:
            name = _normalize(name)
            if name:
                catalog.setdefault(name, entry)
    return catalog


@functools.lru_cache(maxsize=1)
def ocean_pole_tide_catalog() -> tuple:
    """Ocean pole tide loading entries (Desai 2002), in file order.

    Each entry has ``name``, ``code`` (short alias column, possibly empty),
    ``latitude_deg``, ``longitude_deg`` (east) and ``coefficients`` [6]
    (u_r, u_n, u_e real/imaginary pairs), matching ``dOPTLIN``.
    """
    path = os.path.join(_DATA_DIR, "ocean_pole_tide.coef")
    with open(path) as f:
        lines = f.readlines()
    start = None
    for k, line in enumerate(lines):
        if line.lstrip().startswith("Ocean Pole Tide Loading Coefficients"):
            start = k + 5  # dOPTLIN skips four lines after the title
            break
    if start is None:
        raise ValueError(f"No ocean pole tide header found in {path}.")
    entries = []
    for line in lines[start:]:
        name = _normalize(line[1:9])
        code = line[10:13]
        if not name:
            continue
        values = line[13:].split()
        if len(values) < 8:
            continue
        numbers = [float(v) for v in values[:8]]
        entries.append(
            {
                "name": name,
                "code": code,
                "latitude_deg": numbers[0],
                "longitude_deg": numbers[1],
                "coefficients": np.ascontiguousarray(numbers[2:8]),
            }
        )
    return tuple(entries)


def _ocean_pole_tide_for(name: str):
    """dOPTLIN matching: the full name, or the 3-character code column
    equal to the first three characters of the station name."""
    for entry in ocean_pole_tide_catalog():
        if entry["name"] == name:
            return entry
        if entry["code"].strip() and entry["code"].ljust(3) == name.ljust(3)[:3]:
            return entry
    return None


def axis_tilt_catalog(epoch_year: float) -> dict:
    """Antenna fixed-axis tilts (arc-minutes) at a decimal-year epoch.

    Returns ``{name: (east_tilt, north_tilt)}`` for both name columns of
    ``tilt.dat``, reproducing ``dANTILT``: single-epoch entries are
    constant; multi-epoch entries are linearly interpolated between the
    bracketing epochs (clamped before the first, extrapolated from the last
    two beyond the end).
    """
    path = os.path.join(_DATA_DIR, "tilt.dat")
    with open(path) as f:
        lines = f.readlines()[4:]  # dANTILT skips four header lines

    groups = []  # (names, [(epoch, east, north)])
    for line in lines:
        tokens = line.split()
        if not tokens:
            continue
        try:
            float(tokens[0])
            is_header = False
        except ValueError:
            is_header = True
        if is_header:
            if len(tokens) < 7:
                continue
            names = (_normalize(tokens[0]), _normalize(tokens[1]))
            values = tokens[2:]
            groups.append((names, []))
        else:
            if not groups or len(tokens) < 5:
                continue
            values = tokens
        epoch, east, north = float(values[0]), float(values[2]), float(values[4])
        groups[-1][1].append((epoch, east, north))

    catalog = {}
    for names, entries in groups:
        if not entries:
            continue
        a = entries[0]
        if epoch_year <= a[0]:
            b = a
        else:
            b = None
            for entry in entries[1:]:
                if epoch_year <= entry[0] or entry is entries[-1]:
                    b = entry
                    break
                a = entry
            if b is None:  # single-entry group
                b = a
        if abs(b[0] - a[0]) < 0.002:
            east, north = a[1], a[2]
        else:
            fraction = (epoch_year - a[0]) / (b[0] - a[0])
            east = a[1] + fraction * (b[1] - a[1])
            north = a[2] + fraction * (b[2] - a[2])
        for name in names:
            if name:
                catalog.setdefault(name, (east, north))
    return catalog


def telescope_names() -> list:
    """Primary names of every telescope in the ocean loading catalog."""
    return sorted({entry["name"] for entry in ocean_loading_catalog().values()})


def station_geophysics(station_name) -> dict:
    """Per-antenna station arrays for ``almacalc`` from the catalogs.

    Parameters
    ----------
    station_name : sequence of str, [n_antenna]

    Returns
    -------
    dict
        ``ocean_vertical_amplitude`` [n, 11], ``ocean_vertical_phase``
        [n, 11], ``ocean_horizontal_amplitude`` [n, 2, 11],
        ``ocean_horizontal_phase`` [n, 2, 11] and
        ``ocean_pole_tide_coefficients`` [n, 6]. Stations missing from a
        catalog contribute zeros, with a warning naming them (difxcalc11
        warns and continues identically).
    """
    names = [_normalize(name) for name in station_name]
    n = len(names)
    arrays = {
        "ocean_vertical_amplitude": np.zeros((n, 11)),
        "ocean_vertical_phase": np.zeros((n, 11)),
        "ocean_horizontal_amplitude": np.zeros((n, 2, 11)),
        "ocean_horizontal_phase": np.zeros((n, 2, 11)),
        "ocean_pole_tide_coefficients": np.zeros((n, 6)),
    }
    ocean = ocean_loading_catalog()
    missing_ocean = []
    missing_pole_tide = []
    for k, name in enumerate(names):
        entry = ocean.get(name)
        if entry is None:
            missing_ocean.append(name)
        else:
            arrays["ocean_vertical_amplitude"][k] = entry["vertical_amplitude"]
            arrays["ocean_vertical_phase"][k] = entry["vertical_phase"]
            arrays["ocean_horizontal_amplitude"][k] = entry["horizontal_amplitude"]
            arrays["ocean_horizontal_phase"][k] = entry["horizontal_phase"]
        pole_tide = _ocean_pole_tide_for(name)
        if pole_tide is None:
            missing_pole_tide.append(name)
        else:
            arrays["ocean_pole_tide_coefficients"][k] = pole_tide["coefficients"]
    if missing_ocean:
        warnings.warn(
            "No ocean loading coefficients for "
            f"{sorted(set(missing_ocean))}; using zeros (difxcalc11 warns "
            "and continues identically).",
            stacklevel=2,
        )
    if missing_pole_tide:
        warnings.warn(
            "No ocean pole tide loading coefficients for "
            f"{sorted(set(missing_pole_tide))}; using zeros (difxcalc11 "
            "warns and continues identically).",
            stacklevel=2,
        )
    return arrays

"""Extract uvw-comparison inputs from FITS-IDI files.

``read_fits_idi`` pulls everything a uvw computation needs out of a
FITS-IDI file -- antenna ITRF positions, station names, mounts, axis
offsets, sources, per-row epochs/baselines and the stored uvw -- and
``compare_uvw`` runs any of this package's methods over those rows and
returns the residuals (used by the example notebooks and the unit tests;
see ``example_notebooks/EHT_Fits_IDI_UVW_Comparison.ipynb``).

``read_fits_idi`` also reads the compact ``.npz`` extracts this package
uses for unit-test data (``tests/unit/data/eht_e21d15_3c273_uvw.npz``,
sampled from an EHT 2021 DiFX job), so tests and notebooks share one code
path without carrying multi-GB archives.
"""

from __future__ import annotations

import numpy as np

__all__ = ["read_fits_idi", "compare_uvw"]

SPEED_OF_LIGHT = 299792458.0

# FITS-IDI MNTSTA codes -> difxcalc11 mount strings. 4/5 are alt-az with a
# Nasmyth focus (right/left): identical to AZEL for delay/uvw purposes.
_MNTSTA_TO_MOUNT = {0: "AZEL", 1: "EQUA", 3: "XYEW", 4: "AZEL", 5: "AZEL"}


def read_fits_idi(path, max_rows: int | None = None) -> dict:
    """Read the uvw-relevant content of a FITS-IDI file (or ``.npz`` extract).

    Returns a dict with:

    * ``station_name`` [n_antenna] str, ``antenna_position`` [n_antenna, 3]
      ITRF metres (``ARRAYX/Y/Z`` offsets applied), ``mount_type``
      [n_antenna] difxcalc11 mount strings (from ``MNTSTA``),
      ``axis_offset`` [n_antenna] metres;
    * ``source_name`` [n_source] str, ``phase_center_ra_dec``
      [n_source, 2] ICRS/J2000 radians;
    * per cross-correlation row (autocorrelations dropped, optionally
      capped at ``max_rows``): ``time`` (astropy UTC ``Time``),
      ``antenna1``/``antenna2`` (0-based indices into the antenna axis, in
      the file's own order -- a row may have antenna1 > antenna2),
      ``source_index`` (0-based into the source axis) and ``uvw``
      [n_row, 3] metres (the file's light-second values times c).
    """
    from astropy.time import Time

    path = str(path)
    if path.endswith(".npz"):
        data = np.load(path)
        station_name = [str(s) for s in data["station_name"]]
        antenna_position = np.asarray(data["antenna_position"], dtype=np.float64)
        axis_offset = np.asarray(data["axis_offset"], dtype=np.float64)
        mntsta = np.asarray(data["mntsta"], dtype=int)
        source_name = [str(s) for s in data["source_name"]]
        ra_deg = np.asarray(data["source_ra_deg"], dtype=np.float64)
        dec_deg = np.asarray(data["source_dec_deg"], dtype=np.float64)
        date = np.asarray(data["DATE"], dtype=np.float64)
        day_fraction = np.asarray(data["TIME"], dtype=np.float64)
        baseline = np.asarray(data["BASELINE"], dtype=int)
        source_id = np.asarray(data["SOURCE_ID"], dtype=int)
        uvw_seconds = np.stack([data["UU"], data["VV"], data["WW"]], axis=1).astype(
            np.float64
        )
    else:
        from astropy.io import fits

        with fits.open(path, memmap=True) as hdus:
            geometry = hdus["ARRAY_GEOMETRY"]
            station_name = [str(n) for n in geometry.data["ANNAME"]]
            antenna_position = np.asarray(geometry.data["STABXYZ"], dtype=np.float64)
            center = np.array(
                [
                    geometry.header.get(key, 0.0) or 0.0
                    for key in ("ARRAYX", "ARRAYY", "ARRAYZ")
                ]
            )
            antenna_position = antenna_position + center
            axis_offset = np.atleast_2d(
                np.asarray(geometry.data["STAXOF"], dtype=np.float64)
            )[:, 0]
            mntsta = np.asarray(geometry.data["MNTSTA"], dtype=int)

            source = hdus["SOURCE"].data
            source_name = [str(s) for s in source["SOURCE"]]
            ra_deg = np.asarray(source["RAEPO"], dtype=np.float64)
            dec_deg = np.asarray(source["DECEPO"], dtype=np.float64)

            uv = hdus["UV_DATA"]
            columns = uv.columns.names

            def uvw_column(prefix):
                matches = [c for c in columns if c.upper().startswith(prefix)]
                if not matches:
                    raise KeyError(f"No {prefix}* column in UV_DATA.")
                return matches[0]

            rows = uv.header["NAXIS2"]
            index = np.arange(rows)
            if max_rows is not None and rows > max_rows:
                index = index[:: max(1, rows // max_rows)][:max_rows]
            table = uv.data
            date = np.asarray(table["DATE"][index], dtype=np.float64)
            day_fraction = np.asarray(table["TIME"][index], dtype=np.float64)
            baseline = np.asarray(table["BASELINE"][index], dtype=int)
            source_key = "SOURCE" if "SOURCE" in columns else "SOURCE_ID"
            source_id = np.asarray(table[source_key][index], dtype=int)
            uvw_seconds = np.stack(
                [
                    np.asarray(table[uvw_column(p)][index], dtype=np.float64)
                    for p in ("UU", "VV", "WW")
                ],
                axis=1,
            )

    unknown = sorted({int(m) for m in mntsta if int(m) not in _MNTSTA_TO_MOUNT})
    if unknown:
        raise ValueError(f"Unsupported FITS-IDI MNTSTA codes: {unknown}.")
    antenna1 = baseline // 256 - 1
    antenna2 = baseline % 256 - 1
    cross = antenna1 != antenna2
    return {
        "station_name": station_name,
        "antenna_position": antenna_position,
        "mount_type": [_MNTSTA_TO_MOUNT[int(m)] for m in mntsta],
        "axis_offset": axis_offset,
        "source_name": source_name,
        "phase_center_ra_dec": np.deg2rad(np.stack([ra_deg, dec_deg], axis=1)),
        "time": Time(date[cross], day_fraction[cross], format="jd", scale="utc"),
        "antenna1": antenna1[cross],
        "antenna2": antenna2[cross],
        "source_index": source_id[cross] - 1,
        "uvw": uvw_seconds[cross] * SPEED_OF_LIGHT,
    }


def compare_uvw(data: dict, method: str = "difxcalc11") -> dict:
    """Compute uvw for every row of a ``read_fits_idi`` result and compare.

    ``method`` is ``"difxcalc11"`` (the embedded difxcalc11 pipeline -- the
    convention DiFX-era archives store), ``"calc"`` (geometric CALC recipe,
    geocenter reference, with the file's axis offsets) or ``"astropy"``
    (pure geometric projection -- the convention pre-DiFX archives store).

    Rows whose baseline code lists the higher antenna first are compared
    against the negated model baseline (``uvw(i, j) = -uvw(j, i)``).

    Returns ``{"model_uvw", "file_uvw", "residual"}`` with the residual as
    ``model - file`` per row [n_row, 3] metres.
    """
    from radio_telescope_delay_model.calculate_uvw import (
        calculate_uvw_astropy,
        calculate_uvw_calc,
    )

    if method not in ("difxcalc11", "calc", "astropy"):
        raise ValueError(f"Unknown method {method!r}.")

    model_rows = np.empty_like(data["uvw"])
    for source in np.unique(data["source_index"]):
        rows = np.flatnonzero(data["source_index"] == source)
        phase_center = data["phase_center_ra_dec"][source : source + 1]
        unique_unix, inverse = np.unique(data["time"].unix[rows], return_inverse=True)
        from astropy.time import Time

        epochs = Time(unique_unix, format="unix", scale="utc")
        if method == "difxcalc11":
            model, b1, b2 = calculate_uvw_calc(
                data["antenna_position"],
                epochs,
                phase_center,
                mode="difxcalc11",
                station_name=data["station_name"],
                mount_type=data["mount_type"],
            )
        elif method == "calc":
            model, b1, b2 = calculate_uvw_calc(
                data["antenna_position"],
                epochs,
                phase_center,
                reference_position=np.zeros(3),
                station_name=data["station_name"],
                mount_type=data["mount_type"],
                axis_offset_metres=data["axis_offset"],
            )
        else:
            model, b1, b2 = calculate_uvw_astropy(
                data["antenna_position"], epochs, phase_center
            )
        pair_index = {
            (int(x), int(y)): k for k, (x, y) in enumerate(zip(b1, b2, strict=True))
        }
        for where, row in enumerate(rows):
            i = int(data["antenna1"][row])
            j = int(data["antenna2"][row])
            if (i, j) in pair_index:
                model_rows[row] = model[inverse[where], pair_index[(i, j)]]
            else:
                model_rows[row] = -model[inverse[where], pair_index[(j, i)]]
    return {
        "model_uvw": model_rows,
        "file_uvw": data["uvw"],
        "residual": model_rows - data["uvw"],
    }

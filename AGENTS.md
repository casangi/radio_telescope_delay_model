# radio_telescope_delay_model — agent guide

Radio-telescope delay model: baseline uvw and delays via **CALC11** (the
geodetic VLBI model used by correlators) and via **astropy**. The package is a
VIPER family member and follows the AstroVIPER conventions (see
`astroviper/AGENTS.md` for the full text); the deltas and specifics are here.

## Layout

```
src/radio_telescope_delay_model/
  calculate_uvw.py        # calculate_uvw_astropy / calculate_uvw_calc / delays
  earth_orientation.py    # IERS EOP via astropy
  data/                   # JPL DE421 ephemeris + CALC coefficient files
  delay_model_cpp/        # pybind11 extension (C++ driver + Fortran core)
    include/ src/ python/ # C++ per the AstroVIPER §6 pybind11 contract
    fortran/              # CALC11 reference core, ported per PORT_PLAN.md
tests/unit/
cmake/RtdmPybind.cmake    # single source of truth for build flags
```

## Conventions

- **UVW**: archival / VLBI convention adopted by MSv4 —
  `uvw = P(antenna1) - P(antenna2)`, right-handed, W towards the source
  (see `astroviper/experiments/uvw_convention_investigation`).
- **pybind11 (AstroVIPER §6)**: typed `py::array_t` without forcecast,
  outputs written in place, GIL released around native work, errors as
  exceptions. The CALC core keeps COMMON-block state, so the driver holds a
  mutex: concurrent calls are safe (serialized); true concurrency comes with
  the C++ port (PORT_PLAN.md).
- **Numerical reproducibility**: `-ffp-contract=off` everywhere (C++ and
  Fortran) and Fortran pinned to `-O2` — the flags under which the C++ driver
  is verified **bit-identical** to the reference CALC build.
- **Fortran → C++ porting**: mechanism, verification recipe and module order
  in PORT_PLAN.md. Never remove a `.f` without `nm`-checking the extension
  still defines every symbol (macOS `-undefined dynamic_lookup` masks
  omissions at link time).
- **ruff** (`E4/E7/E9, F, B, UP, I` + format) via pre-commit; the GitHub
  workflows mirror the other VIPER packages (NRAO shared templates).

## Testing

`python -m pytest tests/` — the uvw tests mirror the example notebook and
cross-check CALC vs astropy (~2 cm on 230 m baselines, earth-model level) and
vs pycalc11 (independent wrapper of the same Fortran; skipped if not
installed). Thread-safety and determinism of the extension are asserted.

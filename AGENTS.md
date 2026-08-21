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
- **Array-general**: every telescope property is a per-antenna input
  (ITRF position, station name for catalog geophysics, mount type, axis
  offset, weather), so the package serves any array — ALMA/VLA-style
  connected arrays through VLBA/EVN/EHT-scale VLBI
  (`tests/unit/test_telescope_arrays.py` exercises each class). The
  Python-exposed driver is `calc_delay_model` (heritage: ALMA's almacalc.f;
  `almacalc` remains as an alias). For arrays whose mean position lies deep
  inside the Earth (global VLBI), the default array-mean reference is
  physically meaningless and a warning recommends the geocenter
  (`reference_position=np.zeros(3)`), the VLBI convention; the
  calc-vs-astropy methods differ at the annual-aberration level
  (~1e-4 relative) by uvw convention, which dominates at VLBI scales.
- **difxcalc11 bit parity**: `mode="difxcalc11"` runs difxcalc11's OWN
  vendored source (`delay_model_cpp/fortran_difxcalc11/`, a second
  self-contained extension `_difxcalc11_ext`; byte-identical to
  `difx/applications/difxcalc11/src` except five runtime data-path OPEN
  sites, see `rtdm_paths.i`) through its own dSTART/dINITL/dSCAN/dDRIVR
  pipeline on a generated `.calc` job (`difxcalc11_core.py`). Raw samples
  verified **bit-identical (6480/6480)** to an instrumented difxcalc
  binary built with the same toolchain and `-ffp-contract=off` (now in
  `RTDM_FORTRAN_FLAGS` — required; contraction produces 1-ULP drift).
  Arbitrary epochs: barycentric quintic through each 2-minute interval's
  6 samples (the `.im` polynomial without its encoding loss; `.im` files
  themselves differ by their GSL-fit/16-digit rounding, u/v ~1e-5 m).
  **Never link the two cores into one binary** — same symbol names,
  different (ALMA vs difxcalc11) implementations; keep each extension's
  symbol set complete (nm-verify, PORT_PLAN.md). The geocenter reference
  is part of the convention: the consensus formula is nonlinear in the
  baseline, so geocentric-vs-array-referenced per-antenna differencing
  differs at (K·b_geo/c)(K·(ω×B)/c) ~ 1e-4 m per 200 m baseline;
  `reference_position=(0,0,0)` selects CALC's geocenter mode in the
  almacalc driver too (`Zero_site=1`). The default `mode="geometric"`
  keeps uvw purely geometric with the array reference. In the almacalc
  driver, unspecified surface weather uses CALC's height-based
  standard-atmosphere defaults (−999 sentinels), exactly as difxcalc11
  with a weatherless `.calc`; `mode="difxcalc11"` rejects weather/axis
  offsets outright (a `.calc` job carries none).
- **Named telescopes**: `station_name=` applies ocean loading and ocean
  pole tide loading from the packaged difxcalc11 catalogs
  (`station_data.py` parses `data/ocean_load.coef`,
  `ocean_pole_tide.coef`, `tilt.dat` exactly as `dinit.f` does; unknown
  names warn and get zeros). `mount_type=` maps difxcalc11 mount strings
  (AZEL/EQUA/XYNS/XYEW/RICH) to CALC axis types. The `OCEG`/`OPTLG` calls
  in `adrvr.f` were reactivated to match difxcalc11's `ddrvr.f` (zero
  coefficients contribute exactly zero — the pinned reference values
  verify this). Axis tilts are parsed but not applied: this CALC 11
  version's axis module has the tilt rotation disabled (same in
  difxcalc11).
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

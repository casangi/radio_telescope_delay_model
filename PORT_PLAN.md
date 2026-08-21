# CALC11 Fortran → C++ port plan

Goal: rewrite the CALC11 core in C++ (AGENTS.md conventions), numerically
identical to the reference Fortran and thread safe. The port proceeds
**module by module with bit-verification at every step**, because a wholesale
rewrite of a 40 000-line geodetic model cannot be verified.

## Status

| layer | files | lines | status |
|---|---|---|---|
| Interface (driver + wrappers) | `almacalc.f`, `calc11_wrapper.f90`, `almastub.f`, `almaout.c` | ~300 | **PORTED to C++** (`delay_model_cpp/src/delay_model.cpp` + `python/bindings.cpp`), **bit-identical** to the reference (geometric, dry and wet delays; verified against the legacy build in clean processes) |
| Computational core | 28 `.f` files below | ~40 000 | Reference Fortran, compiled into the extension by CMake with the reference flags + `-ffp-contract=off` |

## Porting mechanism (per module)

The core is plain Fortran with COMMON-block state; gfortran mangles
`SUBROUTINE FOO` to the C symbol `foo_` with pass-by-reference arguments.
A module is ported by:

1. Implementing its subroutines in C++ with `extern "C"` symbols of the same
   names/signatures (COMMON blocks are declared as `extern "C"` structs
   mirroring the include files exactly — see `delay_model.cpp` for the
   pattern), so the remaining Fortran calls the C++ transparently.
2. Removing the `.f` from `_calc11_core` in
   `src/radio_telescope_delay_model/delay_model_cpp/CMakeLists.txt`.
3. Verifying bit-identity: build with and without the ported module and
   compare `almacalc` outputs on the test scenarios (the same procedure used
   to verify the driver port; the tests carry the reference values).

Caution learned the hard way: on macOS, pybind11 modules link with
`-undefined dynamic_lookup`, so a missing symbol is NOT a link error — it
resolves at runtime against whatever library defines it (verify with
`nm` that every ported/removed symbol is provided by the extension itself).

## Core inventory and proposed order

Port order favours leaf modules (no calls into unported Fortran) and small
files first; sizes are line counts.

| order | file | lines | role | depends on |
|---|---|---|---|---|
| 1 | `cmabd.f` | 59 | matrix add | — |
| 2 | `cvecu.f` | 656 | vector utilities (dot, cross, rotate) | — |
| 3 | `cmatu.f` | 1480 | matrix utilities / rotations | cvecu |
| 4 | `dkill.f` | ~40 | error termination (`terminate_calc_`, `consen_` helpers live in ctheu) | — |
| 5 | `cdtdb.f` | ~700 | TDB time ephemeris | cvecu |
| 6 | `catiu.f` | ~300 | atomic time (TAI–UTC) | — |
| 7 | `cut1m.f` | 1985 | UT1 module | catiu |
| 8 | `cwobm.f` | 1783 | polar motion (wobble) | — |
| 9 | `cnutu.f` + `cnutu6.f` + `cnutm.f` | 15 000 | IAU nutation series (mostly data tables — mechanical) | cmatu |
| 10 | `cdiuu.f` | ~600 | diurnal rotation | cut1m, cwobm |
| 11 | `cm20u.f`, `crosu.f` | ~700 | precession / rotation composition | cmatu |
| 12 | `cpepu.f` | ~900 | JPL ephemeris reader (binary DE421 I/O) | — |
| 13 | `csitm.f`, `dsitu.f` | ~1100 | site geometry | cmatu |
| 14 | `cetdm.f`, `cptdm.f`, `cocem.f`, `hardisp.f` | ~4200 | solid-earth / pole / ocean tides | csitm |
| 15 | `catmm.f` | 2424 | Niell atmosphere (dry/wet mapping) | csitm |
| 16 | `caxom.f` | 1392 | axis offset | csitm |
| 17 | `astrm.f`, `cuvm.f` | ~1000 | aberration / uvw | ctheu |
| 18 | `ctheu.f` | 1388 | Consensus relativistic delay (`thery_`, `consen_`) | everything above |
| 19 | `adrvr.f`, `ainit.f` | ~2600 | per-observation driver + init | all |

Thread safety end-state: once the COMMON blocks are gone (each ported module
takes a context struct), the driver mutex in `delay_model.cpp` can be dropped
and calls can run concurrently. Until then, thread safety is provided by
serialization (concurrent Python calls are safe and verified in the tests).

## Verified reference

The bit-identity reference is the original `Makefile` build (gfortran
`-ffree-form -ffree-line-length-none -O2`) with `-ffp-contract=off` added to
match the VIPER reproducibility policy (contraction differences otherwise
produce few-ULP drift). The CMake build pins `CMAKE_Fortran_FLAGS_RELEASE=-O2`
for the same reason (see pyproject.toml).

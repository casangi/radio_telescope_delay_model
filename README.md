# radio_telescope_delay_model

Baseline **uvw coordinates and delays** for radio interferometry, by two
methods:

- **CALC11** — the geodetic VLBI delay model used by correlators (delays +
  dry/wet troposphere; uvw via the DiFX finite-difference recipe), wrapped as
  a C++/pybind11 extension.
- **astropy** — geometric ITRF→GCRS projections (the SIRIUS/AstroVIPER
  algorithm).

Both return uvw in the archival / VLBI convention adopted by MSv4:
`uvw = P(antenna1) − P(antenna2)`, W towards the source.

```python
from radio_telescope_delay_model import calculate_uvw_calc, calculate_uvw_astropy

uvw, antenna1, antenna2 = calculate_uvw_calc(
    antenna_position_itrf, times, phase_center_ra_dec
)
uvw, antenna1, antenna2 = calculate_uvw_astropy(
    antenna_position_itrf, times, phase_center_ra_dec
)
```

## Install

```
pip install .          # needs a C++23 compiler and gfortran (CALC core)
python -m pytest tests/
```

The JPL DE421 ephemeris and CALC coefficient files ship with the package.

## Status

The CALC11 interface layer is C++ (bit-identical to the reference Fortran
build); the computational core is the reference Fortran, being ported to C++
module by module — see `PORT_PLAN.md`. Conventions and contributor rules:
`AGENTS.md`.

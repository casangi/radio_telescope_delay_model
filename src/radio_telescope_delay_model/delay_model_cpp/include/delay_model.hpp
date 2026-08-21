// C++ driver of the CALC11 delay model -- the port of the ALMA interface
// layer (almacalc.f + calc11_wrapper.f90; see PORT_PLAN.md).
//
// The computational core (aINITL / aDRIVR and everything below them) is still
// the reference Fortran, compiled into the same extension; this driver fills
// the core's COMMON blocks exactly as almacalc.f did and reproduces its
// results bit for bit (verified in the tests). The COMMON blocks make the
// core stateful and therefore NOT thread safe: all entry points serialize on
// a module mutex, so concurrent calls from Python threads are safe (they run
// one at a time). True concurrency arrives with the module-by-module C++
// port.
#pragma once

#include <cstddef>
#include <string>

namespace rtdm::delay_model {

// Geometric (vacuum), dry and wet delays of every antenna relative to the
// array reference position, for every time. Outputs are written in place
// into caller-owned [n_time, n_antenna] row-major arrays (seconds).
//
// Inputs mirror almacalc.f: ITRF metres for positions, MJD (UTC) times, ICRS
// radians for the source per time, EOP (arcsec, arcsec, seconds), leap
// seconds (TAI - UTC), axis offsets in metres and the JPL DE421 ephemeris
// path. Thread safe via internal serialization.
void alma_delay_model(
    const double reference_position[3],
    std::size_t n_antenna,
    const double* antenna_x,
    const double* antenna_y,
    const double* antenna_z,
    const double* temperature_celsius,
    const double* pressure_hpa,
    const double* humidity_fraction,
    std::size_t n_time,
    const double* mjd_utc,
    const double* right_ascension,
    const double* declination,
    const double* polar_motion_x_arcsec,
    const double* polar_motion_y_arcsec,
    const double* ut1_minus_utc_seconds,
    double leap_seconds,
    const double* axis_offset_metres,
    const std::string& ephemeris_path,
    double* geometric_delay,
    double* dry_delay,
    double* wet_delay);

}  // namespace rtdm::delay_model

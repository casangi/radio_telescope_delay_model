// C++ driver over the vendored difxcalc11 pipeline -- see
// fortran_difxcalc11/rtdm_paths.i for the (only) vendored-source patches.
//
// Unlike delay_model.hpp's almacalc driver, this runs difxcalc11's OWN
// dSTART -> dINITL -> dSCAN -> dDRIVR sequence on a .calc file, so every
// input path (EOP tables and their interpolation, leap seconds, station
// catalogs, met defaults, the sample-epoch grid) is difxcalc11's own code:
// the sample values are bit-identical to what the difxcalc binary computes
// internally before its .im polynomial fit (verified in the tests against
// an instrumented reference build). Only the .im writing (d_out + GSL
// fitpoly) is omitted; samples are read back from the /OUT_C/ common.
#pragma once

#include <cstddef>
#include <string>

namespace rtdm::difxcalc {

// Samples per 2-minute interval and their spacing come from difxcalc11
// (c2poly.i: Max_Epochs = 6, d_interval = 24 s; epochs 0..120 s inclusive,
// so consecutive intervals share their boundary sample).
constexpr int kMaxEpochsPerInterval = 6;

// Run the difxcalc11 pipeline on `calc_file`. Outputs are written in place
// into caller-owned row-major arrays sized [max_samples, ...]:
//   ymdhms       [max_samples, 6]          int32 (Y, M, D, h, m, s UTC)
//   delay        [max_samples, n_antenna]  total delay, CALC sign, seconds
//                                          (the .im DELAY equals -1e6 x this)
//   dry, wet     [max_samples, n_antenna]  atmosphere, seconds
//   u, v, w      [max_samples, n_antenna]  per-antenna geocentric model, m
// Returns the number of samples written. Thread safe via serialization on
// the same mutex as the almacalc driver (the cores are separate extensions
// but share process-global state like Fortran unit numbers).
std::size_t run_difxcalc11(
    const std::string& calc_file,
    const std::string& ocean_loading_file,
    const std::string& tilt_file,
    const std::string& ocean_pole_tide_file,
    const std::string& de421_file,
    const std::string& leap_second_file,
    std::size_t n_antenna,
    std::size_t max_samples,
    int* ymdhms,
    double* delay,
    double* dry,
    double* wet,
    double* u,
    double* v,
    double* w);

}  // namespace rtdm::difxcalc

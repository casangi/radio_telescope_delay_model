// pybind11 bindings for the vendored difxcalc11 core (difxcalc_core.hpp).
//
// Same memory contract as bindings.cpp: typed py::array_t without forcecast,
// outputs written in place, GIL released around the Fortran work.

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <stdexcept>
#include <string>

#include "difxcalc_core.hpp"

namespace py = pybind11;

namespace {

std::size_t difxcalc11_samples(
    const std::string& calc_file, const std::string& ocean_loading_file,
    const std::string& tilt_file, const std::string& ocean_pole_tide_file,
    const std::string& de421_file, const std::string& leap_second_file,
    py::array_t<int, py::array::c_style> ymdhms,
    py::array_t<double, py::array::c_style> delay,
    py::array_t<double, py::array::c_style> dry,
    py::array_t<double, py::array::c_style> wet,
    py::array_t<double, py::array::c_style> u,
    py::array_t<double, py::array::c_style> v,
    py::array_t<double, py::array::c_style> w) {
    if (ymdhms.ndim() != 2 || ymdhms.shape(1) != 6)
        throw std::invalid_argument("ymdhms must be [max_samples, 6] int32.");
    const std::size_t max_samples = ymdhms.shape(0);
    if (delay.ndim() != 2 ||
        static_cast<std::size_t>(delay.shape(0)) != max_samples)
        throw std::invalid_argument("delay must be [max_samples, n_antenna].");
    const std::size_t n_antenna = delay.shape(1);
    for (const auto& [array, name] :
         {std::pair{&dry, "dry"}, std::pair{&wet, "wet"}, std::pair{&u, "u"},
          std::pair{&v, "v"}, std::pair{&w, "w"}}) {
        if (array->ndim() != 2 ||
            static_cast<std::size_t>(array->shape(0)) != max_samples ||
            static_cast<std::size_t>(array->shape(1)) != n_antenna)
            throw std::invalid_argument(std::string(name) +
                                        " must be [max_samples, n_antenna].");
    }

    int* ymdhms_data = ymdhms.mutable_data();
    double* delay_data = delay.mutable_data();
    double* dry_data = dry.mutable_data();
    double* wet_data = wet.mutable_data();
    double* u_data = u.mutable_data();
    double* v_data = v.mutable_data();
    double* w_data = w.mutable_data();

    std::size_t n_samples = 0;
    {
        py::gil_scoped_release release;
        n_samples = rtdm::difxcalc::run_difxcalc11(
            calc_file, ocean_loading_file, tilt_file, ocean_pole_tide_file,
            de421_file, leap_second_file, n_antenna, max_samples, ymdhms_data,
            delay_data, dry_data, wet_data, u_data, v_data, w_data);
    }
    return n_samples;
}

}  // namespace

PYBIND11_MODULE(_difxcalc11_ext, m) {
    m.doc() =
        "Vendored difxcalc11 core (dSTART/dINITL/dSCAN/dDRIVR pipeline over "
        "a .calc file; bit-identical to the difxcalc binary's internal "
        "samples).";
    m.def("difxcalc11_samples", &difxcalc11_samples, py::arg("calc_file"),
          py::arg("ocean_loading_file"), py::arg("tilt_file"),
          py::arg("ocean_pole_tide_file"), py::arg("de421_file"),
          py::arg("leap_second_file"), py::arg("ymdhms"), py::arg("delay"),
          py::arg("dry"), py::arg("wet"), py::arg("u"), py::arg("v"),
          py::arg("w"),
          "Run difxcalc11 on a .calc file; fills the sample arrays in place "
          "and returns the number of samples (6 per 2-minute interval, "
          "boundary samples repeated between intervals). Delays are "
          "CALC-sign seconds; u/v/w are the geocentric per-antenna model in "
          "metres.");
}

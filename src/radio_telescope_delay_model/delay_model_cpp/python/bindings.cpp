// pybind11 bindings for the CALC11 delay model.
//
// Follows the AstroVIPER pybind11 memory contract: typed py::array_t without
// forcecast (dtype/contiguity errors instead of silent copies), outputs
// written in place into caller-owned arrays, GIL released around the C++ /
// Fortran work (the core itself serializes on an internal mutex; see
// delay_model.hpp).

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <stdexcept>
#include <string>

#include "delay_model.hpp"

namespace py = pybind11;

namespace {

void almacalc(py::array_t<double, py::array::c_style> reference_position,
              py::array_t<double, py::array::c_style> antenna_position,
              py::array_t<double, py::array::c_style> temperature_celsius,
              py::array_t<double, py::array::c_style> pressure_hpa,
              py::array_t<double, py::array::c_style> humidity_fraction,
              py::array_t<double, py::array::c_style> mjd_utc,
              py::array_t<double, py::array::c_style> right_ascension,
              py::array_t<double, py::array::c_style> declination,
              py::array_t<double, py::array::c_style> polar_motion_x_arcsec,
              py::array_t<double, py::array::c_style> polar_motion_y_arcsec,
              py::array_t<double, py::array::c_style> ut1_minus_utc_seconds,
              double leap_seconds,
              py::array_t<double, py::array::c_style> axis_offset_metres,
              const std::string& ephemeris_path,
              py::array_t<double, py::array::c_style> geometric_delay,
              py::array_t<double, py::array::c_style> dry_delay,
              py::array_t<double, py::array::c_style> wet_delay) {
    if (reference_position.ndim() != 1 || reference_position.shape(0) != 3)
        throw std::invalid_argument("reference_position must be [3].");
    if (antenna_position.ndim() != 2 || antenna_position.shape(1) != 3)
        throw std::invalid_argument("antenna_position must be [n_antenna, 3].");
    const std::size_t n_antenna = antenna_position.shape(0);
    const std::size_t n_time = mjd_utc.shape(0);
    for (const auto& [array, name] :
         {std::pair{&temperature_celsius, "temperature_celsius"},
          std::pair{&pressure_hpa, "pressure_hpa"},
          std::pair{&humidity_fraction, "humidity_fraction"},
          std::pair{&axis_offset_metres, "axis_offset_metres"}}) {
        if (array->ndim() != 1 ||
            static_cast<std::size_t>(array->shape(0)) != n_antenna)
            throw std::invalid_argument(std::string(name) +
                                        " must be [n_antenna].");
    }
    for (const auto& [array, name] :
         {std::pair{&right_ascension, "right_ascension"},
          std::pair{&declination, "declination"},
          std::pair{&polar_motion_x_arcsec, "polar_motion_x_arcsec"},
          std::pair{&polar_motion_y_arcsec, "polar_motion_y_arcsec"},
          std::pair{&ut1_minus_utc_seconds, "ut1_minus_utc_seconds"}}) {
        if (array->ndim() != 1 ||
            static_cast<std::size_t>(array->shape(0)) != n_time)
            throw std::invalid_argument(std::string(name) +
                                        " must be [n_time].");
    }
    for (const auto& [array, name] :
         {std::pair{&geometric_delay, "geometric_delay"},
          std::pair{&dry_delay, "dry_delay"}, std::pair{&wet_delay, "wet_delay"}}) {
        if (array->ndim() != 2 ||
            static_cast<std::size_t>(array->shape(0)) != n_time ||
            static_cast<std::size_t>(array->shape(1)) != n_antenna)
            throw std::invalid_argument(std::string(name) +
                                        " must be [n_time, n_antenna].");
    }

    // Split the antenna positions into the x/y/z arrays the driver takes.
    std::vector<double> x(n_antenna), y(n_antenna), z(n_antenna);
    const double* positions = antenna_position.data();
    for (std::size_t i = 0; i < n_antenna; ++i) {
        x[i] = positions[3 * i];
        y[i] = positions[3 * i + 1];
        z[i] = positions[3 * i + 2];
    }

    const double* reference = reference_position.data();
    double reference_array[3] = {reference[0], reference[1], reference[2]};

    double* geometric = geometric_delay.mutable_data();
    double* dry = dry_delay.mutable_data();
    double* wet = wet_delay.mutable_data();

    {
        py::gil_scoped_release release;
        rtdm::delay_model::alma_delay_model(
            reference_array, n_antenna, x.data(), y.data(), z.data(),
            temperature_celsius.data(), pressure_hpa.data(),
            humidity_fraction.data(), n_time, mjd_utc.data(),
            right_ascension.data(), declination.data(),
            polar_motion_x_arcsec.data(), polar_motion_y_arcsec.data(),
            ut1_minus_utc_seconds.data(), leap_seconds,
            axis_offset_metres.data(), ephemeris_path, geometric, dry, wet);
    }
}

}  // namespace

PYBIND11_MODULE(_delay_model_ext, m) {
    m.doc() =
        "CALC11 delay model (C++ driver over the reference Fortran core; "
        "thread safe via internal serialization).";
    m.def("almacalc", &almacalc, py::arg("reference_position"),
          py::arg("antenna_position"), py::arg("temperature_celsius"),
          py::arg("pressure_hpa"), py::arg("humidity_fraction"),
          py::arg("mjd_utc"), py::arg("right_ascension"), py::arg("declination"),
          py::arg("polar_motion_x_arcsec"), py::arg("polar_motion_y_arcsec"),
          py::arg("ut1_minus_utc_seconds"), py::arg("leap_seconds"),
          py::arg("axis_offset_metres"), py::arg("ephemeris_path"),
          py::arg("geometric_delay"), py::arg("dry_delay"), py::arg("wet_delay"),
          "Geometric / dry / wet delays per antenna relative to the array "
          "reference, written in place into [n_time, n_antenna] float64 "
          "arrays (seconds).");
}

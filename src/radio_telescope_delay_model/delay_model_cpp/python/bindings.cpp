// pybind11 bindings for the CALC11 delay model.
//
// Follows the AstroVIPER pybind11 memory contract: typed py::array_t without
// forcecast (dtype/contiguity errors instead of silent copies), outputs
// written in place into caller-owned arrays, GIL released around the C++ /
// Fortran work (the core itself serializes on an internal mutex; see
// delay_model.hpp).

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "delay_model.hpp"

namespace py = pybind11;

namespace {

using OptionalDoubles =
    std::optional<py::array_t<double, py::array::c_style>>;

// Validate an optional per-antenna station array against the given trailing
// shape and return its data pointer (nullptr when absent).
const double* station_array(const OptionalDoubles& array, const char* name,
                            std::size_t n_antenna,
                            std::initializer_list<std::size_t> trailing) {
    if (!array.has_value()) return nullptr;
    bool ok = static_cast<std::size_t>(array->ndim()) == 1 + trailing.size() &&
              static_cast<std::size_t>(array->shape(0)) == n_antenna;
    std::size_t dim = 1;
    for (std::size_t extent : trailing) {
        ok = ok && static_cast<std::size_t>(array->shape(dim)) == extent;
        ++dim;
    }
    if (!ok) {
        std::string shape = "[n_antenna";
        for (std::size_t extent : trailing)
            shape += ", " + std::to_string(extent);
        throw std::invalid_argument(std::string(name) + " must be " + shape +
                                    "].");
    }
    return array->data();
}

void calc_delay_model(py::array_t<double, py::array::c_style> reference_position,
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
              py::array_t<double, py::array::c_style> wet_delay,
              std::optional<py::array_t<int, py::array::c_style>> axis_type,
              OptionalDoubles ocean_vertical_amplitude,
              OptionalDoubles ocean_vertical_phase,
              OptionalDoubles ocean_horizontal_amplitude,
              OptionalDoubles ocean_horizontal_phase,
              OptionalDoubles ocean_pole_tide_coefficients) {
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

    const int* axis = nullptr;
    if (axis_type.has_value()) {
        if (axis_type->ndim() != 1 ||
            static_cast<std::size_t>(axis_type->shape(0)) != n_antenna)
            throw std::invalid_argument("axis_type must be [n_antenna] int32.");
        axis = axis_type->data();
    }
    const double* ov_amp = station_array(
        ocean_vertical_amplitude, "ocean_vertical_amplitude", n_antenna, {11});
    const double* ov_phs = station_array(
        ocean_vertical_phase, "ocean_vertical_phase", n_antenna, {11});
    const double* oh_amp =
        station_array(ocean_horizontal_amplitude, "ocean_horizontal_amplitude",
                      n_antenna, {2, 11});
    const double* oh_phs =
        station_array(ocean_horizontal_phase, "ocean_horizontal_phase",
                      n_antenna, {2, 11});
    const double* optl =
        station_array(ocean_pole_tide_coefficients,
                      "ocean_pole_tide_coefficients", n_antenna, {6});

    double* geometric = geometric_delay.mutable_data();
    double* dry = dry_delay.mutable_data();
    double* wet = wet_delay.mutable_data();

    {
        py::gil_scoped_release release;
        rtdm::delay_model::calc_delay_model(
            reference_array, n_antenna, x.data(), y.data(), z.data(),
            temperature_celsius.data(), pressure_hpa.data(),
            humidity_fraction.data(), n_time, mjd_utc.data(),
            right_ascension.data(), declination.data(),
            polar_motion_x_arcsec.data(), polar_motion_y_arcsec.data(),
            ut1_minus_utc_seconds.data(), leap_seconds,
            axis_offset_metres.data(), ephemeris_path, geometric, dry, wet,
            axis, ov_amp, ov_phs, oh_amp, oh_phs, optl);
    }
}

}  // namespace

PYBIND11_MODULE(_delay_model_ext, m) {
    m.doc() =
        "CALC11 delay model (C++ driver over the reference Fortran core; "
        "thread safe via internal serialization).";
    m.def("calc_delay_model", &calc_delay_model, py::arg("reference_position"),
          py::arg("antenna_position"), py::arg("temperature_celsius"),
          py::arg("pressure_hpa"), py::arg("humidity_fraction"),
          py::arg("mjd_utc"), py::arg("right_ascension"), py::arg("declination"),
          py::arg("polar_motion_x_arcsec"), py::arg("polar_motion_y_arcsec"),
          py::arg("ut1_minus_utc_seconds"), py::arg("leap_seconds"),
          py::arg("axis_offset_metres"), py::arg("ephemeris_path"),
          py::arg("geometric_delay"), py::arg("dry_delay"), py::arg("wet_delay"),
          py::arg("axis_type") = py::none(),
          py::arg("ocean_vertical_amplitude") = py::none(),
          py::arg("ocean_vertical_phase") = py::none(),
          py::arg("ocean_horizontal_amplitude") = py::none(),
          py::arg("ocean_horizontal_phase") = py::none(),
          py::arg("ocean_pole_tide_coefficients") = py::none(),
          "Geometric / dry / wet delays per antenna relative to the array "
          "reference, written in place into [n_time, n_antenna] float64 "
          "arrays (seconds). The optional trailing arrays carry per-antenna "
          "station data (axis types, ocean loading, ocean pole tide "
          "loading) as difxcalc11's dinit.f loads them.");
}

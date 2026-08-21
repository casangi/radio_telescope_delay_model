cmake_minimum_required(VERSION 3.15...3.27)
project(radio_telescope_delay_model VERSION ${SKBUILD_PROJECT_VERSION}
        LANGUAGES CXX Fortran)

# All shared C++ / Fortran / pybind11 build settings live in
# cmake/RtdmPybind.cmake -- edit THAT file to change a flag for every
# extension at once (the AstroVIPER pattern).
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(RtdmPybind)

add_subdirectory(src/radio_telescope_delay_model/delay_model_cpp)

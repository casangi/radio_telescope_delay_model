# ---------------------------------------------------------------------------
# RtdmPybind.cmake
#
# Single source of truth for how every radio_telescope_delay_model C++ /
# pybind11 extension is compiled -- the direct analogue of AstroVIPER's
# AstroviperPybind.cmake (same C++23 dialect, reproducibility flags and
# helper-function pattern), extended with Fortran support for the CALC11
# computational core that is still awaiting its C++ port (see PORT_PLAN.md).
#
# Two ways this file is used:
#   * From the repo-root CMakeLists.txt (the normal scikit-build-core package
#     build): included once, then each module's CMakeLists.txt calls
#     rtdm_add_pybind_module().
#   * Standalone (cmake -S <module_dir> -B build): the module's CMakeLists.txt
#     finds and includes this same file, so a standalone build uses
#     byte-for-byte the same flags as the package build.
# ---------------------------------------------------------------------------
include_guard(GLOBAL)

# Toolchain floor: complete C++23 implementations only (see AstroVIPER's
# AGENTS.md for the reasoning; the floors are identical).
set(RTDM_MIN_GNU_VERSION        14)
set(RTDM_MIN_CLANG_VERSION      17)
set(RTDM_MIN_APPLECLANG_VERSION 16)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS RTDM_MIN_GNU_VERSION)
        message(FATAL_ERROR
            "radio_telescope_delay_model requires GCC >= "
            "${RTDM_MIN_GNU_VERSION} for C++23, but found GCC "
            "${CMAKE_CXX_COMPILER_VERSION} (${CMAKE_CXX_COMPILER}).")
    endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS RTDM_MIN_CLANG_VERSION)
        message(FATAL_ERROR
            "radio_telescope_delay_model requires Clang >= "
            "${RTDM_MIN_CLANG_VERSION} for C++23, but found Clang "
            "${CMAKE_CXX_COMPILER_VERSION}.")
    endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS RTDM_MIN_APPLECLANG_VERSION)
        message(FATAL_ERROR
            "radio_telescope_delay_model requires Apple Clang >= "
            "${RTDM_MIN_APPLECLANG_VERSION} for C++23, but found Apple Clang "
            "${CMAKE_CXX_COMPILER_VERSION}.")
    endif()
endif()

# The CALC11 core is (for now) compiled from its reference Fortran; the
# gfortran runtime is linked into the extension automatically because the
# module target carries Fortran sources.
enable_language(Fortran)

set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG QUIET)
if(NOT pybind11_FOUND)
    find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
    execute_process(
        COMMAND "${Python_EXECUTABLE}" -c
                "import pybind11; print(pybind11.get_cmake_dir())"
        OUTPUT_VARIABLE _rtdm_pybind11_cmakedir
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET)
    if(_rtdm_pybind11_cmakedir)
        list(APPEND CMAKE_PREFIX_PATH "${_rtdm_pybind11_cmakedir}")
    endif()
    find_package(pybind11 CONFIG REQUIRED)
endif()

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

add_library(rtdm_cpp_flags INTERFACE)

set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
target_compile_features(rtdm_cpp_flags INTERFACE cxx_std_23)

# Scientific reproducibility: no floating-point expression contraction, as in
# AstroVIPER -- results must be bit-identical across IEEE-754 platforms.
target_compile_options(rtdm_cpp_flags INTERFACE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-ffp-contract=off>
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-pipe>
)
target_compile_options(rtdm_cpp_flags INTERFACE
    $<$<CXX_COMPILER_ID:Clang,AppleClang>:-stdlib=libc++>)
target_link_options(rtdm_cpp_flags INTERFACE
    $<$<CXX_COMPILER_ID:Clang,AppleClang>:-stdlib=libc++>)

if(DEFINED SKBUILD_PROJECT_VERSION)
    set(_rtdm_version "${SKBUILD_PROJECT_VERSION}")
else()
    set(_rtdm_version "0.0.0+local")
endif()
target_compile_definitions(rtdm_cpp_flags INTERFACE
    VERSION_INFO="${_rtdm_version}")

target_link_libraries(rtdm_cpp_flags INTERFACE Threads::Threads)

# Flags for the legacy CALC11 Fortran core: exactly the reference build's
# (-ffree-form -ffree-line-length-none -O2), plus -fPIC for the shared module.
# -frecursive puts local variables on the stack (no static locals), one of the
# prerequisites for eventually calling the core from multiple threads; the
# COMMON blocks still make the core single-threaded today (see PORT_PLAN.md).
set(RTDM_FORTRAN_FLAGS
    -ffree-form -ffree-line-length-none -fPIC -O2)

# ---------------------------------------------------------------------------
# rtdm_add_pybind_module(NAME <name> SOURCES <src>...
#                        [FORTRAN_SOURCES <f>...]
#                        [INCLUDE_DIRS <dir>...] [DESTINATION <dir>])
#
# Declare a pybind11 extension with all shared flags applied; Fortran sources
# (the unported CALC11 core) are compiled into the same module with the
# reference flags and the gfortran runtime linked automatically.
# ---------------------------------------------------------------------------
function(rtdm_add_pybind_module)
    cmake_parse_arguments(PB "" "NAME;DESTINATION"
                          "SOURCES;FORTRAN_SOURCES;INCLUDE_DIRS" ${ARGN})

    if(NOT PB_NAME)
        message(FATAL_ERROR "rtdm_add_pybind_module: NAME is required")
    endif()
    if(NOT PB_SOURCES)
        message(FATAL_ERROR "rtdm_add_pybind_module(${PB_NAME}): SOURCES is required")
    endif()
    if(NOT PB_INCLUDE_DIRS)
        set(PB_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include")
    endif()
    if(NOT PB_DESTINATION)
        if(EXISTS "${CMAKE_SOURCE_DIR}/src")
            file(RELATIVE_PATH PB_DESTINATION
                 "${CMAKE_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}")
        else()
            set(PB_DESTINATION ".")
        endif()
    endif()

    pybind11_add_module(${PB_NAME} ${PB_SOURCES} ${PB_FORTRAN_SOURCES})
    if(PB_FORTRAN_SOURCES)
        set_source_files_properties(${PB_FORTRAN_SOURCES} PROPERTIES
            COMPILE_OPTIONS "${RTDM_FORTRAN_FLAGS}")
    endif()
    target_include_directories(${PB_NAME} PRIVATE ${PB_INCLUDE_DIRS})
    target_link_libraries(${PB_NAME} PRIVATE rtdm_cpp_flags)
    install(TARGETS ${PB_NAME} DESTINATION ${PB_DESTINATION})
endfunction()

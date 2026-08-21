cmake_minimum_required(VERSION 3.15...3.27)

# Standalone build support, as in AstroVIPER's kernels: when configured
# directly, bootstrap a project and locate the shared cmake/ directory.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(rtdm_delay_model LANGUAGES CXX Fortran)
    get_filename_component(_rtdm_repo_root
        "${CMAKE_CURRENT_SOURCE_DIR}/../../.." ABSOLUTE)
    list(APPEND CMAKE_MODULE_PATH "${_rtdm_repo_root}/cmake")
endif()

include(RtdmPybind)

# The CALC11 computational core, compiled from the reference Fortran until its
# C++ port lands module by module (status: PORT_PLAN.md). The interface layer
# (almacalc.f, calc11_wrapper.f90, almastub.f, almaout.c of the original
# package) is already replaced by src/delay_model.cpp + python/bindings.cpp.
set(_calc11_core
    fortran/adrvr.f fortran/ainit.f fortran/astrm.f fortran/catiu.f
    fortran/catmm.f fortran/caxom.f fortran/cctiu.f fortran/cdiuu.f
    fortran/cdtdb.f fortran/cetdm.f fortran/cm20u.f fortran/cmabd.f
    fortran/cmatu.f fortran/cnutm.f fortran/cnutu.f fortran/cnutu6.f
    fortran/cocem.f fortran/cpepu.f fortran/cptdm.f fortran/crosu.f
    fortran/csitm.f fortran/ctheu.f fortran/cut1m.f fortran/cuvm.f fortran/cvecu.f
    fortran/cwobm.f fortran/dkill.f fortran/dsitu.f fortran/hardisp.f
)

rtdm_add_pybind_module(
    NAME _delay_model_ext
    SOURCES src/delay_model.cpp python/bindings.cpp
    FORTRAN_SOURCES ${_calc11_core}
)
target_include_directories(_delay_model_ext PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/fortran")

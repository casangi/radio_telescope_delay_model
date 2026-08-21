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

# The vendored difxcalc11 pipeline (bit parity with the difxcalc binary's
# internal samples). A separate, self-contained extension: several modules
# (catmm, csitm, cut1m, ...) exist in both cores with the same symbol names
# but different (ALMA vs difxcalc11) implementations, so they must never be
# linked together -- keep each module's symbol set complete (nm-verify; see
# PORT_PLAN.md on macOS -undefined dynamic_lookup masking omissions).
# Sources are byte-identical to difx/applications/difxcalc11/src except the
# five runtime-path OPEN sites (see fortran_difxcalc11/rtdm_paths.i);
# dmain.f / d_out.f / fitpoly_.c / difxiohelper.c (the .im writing stage and
# its GSL dependency) are replaced by src/difxcalc_core.cpp.
set(_difxcalc11_core
    fortran_difxcalc11/catiu.f fortran_difxcalc11/catmm.f
    fortran_difxcalc11/caxom.f fortran_difxcalc11/cctiu.f
    fortran_difxcalc11/cdiuu.f fortran_difxcalc11/cdtdb.f
    fortran_difxcalc11/cetdm.f fortran_difxcalc11/cm20u.f
    fortran_difxcalc11/cmabd.f fortran_difxcalc11/cmatu.f
    fortran_difxcalc11/cnutm.f fortran_difxcalc11/cnutu.f
    fortran_difxcalc11/cnutu6.f fortran_difxcalc11/cocem.f
    fortran_difxcalc11/cpepu.f fortran_difxcalc11/cptdm.f
    fortran_difxcalc11/crosu.f fortran_difxcalc11/csitm.f
    fortran_difxcalc11/cstrm.f fortran_difxcalc11/ctheu.f
    fortran_difxcalc11/cut1m.f fortran_difxcalc11/cuvm.f
    fortran_difxcalc11/cvecu.f fortran_difxcalc11/cwobm.f
    fortran_difxcalc11/ddrvr.f fortran_difxcalc11/dinit.f
    fortran_difxcalc11/dkill.f fortran_difxcalc11/dscan.f
    fortran_difxcalc11/dsitu.f fortran_difxcalc11/dstrt.f
    fortran_difxcalc11/hardisp.f fortran_difxcalc11/nfmodel.f
    fortran_difxcalc11/nfstrm.f
)

rtdm_add_pybind_module(
    NAME _difxcalc11_ext
    SOURCES src/difxcalc_core.cpp python/difxcalc_bindings.cpp
    FORTRAN_SOURCES ${_difxcalc11_core}
)
target_include_directories(_difxcalc11_ext PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/fortran_difxcalc11")

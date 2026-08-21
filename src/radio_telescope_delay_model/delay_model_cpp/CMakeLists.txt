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

# The difxcalc11 pipeline (bit parity with the difxcalc binary's internal
# samples). A separate, self-contained extension: several modules
# (catmm, csitm, cut1m, ...) exist in both cores with the same symbol names
# but different (ALMA vs difxcalc11) implementations, so they must never be
# linked together -- keep each module's symbol set complete (nm-verify; see
# PORT_PLAN.md on macOS -undefined dynamic_lookup masking omissions).
#
# Sources come from the pinned DiFX submodule (extern/difx; init with
# scripts/init_difx_submodule.sh) and are staged into the build tree at
# configure time with exactly one rtdm patch applied
# (patches/difxcalc11_runtime_paths.patch): the five data-file OPEN sites
# redirected from param11.i's configure-time paths to the runtime
# /RTDMPATHS/ common that src/difxcalc_core.cpp fills; the patch also adds
# rtdm_paths.i, that common's declaration. dmain.f (CLI main + GETCL) and
# d_out.f (the .im writing stage, whose GSL fitpoly_.c / difxiohelper.c we
# also never compile) are replaced by src/difxcalc_core.cpp. Bumping the
# submodule SHA is a deliberate action gated by the tests -- see AGENTS.md.
get_filename_component(_rtdm_repo_root
    "${CMAKE_CURRENT_SOURCE_DIR}/../../.." ABSOLUTE)
set(_difx_src "${_rtdm_repo_root}/extern/difx/applications/difxcalc11/src")
if(NOT EXISTS "${_difx_src}/dinit.f")
    message(FATAL_ERROR
        "difxcalc11 sources not found at ${_difx_src}. Initialize the "
        "pinned DiFX submodule first: scripts/init_difx_submodule.sh "
        "(shallow + sparse; a plain 'git submodule update --init "
        "extern/difx' also works but materializes the full DiFX worktree).")
endif()
find_program(RTDM_PATCH_EXECUTABLE patch)
if(NOT RTDM_PATCH_EXECUTABLE)
    message(FATAL_ERROR
        "The 'patch' tool is required to stage the difxcalc11 sources.")
endif()

set(_difxcalc11_staged "${CMAKE_CURRENT_BINARY_DIR}/difxcalc11_staged")
file(REMOVE_RECURSE "${_difxcalc11_staged}")
file(MAKE_DIRECTORY "${_difxcalc11_staged}")
file(GLOB _difx_upstream "${_difx_src}/*.f" "${_difx_src}/*.i")
file(COPY ${_difx_upstream} DESTINATION "${_difxcalc11_staged}")
file(REMOVE "${_difxcalc11_staged}/dmain.f" "${_difxcalc11_staged}/d_out.f")

# param11.i: upstream's configure generates it from param11.i.in. Its baked
# paths are dead code here (every OPEN that read them is patched to
# /RTDMPATHS/), so substitute inert placeholders and keep the real constants
# (JPL record sizes etc.) exactly as upstream defines them.
file(READ "${_difx_src}/param11.i.in" _param11)
string(REGEX REPLACE "@[A-Za-z0-9_]+@" "/nonexistent" _param11 "${_param11}")
file(WRITE "${_difxcalc11_staged}/param11.i" "${_param11}")

execute_process(
    COMMAND "${RTDM_PATCH_EXECUTABLE}" -p1 --batch --forward
            -i "${CMAKE_CURRENT_SOURCE_DIR}/patches/difxcalc11_runtime_paths.patch"
    WORKING_DIRECTORY "${_difxcalc11_staged}"
    RESULT_VARIABLE _patch_status
    OUTPUT_VARIABLE _patch_output
    ERROR_VARIABLE _patch_output)
if(NOT _patch_status EQUAL 0)
    message(FATAL_ERROR
        "Applying difxcalc11_runtime_paths.patch to the DiFX submodule "
        "sources failed -- is extern/difx off its pinned SHA, or did "
        "upstream change the OPEN sites? Regenerate the patch against the "
        "new sources if this was a deliberate bump.\n${_patch_output}")
endif()

file(GLOB _difxcalc11_core "${_difxcalc11_staged}/*.f")

rtdm_add_pybind_module(
    NAME _difxcalc11_ext
    SOURCES src/difxcalc_core.cpp python/difxcalc_bindings.cpp
    FORTRAN_SOURCES ${_difxcalc11_core}
)
target_include_directories(_difxcalc11_ext PRIVATE "${_difxcalc11_staged}")

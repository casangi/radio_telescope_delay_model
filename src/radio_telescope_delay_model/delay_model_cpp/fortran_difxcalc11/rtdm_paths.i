!  rtdm_paths.i
!
!  Runtime data-file paths for the vendored difxcalc11 core.
!
!  difxcalc11 bakes its station-catalog and ephemeris paths into param11.i
!  at configure time, which a relocatable Python package cannot do. The
!  five file-OPEN sites (dOCNIN, dANTILT, dOPTLIN in dinit.f; the ephemeris
!  OPEN in cpepu.f; get_leapsec in cut1m.f) are patched -- the only rtdm
!  modifications to the vendored source -- to use these variables, which the
!  C++ driver (difxcalc_core.cpp) fills with the packaged data locations
!  before dSTART/dINITL run. Everything else is byte-identical to
!  difx/applications/difxcalc11/src.
!
      CHARACTER*128 RTDM_OC_FILE, RTDM_TILT_FILE, RTDM_OPTL_FILE,       &
     &              RTDM_DE421_FILE, RTDM_LEAP_FILE
      COMMON / RTDMPATHS / RTDM_OC_FILE, RTDM_TILT_FILE,                &
     &              RTDM_OPTL_FILE, RTDM_DE421_FILE, RTDM_LEAP_FILE

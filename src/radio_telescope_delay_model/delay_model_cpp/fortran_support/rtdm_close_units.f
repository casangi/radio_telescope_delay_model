!  rtdm_close_units.f  (rtdm's own helper -- NOT part of the DiFX submodule)
!
!  difxcalc11 acquires Fortran unit numbers through get4unit() for its data
!  files (.calc, ephemeris, catalogs, leap seconds) and, driven repeatedly
!  in one process, never releases them: the embedded backend would exhaust
!  the unit space after a few dozen runs (the standalone binary runs once
!  per process and never sees this; get4unit's lifetime cap is lifted by
!  patches/difxcalc11_embedded_use.patch).
!
!  Fortran unit numbers are PROCESS-global (one libgfortran unit table is
!  shared with the almacalc extension), so blanket-closing a range would
!  destroy the other core's open-file state behind its back. Instead the
!  driver snapshots which units are open before a run and closes only the
!  ones that run opened. Units 1-6 are never touched (5/6 are stdin/stdout,
!  where CALC writes its warnings).
      Subroutine rtdm_units_snapshot()
      Implicit None
      Integer*4 I
      Logical Od
      Logical Was_open(7:99)
      Common / RTDMUNITS / Was_open
      Do I = 7, 99
        Inquire (Unit=I, Opened=Od)
        Was_open(I) = Od
      Enddo
      Return
      End
!
      Subroutine rtdm_close_new_units()
      Implicit None
      Integer*4 I
      Logical Od
      Logical Was_open(7:99)
      Common / RTDMUNITS / Was_open
      Do I = 7, 99
        Inquire (Unit=I, Opened=Od)
        If (Od .and. .not. Was_open(I)) Close(I)
      Enddo
      Return
      End

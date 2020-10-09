!     Include file param11.i 
!
!     ******************************************************************
!     *      CALC 11 Version                                           *
!     ******************************************************************
!
!     JPL_DE421 is the path name for the JPL ephemeris file.
!
!xx   Integer*4   N_RECL_JPL
!xx   PARAMETER ( N_RECL_JPL =    4 )
!xx   Integer*4    K_SIZE_JPL
!xx   PARAMETER ( K_SIZE_JPL = 2036 )
!xx   Integer*4   I_RECL_JPL
!
!xx   CHARACTER  JPL_DE421*128
!xx   PARAMETER (JPL_DE421 = '/share/DE421_little_Endian' )
!
      CHARACTER A_TILTS*128
      PARAMETER ( A_TILTS = '/share/tilt.dat' )
!
      CHARACTER OC_file*128
      PARAMETER ( OC_file = '/share/ocean_load.coef' )
!
      CHARACTER OPTL_file*128
      PARAMETER (OPTL_file = '/share/ocean_pole_tide.coef' )
!
!  Leap seconds file (Not needed by difx)
      CHARACTER DFLEAP*128
      PARAMETER ( DFLEAP =  '/share/ut1ls.dat' )
!
! EOP file
      CHARACTER EOP_file*128
      PARAMETER (EOP_file =  '/share/usno_finals.erp')

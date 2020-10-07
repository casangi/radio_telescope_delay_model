      SUBROUTINE almacalc(refx, refy, refz,                             &
     &                    nant, antx, anty, antz,                       &
     &                    temp, pressure, humidity,                     &
     &                    ntimes, mjd, ra, dec, ssobj,                  &
     &                    dx, dy, dut, leapsec, axisoff, sourcename,    &
     &                    JPX_DE421,                                    &
     &                    geodelay, drydelay, wetdelay)
      IMPLICIT None
!
      INCLUDE 'd_input.i'
      INCLUDE 'cmxst11.i'
      INCLUDE 'cmwob11.i'
!
      Real*8           PI, TWOPI, HALFPI, CONVD, CONVDS, CONVHS, SECDAY
      COMMON / CMATH / PI, TWOPI, HALFPI, CONVD, CONVDS, CONVHS, SECDAY
!              1. CONVDS  -  The conversion factor of radians per arcsecond.
      Real*8  ATMUTC(3), ROTEPH(2,20), A1UTC(3), A1DIFF(3)
      COMMON / EOPCM / ATMUTC, ROTEPH, A1UTC, A1DIFF
!            1. ATMUTC(3)   - THE 'TAI MINUS UTC' INFORMATION ARRAY. THIS ARRAY
!                             CONTAINS RESPECTIVELY THE EPOCH, VALUE, AND TIME
!                             RATE OF CHANGE OF 'TAI MINUS UTC'. Used in the
!                             atomic time module. (DAYS, SEC, SEC/SEC)
!
      Real*8 UTC, UT1, DUT1AT, WOBXL, WOBYL, dWOBXL, dWOBYL,            &
     &       deltime, XJD, YJD, fjd, yleap(5),                          &
     &       SITLON(2), SITLAT(2), CFSITE(3,2), right_asc, declination
      Real*8 surPr(2), surTp(2), surHm(2)
      Real*8 delay_vac, dry_atm1, dry_atm2, wet_atm1, wet_atm2
      Integer*4 ILU, I, J, ierr, Iref, Iremot,                          &
     &                           Iscan, ierr4, IJD
      Logical*4 ss_obj
!
!   Subroutine almacalc:
!     2019-Feb-25  D. Gordon - Initial preliminary version
!     2019-Mar-01  D. Gordon - Updated preliminary version, adding leapsec,   
!                              axisoff and sourcename.
!     2019-Mar-24  D. Gordon - Removed reference met data and atmosphere corection,
!                              making JPL ephemeris path and name an input parameter.
!     2019-April   D. Gordon - Some cleanup, removed some include files, etc. 
!
!  INPUTS:
! Geocentric position of the array reference point (ITRF).
      REAL*8 refx, refy, refz
! Number of stations
      INTEGER*4 nant
! Number of time slots to compute delays
      INTEGER*4 ntimes
! Geocentric (ITRF) position of each antenna.
      REAL*8 antx(nant), anty(nant), antz(nant)
! Temperature (deg. C), Pressure (hPa/mbar), and humidity (0-1) at each antenna
      REAL*8 temp(nant), pressure(nant), humidity(nant)
! The times when to compute the delays
      REAL*8 mjd(ntimes)
! ICRS Source position (radians) for each time slot
      REAL*8 ra(ntimes), dec(ntimes)
! Earth orientation parameters at each time (arc-sec, arc-sec, sec)
      REAL*8 dx(ntimes), dy(ntimes), dut(ntimes)
! # of leap seconds
      Real*8 leapsec 
! axis offsets (meters) for each antenna
      Real*8 axisoff(nant)
! # source names, for future use with solar system objects
      Character*8 sourcename(ntimes)
! Path name of the JPL ephemeris
      CHARACTER*128  JPX_DE421
! True if the source is a solar system object.
      LOGICAL*4 ssobj(ntimes) 
!
!  OUTPUTS:
! Geometric delays for each antenna at each timeslot
      REAL*8 geodelay(ntimes,nant)
! Dry delay for each antenna at each timeslot
      REAL*8 drydelay(ntimes,nant)
! Wet delay for each antenna at each timeslot
      REAL*8 wetdelay(ntimes,nant) 
!
!  Fill calc site arrays
       NUMSIT = nant + 1
       SITXYZ(1,1) = refx
       SITXYZ(2,1) = refy
       SITXYZ(3,1) = refz
       SITAXO(1)   = 0.0D0 
       KTYPE(1) = 3
      Do i = 1, nant
       SITXYZ(1,i+1) = antx(i)
       SITXYZ(2,i+1) = anty(i)
       SITXYZ(3,i+1) = antz(i)
       SITAXO(i+1)   = axisoff(i)
       KTYPE(i+1) = 3   ! alt-az mount
      Enddo
!
!  Name of ephemeris file
      JPL_DE421 = JPX_DE421
!
!   Initialize the model modules and the necessary utilities.
      CALL aINITL
!
       fjd = mjd(1) + 2400000.5d0
       deltime = (mjd(ntimes)-mjd(1)) * 86400.d0
       DUT1AT = 1.d0 + ((dut(ntimes)-dut(1))/deltime)
       dWOBXL = (dx(ntimes)-dx(1))*1.D3/deltime     ! radians/sec
       dWOBYL = (dy(ntimes)-dy(1))*1.D3/deltime     ! radians/sec
       dWOBX = dWOBXL * CONVDS * 1.0D-3     ! radians/sec
       dWOBY = dWOBYL * CONVDS * 1.0D-3     ! radians/sec
       ATMUTC(1) = fjd 
       ATMUTC(2) = leapsec
       ATMUTC(3) = 0.d0    
!
!   Compute the delays 
      Iref = 1    ! Array reference position 
      Do Iscan = 1, ntimes      ! all time tags
       IJD = mjd(Iscan)     ! Integer MJD at midnight
       YJD = IJD            ! Real*8 MJD at midnight
       XJD = YJD + 2400000.5d0    ! JD at midnight
       UTC = (mjd(Iscan)-YJD)   ! UTC, fraction of day
       UT1 = UTC*86400.d0 + dut(Iscan)        ! UT1 in seconds
       WOBXL = dx(Iscan) * 1.d3      ! milli-arc-sec
       WOBYL = dy(Iscan) * 1.d3      ! milli-arc-sec
       WOBX =   WOBXL * CONVDS * 1.0D-3     ! radians
       WOBY =   WOBYL * CONVDS * 1.0D-3     ! radians
       right_asc = ra(iscan)
       declination = dec(iscan)
       surTp(Iref) = -99.0d0
       surPr(Iref) = -99.0d0
       surHm(Iref) = -99.0d0
       ss_obj = ssobj(iscan)
!        
       Do Iremot = 2, nant+1
!
        CFSITE(1,1) = SITXYZ(1,1)
        CFSITE(2,1) = SITXYZ(2,1)
        CFSITE(3,1) = SITXYZ(3,1)
        CFSITE(1,2) = SITXYZ(1,Iremot)
        CFSITE(2,2) = SITXYZ(2,Iremot)
        CFSITE(3,2) = SITXYZ(3,Iremot)
        SITLON(1)   = XLON(Iref)
        SITLAT(1)   = XLAT(Iref)
        SITLON(2)   = XLON(Iremot)
        SITLAT(2)   = XLAT(Iremot)
        surTp(2) = temp(Iremot-1)
        surPr(2) = pressure(Iremot-1)
        surHm(2) = humidity(Iremot-1)
!
       CALL aDRIVR (Iref, Iremot, UTC, XJD, UT1, DUT1AT, CFSITE,        &
     &              SITLON, SITLAT, right_asc, declination, surTp,      &
     &              surPr, surHm, ss_obj,                               &
     &              delay_vac, dry_atm1, dry_atm2, wet_atm1, wet_atm2)
!
       geodelay(Iscan,Iremot-1) = delay_vac
       drydelay(Iscan,Iremot-1) = dry_atm2         ! remote station only
       wetdelay(Iscan,Iremot-1) = wet_atm2         ! remote station only
!
       Enddo
      Enddo
!
      Return 
      END

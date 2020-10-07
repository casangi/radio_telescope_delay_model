      SUBROUTINE almastub(refx, refy, refz,                             &
                          nant, antx, anty, antz,                       &
                          temp, pressure, humidity,                     &
                          ntimes, mjd, ra, dec, ssobj,                  & 
                          dx, dy, dut, taimutc, axisoff, sourcename,    &
     &                    JPX_DE421,                                    &
                          geodelay, drydelay, wetdelay)
      IMPLICIT NONE
      REAL*8 refx, refy, refz, reftemp, refpress, refhum
      INTEGER*4 nant
      REAL*8 antx(nant), anty(nant), antz(nant)
      REAL*8 temp(nant), pressure(nant), humidity(nant)
      INTEGER*4 ntimes
      REAL*8 mjd(ntimes), ra(ntimes), dec(ntimes)
      LOGICAL*4 ssobj(ntimes)
      Character*8 sourcename(ntimes)
      CHARACTER*128  JPX_DE421
      Real*8 taimutc
      Real*8 axisoff(nant)
      REAL*8 dx(ntimes), dy(ntimes), dut(ntimes)
      REAL*8 geodelay(ntimes,nant)
      REAL*8 wetdelay(ntimes,nant), drydelay(ntimes,nant)

      INTEGER*4 i, j, k
      WRITE(*,*) 'Reference position', refx, refy, refz
      WRITE(*,*) 'JPL DE421 ephemeris is at ', JPX_DE421
      WRITE(*,*) 'Number of antennas ', nant
      DO i = 1, nant
         WRITE(*,*) 'Antenna pos', antx(i), anty(i), antz(i), axisoff(i)
      ENDDO
      DO i = 1, nant
         WRITE(*,*) 'Antenna weather', temp(i), pressure(i), humidity(i) 
      ENDDO
      WRITE(*,*) 'Number of timeslots', ntimes, 'TAI-UTC', taimutc
      DO i = 1, ntimes
         WRITE(*,*) 'MJD', mjd(i)
      ENDDO
      DO i = 1, ntimes
         WRITE(*,*) 'Source ', sourcename(i), ra(i), dec(i), ssobj(i)
      ENDDO
      DO i = 1, ntimes
         WRITE(*,*) 'EOP ', dx(i), dy(i), dut(i)
      ENDDO
      k = 0
      DO i = 1, ntimes
         DO j = 1, nant
            geodelay(i, j) = j*100.0 + i
            drydelay(i, j) = j*100.0 + i/10.0
            wetdelay(i, j) = j*100.0 + i/100.0
         ENDDO
      ENDDO
      RETURN
      END


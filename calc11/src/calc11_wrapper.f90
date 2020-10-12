subroutine wrapper_almacalc(refx, refy, refz, nant, antx, anty, antz, &
    temp, pressure, humidity, ntimes, mjd, ra, dec, ssobj, dx, dy, &
    dut, leapsec, axisoff, sourcename, jpx_de421, geodelay, &
    drydelay, wetdelay) bind(c)

  use iso_c_binding

  implicit none

  integer(c_int), intent(in), value :: nant, ntimes
  real(c_double), intent(in), value :: refx, refy, refz, leapsec
  real(c_double), dimension(nant), intent(in) :: antx, anty, antz, &
      temp, pressure, humidity, axisoff
  real(c_double), dimension(ntimes), intent(in) :: mjd, ra, dec, dx, dy, dut
  !character(len=8), dimension(ntimes), intent(in) :: sourcename
  character(kind=c_char), dimension(ntimes), intent(in) :: sourcename(8)
  !character(len=128), dimension(1), intent(in) :: jpx_de421
  character(kind=c_char), dimension(1), intent(in) :: jpx_de421(128)
  logical(c_bool), dimension(ntimes), intent(in) :: ssobj
  real(c_double), dimension(ntimes, nant), intent(out) :: geodelay, &
      drydelay, wetdelay

  logical, dimension(ntimes) :: ssobj_tmp

  ssobj_tmp = ssobj

  call almacalc(refx, refy, refz, nant, antx, anty, antz, temp, &
      pressure, humidity, ntimes, mjd, ra, dec, ssobj_tmp, dx, dy, &
      dut, leapsec, axisoff, sourcename, jpx_de421, geodelay, &
      drydelay, wetdelay)

end subroutine wrapper_almacalc


subroutine wrapper_adrivr(iref, iremot, utc, xjd, ut1, dut1at, cfsite, &
    sitlon, sitlat, right_asc, declination, &
    surTp, surPr, surHm, ssobj, delay_vac, &
    dry_atm1, dry_atm2, wet_atm1, wet_atm2) bind(c)

  use iso_c_binding

  implicit none

  integer(c_int), intent(in), value :: iref, iremot
  real(c_double), intent(in), value :: utc, xjd, ut1, dut1at, &
       right_asc, declination
  real(c_double), dimension(2), intent(in) :: sitlon, sitlat, &
       surTp, surPr, surHm
  real(c_double), dimension(3,2), intent(in) :: cfsite(3,2)
  logical(c_bool),  intent(in) :: ssobj
  real(c_double), intent(out) :: delay_vac, dry_atm1, dry_atm2, wet_atm1, wet_atm2
  
  call aDRIVR(iref, iremot, utc, xjd, ut1, dut1at, cfsite, &
    sitlon, sitlat, right_asc, declination, &
    surTp, surPr, surHm, ssobj, delay_vac, &
    dry_atm1, dry_atm2, wet_atm1, wet_atm2)

end subroutine wrapper_adrivr

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
  character(len=8), dimension(ntimes), intent(in) :: sourcename
  character(len=128), dimension(1), intent(in) :: jpx_de421
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

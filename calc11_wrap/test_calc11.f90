program test_calc11

  implicit none

  integer, parameter :: dp = kind(0.0d0)

  integer :: nant, ntimes
  real(dp) :: refx, refy, refz, leapsec
  character(len=128) :: jpx_de421
  logical, dimension(:), allocatable :: ssobj
  real(dp), dimension(:), allocatable :: antx, anty, antz, temp, &
      pressure, humidity, mjd, ra, dec, dx, dy, dut, axisoff
  real(dp), dimension(:, :), allocatable :: geodelay, drydelay, &
      wetdelay
  character(len=8), dimension(:), allocatable :: sourcename

  nant = 1
  ntimes = 60

  allocate(antx(nant))
  allocate(anty(nant))
  allocate(antz(nant))
  allocate(temp(nant))
  allocate(pressure(nant))
  allocate(humidity(nant))
  allocate(axisoff(nant))

  allocate(mjd(ntimes))
  allocate(ra(ntimes))
  allocate(dec(ntimes))
  allocate(ssobj(ntimes))
  allocate(dx(ntimes))
  allocate(dy(ntimes))
  allocate(dut(ntimes))
  allocate(sourcename(ntimes))

  allocate(geodelay(ntimes, nant))
  allocate(drydelay(ntimes, nant))
  allocate(wetdelay(ntimes, nant))

  refx = 1000.0_dp
  refy = 1299.0_dp
  refz = 700.0_dp
  antx = -3101.52_dp
  anty = -11245.77_dp
  antz = 8916.26_dp

  temp = 30.0_dp
  pressure = 1000.0_dp
  humidity = 0.7_dp
  mjd = (/ 58701.45833333334_dp, 58701.45902777778_dp, 58701.45972222222_dp, 58701.46041666667_dp, 58701.46111111111_dp, &
      58701.46180555555_dp, 58701.46250000000_dp, 58701.46319444444_dp, 58701.46388888889_dp, 58701.46458333333_dp, &
      58701.46527777778_dp, 58701.46597222222_dp, 58701.46666666667_dp, 58701.46736111111_dp, 58701.46805555555_dp, &
      58701.46875000000_dp, 58701.46944444445_dp, 58701.47013888889_dp, 58701.47083333333_dp, 58701.47152777778_dp, &
      58701.47222222222_dp, 58701.47291666667_dp, 58701.47361111111_dp, 58701.47430555556_dp, 58701.47500000000_dp, &
      58701.47569444445_dp, 58701.47638888889_dp, 58701.47708333333_dp, 58701.47777777778_dp, 58701.47847222222_dp, &
      58701.47916666666_dp, 58701.47986111111_dp, 58701.48055555556_dp, 58701.48125000000_dp, 58701.48194444444_dp, &
      58701.48263888889_dp, 58701.48333333333_dp, 58701.48402777778_dp, 58701.48472222222_dp, 58701.48541666667_dp, &
      58701.48611111111_dp, 58701.48680555556_dp, 58701.48750000000_dp, 58701.48819444444_dp, 58701.48888888889_dp, &
      58701.48958333334_dp, 58701.49027777778_dp, 58701.49097222222_dp, 58701.49166666667_dp, 58701.49236111111_dp, &
      58701.49305555555_dp, 58701.49375000000_dp, 58701.49444444444_dp, 58701.49513888889_dp, 58701.49583333333_dp, &
      58701.49652777778_dp, 58701.49722222222_dp, 58701.49791666667_dp, 58701.49861111111_dp, 58701.49930555555_dp /)

  ra = 2.0_dp
  dec = 1.0_dp

  ssobj = .false.
  dx = 0.0_dp
  dy = 0.0_dp
  dut = 0.0_dp
  leapsec = 37.0_dp
  axisoff = 0.0_dp
  sourcename = 'P'
  jpx_de421 = 'data/DE421_little_Endian'

  call almacalc(refx, refy, refz, nant, antx, anty, antz, temp, &
      pressure, humidity, ntimes, mjd, ra, dec, ssobj, dx, dy, dut, &
      leapsec, axisoff, sourcename, JPX_DE421, geodelay, drydelay, &
      wetdelay)

  print '("geodelay:")'
  print *, geodelay
  print '("drydelay:")'
  print *, drydelay
  print '("wetdelay:")'
  print *, wetdelay

  deallocate(wetdelay)
  deallocate(drydelay)
  deallocate(geodelay)

  deallocate(dut)
  deallocate(dy)
  deallocate(dx)
  deallocate(ssobj)
  deallocate(dec)
  deallocate(ra)
  deallocate(mjd)
  
  deallocate(axisoff)
  deallocate(humidity)
  deallocate(pressure)
  deallocate(temp)
  deallocate(antz)
  deallocate(anty)
  deallocate(antx)

end program test_calc11

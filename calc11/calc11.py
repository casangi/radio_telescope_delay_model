import os
import ctypes as ct
import numpy as np
import numpy.ctypeslib as npc

# Define array pointer types.

a_f64_1 = npc.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS')
a_f64_2 = npc.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS')
a_bool_1 = npc.ndpointer(dtype=bool, ndim=1, flags='C_CONTIGUOUS')
a_char_1 = npc.ndpointer(dtype='S1', ndim=1, flags='C_CONTIGUOUS')
a_char_2 = npc.ndpointer(dtype='S1', ndim=2, flags='C_CONTIGUOUS')

# Import the wrapper library.
path = os.path.dirname(os.path.realpath(__file__)) + '/lib'
#os.system('otool -L ' + path + '/calc11_wrapper.so')
#print('******')
#os.system('install_name_tool -change libcalc11.so ' + path + '/libcalc11.so '+ path + '/calc11_wrapper.so')
#os.system('otool -L ' + path + '/calc11_wrapper.so')

#wrapper = npc.load_library('calc11_wrapper.so', path)
wrapper = npc.load_library('libcalc11.so', path)


# Define wrapper_almacalc result and argument types.

wrapper.wrapper_almacalc.restype = None
wrapper.wrapper_almacalc.argtypes = [ct.c_double, ct.c_double,
                                     ct.c_double, ct.c_int, a_f64_1, a_f64_1,
                                     a_f64_1, a_f64_1, a_f64_1, a_f64_1,
                                     ct.c_int, a_f64_1, a_f64_1, a_f64_1,
                                     a_bool_1, a_f64_1, a_f64_1, a_f64_1,
                                     ct.c_double, a_f64_1, a_char_2, a_char_1,
                                     a_f64_2, a_f64_2, a_f64_2]

def almacalc(refx, refy, refz, antx, anty, antz, temp, pressure,
             humidity, mjd, ra, dec, ssobj, dx, dy, dut, leapsec, axisoff,
             sourcename, jpx_de421):

    # Check arguments.

    nant = antx.shape[0]
    ntimes = mjd.shape[0]

    assert anty.shape[0] == nant
    assert antz.shape[0] == nant
    assert temp.shape[0] == nant
    assert pressure.shape[0] == nant
    assert humidity.shape[0] == nant
    assert axisoff.shape[0] == nant

    assert ra.shape[0] == ntimes
    assert dec.shape[0] == ntimes
    assert ssobj.shape[0] == ntimes
    assert dx.shape[0] == ntimes
    assert dy.shape[0] == ntimes
    assert dut.shape[0] == ntimes
    assert len(sourcename) == ntimes

    # Copy strings into arrays.

    nchar = 8
    sourcename_arr = np.zeros((ntimes, nchar), dtype='S1')

    for i in range(ntimes):
        for j, c in enumerate(sourcename[i]):
            if j == nchar:
                break
            sourcename_arr[i][j] = c

    nchar = 128
    jpx_de421_arr = np.zeros(nchar, dtype='S1')

    for i, c in enumerate(jpx_de421):
        if i == nchar:
            break
        jpx_de421_arr[i] = c

    # Create arrays for output - note dimensions are reversed from
    # Fortran code!

    geodelay = np.zeros((nant, ntimes), dtype=np.float64)
    drydelay = np.zeros((nant, ntimes), dtype=np.float64)
    wetdelay = np.zeros((nant, ntimes), dtype=np.float64)

    # Call the wrapper library function.

    wrapper.wrapper_almacalc(refx, refy, refz, nant, antx, anty, antz,
                             temp, pressure, humidity, ntimes, mjd, ra,
                             dec, ssobj, dx, dy, dut, leapsec, axisoff,
                             sourcename_arr, jpx_de421_arr, geodelay,
                             drydelay, wetdelay)

    return geodelay, drydelay, wetdelay

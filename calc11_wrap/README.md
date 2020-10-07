# Calc11

## How to obtain and compile the Calc11 library

The latest version of the Calc11 library is available in the ALMA
repository:

<https://bitbucket.sco.alma.cl/scm/asw/control>

The Calc11 code is in `thirdPartyPackages/AlmaCalc11`. It is written
in Fortran. The only subroutine that seems to be called from outside
is `almacalc`.

The source code has many `.i` files containing common blocks and
parameters which are incorporated into the code using include
statements. One called `params11.i` contains hard-coded paths to data
files. It is created from a template when building the library, but it
seems not to be used because the corresponding include statements are
commented out.

The build instructions in the `README` file do not work, but it is
possible to do it manually. You need gcc and gfortran to compile it.
The steps are:
```
$ git clone https://bitbucket.sco.alma.cl/scm/asw/control.git
$ cd control/thirdPartyPackages/AlmaCalc11/src
$ sed -i '' '21,25s/^\(.*\)$/#\1/' Makefile
$ make
```
(The third line comments out lines 21-25 of the `Makefile`, the ones
concerning standard ALMA definitions.)

To use the library, you need the `libcalc11.so` shared library file
and the planetary ephemeris data file in `../data/DE421_little_Endian`
(see below).

## Python wrapper and test programs

The test programs are:
* `test_calc11.f90`: Fortran code to call the `almacalc` subroutine
  directly
* `test_calc11.py`: Python code to call `almacalc` using the wrapper

To build the Fortran test program and the Python wrapper, first copy
the Calc11 shared library file into this directory and the planetary
ephemeris data file into the `data` subdirectory. Then do:
```
$ make
```

To run the test programs, do:
```
./test_calc11
./test_calc11.py
```

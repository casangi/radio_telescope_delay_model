# ALMA - Atacama Large Millimeter Array
# (c) Associated Universities Inc., 2019
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA

#
# Standard ALMA definitions like INSTALL_ROOT
#
#MAKEDIRTMP := $(shell searchFile include/acsMakefile)
#ifneq ($(MAKEDIRTMP),\#error\#)
#   MAKEDIR := $(MAKEDIRTMP)/include
#    include $(MAKEDIR)/acsMakefile
#endif

############################################################################
#
#  variables

#
# Fortran source files
#
SOURCES.f = almacalc.f cut1m.f cwobm.f ainit.f csitm.f adrvr.f \
	cctiu.f catiu.f cdtdb.f cpepu.f cnutm.f cnutu.f cnutu6.f cdiuu.f \
	cm20u.f crosu.f cetdm.f cptdm.f cocem.f hardisp.f dsitu.f astrm.f \
	catmm.f caxom.f cuvm.f ctheu.f cmatu.f cmabd.f cvecu.f dkill.f

SOURCES.f += almastub.f

#
# C source files
#
SOURCES.c = almaout.c

#
# A variable listing the object files (C and Fortran)
#
OBJECTS = $(SOURCES.f:%.f=%.o)
OBJECTS += $(SOURCES.c:%.c=%.o)

###########################################################################
#
#  targets
all: ./libcalc11.so

install: all
	cp libcalc11.so $(INSTALL_ROOT)/lib
	cp ../data/DE421_little_Endian $(INSTALL_ROOT)/config

clean:
	rm -f *.o
	rm -f ./libcalc11.so
	rm -f param11.i

#
# Rule to turn Fortran source to object files
#
.f.o: param11.i
	gfortran -ffree-form -ffree-line-length-none -fPIC -I. -g -O2 -c $<

.c.o:
	gcc  -fPIC -c $<

./libcalc11.so: param11.i $(OBJECTS)
	gfortran -shared -o $@ $(OBJECTS)

param11.i: param11.i.in
	sed 's|@prefix@|$(INSTALL_ROOT)|' param11.i.in > param11.i


calc11_PATH := $(MODPATH)

calc11_SOURCES = \
	adrvr \
	ainit \
	almacalc \
	almaout \
	astrm \
	catiu \
	catmm \
	caxom \
	cctiu \
	cdiuu \
	cdtdb \
	cetdm \
	cm20u \
	cmabd \
	cmatu \
	cnutm \
	cnutu \
	cnutu6 \
	cocem \
	cpepu \
	cptdm \
	crosu \
	csitm \
	ctheu \
	cut1m \
	cuvm \
	cvecu \
	cwobm \
	dkill \
	dsitu \
	hardisp

INSTALL_FILES = ../config/DE421_little_Endian ../lib/libcalc11.so

# Standard rules
$(MODRULE)all: $(MODPATH) $(MODDEP) $(MODPATH)/lib/libcalc11.so
	$(AT)echo " . . . $@ done"

$(MODRULE)clean: $(MODPATH) clean_$(MODDEP)
	$(AT)rm -f $(patsubst %, $(calc11_PATH)/object/%.o, $(calc11_SOURCES))
	$(AT)rm -f $(calc11_PATH)/lib/libcalc11.so
	$(AT)rm -f $(calc11_PATH)/src/param11.i
	$(AT)echo " . . . $@ done"

$(MODRULE)clean_dist: $(MODPATH) clean_dist_$(MODDEP)
	$(AT)echo " . . . $@ done"

$(MODRULE)install: $(MODPATH) install_$(MODDEP)
	$(AT)echo " . . . $@ done"

# Custom rules
$(MODPATH)/src/param11.i: $(MODPATH)/src/param11.i.in
	$(AT)sed 's|@prefix@|$(INSTDIR)|' $< > $@

$(MODPATH)/object/%.o: $(MODPATH)/src/%.f $(MODPATH)/src/param11.i | $(MODPATH)/object
	$(AT)gfortran -ffree-form -ffree-line-length-none -fPIC -I. -g -O2 -c -o $@ $<

$(MODPATH)/object/almaout.o: $(MODPATH)/src/almaout.c | $(MODPATH)/object
	$(AT)gcc -fPIC -c -o $@ $<

$(MODPATH)/lib/libcalc11.so: $(patsubst %, $(MODPATH)/object/%.o, $(calc11_SOURCES)) | $(MODPATH)/lib
	$(AT)gfortran -shared -o $@ $^

$(MODPATH)/config/DE421_little_Endian: $(MODPATH)/data/DE421_little_Endian | $(MODPATH)/config
	$(AT)cp $< $@

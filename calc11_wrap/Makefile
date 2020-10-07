FC = gfortran
CALC11_LINK = -L. -lcalc11

.PHONY: all clean

all: test_calc11 calc11_wrapper.so

test_calc11: test_calc11.f90
	$(FC) -o $@ $< $(CALC11_LINK)

calc11_wrapper.so: calc11_wrapper.f90
	$(FC) -fPIC -shared -o $@ $< $(CALC11_LINK)

clean:
	$(RM) test_calc11 calc11_wrapper.so

.PHONY: all

SOURCES=$(shell ls *.f90)
MODULES=$(shell ls module_*)
OBJECTS=$(SOURCES:.f90=.o)
MODS=$(MODULES:.f90=.mod)
COPTS=-O3
all: bryan.x

bryan.x: $(MODS) $(OBJECTS) bryan.f90
	gfortran -g $(COPTS) -o bryan.x  $(OBJECTS) -L/home/titan/rossi/library_links/ -llapack -lfftw3 
module_path.mod:   module_prng.mod

%.o %.mod: %.f90
	gfortran -g  $(COPTS) -c $<

clean:
	rm -f *.o *.mod *.x

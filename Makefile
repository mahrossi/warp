.PHONY: all

SOURCES=$(shell ls *.f90)
MODULES=$(shell ls module_*.f90)
OBJECTS=$(SOURCES:.f90=.o)
MODS=$(MODULES:.f90=.mod)
COPTS=-O3
all: moments.x

moments.x: $(MODS) $(OBJECTS) moments.f90
	gfortran -g $(COPTS) -o moments.x $(OBJECTS)  -llapack

%.o %.mod: %.f90
	gfortran -g $(COPTS) -c $<

clean:
	rm -f *.o *.mod *.x

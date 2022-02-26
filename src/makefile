

SRCDIR = modules
BUILDDIR = build


fc = gfortran
OPT= -O3

program.exe: build main.o boundary.o forces.o initialization.o integrators.o statistics.o thermostats.o tools.o
	$(fc) -o program.exe $(OPT) $(BUILDDIR)/boundary.o $(BUILDDIR)/forces.o $(BUILDDIR)/initialization.o  $(BUILDDIR)/integrators.o $(BUILDDIR)/statistics.o $(BUILDDIR)/thermostats.o  $(BUILDDIR)/tools.o   $(BUILDDIR)/main.o -J$(BUILDDIR)

#Boundary module	
boundary.mod: boundary.o $(SRCDIR)/boundary.f90
	$(fc) -o $(BUILDDIR)/boundary.o -c $(OPT) $(SRCDIR)/boundary.f90 -J$(BUILDDIR)

boundary.o: $(SRCDIR)/boundary.f90
	$(fc) -o $(BUILDDIR)/boundary.o -c $(OPT) $(SRCDIR)/boundary.f90 -J$(BUILDDIR)

#Forces module	
forces.mod: forces.o $(SRCDIR)/forces.f90
	$(fc) -o $(BUILDDIR)/forces.o -c $(OPT) $(SRCDIR)/forces.f90 -J$(BUILDDIR)

forces.o: $(SRCDIR)/forces.f90
	$(fc) -o $(BUILDDIR)/forces.o -c $(OPT) $(SRCDIR)/forces.f90 -J$(BUILDDIR)

#Initialization module
initialization.mod: initialization.o $(SRCDIR)/initialization.f90
	$(fc) -o $(BUILDDIR)/initialization.o -c $(OPT) $(SRCDIR)/initialization.f90 -J$(BUILDDIR)

initialization.o: $(SRCDIR)/initialization.f90
	$(fc) -o $(BUILDDIR)/initialization.o -c $(OPT) $(SRCDIR)/initialization.f90 -J$(BUILDDIR)
	
#Integrators module
integrators.mod: integrators.o $(SRCDIR)/integrators.f90
	$(fc) -o $(BUILDDIR)/integrators.o -c $(OPT) $(SRCDIR)/integrators.f90 -J$(BUILDDIR)

integrators.o: $(SRCDIR)/integrators.f90
	$(fc) -o $(BUILDDIR)/integrators.o -c $(OPT) $(SRCDIR)/integrators.f90 -J$(BUILDDIR)
	
#Statistics module
statistics.mod: statistics.o $(SRCDIR)/statistics.f90
	$(fc) -o $(BUILDDIR)/statistics.o -c $(OPT) $(SRCDIR)/statistics.f90 -J$(BUILDDIR)

statistics.o: $(SRCDIR)/statistics.f90
	$(fc) -o $(BUILDDIR)/statistics.o -c $(OPT) $(SRCDIR)/statistics.f90 -J$(BUILDDIR)

#Thermostat module
thermostats.mod: thermostats.o $(SRCDIR)/thermostats.f90
	$(fc) -o $(BUILDDIR)/thermostats.o -c $(OPT) $(SRCDIR)/thermostats.f90 -J$(BUILDDIR)

thermostats.o: $(SRCDIR)/thermostats.f90
	$(fc) -o $(BUILDDIR)/thermostats.o -c $(OPT) $(SRCDIR)/thermostats.f90 -J$(BUILDDIR)

#Tools module
tools.mod: tools.o $(SRCDIR)/tools.f90
	$(fc) -o $(BUILDDIR)/tools.o -c $(OPT) $(SRCDIR)/tools.f90 -J$(BUILDDIR)

tools.o: $(SRCDIR)/tools.f90
	$(fc) -o $(BUILDDIR)/tools.o -c $(OPT) $(SRCDIR)/tools.f90 -J$(BUILDDIR)


#Main file
main.o: boundary.mod forces.mod initialization.mod integrators.mod statistics.mod thermostats.mod tools.mod main.f90
	$(fc) -o $(BUILDDIR)/main.o -c $(OPT) main.f90 -J$(BUILDDIR)


build:
	mkdir $(BUILDDIR)

run:
	./program.exe

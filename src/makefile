

SRCDIR = modules
BUILDDIR = build


fc = gfortran
OPT= -O1

program.exe: build main.o boundary.o forces.o initialization.o integrators.o thermostats.o
	$(fc) -o program.exe $(OPT) $(BUILDDIR)/boundary.o $(BUILDDIR)/forces.o $(BUILDDIR)/initialization.o  $(BUILDDIR)/integrators.o $(BUILDDIR)/thermostats.o  $(BUILDDIR)/main.o -J$(BUILDDIR)

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

#Thermostat module
thermostats.mod: thermostats.o $(SRCDIR)/thermostats.f90
	$(fc) -o $(BUILDDIR)/thermostats.o -c $(OPT) $(SRCDIR)/thermostats.f90 -J$(BUILDDIR)

thermostats.o: $(SRCDIR)/thermostats.f90
	$(fc) -o $(BUILDDIR)/thermostats.o -c $(OPT) $(SRCDIR)/thermostats.f90 -J$(BUILDDIR)



#Main file
main.o: initialization.mod boundary.mod thermostats.mod forces.mod integrators.mod main.f90
	$(fc) -o $(BUILDDIR)/main.o -c $(OPT) main.f90 -J$(BUILDDIR)


build:
	mkdir $(BUILDDIR)

run:
	./program.exe

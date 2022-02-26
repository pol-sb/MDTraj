

SRCDIR = modules
BUILDDIR = build


fc = gfortran
OPT= -O3

program.exe: build main.o function_module.o 
	$(fc) -o program.exe $(OPT) $(BUILDDIR)/function_module.o  $(BUILDDIR)/main.o -J$(BUILDDIR)

	
function_module.mod: function_module.o $(SRCDIR)/function_module.f90
	$(fc) -o $(BUILDDIR)/function_module.o -c $(OPT) $(SRCDIR)/function_module.f90 -J$(BUILDDIR)

function_module.o: $(SRCDIR)/function_module.f90
	$(fc) -o $(BUILDDIR)/function_module.o -c $(OPT) $(SRCDIR)/function_module.f90 -J$(BUILDDIR)


main.o: function_module.mod main.f90
	$(fc) -o $(BUILDDIR)/main.o -c $(OPT) main.f90 -J$(BUILDDIR)


build:
	mkdir $(BUILDDIR)

clean:
	rm -f $(BUILDDIR)/main.o $(BUILDDIR)/function_module.o $(BUILDDIR)/function_module.mod $(BUILDDIR)/statistic_module.o $(BUILDDIR)/statistic_module.mod program.exe

run:
	./program.exe

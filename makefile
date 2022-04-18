PY_DIST = $(shell which python3)
MAKEFLAGS += --no-print-directory

## By default, it is all (just write make)

## 1. make all - to compile, execute and plot.
all:
	$(MAKE) -C src/
	@make run
	@$(MAKE) plot

## 2. make compile - to compile.
compile:
	$(MAKE) -C src/

## 3. make run - to run.
run:
	 mpirun -np 8 ./src/program.exe

## 4. make plot - to plot the graphics
plot:
	@echo "Python found in: ${PY_DIST}. Plotting now."
	@${PY_DIST} ./src/modules/plot_statistics.py

## 5. make clean - to remove build and program.exe
clean:
	rm -r src/build
	rm -r src/program.exe


## [!] WARNING: If you change the parameters, you need to recompile.

.PHONY : help
help:
	@sed -n 's/^##//p' makefile

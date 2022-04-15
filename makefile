PY_DIST = $(shell which python3)
MAKEFLAGS += --no-print-directory
MPI_NPROC = 4
MPI_minP=4
MPI_maxP=8
MPI_stepP=2

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
	 mpirun -np $(MPI_NPROC)  ./src/program.exe

## 4. make plot - to plot the graphics
plot:
	@echo "Python found in: ${PY_DIST}. Plotting now."
	@${PY_DIST} ./src/modules/plot_statistics.py

## 5. make clean - to remove build and program.exe
clean:
	rm -r src/build
	rm -r src/program.exe

## 6. check paralel#!/usr/bin/bash
## send as nohup make check_paralel &
check_paralel:
	./src/modules/check_paralel.sh -m $($MPI_maxP) -d $($MPI_minP) -s $($MPI_stepP)
## [!] WARNING: If you change the parameters, you need to recompile. 

.PHONY : help
help:
	@sed -n 's/^##//p' makefile

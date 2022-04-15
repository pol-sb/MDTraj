PY_DIST = $(shell which python3)
MAKEFLAGS += --no-print-directory
MPI_NPROC = 4
MPI_minP=4
MPI_maxP=8
MPI_stepP=2
## By default, it is all (just write make)

## 1. make all - to compile, execute and plot.
all:
	@$(MAKE) -C src/
	@make run
	@$(MAKE) plot

## 2. make compile - to compile.
compile:
	@$(MAKE) -C src/
	
## 3. make run - to run.
run:
	@echo "Running the simulation with $(MPI_NPROC) processors."
	mpirun -np $(MPI_NPROC) ./src/program.exe

## 4. make plot - to plot the graphics
plot:
	@echo "Python found in: ${PY_DIST}. Plotting now."
	@${PY_DIST} ./src/modules/plot_statistics.py

## 5. make clean - to remove build and program.exe
clean:
	rm -r src/build
	rm -r src/program.exe

## 6. check paralel#!/usr/bin/bash
while getopts "m:d:s:" flag
do
        case ${flag} in
                m) maxP=${OPTARG};;
                d) minP=${OPTARG};;
                s) step=${OPTARG};;
        esac
done

SOURCE=${BASH_SOURCE[0]}
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR=$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )
  SOURCE=$(readlink "$SOURCE")
  [[ $SOURCE != /* ]] && SOURCE=$DIR/$SOURCE # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR=$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )

#echo $(($(($maxP))*$(($real_step))))

step=$(($step))
maxP=$(($maxP))
minP=$(($minP))
# Data fixing
if [ $step -le 0 ]
then
        echo "Wrong step size in -s option, $step is no a good value"
        exit 1
fi

if [ $maxP -lt $step ]
then
        echo "Wrong step size in -m option, $maxP is not a good value"
        exit 1
fi
touch ./output/performance.dat
echo "" > ./output/performance.dat
# Change in elongation of the simulation
sed -i "s/integer::ntimes=.*\!/integer::ntimes=1000  !/" ./input/parameter.h

for Length in 5 10 22 ;
do

for  (( c=$(($minP)); c<=$(($maxP)); c+=$(($step))));
do

        sed -i "s/integer::nc=.*\!/integer::nc="$Length"  \!/" ./input/parameter.h
        echo "running with $(($Length**3)) particles and $c process"
        grep "integer::nc" ./input/parameter.h
        grep "integer::ntimes" ./input/parameter.h
        touch ./src/main.f90
        make -C ./src/
        mpirun -np $c ./src/program.exe &
        echo "PID $!" 
        wait
#       cat ./output/performance.dat
done
done
echo "----------"
sed -i "s/integer::ntimes=1000/integer::ntimes=200000/" ./input/parameter.h
sed -i "s/integer::nc=./integer::nc=5/" ./input/parameter.h

## send as nohup make check_paralel &
check_paralel:
	./src/modules/check_paralel.sh -m $($MPI_maxP) -d $($MPI_minP) -s $($MPI_stepP)
## [!] WARNING: If you change the parameters, you need to recompile. 

.PHONY : help
help:
	@sed -n 's/^##//p' makefile

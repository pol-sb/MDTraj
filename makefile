PY_DIST != which python3

## By default, it is all (just write make)

## 1. make all - to compile, execute and plot.
all:
	$(MAKE) -C src/
	./src/program.exe
	${PY_DIST} ./src/modules/plot_statistics.py

## 2. make compile - to compile.
compile:
	$(MAKE) -C src/
	
## 3. make run - to run.
run:
	 ./src/program.exe

## 4. make plot - to plot the graphics
plot:
	@echo "Python found in: ${PY_DIST} "
	@${PY_DIST} ./src/modules/plot_statistics.py

## 5. make clean - to remove build and program.exe
clean:
	rm -r src/build
	rm -r src/program.exe


## ADVERTISING: If you change the parameters, you need to recompile. 

.PHONY : help
help:
	@sed -n 's/^##//p' makefile

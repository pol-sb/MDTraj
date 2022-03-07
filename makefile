PY_DIST != which python3

all:src/program.exe
	$(MAKE) -C src/
	./src/program.exe
	${PY_DIST} ./src/modules/plot_statistics.py

compile:
	$(MAKE) -C src/
	
run:
	 ./src/program.exe

src/program.exe:
	$(MAKE) -C src/

clear:
	rm -r src/build
	rm -r src/program.exe

plot:
	@echo "Python found in: ${PY_DIST} "
	@${PY_DIST} ./src/modules/plot_statistics.py
	


all:src/program.exe
	$(MAKE) -C src/
	./src/program.exe

compile:
	$(MAKE) -C src/
	
run:
	 ./src/program.exe

src/program.exe:
	$(MAKE) -C src/

clear:
	rm -r src/build
	rm -r src/program.exe
	

SRC_DIR := ./modules
BUILD_DIR := ./build
# Si no tuvienes dependencias internas los modulos prodriamos usar esta linea para tener todos los modulos
#SRC_FILES := $(wildcard $(SRC_DIR)/* *.f90)

SRC_FILES := $(SRC_DIR)/mtfort.f90 $(SRC_DIR)/observable_measures.f90 $(SRC_DIR)/monte_carlo.f90 $(SRC_DIR)/input.f90 $(SRC_DIR)/initialize_geometry.f90

OBJ_FILES := $(patsubst $(SRC_DIR)/%.f90,$(BUILD_DIR)/%.o,$(SRC_FILES))
FC := mpif90 -o
FFLAGS := -O1

Sources = main.f90 input.f90 initialize_geometry.f90 monte_carlo.f90 mtfort.f90 input.f90

MC.exe: build $(OBJ_FILES)
	@$(FC) $(FFLAGS) -o $@ main.f90 $(OBJ_FILES)
	@echo "Compiled"
	@mv *.mod $(BUILD_DIR)



#%.mod: $(SRC_DIR)/%.f90
#	$(FC) $(FFLAGS) -c -o $(SRC)/$@  $^
#
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90 build
	@$(FC) $(FFLAGS) -c -o 	$@ $<
# si no hubiese interdepencias tambien estaría bien
#rm *.mod


build:
	@mkdir $(BUILD_DIR)

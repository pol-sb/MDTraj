SRC_DIR := ./modules
BUILD_DIR := ./build

FC := mpifort
FFLAGS := -O1

SRC_FILES := $(SRC_DIR)/boundary.f90  $(SRC_DIR)/thermostats.f90   $(SRC_DIR)/initialization.f90   $(SRC_DIR)/forces.f90   $(SRC_DIR)/integrators.f90

OBJ_FILES := $(patsubst $(SRC_DIR)/%.f90,$(BUILD_DIR)/%.o,$(SRC_FILES))

program.exe: build $(OBJ_FILES)
	@echo "Starting compilation..."
	@$(FC) $(FFLAGS) -o $@ main.f90 $(OBJ_FILES)
	@echo "Program compiled successfully!"
	@mv *.mod $(BUILD_DIR)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90 build
	@$(FC) $(FFLAGS) -c -o 	$@ $<

build:
	@mkdir $(BUILD_DIR)
	@echo "'$(BUILD_DIR)' directory created."

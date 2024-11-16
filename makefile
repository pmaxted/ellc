FC = gfortran
FFLAGS = -O2 -fPIC
LDFLAGS = -shared

SRC_DIR = f_src
SRC_FILES = $(SRC_DIR)/constants.f90 $(SRC_DIR)/utils.f90 $(SRC_DIR)/coords.f90 $(SRC_DIR)/gauss_legendre.f90 $(SRC_DIR)/solve_real_poly.f90 $(SRC_DIR)/ellipse.f90 $(SRC_DIR)/stellar.f90 $(SRC_DIR)/spots.f90 $(SRC_DIR)/ellc.f90
TARGET = libellc.so
INSTALL_DIR = ./ellc

$(TARGET): $(SRC_FILES)
	$(FC) $(FFLAGS) $(LDFLAGS) $(SRC_FILES) -o $(TARGET)  

install:
	cp $(TARGET) $(INSTALL_DIR)/

clean:
	rm -f *.mod $(TARGET)
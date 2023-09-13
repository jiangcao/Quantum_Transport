# Fortran 90 compiler
FC90 = ifort
#F90_FLAGS =  -r8 -check bounds -traceback -fpp 
MACFLAGS = -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib
F90_FLAGS = -i8 -r8 -O2 -fpp -mkl -traceback  -qopenmp -qopt-matmul -check bounds  
LIBS = -lpthread -lm -ldl

# Modules directory
MODDIR = compiled

# Source directory
SRCDIR = src/

# Search directories
vpath %.f90 $(SRCDIR)
vpath %.o $(MODDIR)

# Targets.

all: qt.x 

qt.x : main.f90 mkl_dfti.o fft_mod.o static.o linalg.o matrix_c.o graph_partition.o wannierHam3d.o deviceHam_mod.o rgf_mod.o negf_mod.o

	$(FC90) -o $@ $< $(MODDIR)/*.o $(F90_FLAGS) $(LIBS) -module $(MODDIR) $(MACFLAGS)

.PHONY : clean;

clean :
	rm -f $(MODDIR)/*.o $(MODDIR)/*.mod

# implicit rules

%.o : %.f90
	$(FC90) -o $(MODDIR)/$@ $< -c $(F90_FLAGS) -module $(MODDIR)

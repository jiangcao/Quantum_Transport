# Fortran 90 compiler
MKLROOT = /usr/pack/intel_compiler-2020-af/x64/compilers_and_libraries_2020.0.166/linux/mkl/
FC90 = /usr/pack/mpich-3.2.1-af/linux-x64/bin/mpif90
MACFLAGS = -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib
F90_FLAGS = -Wall -Wextra -O0 -march=native -ffast-math -ffree-line-length-none -fopenmp -fbacktrace  -Wno-unused-variable -Wno-unused-parameter -Wno-unused-function
LIBS = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl 

# Modules directory
MODDIR = compiled

# Source directory
SRCDIR = src/

# Search directories
vpath %.f90 $(SRCDIR)
vpath %.o $(MODDIR)

# Targets.

all: qt.x 

qt.x : main.f90  static.o  mkl_dfti.o linalg.o matrix_c.o graph_partition.o wannierHam3d.o deviceHam_mod.o rgf_mod.o negf_mod.o

	$(FC90) -o $@ $< $(MODDIR)/*.o $(F90_FLAGS) $(LIBS) -I$(MODDIR) -J$(MODDIR)

.PHONY : clean;

clean :
	rm -f $(MODDIR)/*.o $(MODDIR)/*.mod

# implicit rules

%.o : %.f90
	$(FC90) -o $(MODDIR)/$@ $< -c $(F90_FLAGS) -I$(MODDIR) -J$(MODDIR) 

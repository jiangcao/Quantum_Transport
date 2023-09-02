# Fortran 90 compiler
MKLROOT = /usr/pack/intel_compiler-2020-af/x64/compilers_and_libraries_2019.0.117/linux/mkl/
FC90 = /usr/sepp/bin/ifort-2020-af
#F90_FLAGS =  -r8 -check bounds -traceback -fpp 
MACFLAGS = -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib
F90_FLAGS = -i8  -I"${MKLROOT}/include" -r8 -O2 -fpp -mkl -traceback  -qopenmp -qopt-matmul -check bounds  
LIBS = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl

# Modules directory
MODDIR = compiled

# Source directory
SRCDIR = src/

# Search directories
vpath %.f90 $(SRCDIR)
vpath %.o $(MODDIR)

# Targets.

all: qt.x 

qt.x : main.f90 static.o matrix_mod.o graph_partition.o deviceHam_mod.o negf_mod.o

	$(FC90) -o $@ $< $(MODDIR)/*.o $(F90_FLAGS) $(LIBS) -module $(MODDIR)

.PHONY : clean;

clean :
	rm -f $(MODDIR)/*.o $(MODDIR)/*.mod

# implicit rules

%.o : %.f90
	$(FC90) -o $(MODDIR)/$@ $< -c $(F90_FLAGS) -module $(MODDIR)

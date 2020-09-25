EXEC = int2-surflap
#HOST = gcc
HOST = gcc-openmp
#HOST = intel
#HOST = intel-ompenmp

#
# For linux systems, it is assumed that the environment
# variable LD_LIBRARY_PATH contains the locations to libfmm3d.so
# and libsolvers3d.so, for Macosx, these .so files also need to be
# copied over /usr/local/lib


ifneq ($(OS),Windows_NT) 
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Darwin)
        LDF = /usr/local/lib
    endif
    ifeq ($(UNAME_S),Linux)
        LDF = ${HOME}/lib
    endif
endif
    

LIBS = -lfmm3dbie -lfmm3d  -lopenblas ${LDFLAGS}
ifeq ($(HOST),gcc)
    FC=gfortran -L${LDF} 
    FFLAGS=-fPIC -O3 -funroll-loops -march=native -std=legacy 
endif

ifeq ($(HOST),gcc-openmp)
    FC = gfortran -L${LDF} 
    FFLAGS=-fPIC -O3 -funroll-loops -march=native -fopenmp -std=legacy
endif

ifeq ($(HOST),intel)
    FC=ifort -L${LDF} 
    FFLAGS= -O3 -fPIC -march=native
endif

ifeq ($(HOST),intel-openmp)
    FC = ifort -L${LDF} 
    FFLAGS= -O3 -fPIC -march=native -qopenmp
endif

SURF=../src/surface_routs
LBW=../src/lapbel_wrap
COM=../src/common

.PHONY: all clean 

OBJECTS =  test_lapbel_hodge_hgfast_mr_exp.o \
    $(SURF)/vtk_routs.o \
    $(SURF)/surf_routs.o \
    $(LBW)/lapbel2fast.o \
    $(LBW)/lapbel_kernels.o \
    $(COM)/legetens.o


#
# use only the file part of the filename, then manually specify
# the build location
#

%.o : %.f  
	$(FC) -c $(FFLAGS) $< -o $@

%.o : %.f90  
	$(FC) -c $(FFLAGS) $< -o $@

all: $(OBJECTS)
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJECTS) -L${LDF} $(LIBS)
	./$(EXEC)  

clean:
	rm -f $(OBJECTS)
	rm -f $(EXEC)

list: $(SOURCES)
	$(warning Requires:  $^)




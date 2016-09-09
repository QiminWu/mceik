OBJ = ./Obj

CC = gcc-5
MPICC = /home/bakerb3/Fortran/openmpi-1.8.4/bin/mpicc
MPIF90 = /home/bakerb3/Fortran/openmpi-1.8.4/bin/mpif90
F90 = /home/bakerb3/Fortran/openmpi-1.8.4/bin/mpif90
CFLAGS = -Wall -O2 -fopenmp -g3
FFLAGS = -Wall -O2 -fopenmp -g3 -fbounds-check
MODDIR = ./modules
MODULE = -J$(MODDIR)
INCF =-I$(MODDIR)

INC_LOC = -I./include

Z_LIB = /home/bakerb3/C/zlib-1.2.8/libz.a
H5_INC = -I/home/bakerb3/Fortran/hdf5-1.8.16/include
H5_LIB = /home/bakerb3/Fortran/hdf5-1.8.16/lib/libhdf5.a $(Z_LIB)
LLAPACK = /usr/lib/liblapack.a /usr/lib/libcblas.so /usr/lib/libblas.a
INI_INC = -I/home/bakerb3/C/iniparser/src
INI_LIB = /home/bakerb3/C/iniparser/libiniparser.a

MKLROOT = /opt/intel/mkl
MKL_INC = -I$(MKLROOT)/include
MKL_LIB = -Wl,--start-group \
	  ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a \
	  ${MKLROOT}/lib/intel64/libmkl_core.a \
	  ${MKLROOT}/lib/intel64/libmkl_sequential.a \
	 -Wl,--end-group -lpthread -lm -ldl

INC_ALL = $(INC_LOC) $(INI_INC) $(MKL_INC)

LIBALL = $(H5_LIB) $(LLAPACK) $(INI_LIB) $(MKL_LIB) -lm -ldl


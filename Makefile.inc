OBJ = ./Obj

CC = icc
MPICC = /opt/intel/impi/2017.0.098/bin64/mpicc
MPIF90 = /opt/intel/impi/2017.0.098/bin64/mpif90
CC = $(MPICC)
F90 = $(MPIF90)
#CC = gcc
#MPICC = /usr/local/bin/mpicc
#MPIF90 = /usr/local/bin/mpif90
#F90 = gfortran
CFLAGS = -Wall -O2 -fopenmp -g3
FFLAGS = -Wall -O2 -fopenmp -g3 -fbounds-check
MODDIR = ./modules
MODULE = -J$(MODDIR)

INC_LOC = -I./include

Z_LIB = /home/bakerb3/C/zlib-1.2.8/libz.a
#H5_INC = -I/home/bakerb3/C/hdf5_mpi/include
#H5_LIB = /home/bakerb3/C/hdf5_mpi/lib/libhdf5.a $(Z_LIB)
H5_INC = -I/home/bakerb3/C/hdf5-1.10.0-patch1_intel/include
H5_LIB = /home/bakerb3/C/hdf5-1.10.0-patch1_intel/lib/libhdf5.a $(Z_LIB)
LLAPACK = /usr/lib/liblapack.a /usr/lib/libcblas.so /usr/lib/libblas.a
INI_INC = -I/home/bakerb3/C/iniparser/src
INI_LIB = /home/bakerb3/C/iniparser/libiniparser.a

MKLROOT = /opt/intel/compilers_and_libraries_2017.0.098/linux/mkl
MKL_INC = -I$(MKLROOT)/include
MKL_LIB = -Wl,--start-group \
	  ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a \
	  ${MKLROOT}/lib/intel64/libmkl_core.a \
	  ${MKLROOT}/lib/intel64/libmkl_sequential.a \
	 -Wl,--end-group -lpthread -lm -ldl

INCF =-I$(MODDIR) $(H5_INC)
INC_ALL = $(INC_LOC) $(INI_INC) $(H5_INC) $(MKL_INC)

LIBALL = $(H5_LIB) $(INI_LIB) $(MKL_LIB) -lm -ldl


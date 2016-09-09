include Makefile.inc

ifeq "$(wildcard $(OBJ) )" ""
-include $(shell mkdir $(OBJ)) $(wildcard $(OBJ)/*)
endif
ifeq "$(wildcard $(MODDIR) )" ""
-include $(shell mkdir $(MODDIR)) $(wildcard $(MODDIR)/*)
endif

XFSM2D = xfsm2d
XFSM3D = xfsm3d
XLOC = xloc
XHOMOG = xhomog

EXECS = $(XFSM3D) $(XLOC) $(XHOMOG)

#OBJ_COMMON = $(OBJ)/eikonal_utils.o
OBJ_REQD = $(OBJ)/module.o
OBJ_COMMON = $(OBJ)/mpiutils.o
OBJHOMOG = $(OBJ_REQD) $(OBJ_COMMON) $(OBJ)/homog.o
OBJLOC = $(OBJ_REQD) $(OBJ_COMMON) $(OBJ)/locate.o
OBJ2D = $(OBJ_REQD) $(OBJ_COMMON) $(OBJ)/sweep2d.o
OBJ3D = $(OBJ_REQD) $(OBJ_COMMON) $(OBJ)/fsm3d.o

all: $(EXECS)

$(XHOMOG): $(OBJHOMOG)
	$(MPIF90) $(CFLAGS) -o $(XHOMOG) $(OBJHOMOG) $(LIBALL)

$(XFSM2D): $(OBJ2D)
	$(MPICC) $(CFLAGS) -o $(XFSM2D) $(OBJ2D) $(LIBALL)

$(XFSM3D): $(OBJ3D)
	$(MPIF90) $(FFLAGS) -o $(XFSM3D) $(OBJ3D) $(LIBALL)

$(XLOC): $(OBJLOC)
	$(MPIF90) $(FFLAGS) -o $(XLOC) $(OBJLOC) $(LIBALL)

$(OBJ)/fsm3d.o: fsm3d.f90
	$(MPIF90) $(FFLAGS) $(INCF) -c fsm3d.f90 -o $(OBJ)/fsm3d.o

$(OBJ)/homog.o: homog.c
	$(MPICC) $(CFLAGS) $(INC_ALL) -c homog.c -o $(OBJ)/homog.o

$(OBJ)/locate.o: locate.f90
	$(MPIF90) $(FFLAGS) $(INCF) -c locate.f90 -o $(OBJ)/locate.o

$(OBJ)/module.o: module.F90
	$(F90) $(FFLAGS) $(INCF) -c module.F90 -o $(OBJ)/module.o
 
$(OBJ)/mpiutils.o: mpiutils.f90
	$(MPIF90) $(FFLAGS) $(INCF) -c mpiutils.f90 -o $(OBJ)/mpiutils.o

$(OBJ)/sweep2d.o: sweep2d.c
	$(MPICC) $(CFLAGS) $(INC_ALL) -c sweep2d.c -o $(OBJ)/sweep2d.o

$(OBJ)/sweep3d.o: sweep3d.c
	$(MPICC) $(CFLAGS) $(INC_ALL) -c sweep3d.c -o $(OBJ)/sweep3d.o

clean:
	@$(RM) $(OBJ)/*.o $(EXECS) $(MODDIR)/*.mod

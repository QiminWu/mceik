      MODULE CONSTANTS_MODULE
         DOUBLE PRECISION, PARAMETER :: zero  = 0.d0
         DOUBLE PRECISION, PARAMETER :: one   = 1.d0
         DOUBLE PRECISION, PARAMETER :: two   = 2.d0
         DOUBLE PRECISION, PARAMETER :: three = 3.d0
         DOUBLE PRECISION, PARAMETER :: four  = 4.d0
         DOUBLE PRECISION, PARAMETER :: third = one/three
         DOUBLE PRECISION, PARAMETER :: two_third = two/three
         DOUBLE PRECISION, PARAMETER :: half = one/two
         DOUBLE PRECISION, PARAMETER :: sqrt2i = one/SQRT(two)
      END MODULE CONSTANTS_MODULE

      MODULE MPIUTILS_MODULE
         INTERFACE
            SUBROUTINE MPIUTILS_SPLITCOMM3D(comm,                             &
                                            lreord, iwt, ndivx, ndivy, ndivz, &
                                            global_comm, intra_table_comm,    &   
                                            inter_table_comm, ierr)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: comm, iwt, ndivx, ndivy, ndivz
            LOGICAL, INTENT(IN) :: lreord
            INTEGER, INTENT(OUT) :: global_comm, intra_table_comm, &
                                    inter_table_comm, ierr
            END SUBROUTINE MPIUTILS_SPLITCOMM3D 

            SUBROUTINE MPIUTILS_SPLITCOMM2D(comm,                            &
                                            lreord, iwt, ndivx, ndivz,       &
                                            global_comm, intra_table_comm,   &   
                                            inter_table_comm, ierr)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: comm, iwt, ndivx, ndivz
            LOGICAL, INTENT(IN) :: lreord
            INTEGER, INTENT(OUT) :: global_comm, intra_table_comm, &
                                    inter_table_comm, ierr
            END SUBROUTINE MPIUTILS_SPLITCOMM2D

            SUBROUTINE MPIUTILS_INITIALIZE3D(comm, ireord, iwt,         &
                                             ndivx, ndivy, ndivz, ierr) &
                       BIND(C, NAME='mpiutils_initialize3d')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN) :: comm, ireord, iwt, ndivx, ndivy, ndivz
            INTEGER(C_INT), INTENT(OUT) :: ierr
            END SUBROUTINE MPIUTILS_INITIALIZE3D

            SUBROUTINE MPIUTILS_INITIALIZE2D(comm, ireord, iwt, ndivx, ndivz, ierr) &
                       BIND(C, NAME='mpiutils_initialize2d')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN) :: comm, ireord, iwt, ndivx, ndivz
            INTEGER(C_INT), INTENT(OUT) :: ierr
            END SUBROUTINE MPIUTILS_INITIALIZE2D

            SUBROUTINE MPIUTILS_GRD2IJK(igrd, nx, ny, nz, i, j, k, ierr) &
                       BIND(C, NAME='mpiutils_grd2ijk')
            USE ISO_C_BINDING
            IMPLICIT NONE 
            INTEGER(C_INT), INTENT(IN) :: igrd, nx, ny, nz
            INTEGER(C_INT), INTENT(OUT) :: i, j, k, ierr
            END SUBROUTINE MPIUTILS_GRD2IJK

         END INTERFACE
         INTEGER, SAVE :: global_comm      ! overall communication
         INTEGER, SAVE :: intra_table_comm ! communicate within the traveltime table
         INTEGER, SAVE :: inter_table_comm ! communicate between the traveltime tables
      END MODULE MPIUTILS_MODULE

      MODULE EIKONAL3D_TYPES
         TYPE levelType
            INTEGER, ALLOCATABLE :: ixyz_level(:) !> (ix,iy,iz) indices at each node
                                                  !> in the level 
            INTEGER, ALLOCATABLE :: level_ptr(:)  !> start index of ixyz_level for
                                                  !> i'th level 
            INTEGER, ALLOCATABLE :: nnl(:)        !> number of nodes in level
            INTEGER nlevels                       !> number of levels
         END TYPE levelType

         TYPE localModelType
            TYPE (levelType) lstruct                 !> local level structure
            DOUBLE PRECISION, ALLOCATABLE :: slow(:) !> local slowness model (s/m)
            DOUBLE PRECISION, ALLOCATABLE :: u(:)    !> local traveltimes (s)
            INTEGER, ALLOCATABLE :: l2g_node(:)      !> maps from the local to the global
                                                     !> node
            LOGICAL, ALLOCATABLE :: lghost_node(:)   !> if true then this is a ghost
                                                     !> node that i must get from a
                                                     !> neighbor
            LOGICAL, ALLOCATABLE :: lupd(:)          !> if true then i must update this
                                                     !> node
            INTEGER nx   !> number of local x locations
            INTEGER ny   !> number of local y locations
            INTEGER nz   !> number of local z locations
            INTEGER nxyz !> number of local grid points
            INTEGER ix1  !> first global x node index of my local model
            INTEGER iy1  !> first global y node index of my local model
            INTEGER iz1  !> first global z node index of my local model
            INTEGER myblock  !> my block number
         END TYPE localModelType

         TYPE solverParametersType
            DOUBLE PRECISION x0   !> x0 of global model (m)
            DOUBLE PRECISION y0   !> y0 of global model (m)
            DOUBLE PRECISION z0   !> z0 of global model (m)
            DOUBLE PRECISION h    !> grid spacing of model (m)
            DOUBLE PRECISION tol  !> tolerance (s) to declare convergence
            INTEGER nx            !> number of x grid points in global model
            INTEGER ny            !> number of y grid points in global model
            INTEGER nz            !> number of z grid points in global model
            INTEGER ndivx         !> number of blocks in x direction
            INTEGER ndivy         !> number of blocks in y direction
            INTEGER ndivz         !> number of blocks in z direction
            INTEGER noverlap      !> number of overlapping nodes in ghost communication
            INTEGER maxit         !> max number of iterations in fast sweeping method
            INTEGER iverb         !> controls verbosity
         END TYPE solverParametersType

         TYPE ghostType
            DOUBLE PRECISION, ALLOCATABLE :: send_buff(:) !> Workspace for sending to 
                                                          !> mydest
            DOUBLE PRECISION, ALLOCATABLE :: recv_buff(:) !> Workspace for receiving from
                                                          !> mysrc
            INTEGER, ALLOCATABLE :: isend_dest(:)         !> indices of my local model
                                                          !> that i'm sending to mydest
            INTEGER, ALLOCATABLE :: irecv_dest(:)         !> indices of my local model
                                                          !> that i'm getting from mysrc
            INTEGER nsend          !> number of nodes to send to mydest
            INTEGER nrecv          !> number of nodes to receive from mydest
            INTEGER mysrc          !> process ID i'm receiving from
            INTEGER mydest         !> process ID i'm sending to
            INTEGER send_request   !> MPI status of non-blocking send request
            INTEGER recv_request   !> MPI status of non-blocking receive request
            LOGICAL lrecv          !> if true then i've received this request
         END TYPE ghostType

      END MODULE EIKONAL3D_TYPES

      MODULE EIKONAL3D_MODULE
         USE CONSTANTS_MODULE, ONLY : zero, one, two, three, four, &
                                      third, two_third, half
         USE EIKONAL3D_TYPES
         INTERFACE
            DOUBLE PRECISION FUNCTION GET_UXMIN3D(nx, ny, ix, iy, iz, u)
            IMPLICIT NONE
            DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: u
            INTEGER, INTENT(IN) :: ix, iy, iz, nx, ny
            END FUNCTION GET_UXMIN3D

            DOUBLE PRECISION FUNCTION GET_UYMIN3D(nx, ny, ix, iy, iz, u)
            IMPLICIT NONE
            DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: u
            INTEGER, INTENT(IN) :: ix, iy, iz, nx, ny
            END FUNCTION GET_UYMIN3D

            DOUBLE PRECISION FUNCTION GET_UZMIN3D(nx, ny, nz, ix, iy, iz, u)
            IMPLICIT NONE
            DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: u
            INTEGER, INTENT(IN) :: ix, iy, iz, nx, ny, nz
            END FUNCTION GET_UZMIN3D

            SUBROUTINE SORT3(a, b, c, a1, a2, a3)
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT(IN) :: a, b, c
            DOUBLE PRECISION, INTENT(OUT) :: a1, a2, a3
            END SUBROUTINE SORT3

            DOUBLE PRECISION FUNCTION SOLVE_HAMILTONIAN2D(a, b, fijh)
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT(IN) :: a, b, fijh
            END FUNCTION SOLVE_HAMILTONIAN2D

            DOUBLE PRECISION FUNCTION SOLVE_HAMILTONIAN3D(a, b, c, fijkh, ierr)
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT(IN) :: a, b, c, fijkh
            INTEGER, INTENT(OUT) :: ierr
            END FUNCTION SOLVE_HAMILTONIAN3D

            SUBROUTINE UPDATE3D(nx, ny, nz, ix, iy, iz, h, slow, u, ierr)
            DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: slow
            INTEGER, INTENT(IN) :: ix, iy, iz, nx, ny, nz
            DOUBLE PRECISION, INTENT(IN) :: h
            DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: u
            INTEGER, INTENT(OUT) :: ierr
            END SUBROUTINE UPDATE3D

            SUBROUTINE EIKONAL3D_FSM(iverb, maxit, nx, ny, nz, h, tol,  &
                                     lstruct, lupd, slow,        &
                                     u, ierr)
            USE EIKONAL3D_TYPES, ONLY : levelType
            IMPLICIT NONE
            TYPE(levelType), INTENT(IN) :: lstruct
            DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: slow
            DOUBLE PRECISION, INTENT(IN) :: h, tol
            INTEGER, INTENT(IN) :: iverb, maxit, nx, ny, nz
            LOGICAL, DIMENSION(:), INTENT(IN) :: lupd
            INTEGER, INTENT(OUT) :: ierr
            DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: u
            END SUBROUTINE EIKONAL3D_FSM

            SUBROUTINE EVAL_UPDATE3D(nx, ny, nz,          &
                                     lrevx, lrevy, lrevz, &
                                     h, nnl, ixyz_level,  &
                                     lupd, slow,          &
                                     u, ierr)
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT(IN) :: h
            INTEGER, INTENT(IN) :: nnl, nx, ny, nz
            LOGICAL, INTENT(IN) :: lrevx, lrevy, lrevz
            DOUBLE PRECISION, INTENT(IN), DIMENSION(:) :: slow
            INTEGER, INTENT(IN), DIMENSION(:) :: ixyz_level
            LOGICAL, INTENT(IN), DIMENSION(:) :: lupd
            DOUBLE PRECISION, INTENT(INOUT), DIMENSION(:) :: u
            INTEGER, INTENT(OUT) :: ierr
            END SUBROUTINE EVAL_UPDATE3D

            SUBROUTINE EIKONAL_SCATTER_MODEL(comm, l2g_node, slow, slow_loc)
            IMPLICIT NONE
            DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: slow
            INTEGER, DIMENSION(:), INTENT(IN) :: l2g_node
            INTEGER, INTENT(IN) :: comm
            DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: slow_loc
            END SUBROUTINE EIKONAL_SCATTER_MODEL

            SUBROUTINE EIKONAL_SCATTER_BCS(comm, l2g_node, lisbc, u, lisbc_loc, u_loc)
            IMPLICIT NONE
            DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: u
            INTEGER, DIMENSION(:), INTENT(IN) :: l2g_node
            LOGICAL, DIMENSION(:), INTENT(IN) :: lisbc
            INTEGER, INTENT(IN) :: comm
            DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: u_loc
            LOGICAL, DIMENSION(:), INTENT(OUT) :: lisbc_loc
            END SUBROUTINE EIKONAL_SCATTER_BCS

            SUBROUTINE EIKONAL_GATHER_TRAVELTIMES(comm, l2g_node, u_loc, u)
            IMPLICIT NONE
            DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: u_loc
            INTEGER, DIMENSION(:), INTENT(IN) :: l2g_node
            INTEGER, INTENT(IN) :: comm
            DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: u
            END SUBROUTINE EIKONAL_GATHER_TRAVELTIMES

            SUBROUTINE EIKONAL3D_SETBCS(nx, ny, nz, nsrc,       &
                                        dx, dy, dz, x0, y0, z0, &
                                        ts, xs, ys, zs, slow,       &
                                        lisbc, u, ierr)
            DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: slow, ts, xs, ys, zs
            DOUBLE PRECISION, INTENT(IN) :: dx, dy, dz, x0, y0, z0
            INTEGER, INTENT(IN) :: nx, ny, nz, nsrc
            DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: u
            LOGICAL, DIMENSION(:), INTENT(OUT) :: lisbc
            INTEGER, INTENT(OUT) :: ierr
            END SUBROUTINE EIKONAL3D_SETBCS

            SUBROUTINE EIKONAL3D_FSM_MPI(sweepComm, blockComm,      &
                                         iverb, maxit,              &
                                         nx_loc, ny_loc, nz_loc,    &
                                         h, tol,                    &
                                         lstruct, ghosts,           &
                                         lupd_loc, lghost_node,     &
                                         slow_loc, u, ierr)
            USE EIKONAL3D_TYPES, ONLY : ghostType, levelType
            IMPLICIT NONE
            TYPE(levelType), INTENT(IN) :: lstruct
            TYPE(ghostType), DIMENSION(:), INTENT(INOUT) :: ghosts
            DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: slow_loc
            DOUBLE PRECISION, INTENT(IN) :: h, tol
            LOGICAL, DIMENSION(:), INTENT(IN) :: lghost_node, lupd_loc
            INTEGER, INTENT(IN) :: blockComm, sweepComm, iverb, maxit, &
                                   nx_loc, ny_loc, nz_loc
            DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: u
            INTEGER, INTENT(OUT) :: ierr
            END SUBROUTINE EIKONAL3D_FSM_MPI

            SUBROUTINE EIKONAL3D_GET_LOCAL_TTIMES4(nloc, ttimes4, ierr)
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN) :: nloc
            REAL(C_FLOAT), INTENT(OUT) :: ttimes4(nloc)
            INTEGER(C_INT), INTENT(OUT) :: ierr
            END SUBROUTINE EIKONAL3D_GET_LOCAL_TTIMES4

            SUBROUTINE EIKONAL_INIT_GRID(nx, isx, x0, dx, xs, ixloc, ierr)
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT(IN) :: x0, dx, xs
            INTEGER, INTENT(IN) :: nx, isx
            INTEGER, INTENT(OUT) :: ixloc(3), ierr
            END SUBROUTINE EIKONAL_INIT_GRID

            SUBROUTINE EIKONAL_SOURCE_INDEX(nx, x0, dx, xs, isx)
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT(IN) :: x0, dx, xs
            INTEGER, INTENT(IN) :: nx
            INTEGER, INTENT(OUT) :: isx
            END SUBROUTINE EIKONAL_SOURCE_INDEX

            SUBROUTINE EIKONAL_GRD2IJK(igrd, nx, ny, nz, i, j, k, ierr)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: igrd, nx, ny, nz
            INTEGER, INTENT(OUT) :: i, j, k, ierr
            END SUBROUTINE EIKONAL_GRD2IJK

            SUBROUTINE EIKONAL_EXCHANGE_BLOCKING(comm, ghosts, u)
            USE EIKONAL3D_TYPES, ONLY : ghostType
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: comm
            TYPE(ghostType), DIMENSION(:), INTENT(INOUT) :: ghosts
            DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: u
            END SUBROUTINE EIKONAL_EXCHANGE_BLOCKING

            SUBROUTINE EIKONAL_EXCHANGE(comm, ghosts, u)
            USE EIKONAL3D_TYPES, ONLY : ghostType
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: comm
            TYPE(ghostType), DIMENSION(:), INTENT(INOUT) :: ghosts
            DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: u
            END SUBROUTINE EIKONAL_EXCHANGE

            SUBROUTINE EIKONAL3D_GHOST_COMM(comm, ndivx, ndivy, ndivz, &
                                            noverlap, nx, ny, nz,      &
                                            ghosts, model_loc, ierr)
            USE EIKONAL3D_TYPES, ONLY : ghostType
            USE EIKONAL3D_TYPES, ONLY : localModelType
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: comm, ndivx, ndivy, ndivz, noverlap, nx, ny, nz
            TYPE(ghostType), INTENT(OUT), ALLOCATABLE :: ghosts(:)
            TYPE (localModelType), INTENT(OUT) :: model_loc
            INTEGER, INTENT(OUT) :: ierr
            END SUBROUTINE EIKONAL3D_GHOST_COMM

            SUBROUTINE EIKONAL3D_INITIALIZE(comm, iverb, nx, ny, nz,  &
                                            ndivx, ndivy, ndivz,      &
                                            noverlap, maxit,          &
                                            x0, y0, z0, h, tol,       &
                                            ierr) &
                       BIND(C, NAME='eikonal3d_initialize')
            USE ISO_C_BINDING
            IMPLICIT NONE
            REAL(C_DOUBLE), INTENT(IN) :: x0, y0, z0, h, tol
            INTEGER(C_INT), INTENT(IN) :: comm, iverb, maxit, ndivx, ndivy, ndivz, &
                                          noverlap, nx, ny, nz
            INTEGER(C_INT), INTENT(OUT) :: ierr
            END SUBROUTINE EIKONAL3D_INITIALIZE

            SUBROUTINE EIKONAL3D_SOLVE(comm, nsrc, n,         &
                                       ts, xs, ys, zs, slow,  &
                                       u, ierr)               &
                       BIND(C, NAME='eikonal3d_solve')
            USE ISO_C_BINDING
            INTEGER(C_INT), INTENT(IN) :: comm, n, nsrc
            REAL(C_DOUBLE), INTENT(IN) :: slow(n), ts(nsrc), xs(nsrc), ys(nsrc), zs(nsrc)
            REAL(C_DOUBLE), INTENT(OUT) ::  u(n)
            INTEGER(C_INT), INTENT(OUT) :: ierr
            END SUBROUTINE EIKONAL3D_SOLVE

            SUBROUTINE EIKONAL3D_SERIAL_DRIVER(job, iverb, maxit, nsrc, &
                                               nx, ny, nz,              &
                                               tol, h, x0, y0, z0,      &
                                               ts, xs, ys, zs, slow,    &
                                               u, ierr)                 &
                       BIND(C, NAME='eikonal3d_serial_driver')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN) :: job, iverb, maxit, nsrc, nx, ny, nz
            REAL(C_DOUBLE), INTENT(IN) :: slow(nx*ny*nz), ts(nsrc), xs(nsrc),  &
                                          ys(nsrc), zs(nsrc), h, tol, x0, y0, z0
            REAL(C_DOUBLE), INTENT(OUT) :: u(nx*ny*nz)
            INTEGER(C_INT), INTENT(OUT) :: ierr
            END SUBROUTINE EIKONAL3D_SERIAL_DRIVER

            SUBROUTINE EIKONAL3D_FINALIZE(comm, ierr) BIND(C, NAME='eikonal3d_finalize')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN) :: comm
            INTEGER(C_INT), INTENT(OUT) :: ierr
            END SUBROUTINE EIKONAL3D_FINALIZE

            SUBROUTINE FREE_LEVEL_STRUCT(level)
            USE EIKONAL3D_TYPES, ONLY : levelType
            IMPLICIT NONE
            TYPE(levelType), INTENT(INOUT) :: level
            END SUBROUTINE FREE_LEVEL_STRUCT

            SUBROUTINE FREE_LOCAL_MODEL_STRUCT(model_loc)
            USE EIKONAL3D_TYPES, ONLY : localModelType
            IMPLICIT NONE
            TYPE(localModelType), INTENT(INOUT) :: model_loc
            END SUBROUTINE FREE_LOCAL_MODEL_STRUCT

            SUBROUTINE FREE_GHOST_STRUCT(ghosts)
            USE EIKONAL3D_TYPES, ONLY : ghostType
            IMPLICIT NONE
            TYPE(ghostType), INTENT(INOUT), ALLOCATABLE :: ghosts(:)
            END SUBROUTINE FREE_GHOST_STRUCT

            SUBROUTINE MAKE_LEVEL_STRUCT(nx, ny, nz, lstruct, ierr)
            USE EIKONAL3D_TYPES, ONLY : levelType
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: nx, ny, nz
            TYPE(levelType), INTENT(OUT) :: lstruct
            INTEGER, INTENT(OUT) :: ierr
            END SUBROUTINE

         END INTERFACE

         DOUBLE PRECISION, PARAMETER :: u_nan = HUGE(one)
         TYPE(solverParametersType), SAVE :: parms
         TYPE(ghostType), ALLOCATABLE, SAVE :: ghosts(:)
         TYPE(localModelType), SAVE :: model_loc
         LOGICAL, SAVE :: leik_init
      END MODULE EIKONAL3D_MODULE

      MODULE LOCATE_TYPES
         TYPE locateType
            INTEGER, ALLOCATABLE :: l2g_node(:)
            DOUBLE PRECISION, ALLOCATABLE :: xlocs(:)
            DOUBLE PRECISION, ALLOCATABLE :: ylocs(:)
            DOUBLE PRECISION, ALLOCATABLE :: zlocs(:)
            INTEGER ix0
            INTEGER iy0
            INTEGER iz0
            INTEGER nx
            INTEGER ny
            INTEGER nz
            INTEGER nxLoc
            INTEGER nyLoc
            INTEGER nzLoc
            INTEGER ngrd
         END TYPE locateType
      END MODULE LOCATE_TYPES

      MODULE H5IO_MODULE
         INTERFACE
            SUBROUTINE H5IO_READ_TRAVELTIMESF(comm, tttFileID,        &
                                              station, model, iphase, &
                                              ix0f, iy0f, iz0f,       &
                                              nxLoc, nyLoc, nzLoc,    &
                                              ttimes, ierr)           &
                       BIND(C, NAME='eikonal_h5io_readTraveltimesF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN) :: comm, ix0f, iy0f, iz0f, &
                                          iphase, model,          &
                                          nxLoc, nyLoc, nzLoc,    &
                                          station, tttFileID
            REAL(C_FLOAT), INTENT(OUT) :: ttimes(nxLoc*nyLoc*nzLoc)
            INTEGER(C_INT), INTENT(OUT) :: ierr
            END SUBROUTINE H5IO_READ_TRAVELTIMESF
         END INTERFACE
      END MODULE H5IO_MODULE

      MODULE LOCATE_MODULE
         USE CONSTANTS_MODULE, ONLY : one, sqrt2i, two, zero
         USE LOCATE_TYPES, ONLY : locateType
         INTERFACE
            SUBROUTINE LOCATE_INITIALIZE_PDF(ngrd, pdf)
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN) :: ngrd
            REAL(C_DOUBLE), INTENT(OUT) :: pdf(ngrd)
            END SUBROUTINE LOCATE_INITIALIZE_PDF
         END INTERFACE
         TYPE(locateType), SAVE :: locate_loc
         INTEGER, PARAMETER :: COMPUTE_NONE = 0
         INTEGER, PARAMETER :: COMPUTE_LOCATION_ONLY = 1
         INTEGER, PARAMETER :: COMPUTE_LOCATION_AND_ORIGIN_TIME = 2 
         INTEGER, PARAMETER :: COMPUTE_LOCATION_AND_STATICS = 3 
         INTEGER, PARAMETER :: COMPUTE_ORIGIN_TIME_AND_STATICS = 4
         INTEGER, PARAMETER :: COMPUTE_LOCATION_ALL = 5
      END MODULE LOCATE_MODULE


!>    @brief Splits the communicator and renumbers the ranks for the eikonal solver
!>           and locator while respecting the two primary types of interprocess
!>           communication.  
!>
!>    @param[in] comm               global communicator to split
!>    @param[in] lreord             if this is true then the processes IDs on global_comm
!>                                  will be renumbered.
!>                                  if this is false then the processes IDs on global_comm
!>                                  will retain values identical to that of the IDs
!>                                  on comm.
!>    @param[in] iwt                if this is 1 then weighting based on anticipated
!>                                  volume of communication will be employed.
!>                                  if this is 0 then no weighting will be applied.
!>    @param[in] ndivx              number of blocks to divide domain into in x
!>    @param[in] ndivy              number of blocks to divide domain into in y
!>    @param[in] ndivz              number of blocks to divide domain into in z
!>
!>    @param[out] global_comm       new global communicator to take the place of comm
!>    @param[out] intra_table_comm  the communicator between the blocks comprising
!>                                  the traveltime table. 
!>    @param[out] inter_table_comm  the communicator between block numbers of 
!>                                  each traveltime table.
!>    @param[out] ierr              0 indicates succcess
!>
!>    @author Ben Baker
!>
!>    @copyright Apache 2 license
!>
!>    @details To get a better idea of this, consider the problem where we locate an
!>             earthquake in a 2D physical domain with 36 processors and 3 observations.
!>             Because there are 3 observations, maybe two P times and an S time, we
!>             require 3 traveltime grids be computed.  Each grid then uses 12 processes
!>             to solve the eikonal equation.   This means2D physical domain partitioned
!>             into 12 processes. 
!>
!>             Processes holding travel time table 1:
!>
!>             x-----x-----x-----x-----x
!>             |  8  |  9  | 10  | 11  |
!>             x-----x-----x-----x-----x
!>             |  4  |  5  |  6  |  7  |
!>             x-----x-----x-----x-----x
!>             |  0  |  1  |  2  |  4  |
!>             x-----x-----x-----x-----x
!>
!>             Processes holding travel time table 2:
!>
!>             x-----x-----x-----x-----x
!>             | 20  | 21  | 22  | 23  |
!>             x-----x-----x-----x-----x
!>             | 16  | 17  | 18  | 19  |
!>             x-----x-----x-----x-----x
!>             | 12  | 13  | 14  | 15  |
!>             x-----x-----x-----x-----x
!>
!>             Processes holding travel time table 3:
!>
!>             x-----x-----x-----x-----x
!>             | 32  | 33  | 34  | 35  |
!>             x-----x-----x-----x-----x
!>             | 28  | 29  | 30  | 31  |
!>             x-----x-----x-----x-----x
!>             | 24  | 25  | 26  | 27  |
!>             x-----x-----x-----x-----x
!>
!>             Now, we analyze the communication for process 17.  It must talk to
!>             processes 12, 13, 14, 16, 18, 20, 21, and 22 in the finite differencing.
!>             However, because it shares a face with processes, 13, 16, 18, and 21 it
!>             must communicate more - hence the weights are different.
!>
!>             Additionally, during the location grid search we must reduce the norm
!>             of the residuals across the observations.  Thus, process 17 must talk
!>             to processes 5 and 29 with a very heavy communication volume (it will
!>             send an entire domain) - hence this is weighted the greatest. 
!>   
!>             With this in mind we create a topology for process 17:
!>               ids = {5,   12,  13,   14,  16,  18,   20,  21,  22, 29}
!>               wts = {1,    1,   1,    1,   1,   1,    1,   1,   1,  1}
!>             or 
!>               wts = {3,    1,   2,    1,   2,   2,    3,   2,   3,  1} 
!>             depending if you use the weight model or not.
!>
!>             Finally, the scatters and gathers from the master of each table
!>             to all other elements in the table is considered.  Hence, the table
!>             master ID will communicate with all other blocks in the table and
!>             each block will report back to the master all with weight 1. 
!>
!>             Then, MPI will take this information and renumber the processes on the
!>             input communicator because the job scheduler has given us processors on
!>             nodes that are scattered across the cluster/supercomputer.  Hence, while
!>             the program sees 36 processes, MPI actually sees their physical locations
!>             and will pick a reordering that can allow communication to happen faster.
!>             And this makes sense, two chips connected via a serial bus will exchange
!>             information faster than two chips on different nodes who must send their
!>             data to the bus, which goes to the nic, which puts the data into the tubes
!>             connecting the nodes, back down to the nic on the other node, to the bus,
!>             and finally onto the correct chip.
!> 
      SUBROUTINE MPIUTILS_SPLITCOMM3D(comm,                               &
                                      lreord, iwt, ndivx, ndivy, ndivz,   &
                                      global_comm, intra_table_comm,      &
                                      inter_table_comm, ierr)
      USE MPIUTILS_MODULE, ONLY : MPIUTILS_GRD2IJK, linitComm
      USE MPI
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: comm, iwt, ndivx, ndivy, ndivz
      LOGICAL, INTENT(IN) :: lreord
      INTEGER, INTENT(OUT) :: global_comm, intra_table_comm, &
                              inter_table_comm, ierr
      ! local variables
      INTEGER, ALLOCATABLE :: dstwts(:), dsts(:), iwts(:), srcs(:), srcwts(:)
      INTEGER i, idivx, idivy, idivz, info, ixn, j, jyn, k, kzn, &
              master_block_id, mpierr, &
              myblock_id, myid, mytable_id, nblocks, ncomm, ndst, &
              neigh, nprocs, nsrc, ntables
      INTEGER, PARAMETER :: wts(27) = [1, 2, 1,  2, 3, 2,  1, 2, 1, & ! base
                                       2, 3, 2,  3, 0, 3,  2, 3, 2, & ! middle
                                       1, 2, 1,  2, 3, 2,  1, 2, 1]   ! top
      LOGICAL, PARAMETER :: dupComm = .TRUE.
      IF (linitComm) THEN
         WRITE(*,*) 'mpiutils_splitcomm3d: Warning communicator already initialized'
      ENDIF
      ! get MPI information
      CALL MPI_COMM_SIZE(comm, nprocs, mpierr)
      CALL MPI_COMM_RANK(comm, myid,   mpierr)
      ! error checks
      IF (ndivx < 1 .OR. ndivy < 1 .OR. ndivz < 1) THEN
         IF (ndivx < 1) WRITE(*,*) 'mpiutils_splitcomm3d: ndivx must be positive', ndivx
         IF (ndivy < 1) WRITE(*,*) 'mpiutils_splitcomm3d: ndivy must be positive', ndivy
         IF (ndivz < 1) WRITE(*,*) 'mpiutils_splitcomm3d: ndivz must be positive', ndivz
         ierr = 1
         RETURN
      ENDIF
      nblocks = ndivx*ndivy*ndivz
      IF (nblocks > nprocs) THEN
         WRITE(*,*) 'mpiutils_splitcomm3d: Error insufficient processes for blocks', &
                    nprocs, nblocks, myid
         ierr = 1
         RETURN
      ENDIF
      IF (MOD(nprocs, nblocks) /= 0) THEN
         WRITE(*,*) 'mpiutils_splitcomm3d: Error cant divide groups and blocks'
         ierr = 1
         RETURN
      ENDIF
      ! figure out local IDs
      ntables = nprocs/nblocks
      mytable_id = myid/nblocks
      myblock_id = MOD(myid, nblocks) 
      master_block_id = mytable_id*nblocks
      ! figure out local IDs
      CALL MPIUTILS_GRD2IJK(myblock_id, ndivx, ndivy, ndivz, &
                            idivx, idivy, idivz, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'mpiutils_splitcomm3d: Error computing rank in grid', myid, myblock_id
         RETURN
      ENDIF
      ! set space and intiailize
      ncomm = 0
      nsrc = 0
      ndst = 0
      ALLOCATE(dsts(nprocs)) !27 + ntables - 1))
      ALLOCATE(iwts(nprocs)) !27 + ntables - 1))
      ALLOCATE(srcs(nprocs)) !27 + ntables - 1))
      dsts(:) =-1
      iwts(:) =-1
      srcs(:) =-1
      ! tabulate who in this traveltime table i have to talk to
      DO 1 k=-1,1
         DO 2 j=-1,1
            DO 3 i=-1,1
               IF (i == 0 .AND. j == 0 .AND. k == 0) CYCLE ! this is me
               ixn = idivx + i
               jyn = idivy + j
               kzn = idivz + k
               IF (ixn >= 0 .AND. ixn < ndivx .AND. &
                   jyn >= 0 .AND. jyn < ndivy .AND. &
                   kzn >= 0 .AND. kzn < ndivz) THEN
                  ncomm = ncomm + 1
                  ndst = ndst + 1
                  nsrc = nsrc + 1
                  neigh = mytable_id*nblocks + kzn*ndivx*ndivy + jyn*ndivx + ixn
                  dsts(ndst) = neigh
                  srcs(nsrc) = neigh
                  iwts(ncomm) = wts((k + 1)*9 + (j + 1)*3 + i + 2)
               ENDIF
    3       CONTINUE 
    2    CONTINUE
    1 CONTINUE
      ! tabulate who on the other tables i have to talk to
      DO 4 i=0,ntables-1
         IF (i == mytable_id) CYCLE
         ncomm = ncomm + 1
         iwts(ncomm) = 4
         dsts(ncomm) = nblocks*i + myblock_id
         srcs(ncomm) = nblocks*i + myblock_id
    4 CONTINUE
      ! remeber the scatter and gathers in the domain
      IF (myid == master_block_id) THEN
         DO 5 i=1,nblocks-1
            DO 6 j=1,ncomm
               IF (srcs(j) == mytable_id*nblocks + i) GOTO 600
    6       CONTINUE
            ncomm = ncomm + 1
            srcs(ncomm) = mytable_id*nblocks + i
            dsts(ncomm) = mytable_id*nblocks + i
            iwts(ncomm) = 1
  600       CONTINUE
    5    CONTINUE 
      ELSE
         DO 7 i=1,ncomm
            IF (dsts(i) == master_block_id) GOTO 700 
    7    CONTINUE
         ncomm = ncomm + 1 
         dsts(ncomm) = master_block_id
         srcs(ncomm) = master_block_id
         iwts(ncomm) = 1
  700    CONTINUE
      ENDIF
      ALLOCATE(srcwts(ncomm))
      ALLOCATE(dstwts(ncomm))
      IF (iwt == 0) THEN
         srcwts(1:ncomm) = 1
         dstwts(1:ncomm) = 1
      ELSE
         srcwts(1:ncomm) = iwts(1:ncomm)
         dstwts(1:ncomm) = iwts(1:ncomm)
      ENDIF
   
print *, myid, srcs(1:ncomm), dsts(1:ncomm), srcwts, dstwts 
      IF (dupComm) THEN
print *, 'fix me'
         CALL MPI_COMM_DUP(comm, global_comm, mpierr) 
      ELSE
         CALL MPI_DIST_GRAPH_CREATE_ADJACENT(comm,                      &
                                             ncomm, srcs, srcwts,       &
                                             ncomm, dsts, dstwts,       &
                                             info, lreord, global_comm, &
                                             mpierr)
      ENDIF
      DEALLOCATE(dsts)
      DEALLOCATE(iwts)
      DEALLOCATE(srcs)
      DEALLOCATE(srcwts)
      DEALLOCATE(dstwts)
      ! now split the new communicator 
      CALL MPI_COMM_RANK(global_comm, myid, mpierr)
      CALL MPI_COMM_SPLIT(global_comm, mytable_id, myid, intra_table_comm, mpierr)
      CALL MPI_COMM_SPLIT(global_comm, myblock_id, myid, inter_table_comm, mpierr)
      CALL MPI_COMM_SIZE(intra_table_comm, nprocs, mpierr)
      IF (nprocs /= nblocks) THEN
         WRITE(*,*) 'mpiutils_splitcomm3d: Error splitting intra_table_comm'
         ierr = 1
         RETURN
      ENDIF
      CALL MPI_COMM_SIZE(inter_table_comm, nprocs, mpierr)
      IF (nprocs /= ntables) THEN
         WRITE(*,*) 'mpiutils_splitcomm3d: Error splitting inter_table_comm'
         ierr = 1
         RETURN
      ENDIF
      linitComm = .TRUE.
      RETURN
      END !SUBROUTINE MPIUTILS_SPLITCOMM3
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Convenience function - when given a grid point returns the (i,jk) index in
!>           index in (x, y, z) such that: igrd = (k-1)*nx*ny + (j-1)*nx + i
!>
!>    @param[in] igrd     grid point number [1,nx*ny*nz]
!>    @param[in] nx       number of x grid points in domain
!>    @param[in] ny       number of y grid points in domain
!>    @param[in] nz       number of z grid points in domain
!>
!>    @param[out] i       corresponding i grid point in x
!>    @param[out] j       corresponding j grid point in y
!>    @param[out] k       corresponding k grid point in z
!>    @param[out] ierr    0 indicates success
!>
!>    @author Ben Baker
!>
!>    @copyright Apache 2 license
!>
      SUBROUTINE MPIUTILS_GRD2IJK(igrd, nx, ny, nz, i, j, k, ierr) &
                 BIND(C, NAME='mpiutils_grd2ijk')
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: igrd, nx, ny, nz
      INTEGER(C_INT), INTENT(OUT) :: i, j, k, ierr 
      INTEGER nxy
      ierr = 0
      nxy = nx*ny
      k = igrd/(nxy)
      j = (igrd - k*nxy)/nx
      i =  igrd - k*nxy - j*nx
      IF (i < 0 .OR. i >= nx) ierr = ierr + 1
      IF (j < 0 .OR. j >= ny) ierr = ierr + 1
      IF (k < 0 .OR. k >= nz) ierr = ierr + 1  
      IF (k*nxy + j*nx + i /= igrd) ierr = ierr + 1
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !

      SUBROUTINE MPIUTILS_SPLITCOMM2D(comm,                            &
                                      lreord, iwt, ndivx, ndivz,       &
                                      global_comm, intra_table_comm,   &
                                      inter_table_comm, ierr)
      USE MPIUTILS_MODULE, ONLY : MPIUTILS_SPLITCOMM3D
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: comm, iwt, ndivx, ndivz
      LOGICAL, INTENT(IN) :: lreord
      INTEGER, INTENT(OUT) :: global_comm, intra_table_comm, &
                              inter_table_comm, ierr
      CALL MPIUTILS_SPLITCOMM3D(comm, lreord, iwt, ndivx, 1, ndivz, &
                                global_comm, intra_table_comm,      &
                                inter_table_comm, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'mpiutils_splitcomm2d: Error splitting communicator!'
         RETURN
      ENDIF
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Initializes the communicators for the 3D tomography and location.
!>           Further details can be found in the splitcomm3d routine.
!>
!>    @param[in] comm      global communicator to split and set.
!>    @param[in] ireord    if 1 then the MPI IDs will be reordered based on proximity. 
!>    @param[in] iwt       if 1 then the estimated communication will volume be
!>                         applied to reordering ranks.
!>    @param[in] ndivx     number of divisions of domain in x direction.
!>    @param[in] ndivy     number of divisions of domain in y direction.
!>    @param[in] ndivz     number of divisions of domain in z direction.
!>
!>    @param[out] ierr     0 indicates success
!>
!>    @author Ben Baker
!>
!>    @copyright Apache 2
!> 
      SUBROUTINE MPIUTILS_INITIALIZE3D(comm, ireord, iwt,         &
                                       ndivx, ndivy, ndivz, ierr) &
                 BIND(C, NAME='mpiutils_initialize3d')
      USE MPIUTILS_MODULE, ONLY : global_comm, intra_table_comm, inter_table_comm
      USE MPIUTILS_MODULE, ONLY : MPIUTILS_SPLITCOMM3D
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: comm, ireord, iwt, ndivx, ndivy, ndivz
      INTEGER(C_INT), INTENT(OUT) :: ierr
      LOGICAL lreord
      lreord = .FALSE.
      IF (ireord == 1) lreord = .TRUE.
      CALL MPIUTILS_SPLITCOMM3D(comm,                             &
                                lreord, iwt, ndivx, ndivy, ndivz, &
                                global_comm, intra_table_comm,    &
                                inter_table_comm, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'mpiutils_intialize3d: Error generating communicators'
         RETURN
      ENDIF
      RETURN
      END

      SUBROUTINE MPIUTILS_INITIALIZE2D(comm, ireord, iwt, ndivx, ndivz, ierr) &
                 BIND(C, NAME='mpiutils_initialize2d')
      USE MPIUTILS_MODULE, ONLY : global_comm, intra_table_comm, inter_table_comm
      USE MPIUTILS_MODULE, ONLY : MPIUTILS_SPLITCOMM2D
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: comm, ireord, iwt, ndivx, ndivz
      INTEGER(C_INT), INTENT(OUT) :: ierr
      LOGICAL lreord
      lreord = .FALSE.
      IF (ireord == 1) lreord = .TRUE.
      CALL MPIUTILS_SPLITCOMM2D(comm,                            &
                                lreord, iwt, ndivx, ndivz,       &
                                global_comm, intra_table_comm,   &
                                inter_table_comm, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'mpiutils_intialize2d: Error generating communicators'
         RETURN
      ENDIF
      RETURN
      END

      SUBROUTINE MPIUTILS_FINALIZE() &
                 BIND(C, NAME='mpiutils_finalize')
      USE MPIUTILS_MODULE, ONLY : global_comm, intra_table_comm, &
                                  inter_table_comm, linitComm
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER mpierr
      IF (linitComm) THEN
         CALL MPI_COMM_FREE(global_comm, mpierr)
         CALL MPI_COMM_FREE(intra_table_comm, mpierr)
         CALL MPI_COMM_FREE(inter_table_comm, mpierr)
         linitComm = .FALSE.
      ENDIF
      RETURN
      END

      SUBROUTINE MPIUTILS_GET_COMMUNICATORS(globalComm, intraTableComm, &
                                            interTableComm, ierr)       &
                 BIND(C, NAME='mpiutils_getCommunicators')
      USE MPIUTILS_MODULE, ONLY : global_comm, intra_table_comm, &
                                  inter_table_comm, linitComm
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(OUT) :: globalComm, intraTableComm, &
                                     interTableComm, ierr
      ierr = 0
      IF (.NOT.linitComm) THEN
         WRITE(*,*) 'mpiutils_get_communicators: Never initialized communicators'
         ierr = 1
         RETURN
      ENDIF
      globalComm = global_comm
      intraTableComm = intra_table_comm
      interTableComm = inter_table_comm
      RETURN
      END

      subroutine test_mpi()
      USE MPI
      CALL MPI_INIT(mpierr)
      ndivx = 2
      ndivy = 3
      ndivz = 4
!     CALL MPIUTILS_SPLITCOMM3D(MPI_COMM_WORLD, ndivx, ndivy, ndivz, ierr)

      CALL MPI_FINALIZE(mpierr) 
      STOP
      END

!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Solves the eikonal equation using the fast sweeping method described
!>           by Zhao (2004) and Detrixhe, Gibou, and Min (2012).
!>
!>    @param[in] maxit    max number of iterations in fast sweeping method
!>    @param[in] nx       number of x grid points in domain
!>    @param[in] ny       number of y grid points in domain
!>    @param[in] nz       number of z grid points in domain
!>    @param[in] tol      the iterative method will terminate if the maximum update
!>                        to the traveltimes from iteration k, to iteration k + 1
!>                        is less than tol (tol is specified in seconds)
!>    @param[in] lstruct  level structure corresponding to the Cuthill-Mckee ordering 
!>    @param[in] lupd     if true then update the i'th node.  otherwise, the input
!>                        nodal value will not be altered
!>
!>    @param[inout] u     on input contains the initialized traveltimes at the
!>                        boundary conditions.
!>                        on exit contains the traveltimes to all nodes in the grid
!>
!>    @param[out] ierr    0 indicates success
!>
!>    @author Ben Baker
!>
!>    @copyright Apache 2 license
!>
      SUBROUTINE EIKONAL3D_FSM(iverb, maxit, nx, ny, nz, h, tol, &
                               lstruct, lupd, slow,       &
                               u, ierr)
      USE EIKONAL3D_TYPES, ONLY : levelType
      USE EIKONAL3D_MODULE, ONLY : EVAL_UPDATE3D 
      IMPLICIT NONE
      TYPE(levelType), INTENT(IN) :: lstruct
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: slow
      DOUBLE PRECISION, INTENT(IN) :: h, tol
      INTEGER, INTENT(IN) :: iverb, maxit, nx, ny, nz
      LOGICAL, DIMENSION(:), INTENT(IN) :: lupd
      INTEGER, INTENT(OUT) :: ierr
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: u
      ! local variables
      DOUBLE PRECISION, ALLOCATABLE :: u0(:)
      DOUBLE PRECISION t1, t2
      INTEGER i, i1, i2, isweep, k, lconv, level, nnl, nxyz
      LOGICAL lrevx, lrevy, lrevz
      LOGICAL, PARAMETER :: evalSweep(24) = [.FALSE., .FALSE., .FALSE., &
                                             .TRUE.,  .FALSE., .FALSE., &
                                             .FALSE., .TRUE.,  .FALSE., &
                                             .TRUE.,  .TRUE.,  .FALSE., &
                                             .FALSE., .FALSE., .TRUE.,  &
                                             .TRUE.,  .FALSE., .TRUE.,  &
                                             .FALSE., .TRUE.,  .TRUE.,  &
                                             .TRUE.,  .TRUE.,  .TRUE.]
      !----------------------------------------------------------------------------------!
      !
      ! initialize
      ierr = 0
      nxyz = nx*ny*nz
      ALLOCATE(u0(nxyz))
      u0(:) = u(:)
      ! loop on the iterations
      DO 1000 k=1,maxit
         IF (iverb > 2) THEN
            WRITE(*,*) 'eikonal_fsm: Beginning sweep:', k
         ENDIF 
         CALL CPU_TIME(t1) 
         ! apply the sweeps
         DO 2000 isweep=1,8
            lrevx = evalSweep(3*(isweep-1)+1)
            lrevy = evalSweep(3*(isweep-1)+2)
            lrevz = evalSweep(3*(isweep-1)+3)
            ! for each level perform the parallel update
            DO 3000 level=1,lstruct%nlevels
               nnl = lstruct%nlevels
               i1 = lstruct%level_ptr(level)
               i2 = lstruct%level_ptr(level+1) - 1
               nnl = lstruct%nnl(level) 
               CALL EVAL_UPDATE3D(nx, ny, nz,                        &
                                  lrevx, lrevy, lrevz,               &
                                  h, nnl, lstruct%ixyz_level(i1:i2), &
                                  lupd, slow,                        &
                                  u, ierr)
 3000       CONTINUE 
 2000    CONTINUE
         ! test the convergence of the iteration
         lconv = 0
         DO 1005 i=1,nxyz
            IF (ABS(u0(i) - u(i)) < tol) lconv = lconv + 1
            u0(i) = u(i)
 1005    CONTINUE
         IF (iverb > 2) THEN
            CALL CPU_TIME(t2)
            WRITE(*,*) 'eikononal_fsm: Sweep time:', t2 - t1
         ENDIF 
         IF (lconv == nxyz) EXIT
 1000 CONTINUE
      IF (ALLOCATED(u0)) DEALLOCATE(u0)
      RETURN
      END SUBROUTINE EIKONAL3D_FSM
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE EIKONAL3D_FSM_MPI(sweepComm, blockComm,     &
                                   iverb, maxit,             &
                                   nx_loc, ny_loc, nz_loc,   &
                                   h, tol,                   &
                                   lstruct, ghosts,          &
                                   lghost_node, lupd_loc,    & 
                                   slow_loc, uloc, ierr)
      USE MPI
      USE EIKONAL3D_MODULE, ONLY : EIKONAL_EXCHANGE, EVAL_UPDATE3D
      USE EIKONAL3D_TYPES, ONLY : ghostType, levelType
      IMPLICIT NONE
      TYPE(levelType), INTENT(IN) :: lstruct
      TYPE(ghostType), DIMENSION(:), INTENT(INOUT) :: ghosts
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: slow_loc 
      DOUBLE PRECISION, INTENT(IN) :: h, tol
      LOGICAL, DIMENSION(:), INTENT(IN) :: lghost_node, lupd_loc
      INTEGER, INTENT(IN) :: blockComm, sweepComm, iverb, maxit, nx_loc, ny_loc, nz_loc
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: uloc
      INTEGER, INTENT(OUT) :: ierr
      ! local variables
      DOUBLE PRECISION, ALLOCATABLE :: u0loc(:)
      DOUBLE PRECISION t1, t2
      INTEGER i, i1, i2, isweep, k, ksweep, level, mpierr,       &
              my_blockid, my_sweepid, myid, nblocks, nconv, nnl, &
              nxyz, nxyz_loc, nbuf, nsg
      LOGICAL lconv, lrevx, lrevy, lrevz
      INTEGER, PARAMETER :: master = 0
      LOGICAL, PARAMETER :: evalSweep(24) = [.FALSE., .FALSE., .FALSE., &
                                             .TRUE.,  .FALSE., .FALSE., &
                                             .FALSE., .TRUE.,  .FALSE., &
                                             .TRUE.,  .TRUE.,  .FALSE., &
                                             .FALSE., .FALSE., .TRUE.,  &
                                             .TRUE.,  .FALSE., .TRUE.,  &
                                             .FALSE., .TRUE.,  .TRUE.,  &
                                             .TRUE.,  .TRUE.,  .TRUE.]
      !----------------------------------------------------------------------------------!
      !
      ! get MPI information
      CALL MPI_COMM_SIZE(blockComm, nblocks,    mpierr) ! number of blocks in domain
      CALL MPI_COMM_RANK(blockComm, my_blockid, mpierr) ! my block number
      CALL MPI_COMM_SIZE(sweepComm, nsg, mpierr)        ! number of sweep groups
!     CALL MPI_COMM_RANK(sweepComm, my_sweepid, mpierr) ! my sweep number
nsg = 1
my_sweepid = 0
      myid = my_blockid + my_sweepid
      nxyz_loc = nx_loc*ny_loc*nz_loc 
      nbuf = 0
      DO i=1,nxyz_loc
         IF (.NOT.lghost_node(i)) nbuf = nbuf + 1   
      ENDDO
      CALL MPI_ALLREDUCE(nbuf, nxyz, 1, MPI_INTEGER, MPI_SUM, blockComm, mpierr)
      ALLOCATE(u0loc(MAX(nxyz_loc, 1)))
      u0loc(:) = uloc(:)
      lconv = .FALSE.
      ! loop on max iterations
      DO 1000 k=1,maxit
         IF (myid == master .AND. iverb > 2) THEN
            WRITE(*,*) 'eikonal3d_fsm_mpi: Beginning iteration:', k
         ENDIF
         t1 = MPI_WTIME()
         ! parallel loop on sweeps
         DO 2000 ksweep=1,8,nsg
            isweep = (ksweep - 1)*nsg + my_sweepid + 1
            IF (isweep > 8) GOTO 2500
            ! scatter grid
            IF (nsg > 1) THEN

            ENDIF
            lrevx = evalSweep(3*(isweep-1)+1)
            lrevy = evalSweep(3*(isweep-1)+2)
            lrevz = evalSweep(3*(isweep-1)+3)
            ! for each level perform the parallel update
            DO 3000 level=1,lstruct%nlevels
               ! evaluate the sweep
               nnl = lstruct%nlevels
               i1 = lstruct%level_ptr(level)
               i2 = lstruct%level_ptr(level+1) - 1
               nnl = lstruct%nnl(level)
               CALL EVAL_UPDATE3D(nx_loc, ny_loc, nz_loc,            &
                                  lrevx, lrevy, lrevz,               &
                                  h, nnl, lstruct%ixyz_level(i1:i2), &
                                  lupd_loc, slow_loc,                &
                                  uloc, ierr)
               ! swap with neighbors
 3000       CONTINUE 
            CALL EIKONAL_EXCHANGE(blockComm, ghosts, uloc)
 2500       CONTINUE ! get next sweep
            ! take the minimum from each grid
            IF (nsg > 1) THEN
!              ubuff(:) = uloc(:)
!              MPI_ALLREDUCE(ubuff, u, nx_loc*ny_loc*nz_loc,         &
!                            MPI_DOUBLE_PRECISION, MPI_MIN, sweepComm)
            ENDIF 
 2000    CONTINUE
         ! convergence test
         nbuf = 0
         DO 11 i=1,nxyz_loc
            IF (.NOT.lghost_node(i) .AND. ABS(uloc(i) - u0loc(i)) < tol) nbuf = nbuf + 1
            u0loc(i) = uloc(i) 
   11    CONTINUE 
         CALL MPI_ALLREDUCE(nbuf, nconv, 1, MPI_INTEGER, MPI_SUM, blockComm, mpierr)
         t2 = MPI_WTIME()
         IF (myid == master .AND. iverb > 2) THEN
            WRITE(*,*) 'eikonal3d_fsm_mpi: Sweep time (s)', t2 - t1
         ENDIF
         IF (nconv == nxyz) THEN
            lconv = .TRUE.
            EXIT
         ENDIF
 1000 CONTINUE
      IF (myid == master .AND. iverb > 2) THEN
         IF (lconv) THEN
            WRITE(*,*) 'eikonal3d_fsm_mpi: Number of iterations for convergence:', k
         ELSE
            WRITE(*,*) 'eikonal3d_fsm_mpi: Failed to converge after iterations', maxit
         ENDIF
      ENDIF
      DEALLOCATE(u0loc)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE MAKE_LEVEL_STRUCT(nx, ny, nz, lstruct, ierr)
      USE EIKONAL3D_TYPES, ONLY : levelType
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nx, ny, nz
      TYPE(levelType), INTENT(OUT) :: lstruct
      ! local variables
      INTEGER, ALLOCATABLE :: linit(:)
      TYPE two2d_array
         INTEGER, ALLOCATABLE :: dim2(:)
      END TYPE
      TYPE(two2d_array), ALLOCATABLE :: ixyz_level(:)
      INTEGER, INTENT(OUT) :: ierr
      INTEGER i1, i2, ijk, ix, iy, iz, j1, j2, k1, k2, level, np

      ierr = 0
      lstruct%nlevels = nx + ny + nz - 2
      ALLOCATE(lstruct%level_ptr(lstruct%nlevels+1))
      ALLOCATE(lstruct%nnl(lstruct%nlevels))
      ALLOCATE(ixyz_level(lstruct%nlevels))
      lstruct%level_ptr(:) = 0
      ALLOCATE(linit(nx*ny*nz))
      lstruct%level_ptr(1) = 1
      linit(:) = 0
      DO 1 level=2,nx+ny+nz-1
         k1 = MAX(0, level - 2 - (nx - 1) - (ny - 1))
         k1 = MIN(k1, nz - 1)
         k2 = MIN(nz - 1, level - 2)
         np = 0
         DO 2 iz=k1,k2
            j1 = MAX(0, level - 2 - iz - (nx - 1))
            j1 = MIN(j1, ny - 1)
            j2 = MIN(ny - 1, level - 2 - iz)
            DO 3 iy=j1,j2
               i1 = MAX(0, level - 2 - iz - iy)
               i1 = MIN(i1, nx - 1)
               i2 = MIN(nx - 1, level - 2 - iz - iy)
               DO 4 ix=i1,i2
                  np = np + 1
                  ijk = iz*nx*ny + iy*nx + ix + 1
                  linit(ijk) = linit(ijk) + 1
    4          CONTINUE 
    3       CONTINUE
    2    CONTINUE
         lstruct%nnl(level-1) = np
         lstruct%level_ptr(level) = lstruct%level_ptr(level-1) + 3*np
         IF (.NOT.ALLOCATED(ixyz_level(level-1)%dim2)) THEN
            ALLOCATE(ixyz_level(level-1)%dim2(3*np))
         ENDIF
         np = 0
         DO 11 iz=k1,k2
            j1 = MAX(0, level - 2 - iz - (nx - 1)) 
            j1 = MIN(j1, ny - 1)
            j2 = MIN(ny - 1, level - 2 - iz) 
            DO 12 iy=j1,j2
               i1 = MAX(0, level - 2 - iz - iy) 
               i1 = MIN(i1, nx - 1)
               i2 = MIN(nx - 1, level - 2 - iz - iy) 
               DO 13 ix=i1,i2
                  np = np + 1 
                  ixyz_level(level-1)%dim2(3*(np-1)+1) = ix + 1
                  ixyz_level(level-1)%dim2(3*(np-1)+2) = iy + 1
                  ixyz_level(level-1)%dim2(3*(np-1)+3) = iz + 1
   13          CONTINUE 
   12       CONTINUE
   11    CONTINUE
    1 CONTINUE
      ALLOCATE(lstruct%ixyz_level(lstruct%level_ptr(lstruct%nlevels+1)-1))
      lstruct%ixyz_level(:) = 0
      ! copy the level structure
      DO 21 level=2,nx+ny+nz-1
         i1 = lstruct%level_ptr(level-1)
         i2 = lstruct%level_ptr(level) - 1
         lstruct%ixyz_level(i1:i2) = ixyz_level(level-1)%dim2(:)
         DEALLOCATE(ixyz_level(level-1)%dim2)
   21 CONTINUE
      IF (MINVAL(linit) < 1 .OR. MAXVAL(linit) > 1) THEN
         WRITE(*,*) 'Failed to initialize node'
         ierr = 1
      ENDIF
      DEALLOCATE(ixyz_level) 
      DEALLOCATE(linit)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Frees the level structure
!>
!>    @param[inout] level   on input this is the level structure.
!>                          on output all memory has been freed from level.
!>
!>    @author Ben Baker
!>
!>    @copyright Apache 2 license
!>
      SUBROUTINE FREE_LEVEL_STRUCT(level)
      USE EIKONAL3D_TYPES, ONLY : levelType
      IMPLICIT NONE
      TYPE(levelType), INTENT(INOUT) :: level
      IF (ALLOCATED(level%ixyz_level)) DEALLOCATE(level%ixyz_level)
      IF (ALLOCATED(level%nnl))        DEALLOCATE(level%nnl)
      IF (ALLOCATED(level%level_ptr))  DEALLOCATE(level%level_ptr)
      level%nlevels = 0
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Frees the memory on the local model structure
!>
!>    @param[inout] model_loc    on input this is the local model structure.
!>                               on output all memroy has bene freed from model_loc.
!>
!>    @author Ben Baker
!>
      SUBROUTINE FREE_LOCAL_MODEL_STRUCT(model_loc)
      USE EIKONAL3D_TYPES, ONLY : localModelType
      USE EIKONAL3D_MODULE, ONLY : FREE_LEVEL_STRUCT
      IMPLICIT NONE
      TYPE(localModelType), INTENT(INOUT) :: model_loc
      CALL FREE_LEVEL_STRUCT(model_loc%lstruct)
      IF (ALLOCATED(model_loc%slow))        DEALLOCATE(model_loc%slow)
      IF (ALLOCATED(model_loc%u))           DEALLOCATE(model_loc%u)
      IF (ALLOCATED(model_loc%l2g_node))    DEALLOCATE(model_loc%l2g_node)
      IF (ALLOCATED(model_loc%lghost_node)) DEALLOCATE(model_loc%lghost_node)
      IF (ALLOCATED(model_loc%lupd))        DEALLOCATE(model_loc%lupd)
      model_loc%nx = 0
      model_loc%ny = 0
      model_loc%nz = 0
      model_loc%myblock = 0
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Frees the ghost communication structure
!>
!>    @param[inout] ghosts    on input contains the ghost structure.
!>                            on output all memory freed from the ghost structure.
!>
!>    @author Ben Baker
!>
!>    @copyright Apache 2 license
!>
      SUBROUTINE FREE_GHOST_STRUCT(ghosts)
      USE EIKONAL3D_TYPES, ONLY : ghostType
      IMPLICIT NONE
      TYPE(ghostType), INTENT(INOUT), ALLOCATABLE :: ghosts(:)
      INTEGER i
      IF (ALLOCATED(ghosts)) THEN
         DO 1 i=1,SIZE(ghosts)
            IF (ALLOCATED(ghosts(i)%send_buff))  DEALLOCATE(ghosts(i)%send_buff)
            IF (ALLOCATED(ghosts(i)%recv_buff))  DEALLOCATE(ghosts(i)%recv_buff)
            IF (ALLOCATED(ghosts(i)%isend_dest)) DEALLOCATE(ghosts(i)%isend_dest)
            IF (ALLOCATED(ghosts(i)%irecv_dest)) DEALLOCATE(ghosts(i)%irecv_dest)
            ghosts(i)%nsend = 0
            ghosts(i)%nrecv = 0
            ghosts(i)%mysrc = 0
            ghosts(i)%mydest = 0
            ghosts(i)%send_request = 0 
            ghosts(i)%recv_request = 0
            ghosts(i)%lrecv = .FALSE.
    1    CONTINUE
         DEALLOCATE(ghosts)
      ENDIF
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Updates all nodes in a level
!>
!>    @param[in] nx      number of x grid points in domain
!>    @param[in] ny      number of y grid points in domain
!>    @param[in] nz      number of z grid points in domain
!>    @param[in] lrevx   if true then the x index is to be reversed (ix <- nx + 1 - ix)
!>    @param[in] lrevy   if true then the y index is to be reversed (iy <- ny + 1 - iy)
!>    @param[in] lrevz   if true then the z index is to be reversed (iz <- nz + 1 - iz)
!>    @param[in] lupd    if true then the i'th node will be updated [nx*ny*nz]
!>    @param[in] slow    slowness (s/m) at each grid points [nx*ny*nz]
!>
!>    @param[inout] u    on input these are the traveltimes (s) at each grid point in
!>                       the level.
!>                       on output these are the new traveltimes at each grid point 
!>                       in the level [nx*ny*nz]
!>
!>    @param[out] ierr   0 indicates success
!>
!>    @author Ben Baker
!>
!>    @copyright Apache 2 license
!>
      SUBROUTINE EVAL_UPDATE3D(nx, ny, nz,          &
                               lrevx, lrevy, lrevz, &
                               h, nnl,  ixyz_level, &
                               lupd, slow,          &
                               u, ierr)
      USE EIKONAL3D_MODULE, ONLY : UPDATE3D
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: h
      INTEGER, INTENT(IN) :: nnl, nx, ny, nz
      LOGICAL, INTENT(IN) :: lrevx, lrevy, lrevz
      DOUBLE PRECISION, INTENT(IN), DIMENSION(:) :: slow
      INTEGER, INTENT(IN), DIMENSION(:) :: ixyz_level
      LOGICAL, INTENT(IN), DIMENSION(:) :: lupd
      DOUBLE PRECISION, INTENT(INOUT), DIMENSION(:) :: u
      INTEGER, INTENT(OUT) :: ierr
      ! local variables
      INTEGER ierr_all, ierr1, ijk, ip, ix, iy, iz
      ierr_all = 0
!$OMP PARALLEL DO DEFAULT(none) &
!$OMP PRIVATE(ijk, ierr1, ix, iy, iz) &
!$OMP SHARED(h, ixyz_level, lrevx, lrevy, lrevz, lupd, nnl, nx, ny, nz, slow, u) &
!$OMP REDUCTION(+:ierr_all)
      DO 1 ip=1,nnl
         ix = ixyz_level(3*(ip-1)+1)
         iy = ixyz_level(3*(ip-1)+2)
         iz = ixyz_level(3*(ip-1)+3)
         IF (lrevx) ix = nx + 1 - ix
         IF (lrevy) iy = ny + 1 - iy
         IF (lrevz) iz = nz + 1 - iz
         ijk = (iz - 1)*nx*ny + (iy - 1)*nx + ix
         IF (lupd(ijk)) THEN
            CALL UPDATE3D(nx, ny, nz, ix, iy, iz, h, slow, u, ierr1)
            ierr_all = ierr_all + ierr1
         ENDIF 
    1 CONTINUE ! loop on levels
      ierr = ierr_all
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE UPDATE3D(nx, ny, nz, ix, iy, iz, h, slow, u, ierr)
      USE EIKONAL3D_MODULE, ONLY : GET_UXMIN3D, GET_UYMIN3D, GET_UZMIN3D, &
                                   SOLVE_HAMILTONIAN3D
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: slow
      INTEGER, INTENT(IN) :: ix, iy, iz, nx, ny, nz
      DOUBLE PRECISION, INTENT(IN) :: h
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: u
      INTEGER, INTENT(OUT) :: ierr
      DOUBLE PRECISION fijkh, ubar, uxmin, uymin, uzmin
      INTEGER ijk
      ierr = 0
      ijk = (iz - 1)*nx*ny + (iy - 1)*nx + ix
      fijkh = slow(ijk)*h
      uxmin = GET_UXMIN3D(nx, ny, ix, iy, iz, u)
      uymin = GET_UYMIN3D(nx, ny, ix, iy, iz, u)
      uzmin = GET_UZMIN3D(nx, ny, nz, ix, iy, iz, u)
      ubar = SOLVE_HAMILTONIAN3D(uxmin, uymin, uzmin, fijkh, ierr)
      u(ijk) = MIN(u(ijk), ubar)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      DOUBLE PRECISION FUNCTION GET_UXMIN3D(nx, ny, ix, iy, iz, u)
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: u
      INTEGER, INTENT(IN) :: ix, iy, iz, nx, ny 
      INTEGER ijk, im, ip
      ijk = (iz - 1)*nx*ny + (iy - 1)*nx + ix
      im = ijk - 1
      ip = ijk + 1
      IF (ix > 1 .AND. ix < nx) THEN
         get_uxmin3d = MIN(u(im), u(ip))
      ELSE
         IF (ix == 1) THEN
            get_uxmin3d = MIN(u(ijk), u(ip))
         ELSE
            get_uxmin3d = MIN(u(im), u(ijk))
         ENDIF 
      ENDIF
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      DOUBLE PRECISION FUNCTION GET_UYMIN3D(nx, ny, ix, iy, iz, u)
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: u
      INTEGER, INTENT(IN) :: ix, iy, iz, nx, ny
      INTEGER ijk, jm, jp
      ijk = (iz - 1)*nx*ny + (iy - 1)*nx + ix
      jm = ijk - nx
      jp = ijk + nx
      IF (iy > 1 .AND. iy < ny) THEN
         get_uymin3d = MIN(u(jm), u(jp))
      ELSE
         IF (iy == 1) THEN
            get_uymin3d = MIN(u(ijk), u(jp))
         ELSE
            get_uymin3d = MIN(u(jm), u(ijk))
         ENDIF
      ENDIF 
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      DOUBLE PRECISION FUNCTION GET_UZMIN3D(nx, ny, nz, ix, iy, iz, u)
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: u
      INTEGER, INTENT(IN) :: ix, iy, iz, nx, ny, nz
      INTEGER ijk, km, kp
      ijk = (iz - 1)*nx*ny + (iy - 1)*nx + ix
      km = ijk - nx*ny
      kp = ijk + nx*ny
      IF (iz > 1 .AND. iz < nz) THEN
         get_uzmin3d = MIN(u(km), u(kp))
      ELSE
         IF (iz == 1) THEN
            get_uzmin3d = MIN(u(ijk), u(kp))
         ELSE
            get_uzmin3d = MIN(u(km), u(ijk))
         ENDIF
      ENDIF
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sorts three numbers, a, b, c, into ascending order a1, a2, a3
!>
!>    @param[in] a    first number to sort
!>    @param[in] b    second number to sort
!>    @param[in] c    third number to sort
!>
!>    @param[out] a1  smallest of a, b, and c 
!>    @param[out] a2  intermediate value of a, b, and c
!>    @param[out] a3  largest of a, b, and c
!>
!>    @copyright Apache 2 license
!>
      SUBROUTINE SORT3(a, b, c, a1, a2, a3)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: a, b, c
      DOUBLE PRECISION, INTENT(OUT) :: a1, a2, a3
      LOGICAL lab, lac, lbc
      !----------------------------------------------------------------------------------!
      !
      ! compute comparisons
      lab = .TRUE.
      lac = .TRUE.
      lbc = .TRUE.
      IF (a > b) lab = .FALSE.
      IF (a > c) lac = .FALSE.
      IF (b > c) lbc = .FALSE.
      ! sort input - a is the smallest 
      IF (lab .AND. lac) THEN !a <= b .AND. a <= c) THEN
         ! a <= b <= c
         IF (lbc) THEN !b <= c) THEN
            a1 = a
            a2 = b
            a3 = c
         ELSE ! a <= c <= b
            a1 = a
            a2 = c
            a3 = b
         ENDIF
      ! b is the smallest
      ELSEIF (.NOT. lab .AND. lbc) THEN !b <= a .AND. b <= c) THEN
         ! b <= a <= c
         IF (lac) THEN !a <= c) THEN
            a1 = b
            a2 = a
            a3 = c
         ELSE ! b <=c <= a
            a1 = b
            a2 = c
            a3 = a
         ENDIF
      ! c is the smallest
      ELSE
         ! c <= a <= b
         IF (lab) THEN !a <= b) THEN
            a1 = c
            a2 = a
            a3 = b
         ELSE ! c <= b <= a
            a1 = c
            a2 = b
            a3 = a
         ENDIF 
      ENDIF 
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Solves the 2D Hamiltonian (Equation 2.4 of Zhao 2004)
!>
!>    @author Ben Baker
!>
!>    @copyright Apache 2 license
!>
      DOUBLE PRECISION FUNCTION SOLVE_HAMILTONIAN2D(a, b, fijh)
      USE EIKONAL3D_MODULE, ONLY : half, two 
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: a, b, fijh
      DOUBLE PRECISION xbar, amb, arg
      amb = a - b
      IF (ABS(amb) >= fijh) THEN
         xbar = MIN(a, b) + fijh
      ELSE
         arg  =  two*fijh*fijh - amb*amb
         xbar = half*(a + b + SQRT(arg))
      ENDIF 
      solve_hamiltonian2D = xbar
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        ! 
!>    @brief Solves the 3D Hamiltonian (Equation 2.5 and 2.6 of Zhao 2004.)
!>
!>    @author Ben Baker
!>
!>    @copyright Apache 2 license
!>
      DOUBLE PRECISION FUNCTION SOLVE_HAMILTONIAN3D(a, b, c, fijkh, ierr)
      USE EIKONAL3D_MODULE, ONLY : SOLVE_HAMILTONIAN2D, four, half, u_nan, &
                                   third, two_third, zero
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: a, b, c, fijkh
      INTEGER, INTENT(OUT) :: ierr
      DOUBLE PRECISION xtilde, a1, a2, a3, qb, qc, disc
      DOUBLE PRECISION, PARAMETER :: a4 = u_nan
      ierr = 0
      solve_hamiltonian3d = u_nan
      ! Order a, b, and c
      CALL SORT3(a, b, c, a1, a2, a3)
      ! Check on type of estimate solution we could find before working
      IF (a1 == u_nan) RETURN
      ! p == 1
      xtilde = a1 + fijkh
      IF (xtilde <= a2) THEN
         solve_hamiltonian3d = xtilde
         RETURN
      ELSE
         ! Set p == 2 and solve (x - a1)^2 + (x - a2)^2 = fij^2 h^2
         xtilde = SOLVE_HAMILTONIAN2D(a1, a2, fijkh)
         IF (xtilde <= a3) THEN
            solve_hamiltonian3d = xtilde
            RETURN
         ELSE
            ! Set p == 3 and solve 
            ! (x - a1)^2 + (x - a2)^2 + (x - a3)^2 = fij^2 h^2 
            ! 3x**2 - 2x(a1 + a2 + a3) + a1^2 + a2^2 + a3^3 - fij^2 h^2 = 0
            qb =-two_third*(a1 + a2 + a3)
            qc = (a1*a1 + a2*a2 + a3*a3 - fijkh*fijkh)*third
            disc = qb*qb  - four*qc
            ! Complex root
            IF (disc < zero) ierr = 1 
            xtilde = half*(-qb + SQRT(disc))
            ! Negative traveltime
            IF (xtilde < zero) ierr = 2
            IF (xtilde < a4) THEN
               solve_hamiltonian3d = xtilde
               RETURN
            ENDIF
         ENDIF
      ENDIF
      ierr = 3
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE EIKONAL_SOURCE_INDEX(nx, x0, dx, xs, isx)
      USE EIKONAL3D_MODULE, ONLY : half
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: x0, dx, xs
      INTEGER, INTENT(IN) :: nx
      INTEGER, INTENT(OUT) :: isx
      IF (xs <= x0) THEN
         isx = 1
      ELSEIF (xs >= x0 + FLOAT(nx - 1)*dx) THEN
         isx = nx
      ELSE
         isx = INT( (xs - x0)/dx + half ) + 1 
      ENDIF
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !

      SUBROUTINE EIKONAL_INIT_GRID(nx, isx, x0, dx, xs, ixloc, ierr)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: x0, dx, xs
      INTEGER, INTENT(IN) :: nx, isx
      INTEGER, INTENT(OUT) :: ixloc(3), ierr
      DOUBLE PRECISION xs_est
      INTEGER i, npinit
      ierr = 0
      ixloc(:) =-1
      xs_est = x0 + FLOAT(isx - 1)*dx
      IF (xs_est > xs) THEN 
         npinit = 2
         ixloc(1) = isx - 1
         ixloc(2) = isx
      ELSE IF (xs_est < xs) THEN 
         npinit = 2
         ixloc(1) = isx
         ixloc(2) = isx + 1
      ELSE
         npinit = 0
         IF (isx > 0) THEN
            npinit = npinit + 1
            ixloc(npinit) = isx - 1
         ENDIF
         npinit = npinit + 1
         ixloc(npinit) = isx
         IF (isx < nx - 1) THEN
            npinit = npinit + 1
            ixloc(npinit) = isx + 1
         ENDIF 
      ENDIF
      ! make sure source is in bounds
      DO 1 i=1,npinit
         IF (ixloc(i) < 1 .OR. ixloc(i) > nx) THEN
            WRITE(*,*) 'Source initialization error'
            ierr = 1
         ENDIF 
    1 CONTINUE
      RETURN
      END

!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Sets the boundary conditions (sources) and traveltimes in an around
!>           the source
      SUBROUTINE EIKONAL3D_SETBCS(nx, ny, nz, nsrc,       &
                                  dx, dy, dz, x0, y0, z0, &
                                  ts, xs, ys, zs, slow,   &
                                  lisbc, u, ierr)
      USE EIKONAL3D_MODULE, ONLY : EIKONAL_INIT_GRID, EIKONAL_SOURCE_INDEX, u_nan
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: slow, ts, xs, ys, zs
      DOUBLE PRECISION, INTENT(IN) :: dx, dy, dz, x0, y0, z0
      INTEGER, INTENT(IN) :: nx, ny, nz, nsrc
      DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: u
      LOGICAL, DIMENSION(:), INTENT(OUT) :: lisbc
      INTEGER, INTENT(OUT) :: ierr
      ! local variables
      INTEGER, ALLOCATABLE :: isx(:), isy(:), isz(:)
      DOUBLE PRECISION d, x, y, z
      INTEGER ixloc(3), iyloc(3), izloc(3), isrc, i, ijk, ix, iy, iz, j, k
      ALLOCATE(isx(nsrc))
      ALLOCATE(isy(nsrc))
      ALLOCATE(isz(nsrc))
      ierr = 0
      lisbc(:) = .FALSE.
      u(:) = u_nan
      ! locate closest source in grid
      DO 1 isrc=1,nsrc
         CALL EIKONAL_SOURCE_INDEX(nx, x0, dx, xs(isrc), isx(isrc))
         CALL EIKONAL_SOURCE_INDEX(ny, y0, dy, ys(isrc), isy(isrc))
         CALL EIKONAL_SOURCE_INDEX(nz, z0, dz, zs(isrc), isz(isrc))
    1 CONTINUE
      ! compute the times around the source
      ixloc(:) =-1
      iyloc(:) =-1
      izloc(:) =-1 
      DO 11 isrc=1,nsrc
         CALL EIKONAL_INIT_GRID(nx, isx(isrc), x0, dx, xs(isrc), ixloc, ierr)
         IF (ierr /= 0) THEN
            WRITE(*,*) 'eikonal3d_setbcs: Error setting ixloc'
            GOTO 500
         ENDIF
         CALL EIKONAL_INIT_GRID(ny, isy(isrc), y0, dy, ys(isrc), iyloc, ierr)
         IF (ierr /= 0) THEN
            WRITE(*,*) 'eikonal3d_setbcs: Error setting iyloc'
            GOTO 500
         ENDIF
         CALL EIKONAL_INIT_GRID(nz, isz(isrc), z0, dz, zs(isrc), izloc, ierr)
         IF (ierr /= 0) THEN
            WRITE(*,*) 'eikonal3d_setbcs: Error setting izloc'
            GOTO 500
         ENDIF
         DO 12 i=1,3
            IF (ixloc(i) ==-1) CYCLE
            DO 13 j=1,3
               IF (iyloc(j) ==-1) CYCLE
               DO 14 k=1,3
                  IF (izloc(k) ==-1) CYCLE
                  ix = ixloc(i)
                  iy = iyloc(j)
                  iz = izloc(k)
                  ijk = (iz - 1)*nx*ny + (iy - 1)*nx + ix
                  x = x0 + FLOAT(ix - 1)*dx
                  y = y0 + FLOAT(iy - 1)*dy
                  z = z0 + FLOAT(iz - 1)*dz
                  d = SQRT( (xs(isrc) - x)**2 + (ys(isrc) - y)**2 + (zs(isrc) - z)**2 )
                  ! collocate
                  IF (ABS(d) < 1.d-10) THEN
                     u(ijk) = ts(isrc) + d*slow(ijk)
                  ELSE ! take minimum
                     u(ijk) = MIN(u(ijk), ts(isrc) + d*slow(ijk))
                  ENDIF
                  lisbc(ijk) = .TRUE.
   14          CONTINUE 
   13       CONTINUE
   12    CONTINUE
   11 CONTINUE
  500 CONTINUE
      IF (ALLOCATED(isx)) DEALLOCATE(isx)
      IF (ALLOCATED(isy)) DEALLOCATE(isy)
      IF (ALLOCATED(isz)) DEALLOCATE(isz)
      RETURN
      END
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
      SUBROUTINE EIKONAL_GRD2IJK(igrd, nx, ny, nz, i, j, k, ierr)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: igrd, nx, ny, nz
      INTEGER, INTENT(OUT) :: i, j, k, ierr
      INTEGER igrd1, nxy
      ierr = 0
      igrd1 = igrd - 1
      nxy = nx*ny
      k = igrd1/(nxy)
      j = (igrd1 - k*nxy)/nx
      i =  igrd1 - k*nxy - j*nx
      i = i + 1
      j = j + 1
      k = k + 1
      IF (i < 1 .OR. i > nx) ierr = ierr + 1
      IF (j < 1 .OR. j > ny) ierr = ierr + 1
      IF (k < 1 .OR. k > nz) ierr = ierr + 1 
      IF ((k - 1)*nxy + (j - 1)*nx + i /= igrd) ierr = ierr + 1
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Exchanges values in the ghosts communication region.  This is a blocking
!>           version and should only be used for testing.
!>
!>    @param[in] comm       MPI communicator
!>
!>    @param[inout] ghosts  on input contains the sending and recieving information
!>                          and preallocated workspace for this process to all other
!>                          processes in communicator.
!>                          on output the buffers have been filled but they should
!>                          not be accessed by another subroutine. [mpi communicator size]
!>    @param[inout] u       on input contains the nodes this process has updated.
!>                          on output contains the ndoes this process needs and others
!>                          have updated.
!>
!>    @author Ben Baker
!>
!>    @copyright Apache 2 license
!>
      SUBROUTINE EIKONAL_EXCHANGE_BLOCKING(comm, ghosts, u)
      USE EIKONAL3D_TYPES, ONLY : ghostType
      USE MPI
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: comm
      TYPE(ghostType), DIMENSION(:), INTENT(INOUT) :: ghosts
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: u
      ! local variables
      INTEGER stat(MPI_STATUS_SIZE), i, j, k, mpierr, myid, nblocks 
      !----------------------------------------------------------------------------------!
      !
      ! get mpi information 
      CALL MPI_COMM_SIZE(comm, nblocks, mpierr)
      CALL MPI_COMM_RANK(comm, myid, mpierr)
      DO 11 i=1,nblocks
         IF (i - 1 == myid) THEN
            DO 12 j=1,nblocks
               IF (i == j) CYCLE ! don't talk to myself
               ! require there be something to send 
               IF (ghosts(j)%nsend > 0) THEN
                  ! extract elements of u onto buffer 
                  DO 13 k=1,ghosts(j)%nsend
                     ghosts(j)%send_buff(k) = u(ghosts(j)%isend_dest(k))
   13             CONTINUE
                  ! send it
                  CALL MPI_SEND(ghosts(j)%send_buff, ghosts(j)%nsend,         &
                                MPI_DOUBLE_PRECISION, ghosts(j)%mydest, myid, &
                                comm, mpierr)
               ENDIF
   12       CONTINUE
         ELSE
            ! receive data (if i'm expecting it) 
            IF (ghosts(i)%nrecv > 0) THEN
               CALL MPI_RECV(ghosts(i)%recv_buff, ghosts(i)%nrecv,   &
                             MPI_DOUBLE_PRECISION, ghosts(i)%mysrc,  &
                             MPI_ANY_TAG, comm, stat, mpierr)
            ENDIF
         ENDIF
         CALL MPI_BARRIER(comm, mpierr)
   11 CONTINUE
      DO 21 i=1,nblocks
         IF (ghosts(i)%nrecv > 0) THEN
            DO 22 k=1,ghosts(i)%nrecv
               u(ghosts(i)%irecv_dest(k)) = ghosts(i)%recv_buff(k)
   22       CONTINUE
         ENDIF
   21 CONTINUE
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Exchanges values in the ghosts communication region
!>
!>    @param[in] comm       MPI communicator
!>
!>    @param[inout] ghosts  on input contains the sending and recieving information
!>                          and preallocated workspace for this process to all other
!>                          processes in communicator.
!>                          on output the buffers have been filled but they should
!>                          not be accessed by another subroutine. [mpi communicator size]
!>    @param[inout] u       on input contains the nodes this process has updated.
!>                          on output contains the ndoes this process needs and others
!>                          have updated.
!>
!>    @author Ben Baker
!>
!>    @copyright Apache 2 license
!>
      SUBROUTINE EIKONAL_EXCHANGE(comm, ghosts, u)
      USE EIKONAL3D_TYPES, ONLY : ghostType
      USE MPI
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: comm
      TYPE(ghostType), DIMENSION(:), INTENT(INOUT) :: ghosts
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: u
      ! local variables
      INTEGER stat(MPI_STATUS_SIZE), i, irecv, j, k, mpierr, myid, nblocks, nsend, nrecv
      !----------------------------------------------------------------------------------!
      !
      ! get mpi information 
      CALL MPI_COMM_SIZE(comm, nblocks, mpierr)
      CALL MPI_COMM_RANK(comm, myid, mpierr)
      ! initialize 
      DO 1 i=1,nblocks
         ghosts(i)%send_request =-1
         ghosts(i)%recv_request =-1
         ghosts(i)%lrecv = .FALSE.
    1 CONTINUE
      ! send out buffers for sending and receiver
      nsend = 0
      nrecv = 0
      DO 11 i=1,nblocks
         IF (i - 1 == myid) THEN
            DO 12 j=1,nblocks
               IF (i == j) CYCLE ! don't talk to myself
               ! require there be something to send 
               IF (ghosts(j)%nsend > 0) THEN
                  ! extract elements of u onto buffer 
                  DO 13 k=1,ghosts(j)%nsend
                     ghosts(j)%send_buff(k) = u(ghosts(j)%isend_dest(k))
   13             CONTINUE
                  ! send it
                  nsend = nsend + 1
                  CALL MPI_ISEND(ghosts(j)%send_buff, ghosts(j)%nsend,         &
                                 MPI_DOUBLE_PRECISION, ghosts(j)%mydest, myid, &
                                 comm, ghosts(j)%send_request, mpierr)
               ENDIF
   12       CONTINUE
         ELSE
            ! receive data (if i'm expecting it) 
            IF (ghosts(i)%nrecv > 0) THEN
               nrecv = nrecv + 1
               CALL MPI_IRECV(ghosts(i)%recv_buff, ghosts(i)%nrecv,             &
                              MPI_DOUBLE_PRECISION, ghosts(i)%mysrc,            &
                              MPI_ANY_TAG, comm, ghosts(i)%recv_request, mpierr)
            ENDIF
         ENDIF
   11 CONTINUE
      ! as i receive data unpack it 
      irecv = 0
      DO WHILE (irecv < nrecv)
         DO 21 i=1,nblocks
            IF (ghosts(i)%nrecv > 0 .AND. .NOT.ghosts(i)%lrecv) THEN
               CALL MPI_TEST(ghosts(i)%recv_request, ghosts(i)%lrecv, stat, mpierr)
               IF (ghosts(i)%lrecv) THEN
                  DO 22 k=1,ghosts(i)%nrecv
                     u(ghosts(i)%irecv_dest(k)) = ghosts(i)%recv_buff(k)
   22             CONTINUE
                  ghosts(i)%lrecv = .TRUE.
                  irecv = irecv + 1
               ENDIF 
            ENDIF
   21    CONTINUE 
      ENDDO
      ! wait until i've sent everything before continuiing
      DO 31 i=1,nblocks
         IF (ghosts(i)%nsend > 0) THEN
            CALL MPI_WAIT(ghosts(i)%send_request, stat, mpierr)
         ENDIF
   31 CONTINUE
      CALL MPI_BARRIER(comm, mpierr)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE EIKONAL3D_GHOST_COMM(comm, ndivx, ndivy, ndivz, &
                                      noverlap, nx, ny, nz, &
                                      ghosts, model_loc, ierr )
      USE MPIUTILS_MODULE, ONLY : MPIUTILS_GRD2IJK
      USE EIKONAL3D_MODULE, ONLY : EIKONAL_GRD2IJK
      USE EIKONAL3D_TYPES, ONLY : ghostType
      USE EIKONAL3D_TYPES, ONLY : localModelType
      USE MPI
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: comm, ndivx, ndivy, ndivz, noverlap, nx, ny, nz
      TYPE(ghostType), INTENT(OUT), ALLOCATABLE :: ghosts(:)
      TYPE (localModelType), INTENT(OUT) :: model_loc
      INTEGER, INTENT(OUT) :: ierr
      INTEGER, ALLOCATABLE :: ibuff(:), need_nodes(:), need_work(:), norigin(:), &
                              send_work(:)
      INTEGER stat(MPI_STATUS_SIZE), i, i1, i2, imbx, imby, imbz, ineed, isend, &
              ix, ix1, ix2, iy, iy1, iy2, iz, iz1, iz2, j, j1, j2, k, k1, k2, &
              nbcast, nblocks, ndx, ndy, ndz, nneed, nsend, nupd_loc, nupd_tot, &
              nwork, nx_loc, nxyz_loc, ny_loc, nz_loc, myid, mpierr
      INTEGER, PARAMETER :: master = 0
      CALL MPI_COMM_SIZE(comm, nblocks, mpierr)
      CALL MPI_COMM_RANK(comm, myid, mpierr)
      IF (ALLOCATED(ghosts)) DEALLOCATE(ghosts)
      ALLOCATE(ghosts(nblocks))
      !CALL EIKONAL_GRD2IJK(myid+1, ndivx, ndivy, ndivz, imbx, imby, imbz, ierr)
      CALL MPIUTILS_GRD2IJK(myid, ndivx, ndivy, ndivz, imbx, imby, imbz, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'Error finding process block', myid
         RETURN
      ENDIF
      DO 1 i=1,nblocks
         ghosts(i)%nrecv = 0
         ghosts(i)%nsend = 0
    1 CONTINUE 
!     imbx = imbx - 1
!     imby = imby - 1
!     imbz = imbz - 1
      ! now compute the nodes i need from others
      ndx = MAX(nx/ndivx, 1)
      ndy = MAX(ny/ndivy, 1)
      ndz = MAX(nz/ndivz, 1)
      i1 = ndx*imbx + 1
      i2 = ndx*(imbx + 1)
      j1 = ndy*imby + 1
      j2 = ndy*(imby + 1)
      k1 = ndz*imbz + 1
      k2 = ndz*(imbz + 1)
      IF (imbx + 1 == ndivx) i2 = nx
      IF (imby + 1 == ndivy) j2 = ny
      IF (imbz + 1 == ndivz) k2 = nz
      ix1 = MAX(1,  i1 - noverlap)
      ix2 = MIN(nx, i2 + noverlap)
      iy1 = MAX(1,  j1 - noverlap)
      iy2 = MIN(ny, j2 + noverlap)
      iz1 = MAX(1,  k1 - noverlap)
      iz2 = MIN(nz, k2 + noverlap)
      nupd_loc = (i2 - i1 + 1)*(j2 - j1 + 1)*(k2 - k1 + 1)
!print *, ierr, imbx, imby, imbz, nupd_loc, nx*ny*nz
      CALL MPI_ALLREDUCE(nupd_loc, nupd_tot, 1, MPI_INTEGER, MPI_SUM, comm, mpierr)
      IF (nupd_tot /= nx*ny*nz) THEN
         WRITE(*,*) 'Error computing local domain'
         ierr = 1
         RETURN
      ENDIF 
      nx_loc = ix2 - ix1 + 1
      ny_loc = iy2 - iy1 + 1
      nz_loc = iz2 - iz1 + 1
      nxyz_loc = nx_loc*ny_loc*nz_loc
      nneed = 6*MAX(nx*ny, nx*nz, ny*nz)
      CALL MPI_ALLREDUCE(nneed, nwork, 1, MPI_INTEGER, MPI_MAX, comm, mpierr) 
      ! initialize the local model
      model_loc%nx = nx_loc
      model_loc%ny = ny_loc
      model_loc%nz = nz_loc
      model_loc%ix1 = ix1
      model_loc%iy1 = iy1
      model_loc%iz1 = iz1
      model_loc%myblock = myid + 1
      ALLOCATE(model_loc%l2g_node(nxyz_loc))
      ALLOCATE(model_loc%lghost_node(nxyz_loc))
      model_loc%l2g_node(:) = 0
      model_loc%lghost_node(:) = .FALSE.
      ! now get the global coordinates of the nodes i need
      ALLOCATE(need_nodes(nneed))
      need_nodes(:) = 0
      ineed = 0
      DO 2 iz=iz1,iz2
         DO 3 iy=iy1,iy2
            DO 4 ix=ix1,ix2
               i = (iz - iz1)*nx_loc*ny_loc &
                 + (iy - iy1)*nx_loc        &
                 + (ix - ix1) + 1
               model_loc%l2g_node(i) = (iz - 1)*nx*ny + (iy - 1)*nx + ix
               IF ((ix >= i1 .AND. ix <= i2) .AND. &
                   (iy >= j1 .AND. iy <= j2) .AND. &
                   (iz >= k1 .AND. iz <= k2)) CYCLE
               ineed = ineed + 1
               need_nodes(ineed) = (iz - 1)*nx*ny + (iy - 1)*nx + ix
               model_loc%lghost_node(i) = .TRUE.
    4       CONTINUE 
    3    CONTINUE
    2 CONTINUE
      nneed = ineed
      ! verify the local to global map makes sense
      IF (MINVAL(model_loc%l2g_node) < 1 .OR.       &
          MAXVAL(model_loc%l2g_node) > nx*ny*nz) THEN
         WRITE(*,*) 'local to global node map is wrong'
         ierr = 1
         RETURN
      ENDIF
      ! verify need_nodes is sorted
      DO 5 i=2,nneed
         IF (need_nodes(i-1) >= need_nodes(i)) THEN
            WRITE(*,*) 'nodes not sorted'
            ierr = 1
            RETURN
         ENDIF
    5 CONTINUE
      ALLOCATE(need_work(nwork))
      ALLOCATE(send_work(nwork))
      ALLOCATE(ibuff(nblocks))
      ALLOCATE(norigin(nblocks))
      ! loop on processes and tell others what i need
      DO 21 i=1,nblocks
         need_work(:) = 0
         IF (i - 1 == myid) THEN
            nbcast = nneed
            need_work(1:nbcast) = need_nodes(1:nbcast) 
         ENDIF
         CALL MPI_BCAST(nbcast,         1, MPI_INTEGER, i-1, comm, mpierr)
         CALL MPI_BCAST(need_work, nbcast, MPI_INTEGER, i-1, comm, mpierr)
         ! should i be sending this node?
         isend = 0
         ibuff(:) = 0
         IF (myid /= i - 1) THEN
            DO 22 k=1,nbcast
               CALL EIKONAL_GRD2IJK(need_work(k), nx, ny, nz, ix, iy, iz, ierr)
               IF (ierr /= 0) WRITE(*,*) 'grd2ijk mistake'
               IF (ix >= i1 .AND. ix <= i2 .AND. &
                   iy >= j1 .AND. iy <= j2 .AND. &
                   iz >= k1 .AND. iz <= k2) THEN
                  isend = isend + 1
                  send_work(isend) = need_work(k) !(iz - iz1)*nx_loc*ny_loc &
!                                  + (iy - iy1)*nx_loc &
!                                  + (ix - ix1 + 1)
               ENDIF
  22        CONTINUE 
         ENDIF
         CALL MPI_ALLREDUCE(isend, nsend, 1, MPI_INTEGER, MPI_SUM, comm, mpierr)
         IF (nsend /= nneed) THEN
            WRITE(*,*) 'Some nodes do not overlap!'
            ierr = 1
            RETURN
         ENDIF
         ! let myid know who is sending how much 
         ibuff(myid+1) = isend
         norigin(:) = 0
!print *, ibuff
         CALL MPI_REDUCE(ibuff, norigin, nblocks, MPI_INTEGER, MPI_SUM, &
                         i-1, comm, mpierr) 
         ! if i'm sending update my ghosts
!        ghosts(i)%nsend = isend
!        IF (ghosts(i)%nsend > 0) THEN

         IF (myid == i - 1) THEN
            ! get the send destinations from others 
            DO j=1,nblocks
               IF (norigin(j) < 1) CYCLE
               IF (myid == j - 1) THEN
                  WRITE(*,*) 'Should not be here'
               ENDIF
               ghosts(j)%nrecv = norigin(j)
               ghosts(j)%mysrc = j - 1
               ALLOCATE(ghosts(j)%recv_buff(ghosts(j)%nrecv))
               ALLOCATE(ghosts(j)%irecv_dest(ghosts(j)%nrecv))
               ghosts(j)%recv_buff(:) = 0.d0
               ghosts(j)%irecv_dest(:) = 0
               CALL MPI_RECV(ghosts(j)%irecv_dest, norigin(j), MPI_INTEGER, j - 1, &
                             MPI_ANY_TAG, comm, stat, mpierr)
               ! map the global receive nodes to local nodes
               DO k=1,ghosts(j)%nrecv
                  CALL EIKONAL_GRD2IJK(ghosts(j)%irecv_dest(k), nx, ny, nz, &
                                       ix, iy, iz, ierr)
                  IF (ierr /= 0) THEN
                     WRITE(*,*) 'Error fixing recv nodes' 
                     RETURN
                  ENDIF
                  ghosts(j)%irecv_dest(k) = (iz - iz1)*nx_loc*ny_loc &
                                          + (iy - iy1)*nx_loc        &
                                          + (ix - ix1) + 1
                  IF (ghosts(j)%irecv_dest(k) < 1 .OR. &
                      ghosts(j)%irecv_dest(k) > nxyz_loc) THEN
                     WRITE(*,*) 'Error receive destination is out of bounds'
                     ierr = 1
                     RETURN
                  ENDIF 
                  ! not a critical error but really should be ordered
                  IF (K > 1) THEN
                     IF (ghosts(j)%irecv_dest(k-1) >= &
                         ghosts(j)%irecv_dest(k)) PRINT *, 'out of order!'
                  ENDIF
               ENDDO
            ENDDO
         ELSE
            ghosts(i)%mydest = i - 1
            ghosts(i)%nsend = isend
            IF (ghosts(i)%nsend > 0) THEN
               ALLOCATE(ghosts(i)%send_buff(ghosts(i)%nsend))
               ALLOCATE(ghosts(i)%isend_dest(ghosts(i)%nsend))
               ghosts(i)%send_buff(:) = 0.d0
               ghosts(i)%isend_dest(:) = 0
               CALL MPI_SEND(send_work, isend, MPI_INTEGER, i - 1, &
                             myid, comm, mpierr)
               ! map the send nodes to local nodes 
               DO k=1,isend
                  CALL EIKONAL_GRD2IJK(send_work(k), nx, ny, nz, ix, iy, iz, ierr)
                  IF (ierr /= 0) THEN
                     WRITE(*,*) 'Error fixing recv nodes' 
                     RETURN
                  ENDIF
                  ghosts(i)%isend_dest(k) = (iz - iz1)*nx_loc*ny_loc &
                                          + (iy - iy1)*nx_loc        &
                                          + (ix - ix1) + 1
               ENDDO
            ENDIF
         ENDIF
  21  CONTINUE ! loop on blocks in domain
      DEALLOCATE(need_work)
      DEALLOCATE(send_work)
      DEALLOCATE(ibuff)
      DEALLOCATE(norigin)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Scatters the slowness model to all processes in the sweep group.
!>
!>    @param[in] comm       MPI communicator for sweep group 
!>    @param[in] l2g_node   maps from each process' local node to the global node
!>    @param[in] slow       global slowness model defined on the master node
!>
!>    @param[out] slow_loc  local model on each process
!>
!>    @author Ben Baker
!>
!>    @copyright Apache 2 license
!>
      SUBROUTINE EIKONAL_SCATTER_MODEL(comm, l2g_node, slow, slow_loc)
      USE MPI
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: slow
      INTEGER, DIMENSION(:), INTENT(IN) :: l2g_node
      INTEGER, INTENT(IN) :: comm
      DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: slow_loc
      ! local variables
      DOUBLE PRECISION, ALLOCATABLE :: buf(:)
      INTEGER, ALLOCATABLE :: idest(:)
      INTEGER stat(MPI_STATUS_SIZE), i, k, myid, mpierr, nblocks, &
              nbuf, nwork, nxyz_loc
      INTEGER, PARAMETER :: master = 0
      !----------------------------------------------------------------------------------!
      !
      ! get communicator information
      CALL MPI_COMM_RANK(comm, myid,    mpierr)
      CALL MPI_COMM_SIZE(comm, nblocks, mpierr)
      ! allocate buffer space on master 
      nxyz_loc = SIZE(l2g_node)
      CALL MPI_REDUCE(nxyz_loc, nwork, 1, MPI_INTEGER, MPI_MAX, master, comm, mpierr) 
      IF (myid == master) THEN
         ALLOCATE(idest(nwork))
         ALLOCATE(buf(nwork))
      ELSE
         ALLOCATE(idest(1))
         ALLOCATE(buf(1))
      ENDIF
      slow_loc(:) = 0.d0 ! initialize result
      ! loop on blocks and have master send model
      DO 1 i=1,nblocks
         IF (myid == master) THEN
            ! i'm simply sending to myself
            IF (i - 1 == myid) THEN
               DO 11 k=1,nxyz_loc
                  slow_loc(k) = slow(l2g_node(k))
   11          CONTINUE
            ELSE
               ! get the requisite coordinates
               CALL MPI_RECV(nbuf,  1,    MPI_INTEGER, i - 1, &
                             MPI_ANY_TAG, comm, stat, mpierr)
               CALL MPI_RECV(idest, nbuf, MPI_INTEGER, i - 1, &
                             MPI_ANY_TAG, comm, stat, mpierr)
               ! extract the model
               DO 12 k=1,nbuf
                  buf(k) = slow(idest(k))
   12          CONTINUE
               ! send back
               CALL MPI_SEND(buf, nbuf, MPI_DOUBLE_PRECISION, i - 1, &
                             myid, comm, mpierr) 
            ENDIF
         ELSE
            ! does the boss want to talk to me?
            IF (myid == i - 1) THEN
               ! send the master the indices of the locations i need
               CALL MPI_SEND(nxyz_loc,        1, MPI_INTEGER, master, &
                             myid, comm, mpierr) 
               CALL MPI_SEND(l2g_node, nxyz_loc, MPI_INTEGER, master, &
                             myid, comm, mpierr)
               ! receive the subset of the model i need 
               CALL MPI_RECV(slow_loc, nxyz_loc, MPI_DOUBLE_PRECISION, master, &
                             MPI_ANY_TAG, comm, stat, mpierr)
            ENDIF
         ENDIF
         CALL MPI_BARRIER(comm, mpierr) ! block until communication is done
    1 CONTINUE ! loop on blocks
      ! free space
      IF (ALLOCATED(idest)) DEALLOCATE(idest)
      IF (ALLOCATED(buf))   DEALLOCATE(buf)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Scatters the boundary conditions to all processes in the sweep group.
!>
!>    @param[in] comm       MPI communicator for sweep group 
!>    @param[in] lisbc      if true then the global node is a boundary condition
!>    @param[in] u          initial traveltimes at each node (with boundary conditions)
!>
!>    @param[out] slow_loc  local model on each process
!>    @param[out] u_loc     traveltimes on each process
!>
!>    @author Ben Baker
!>
!>    @copyright Apache 2 license
!>
      SUBROUTINE EIKONAL_SCATTER_BCS(comm, l2g_node, lisbc, u, lisbc_loc, u_loc)
      USE MPI
      USE EIKONAL3D_MODULE, ONLY : u_nan
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: u
      INTEGER, DIMENSION(:), INTENT(IN) :: l2g_node
      LOGICAL, DIMENSION(:), INTENT(IN) :: lisbc 
      INTEGER, INTENT(IN) :: comm
      DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: u_loc 
      LOGICAL, DIMENSION(:), INTENT(OUT) :: lisbc_loc
      ! local variables
      DOUBLE PRECISION, ALLOCATABLE :: buf(:)
      LOGICAL, ALLOCATABLE :: lbuf(:)
      INTEGER, ALLOCATABLE :: idest(:)
      INTEGER stat(MPI_STATUS_SIZE), i, k, myid, mpierr, nblocks, &
              nbuf, nwork, nxyz_loc
      INTEGER, PARAMETER :: master = 0
      !----------------------------------------------------------------------------------!
      !
      ! get communicator information
      CALL MPI_COMM_RANK(comm, myid,    mpierr)
      CALL MPI_COMM_SIZE(comm, nblocks, mpierr)
      ! allocate buffer space on master 
      nxyz_loc = SIZE(l2g_node)
      CALL MPI_REDUCE(nxyz_loc, nwork, 1, MPI_INTEGER, MPI_MAX, master, comm, mpierr)
      IF (myid == master) THEN
         ALLOCATE(idest(nwork))
         ALLOCATE(buf(nwork))
         ALLOCATE(lbuf(nwork))
      ELSE
         ALLOCATE(idest(1))
         ALLOCATE(buf(1))
         ALLOCATE(lbuf(1))
      ENDIF
      lisbc_loc(:) = .FALSE.
      u_loc(:) = u_nan
      ! loop on blocks and have master send model
      DO 1 i=1,nblocks
         IF (myid == master) THEN 
            ! i'm simply sending to myself
            IF (i - 1 == myid) THEN 
               DO 11 k=1,nxyz_loc
                  u_loc(k) = u(l2g_node(k))
                  lisbc_loc(k) = lisbc(l2g_node(k))
   11          CONTINUE
            ELSE
               ! get the requisite coordinates
               CALL MPI_RECV(nbuf,  1,    MPI_INTEGER, i - 1, & 
                             MPI_ANY_TAG, comm, stat, mpierr)
               CALL MPI_RECV(idest, nbuf, MPI_INTEGER, i - 1, & 
                             MPI_ANY_TAG, comm, stat, mpierr)
               ! extract the model
               DO 12 k=1,nbuf
                  buf(k) = u(idest(k))
                  lbuf(k) = lisbc(idest(k))
   12          CONTINUE
               ! send back
               CALL MPI_SEND(buf,  nbuf, MPI_DOUBLE_PRECISION, i - 1, & 
                             myid, comm, mpierr) 
               CALL MPI_SEND(lbuf, nbuf, MPI_LOGICAL,          i - 1, &
                             myid, comm, mpierr)
            ENDIF
         ELSE
            ! does the boss want to talk to me?
            IF (myid == i - 1) THEN 
               ! send the master the indices of the locations i need
               CALL MPI_SEND(nxyz_loc,        1, MPI_INTEGER, master, &
                             myid, comm, mpierr) 
               CALL MPI_SEND(l2g_node, nxyz_loc, MPI_INTEGER, master, &
                             myid, comm, mpierr)
               ! receive the subset of the boundary conditions i need 
               CALL MPI_RECV(u_loc,      nxyz_loc, MPI_DOUBLE_PRECISION, master, &
                             MPI_ANY_TAG, comm, stat, mpierr)
               CALL MPI_RECV(lisbc_loc, nxyz_loc, MPI_LOGICAL,           master, &
                             MPI_ANY_TAG, comm, stat, mpierr)
            ENDIF
         ENDIF
         CALL MPI_BARRIER(comm, mpierr) ! block until communication is done
    1 CONTINUE ! loop on blocks
      ! free space
      IF (ALLOCATED(idest)) DEALLOCATE(idest)
      IF (ALLOCATED(buf))   DEALLOCATE(buf)
      IF (ALLOCATED(lbuf))  DEALLOCATE(lbuf)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Gathers the traveltimes the slaves onto the master.  Note that while
!>           the domains overlap the overlapping segments should return identical
!>           numbers from each process in the group.
!>
!>    @param[in] comm       MPI communicator
!>    @param[in] l2g_node   maps from the local node number to the global node number
!>    @param[in] u_loc      traveltimes computed by process in local domain
!>
!>    @param[out] u         aggregrated traveltimes from all processes and collected
!>                          by the master 
!>
!>    @author Ben Baker
!>
!>    @copyright Apache 2 license
!>
      SUBROUTINE EIKONAL_GATHER_TRAVELTIMES(comm, l2g_node, u_loc, u)
      USE MPI
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: u_loc
      INTEGER, DIMENSION(:), INTENT(IN) :: l2g_node
      INTEGER, INTENT(IN) :: comm
      DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: u
      ! local variables
      DOUBLE PRECISION, ALLOCATABLE :: buf(:)
      INTEGER, ALLOCATABLE :: idest(:)
      INTEGER stat(MPI_STATUS_SIZE), i, k, mpierr, myid, nblocks, nbuf, nxyz_loc, nwork
      INTEGER, PARAMETER :: master = 0
      !----------------------------------------------------------------------------------!
      !
      ! get mpi information
      CALL MPI_COMM_RANK(comm, myid,    mpierr)
      CALL MPI_COMM_SIZE(comm, nblocks, mpierr)
      ! allocate buffer space on master 
      nxyz_loc = SIZE(l2g_node)
      CALL MPI_REDUCE(nxyz_loc, nwork, 1, MPI_INTEGER, MPI_MAX, master, comm, mpierr) 
      IF (myid == master) THEN 
         ALLOCATE(idest(nwork))
         ALLOCATE(buf(nwork))
      ELSE 
         ALLOCATE(idest(1))
         ALLOCATE(buf(1))
      ENDIF
      ! loop on blocks and have master receive traveltimes
      DO 1 i=1,nblocks
         IF (myid == master) THEN
            ! i'm receiving from myself
            IF (i - 1 == myid) THEN
               DO 11 k=1,nxyz_loc
                  u(l2g_node(k)) = u_loc(k)
   11          CONTINUE
            ELSE
               ! get the destination locations and the values from the process
               CALL MPI_RECV(nbuf,     1, MPI_INTEGER,          i - 1, &
                             MPI_ANY_TAG, comm, stat, mpierr)
               CALL MPI_RECV(idest, nbuf, MPI_INTEGER,          i - 1, &
                             MPI_ANY_TAG, comm, stat, mpierr)
               CALL MPI_RECV(buf, nbuf,   MPI_DOUBLE_PRECISION, i - 1, &
                             MPI_ANY_TAG, comm, stat, mpierr) 
               ! copy the result into u
               DO 12 k=1,nbuf
                  u(idest(k)) = buf(k)
   12          CONTINUE
            ENDIF
         ELSE
            ! my turn to send my answer to the boss
            IF (i - 1 == myid) THEN
               CALL MPI_SEND(nxyz_loc,        1, MPI_INTEGER,          master, &
                             myid, comm, mpierr)
               CALL MPI_SEND(l2g_node, nxyz_loc, MPI_INTEGER,          master, &
                             myid, comm, mpierr)
               CALL MPI_SEND(u_loc, nxyz_loc,    MPI_DOUBLE_PRECISION, master, & 
                             myid, comm, mpierr)
            ENDIF
         ENDIF
         CALL MPI_BARRIER(comm, mpierr)
    1 CONTINUE ! loop on blocks
      ! free space
      IF (ALLOCATED(idest)) DEALLOCATE(idest)
      IF (ALLOCATED(buf))   DEALLOCATE(buf)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Initializes the solver parameters from the ini file and the corresponding
!>           solver data structures.  Notice, the coordinate system is right handed 
!>           with z positive up (x positive right, and y positive away)
!>
!>    @param[in] comm       MPI communicator
!>    @param[in] iverb      controls verbosity (0 is quiet)
!>    @param[in] ndivx      number of blocks to divide domain up in x direction
!>    @param[in] ndivy      number of blocks to divide domain up in y direction
!>    @param[in] ndivz      number of blocks to divide domain up in z direction
!>    @param[in] noverlap   controls the ghost communication overlap (should be 1 
!>                          since i only have a first order finite diference)
!>    @param[in] maxit      max number of iterations in fast sweeping method
!>    @param[in] x0         x origin of model (m)
!>    @param[in] y0         y origin of model (m)
!>    @param[in] z0         z origin of model (m)
!>    @param[in] h          grid spacing in x, y, and z (m)
!>    @param[in] tol        convergence criteria in fast sweeping method.  if an 
!>                          update for all points results in a change less than
!>                          the given tolerance (seconds) the iterations will terminate.
!>
!>    @param[out] ierr      0 indicates success
!>
!>    @author Ben Baker
!>
!>    @copyright Apache 2 license
!>
      SUBROUTINE EIKONAL3D_INITIALIZE(comm, iverb, nx, ny, nz,    &
                                      ndivx, ndivy, ndivz,        &
                                      noverlap, maxit,            &
                                      x0, y0, z0, h, tol,         &
                                      ierr)                       &
                 BIND(C, NAME='eikonal3d_initialize')
      USE MPI
      USE EIKONAL3D_MODULE, ONLY : ghosts, leik_init, model_loc, parms
      USE EIKONAL3D_MODULE, ONLY : EIKONAL3D_GHOST_COMM, &
                                   MAKE_LEVEL_STRUCT
      USE ISO_C_BINDING
      IMPLICIT NONE
      REAL(C_DOUBLE), INTENT(IN) :: x0, y0, z0, h, tol
      INTEGER(C_INT), INTENT(IN) :: comm, iverb, maxit, ndivx, ndivy, ndivz, &
                                    noverlap, nx, ny, nz
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ! local variables
      INTEGER mpierr, myid
      INTEGER, PARAMETER :: master = 0
      !----------------------------------------------------------------------------------!
      CALL MPI_COMM_RANK(comm, myid, mpierr)
      IF (myid == master) THEN
         parms%iverb = iverb
         parms%nx = nx
         parms%ny = ny
         parms%nz = nz
         parms%ndivx = ndivx
         parms%ndivy = ndivy
         parms%ndivz = ndivz
         parms%noverlap = noverlap 
         parms%maxit = maxit
         parms%x0 = x0
         parms%y0 = y0
         parms%z0 = z0
         parms%h = h
         parms%tol = tol 
         IF (parms%iverb > 0) THEN
            WRITE(*,*) 'eikonal3d_initialize: Broadcasting parameters...'
         ENDIF
      ENDIF
      CALL MPI_BCAST(parms%iverb,    1, MPI_INTEGER,          master, comm, mpierr)
      CALL MPI_BCAST(parms%nx,       1, MPI_INTEGER,          master, comm, mpierr)
      CALL MPI_BCAST(parms%ny,       1, MPI_INTEGER,          master, comm, mpierr)
      CALL MPI_BCAST(parms%nz,       1, MPI_INTEGER,          master, comm, mpierr)
      CALL MPI_BCAST(parms%ndivx,    1, MPI_INTEGER,          master, comm, mpierr)
      CALL MPI_BCAST(parms%ndivy,    1, MPI_INTEGER,          master, comm, mpierr)
      CALL MPI_BCAST(parms%ndivz,    1, MPI_INTEGER,          master, comm, mpierr)
      CALL MPI_BCAST(parms%noverlap, 1, MPI_INTEGER,          master, comm, mpierr)
      CALL MPI_BCAST(parms%maxit,    1, MPI_INTEGER,          master, comm, mpierr)
      CALL MPI_BCAST(parms%x0,       1, MPI_DOUBLE_PRECISION, master, comm, mpierr)
      CALL MPI_BCAST(parms%y0,       1, MPI_DOUBLE_PRECISION, master, comm, mpierr)
      CALL MPI_BCAST(parms%z0,       1, MPI_DOUBLE_PRECISION, master, comm, mpierr)
      CALL MPI_BCAST(parms%h,        1, MPI_DOUBLE_PRECISION, master, comm, mpierr)
      CALL MPI_BCAST(parms%tol,      1, MPI_DOUBLE_PRECISION, master, comm, mpierr)
      IF (myid == master) THEN
         WRITE(*,*) 'eikonal3d_initialize: Generating ghost communication structure...'
      ENDIF
      CALL EIKONAL3D_GHOST_COMM(comm, parms%ndivx, parms%ndivy, parms%ndivz,  &
                                parms%noverlap, parms%nx, parms%ny, parms%nz, &
                                ghosts, model_loc, ierr)
      IF (ierr /= 0) THEN 
         WRITE(*,*) 'eikonal3d_initialize: Error making ghosts on process', myid
         GOTO 500
      ENDIF
      IF (myid == master) THEN
         WRITE(*,*) 'eikonal3d_initialize: Generating level structure...'
      ENDIF
      CALL MAKE_LEVEL_STRUCT(model_loc%nx, model_loc%ny, model_loc%nz, &
                             model_loc%lstruct, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'eikonal3d_initialize: Error generating local levels on process', myid
         GOTO 500
      ENDIF
      leik_init = .TRUE.
      RETURN
  500 CONTINUE
      WRITE(*,*) 'eikonal3d_initialize: An error occurred on process', myid 
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Function to test swap exchange blocking and non-blocking exchange
      SUBROUTINE TEST_SWAP(comm, ierr)
      USE MPI
      USE EIKONAL3D_MODULE, ONLY : ghosts, model_loc
      USE EIKONAL3D_MODULE, ONLY : EIKONAL_EXCHANGE, EIKONAL_EXCHANGE_BLOCKING
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: comm
      INTEGER, INTENT(OUT) :: ierr
      DOUBLE PRECISION, ALLOCATABLE :: u(:)
      INTEGER i, myid, mpierr, nxyz_loc
      CALL MPI_COMM_RANK(comm, myid, mpierr)
      ierr = 0
      nxyz_loc = model_loc%nx*model_loc%ny*model_loc%nz
      ALLOCATE(u(nxyz_loc))
      u(:) =-1.d0
      DO i=1,nxyz_loc
         IF (.NOT.model_loc%lghost_node(i)) THEN
            u(i) = FLOAT(model_loc%l2g_node(i))
         ENDIF
      ENDDO
      CALL EIKONAL_EXCHANGE_BLOCKING(comm, ghosts, u)
      DO i=1,nxyz_loc
         IF (INT(u(i)) /= model_loc%l2g_node(i)) THEN
            ierr = ierr + 1
         ENDIF
      ENDDO
      IF (ierr /= 0) THEN
         WRITE(*,*) 'Failed exchange blocking', myid
         RETURN
      ENDIF
      CALL MPI_BARRIER(comm, mpierr)

      u(:) =-1.d0
      DO i=1,nxyz_loc
         IF (.NOT.model_loc%lghost_node(i)) THEN 
            u(i) = FLOAT(model_loc%l2g_node(i))
         ENDIF
      ENDDO
      CALL EIKONAL_EXCHANGE(comm, ghosts, u)
      DO i=1,nxyz_loc
         IF (INT(u(i)) /= model_loc%l2g_node(i)) THEN 
            ierr = ierr + 1
         ENDIF
      ENDDO
      IF (ierr /= 0) THEN 
         WRITE(*,*) 'Failed exchange nonblocking', myid 
         RETURN
      ENDIF
 
      DEALLOCATE(u)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Computes the traveltimes in the given slowness model.  Note, the
!>           (ix,iy,iz)'th coordinate of the slowness model and traveltimes is
!>           accessed by: (iz - 1)*nx*ny + (iy - 1)*nx + ix.
!>
!>    @param[in] block_comm   MPI communicator between domain partitions
!>    @param[in] sweep_comm   MPI communicator between sweeps
!>    @param[in] n            number of grid points in domain (=nx*ny*nz)
!>    @param[in] ts           initial times (s) of source [nsrc]
!>    @param[in] xs           x source positions (m) [nsrc]
!>    @param[in] ys           y source positions (m) [nsrc]
!>    @param[in] zs           z source positions (m) [nsrc]
!>    @param[in] slow         slowness model (s/m) at each grid point [n]
!>
!>    @param[out] u           on successful exit this is traveltimes (seconds) from
!>                            the sources to all points in the medium [n] 
!>    @param[out] ierr        0 indicates success
!>
!>    @author Ben Baker
!>
!>    @copyright Apache 2 license
!>
      SUBROUTINE EIKONAL3D_SOLVE(comm, nsrc, n,         &
                                 ts, xs, ys, zs, slow,  &
                                 u, ierr)               &
                 BIND(C, NAME='eikonal3d_solve')
      USE MPI
      USE ISO_C_BINDING
      USE EIKONAL3D_MODULE, ONLY : ghosts, model_loc, parms
      USE EIKONAL3D_MODULE, ONLY : EIKONAL3D_FSM_MPI,          &
                                   EIKONAL_GATHER_TRAVELTIMES, &
                                   EIKONAL_SCATTER_BCS,        &
                                   EIKONAL_SCATTER_MODEL,      &
                                   EIKONAL3D_SETBCS
      INTEGER(C_INT), INTENT(IN) :: comm, n, nsrc
      REAL(C_DOUBLE), INTENT(IN) :: slow(n), ts(nsrc), xs(nsrc), ys(nsrc), zs(nsrc)
      REAL(C_DOUBLE), INTENT(OUT) :: u(n)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ! local variables
      LOGICAL, ALLOCATABLE :: lisbc(:), lisbc_loc(:)
      INTEGER i, nxyz_loc
      INTEGER, PARAMETER :: master = 0
      !----------------------------------------------------------------------------------!
      ierr = 0
      CALL MPI_COMM_RANK(comm, myid, mpierr)
      ! set the boundary conditions
      IF (myid == master) THEN
         IF (parms%iverb > 0) WRITE(*,*) 'eikonal3d_solve: Setting boundary conditions...'
         ALLOCATE(lisbc(n))
         CALL EIKONAL3D_SETBCS(parms%nx, parms%ny, parms%nz, nsrc, &
                               parms%h, parms%h, parms%h,          &
                               parms%x0, parms%y0, parms%z0,       &
                               ts, xs, ys, zs, slow,               &
                               lisbc, u, ierr)
         IF (ierr /= 0) THEN
            WRITE(*,*) 'eikonal3d_solve: Error setting bcs', myid
         ENDIF
      ELSE
         ALLOCATE(lisbc(1))
      ENDIF
      CALL MPI_BCAST(ierr, 1, MPI_INTEGER, master, comm, mpierr) 
      IF (ierr /= 0) RETURN
      ! initialize space
      nxyz_loc = model_loc%nx*model_loc%ny*model_loc%nz
      ALLOCATE(model_loc%slow(nxyz_loc))
      ALLOCATE(model_loc%u(nxyz_loc))
      ALLOCATE(lisbc_loc(nxyz_loc))
      ALLOCATE(model_loc%lupd(nxyz_loc))
      ! distribute the model
      CALL EIKONAL_SCATTER_MODEL(comm, model_loc%l2g_node, slow, model_loc%slow)
      IF (MINVAL(model_loc%slow) == 0.d0) THEN
         WRITE(*,*) 'eikonal3d_model: Error scattering model!'
         ierr = 1
         GOTO 500
      ENDIF
      ! distribute the boundary conditions
      CALL EIKONAL_SCATTER_BCS(comm, model_loc%l2g_node, lisbc, &
                               u, lisbc_loc, model_loc%u)
      IF (ALLOCATED(lisbc)) DEALLOCATE(lisbc)
      ! compute the nodes i am responsible for updating
      DO 1 i=1,nxyz_loc
         model_loc%lupd(i) = .TRUE.
         IF (lisbc_loc(i) .OR. model_loc%lghost_node(i)) THEN
            model_loc%lupd(i) = .FALSE.
         ENDIF
    1 CONTINUE
      ! free any unnecessary memory
      IF (ALLOCATED(lisbc_loc)) DEALLOCATE(lisbc_loc)
      IF (ALLOCATED(lisbc))     DEALLOCATE(lisbc)
      ! call the solver
      CALL EIKONAL3D_FSM_MPI(comm, comm,                               &
                             parms%iverb, parms%maxit,                 &
                             model_loc%nx, model_loc%ny, model_loc%nz, &
                             parms%h, parms%tol,                       &
                             model_loc%lstruct, ghosts,                &
                             model_loc%lghost_node, model_loc%lupd,    &
                             model_loc%slow, model_loc%u, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'eikonal3d_solve: Error calling solver on process', myid
         RETURN
      ENDIF
      ! fetch the solution onto the master
      CALL EIKONAL_GATHER_TRAVELTIMES(comm, model_loc%l2g_node, model_loc%u, u)
      ! free the memory
  500 CONTINUE
      DEALLOCATE(model_loc%slow)
      DEALLOCATE(model_loc%u)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Finalizes the 3D solver and releases any saved memory
!>
!>    @param[in] comm    MPI communicator
!>
!>    @param[out] ierr   0 indicates success
!>
!>    @author Ben Baker
!>
!>    @copyright Apache 2 license
!>
      SUBROUTINE EIKONAL3D_FINALIZE(comm, ierr) &
                 BIND(C, NAME='eikonal3d_finalize')
      USE EIKONAL3D_MODULE, ONLY : FREE_LEVEL_STRUCT, &
                                   FREE_GHOST_STRUCT, FREE_LOCAL_MODEL_STRUCT, &
                                   ghosts, leik_init, model_loc, parms
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: comm
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER mpierr, myid
      INTEGER, PARAMETER :: master = 0
      ierr = 0
      CALL MPI_COMM_RANK(comm, myid, mpierr)
      IF (.NOT.leik_init) THEN
         WRITE(*,*) 'eikonal3d_finalize: Solver was never initialized'
         ierr = 1
      ENDIF
      IF (myid == master .AND. parms%iverb > 0) THEN
         WRITE(*,*) 'eikonal3d_finalize: Freeing memory...'
      ENDIF
      CALL FREE_GHOST_STRUCT(ghosts)
      CALL FREE_LOCAL_MODEL_STRUCT(model_loc)
      parms%x0 = 0.d0
      parms%y0 = 0.d0
      parms%z0 = 0.d0
      parms%h  = 0.d0
      parms%tol = 0.d0
      parms%nx = 0
      parms%ny = 0
      parms%nz = 0
      parms%ndivx = 0
      parms%ndivy = 0
      parms%ndivz = 0
      parms%noverlap = 0
      parms%maxit = 0
      parms%iverb = 0
      leik_init = .FALSE.
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Utility function for solving the eikonal equation serially.  Note, after
!>           the initialization phase where you must only specify nx, ny, and nz
!>           the solver can be called repeatedly with different sources and velocity
!>           models.  The reference frame is right handed (+x right, +y away, +z up)
!>           and the (ix,iy,iz)'th model index is given by:
!>              (iz - 1)*nx*ny + (iy - 1)*nx + ix
!>
!>    @param[in] job       1 -> initialize.
!>                         2 -> solve given the source.
!>                         3 -> finalize.
!>    @param[in] iverb     controls the verbosity (0 is quiet)
!>    @param[in] maxit     max number of iterations in fast sweeping method
!>    @param[in] nsrc      number of sources
!>    @param[in] nx        number of x grid points in domain
!>    @param[in] ny        number of y grid points in domain
!>    @param[in] nz        number of z grid points in domain
!>    @param[in] tol       convergence tolerance (s)
!>    @param[in] h         grid spacing (m)
!>    @param[in] x0        x origin (m)
!>    @param[in] y0        y origin (m)
!>    @param[in] z0        z origin (m)
!>    @param[in] ts        initial time of source (s) [nsrc]
!>    @param[in] xs        x location of sources (m) [nsrc]
!>    @param[in] ys        y location of sources (m) [nsrc]
!>    @param[in] zs        z location of sources (m) [nsrc]
!>    @param[in] slow      slowness model at each grid point (s/m) [nx*ny*nz]
!>
!>    @param[out] u        when job = 2, on successful exit these are the traveltimes (s)
!>                         to each grid point [nx*ny*nz]
!>    @param[out] ierr     0 indicates success
!>  
!>    @author Ben Baker
!>
!>    @copyright Apache 2 license
!>
      SUBROUTINE EIKONAL3D_SERIAL_DRIVER(job, iverb, maxit, nsrc, &
                                         nx, ny, nz,              &
                                         tol, h, x0, y0, z0,      &
                                         ts, xs, ys, zs, slow,    &
                                         u, ierr)                 &
                 BIND(C, NAME='eikonal3d_serial_driver')
      USE EIKONAL3D_TYPES, ONLY : levelType
      USE EIKONAL3D_MODULE, ONLY : EIKONAL3D_FSM, EIKONAL3D_SETBCS,    &
                                   FREE_LEVEL_STRUCT, MAKE_LEVEL_STRUCT
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: job, iverb, maxit, nsrc, nx, ny, nz
      REAL(C_DOUBLE), INTENT(IN) :: slow(nx*ny*nz), ts(nsrc), xs(nsrc),    &
                                    ys(nsrc), zs(nsrc), h, tol, x0, y0, z0
      REAL(C_DOUBLE), INTENT(OUT) :: u(nx*ny*nz)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      ! local variables
      TYPE(levelType), SAVE :: lstruct
      LOGICAL, ALLOCATABLE :: lisbc(:), lupd(:)
      INTEGER i
      LOGICAL, SAVE :: linit = .FALSE.
      !----------------------------------------------------------------------------------!
      !
      ! initialize
      ierr = 0
      IF (job == 1) THEN
         IF (linit) THEN
            WRITE(*,*) 'eikonal3d_serial_driver: Already initialized!'
            ierr = 1
            RETURN
         ENDIF
         IF (iverb > 0) WRITE(*,*) 'eikonal3d_serial_driver: Generating levels...'
         CALL MAKE_LEVEL_STRUCT(nx, ny, nz, lstruct, ierr)
         IF (ierr /= 0) THEN
            WRITE(*,*) 'eikonal3d_serial_driver: Error generating level structure'
            RETURN
         ENDIF
         linit = .TRUE.
      ! run solver
      ELSE IF (job == 2) THEN
         ! check level structure exists
         IF (.NOT.linit) THEN
            WRITE(*,*) 'eikonal3d_serial_driver: Solver not initalized!'
            ierr = 1
            RETURN
         ENDIF
         ! set the boundary conditions 
         IF (iverb > 0) THEN
            WRITE(*,*) 'eikonal3d_serial_driver: Setting boundary conditions...'
         ENDIF
         ALLOCATE(lupd(nx*ny*nz))
         ALLOCATE(lisbc(nx*ny*nz))
         CALL EIKONAL3D_SETBCS(nx, ny, nz, nsrc,     &
                               h, h, h, x0, y0, z0,  &
                               ts, xs, ys, zs, slow, &
                               lisbc, u, ierr)
         IF (ierr /= 0) THEN
            WRITE(*,*) 'eikonal3d_serial_driver: Error setting boundary conditions'
            RETURN
         ENDIF
         ! determine who i am to update
         DO 1 i=1,nx*ny*nz
            lupd(i) = .TRUE.
            IF (lisbc(i)) lupd(i) = .FALSE.
   1     CONTINUE
         DEALLOCATE(lisbc)
         ! run the solver
         IF (iverb > 0) WRITE(*,*) 'eikonal3d_serial_driver: Solving...'
         CALL EIKONAL3D_FSM(iverb, maxit, nx, ny, nz, h, tol, &
                            lstruct, lupd, slow,       &
                            u, ierr)
         IF (ierr /= 0) THEN
            WRITE(*,*) 'eikonal3d_serial_driver: Error solving eikonal equation'
            RETURN
         ENDIF
         DEALLOCATE(lupd)
      ELSE
         IF (.NOT.linit) THEN
            WRITE(*,*) 'eikonal3d_serial_driver: Never initialized!'
         ENDIF 
         CALL FREE_LEVEL_STRUCT(lstruct)
         linit = .FALSE.
      ENDIF
      RETURN
      END
!========================================================================================!

      USE EIKONAL3D_TYPES
      USE EIKONAL3D_MODULE
      USE MPIUTILS_MODULE, ONLY : global_comm, intra_table_comm, inter_table_comm
      USE MPIUTILS_MODULE, ONLY : MPIUTILS_INITIALIZE3D
      USE MPI
!     TYPE (localModelType) model_loc
      TYPE(levelType) lstruct
!     TYPE(ghostType), ALLOCATABLE :: ghosts(:)
      DOUBLE PRECISION, ALLOCATABLE :: slow(:), u(:)
      INTEGER comm
      LOGICAL, ALLOCATABLE :: lupd(:), lisbc(:)
      DOUBLE PRECISION ts(1), xs(1), ys(1), zs(1), tol, h, x0, y0, z0
      DOUBLE PRECISION d, linf, t_est, x, y, z
      DOUBLE PRECISION, PARAMETER :: const_vel = 5.d3
      INTEGER, PARAMETER :: master = 0
      CALL MPI_INIT(mpierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, mpierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, mpierr)
      comm = MPI_COMM_WORLD
      ndivx = 2
      ndivy = 2
      ndivz = 2
      ! split the communicator
      CALL MPIUTILS_INITIALIZE3D(MPI_COMM_WORLD, 1, 1,         &
                                 ndivx, ndivy, ndivz, ierr) 
      IF (ierr /= 0) THEN
         WRITE(*,*) 'Error dividing up MPI groups'
         GOTO 500
      ENDIF
      comm = global_comm 
      nx = 70
      ny = 80
      nz = 90
      nsrc = 1
      maxit = 5
      x0 = 0.d0
      y0 = 0.d0
      z0 = 0.d0
      noverlap = 1
      tol = 1.d-7
      h = 100.d0
      iverb = 4
      ts(1) = 0.d0 
      xs(1) = x0 + h*DBLE(nx)/2.d0
      ys(1) = y0 + h*DBLE(ny)/2.d0
      zs(1) = z0 + h*DBLE(nz)/2.d0
      ! call the solver
      IF (myid /= master) THEN 
         n = 1
      ELSE 
         n = nx*ny*nz
      ENDIF
      ALLOCATE(u(n))
      ALLOCATE(slow(n))
      IF (myid == master) slow(:) = 1.d0/const_vel

      ! initialize the solver
      CALL EIKONAL3D_INITIALIZE(comm, iverb, nx, ny, nz,    &    
                                ndivx, ndivy, ndivz,        &
                                noverlap, maxit,            &
                                x0, y0, z0, h, tol,         &
                                ierr)
      CALL MPI_BARRIER(comm, mpierr)
!call TEST_SWAP(comm, ierr)
      ! call the solver
      CALL EIKONAL3D_SOLVE(comm, nsrc, n,         &
                           ts, xs, ys, zs, slow,  &
                           u, ierr)
if (myid == master) print *, minval(u), maxval(u)
      CALL MPI_BARRIER(comm, mpierr)
      ! finalize the solver
      CALL EIKONAL3D_FINALIZE(comm, ierr)
      IF (myid == master) THEN
         ! intialize
         CALL EIKONAL3D_SERIAL_DRIVER(1, iverb, maxit, nsrc,   &
                                      nx, ny, nz,              &
                                      tol, h, x0, y0, z0,      &
                                      ts, xs, ys, zs, slow,    &
                                      u, ierr)
         ! run
         CALL EIKONAL3D_SERIAL_DRIVER(2, iverb, maxit, nsrc,   &
                                      nx, ny, nz,              &    
                                      tol, h, x0, y0, z0,      &    
                                      ts, xs, ys, zs, slow,    &    
                                      u, ierr)
         ! deallocate
         CALL EIKONAL3D_SERIAL_DRIVER(3, iverb, maxit, nsrc,   &
                                      nx, ny, nz,              &    
                                      tol, h, x0, y0, z0,      &    
                                      ts, xs, ys, zs, slow,    &    
                                      u, ierr)
         print *, minval(u), maxval(u)
      ENDIF
goto 500
      linf = 0.d0
      ! compute error
      DO iz=1,nz
         DO iy=1,ny
            DO ix=1,nx
               x = x0 + h*FLOAT(ix - 1)
               y = y0 + h*FLOAT(iy - 1)
               z = z0 + h*FLOAT(iz - 1)
               d = SQRT( (x - xs(1))**2 + (y - ys(1))**2 + (z - zs(1))**2 )
               indx = (iz - 1)*nx*ny + (iy - 1)*nx + ix
               t_est = d*slow(indx)
               linf = MAX(ABS(u(indx) - t_est), linf)
            ENDDO
         ENDDO
      ENDDO 
      print *, linf
      ! free memory
      CALL FREE_LEVEL_STRUCT(lstruct)
      CALL FREE_GHOST_STRUCT(ghosts)
      CALL FREE_LOCAL_MODEL_STRUCT(model_loc)
      IF (ALLOCATED(u))     DEALLOCATE(u)
      IF (ALLOCATED(slow))  DEALLOCATE(slow)
      IF (ALLOCATED(lisbc)) DEALLOCATE(lisbc)
      IF (ALLOCATED(lupd))  DEALLOCATE(lupd)
500   CONTINUE
      IF (ierr /= 0) THEN
         WRITE(*,*) 'Error on process:', myid
         CALL MPI_ABORT(MPI_COMM_WORLD, myid, mpierr)
      ENDIF
      CALL MPI_COMM_FREE(intra_table_comm, mpierr)
      CALL MPI_COMM_FREE(inter_table_comm, mpierr)
      CALL MPI_COMM_FREE(global_comm, mpierr)
      CALL MPI_FINALIZE(mpierr)
      STOP
      END

!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Zeros out the PDF
!>
!>    @param[in] ngrd   number of grid points in grid search
!>    @param[out] pdf   zero-ed out PDF
!>
!>    @author Ben Baker
!>
!>    @copyright Apache 2 license
!>
      SUBROUTINE LOCATE_INITIALIZE_PDF(ngrd, pdf)
      USE LOCATE_MODULE, ONLY : zero
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: ngrd
      REAL(C_DOUBLE), INTENT(OUT) :: pdf(ngrd)
      INTEGER i
!$OMP DO SIMD
      DO 1 i=1,ngrd
         pdf(i) = zero 
    1 CONTINUE 
!$OMP END DO SIMD NOWAIT
      RETURN
      END


!     SUBROUTINE GRIDSEARCH( )
!     IF (nobs == 1) THEN
!        DO 1 i=1,ngrd
!           t0 = t0 + wtmtx
!   1    CONTINUE
!     ELSE
!        DO 1 i=1,ngrd
!           t0 = tori(i)
!           DO 2 iobs=1,nobs
!              t0 = t0 + wtmtx(iobs)*(tobs(iobs) -  
!   2       CONTINUE
!   1    CONTINUE
!     ENDIF

      SUBROUTINE LOCATE_NORMALIZE_PDF(comm, ngrd_loc, pdfloc, ierr)
      USE MPI
      USE LOCATE_MODULE, ONLY : one, zero
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: comm, ngrd_loc
      DOUBLE PRECISION, INTENT(INOUT) :: pdfloc(ngrd_loc)
      INTEGER, INTENT(OUT) :: ierr
      DOUBLE PRECISION xsumbuf, xsum, xsumi
      INTEGER mpierr, myid
      CALL MPI_COMM_RANK(comm, myid, mpierr)
      xsumbuf = SUM(pdfloc)
      CALL MPI_ALLREDUCE(xsumbuf, xsum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                         comm, mpierr)
      IF (xsum == zero) THEN
         WRITE(*,*) 'locate_optloc: Division by zero on rank', myid
         ierr = 1 
         RETURN
      ENDIF
      xsumi = one/xsum
      CALL DSCAL(ngrd_loc, xsumi, pdfloc, 1)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Finds the optimal global nodal index from the distributed PDF
!>
!>    @author Ben Baker
!>
!>    @copyright Apache 2 license
!>
      SUBROUTINE LOCATE_OPTNODE(comm, ngrd_loc, l2gnode, pdfloc, & 
                                nodeopt, ierr)
      USE MPI
      USE LOCATE_MODULE, ONLY : zero
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: pdfloc(ngrd_loc)
      INTEGER, INTENT(IN) :: l2gnode(ngrd_loc), comm, ngrd_loc
      INTEGER, INTENT(OUT) :: ierr, nodeopt
      ! local variables
      DOUBLE PRECISION, ALLOCATABLE :: xopt(:), xoptbuf(:)
      DOUBLE PRECISION xmax
      INTEGER igrdopt(1), i, ioptid, mpierr, myid, nprocs
      !----------------------------------------------------------------------------------!
      !
      ! initialize and get the MPI information
      nodeopt = 0
      CALL MPI_COMM_SIZE(comm, nprocs, mpierr)
      CALL MPI_COMM_RANK(comm, myid,   mpierr)
      ALLOCATE(xoptbuf(nprocs))
      ALLOCATE(xopt(nprocs))
      ! compute the optimal location 
      xoptbuf(:) = zero
      igrdopt = MAXLOC(pdfloc)
      xoptbuf(myid+1) = pdfloc(igrdopt(1))
      CALL MPI_ALLREDUCE(xoptbuf, xopt, nprocs, MPI_DOUBLE_PRECISION, MPI_SUM, &
                         comm, mpierr) 
      ! which process has the most likely origin location? 
      xmax = MAXVAL(xopt)
      ioptid =-1
      DO 1 i=1,nprocs
         IF (xmax == xopt(i)) ioptid = i - 1
    1 CONTINUE 
      IF (ioptid ==-1) THEN
         WRITE(*,*) 'locate_optloc: Failed to find optimum process!', myid
         ierr = 1
         RETURN
      ENDIF
      IF (ioptid == myid) nodeopt = l2gnode(igrdopt(1))
      CALL MPI_BCAST(nodeopt, 1, MPI_INTEGER, ioptid, comm, mpierr)
      ! free space
      DEALLOCATE(xoptbuf)
      DEALLOCATE(xopt) 
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE LOCATE3D_GRIDSEARCH_MPI(tttFileID, model, &
                                         job, ngrd, nobs, nsrc,  &
                                         statCor, tori, wtobs, tobs,          &
                                       hypo, ierr)  &
                 BIND(C,NAME='locate3d_gridSearch_mpi')
      USE MPI
      USE H5IO_MODULE, ONLY : H5IO_READ_TRAVELTIMESF
      USE LOCATE_MODULE, ONLY : COMPUTE_LOCATION_ONLY,            &
                                COMPUTE_LOCATION_AND_ORIGIN_TIME, &
                                COMPUTE_LOCATION_AND_STATICS,     &
                                COMPUTE_LOCATION_ALL, zero, one, sqrt2i
      USE LOCATE_TYPES, ONLY : locateType
      USE ISO_C_BINDING
      INTEGER(C_INT), INTENT(IN) :: job, model, nobs, ngrd, nsrc, tttFileID
      REAL(C_DOUBLE), INTENT(IN) :: statCor(nobs), tori(nsrc), wtobs(nsrc*nobs), tobs(nsrc*nobs)
      REAL(C_DOUBLE), INTENT(OUT) :: hypo(4*nsrc)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      DOUBLE PRECISION, ALLOCATABLE :: logPDF(:), logPDFbuf(:), t0(:), t0buf(:)
      REAL, ALLOCATABLE :: test4(:)
      DOUBLE PRECISION res, tobs_i, xnorm, wt_i, wt_i_sqrt2i
      INTEGER igrd, isrc
      INTEGER, PARAMETER :: master = 0
TYPE(locateType) :: locate
      ierr = 0
      hypo(1:4*nsrc) = zero
      nobs_groups = 1
      npdomain = 1
      mydomain_id = 0
      myobs_id = 0
      !CALL MPI_COMM_RANK(obs_comm,    myobs_id,    mpierr) ! 
      !CALL MPI_COMM_RANK(domain_comm, mydomain_id, mpierr) ! block of domain 
      !CALL MPI_COMM_SIZE(obs_comm,    nobs_groups, mpierr) ! observations groups
      !CALL MPI_COMM_SIZE(domain_comm, npdomain,    mpierr) ! processors in domain
      ALLOCATE(logPDF(ngrd))
      ALLOCATE(logPDFbuf(ngrd))
      ALLOCATE(test4(ngrd))
      ALLOCATE(t0(ngrd))
      logPDF(:) = zero
      logPDFbuf(:) = zero
      t0(:) = zero
      ! classify the job
      IF (job == COMPUTE_LOCATION_ONLY        .OR.    &
          job == COMPUTE_LOCATION_AND_ORIGIN_TIME) THEN
         DO 1 isrc=1,nsrc
            logPDF(:) = zero
            logPDFbuf(:) = zero
            hypo(4*(isrc - 1)+1) = tori(isrc)
            ! compute the analytic origin time (e.g. moser 1992 eqn 19)
            IF (job == COMPUTE_LOCATION_AND_ORIGIN_TIME) THEN
               ALLOCATE(t0buf(ngrd))
               t0buf(:) = zero 
               t0(:) = zero
               ! load test4 from disk
               CALL H5IO_READ_TRAVELTIMESF(MPI_COMM_WORLD, tttFileID,                &   
                                           iobs, model, iphase,                      &   
                                           locate%ix0, locate%iy0, locate%iz0,       &   
                                           locate%nxLoc, locate%nyLoc, locate%nzLoc, &
                                           test4, ierr)
               IF (ierr /= 0) THEN
                  WRITE(*,*) 'Error reading observed traveltimes 1'
                  GOTO 500 
               ENDIF
               ! tabulate the common residual
               DO 2 iobs=1,nobs
                  tobs_i =  tobs((isrc-1)*nobs+iobs) - statCor(iobs)
                  wt_i   = wtobs((isrc-1)*nobs+iobs)
                  IF (wt_i > zero) THEN
                     !$OMP DO SIMD
                     DO 3 igrd=1,ngrd
                        t0buf(igrd) = t0buf(igrd) + wt_i*(tobs_i - DBLE(test4(igrd)))
    3                CONTINUE
                     !$OMP END DO SIMD
                  ENDIF
    2          CONTINUE
               xnorm = SUM(wtobs((isrc-1)*nobs+1:isrc*nobs))
               IF (xnorm > zero) xnorm = one/xnorm
               CALL DSCAL(ngrd, xnorm, t0buff, 1) ! divide by the sum of the weights
               ! add up the common residual (which is the origin time at each grid point)
               !CALL MPI_ALLREDUCE(t0buf, t0, MPI_DOUBLE_PRECISION, MPI_SUM, &
               !                   domain_comm, mpierr)
            ELSE
               t0(:) = tori(isrc)
            ENDIF
            ! stack the residuals for this event
            DO 4 iobs=1,nobs,nobs_groups
               ! load test4 from disk
!              CALL H5IO_READ_TRAVELTIMESF(myobs_comm, tttFileID,                &
!                                          iobs, model, iphase,                      &
!                                          locate%ix0, locate%iy0, locate%iz0,       &
!                                          locate%nxLoc, locate%nyLoc, locate%nzLoc, &
!                                          test4, ierr)
               IF (ierr /= 0) THEN
                  WRITE(*,*) 'Error reading observed traveltimes 2'
                  GOTO 500
               ENDIF
               ! set the observations and weights 
               tobs_i =  tobs((isrc-1)*nobs+iobs) - statCor(iobs)
               wt_i   = wtobs((isrc-1)*nobs+iobs)
               wt_i_sqrt2i = wt_i*sqrt2i ! accounts for 1/2 scaling in L2 norm
               ! sum the residuals on the grid [log(e^{-(L_1+L2+...)}) =-L_1 - L_2 ...
               !$OMP DO SIMD
               DO 5 igrd=1,ngrd
                  res = wt_i_sqrt2i*(tobs_i - (DBLE(test4(igrd)) + t0(igrd)))
                  logPDFbuf(igrd) = logPDFbuf(igrd) - res*res
    5          CONTINUE
               !$OMP END DO SIMD
    4       CONTINUE
            ! reduce the log of PDFs onto the first observation group
            !CALL MPI_REDUCE(logPDFbuf, logPDF, ngrd, MPI_DOUBLE_PRECISION, MPI_SUM, master, &
            !                cross_comm, mpierr)
            ! now the head group can locate the event
            IF (mydomain_id == master) THEN

            ENDIF
    1    CONTINUE ! loop on sources
      ! locate many sources
      ELSEIF (job == COMPUTE_LOCATION_AND_STATICS) THEN
         WRITE(*,*) 'Not yet done'
         ierr = 1
!        DO 100 iobs=1,nobs
!           DO 200 isrc=1,nsrc

! 200       CONTINUE
! 100    CONTINUE
      ELSEIF (job == COMPUTE_LOCATION_ALL) THEN
         WRITE(*,*) 'Not yet done'
         ierr = 1
      ELSE
         WRITE(*,*) 'locate_gridsearch: Invalid job'
         ierr = 1
      ENDIF
  500 CONTINUE
      RETURN
      END

!     SUBROUTINE LOCATE_MAKE_EIK2TTABLE_MAP(model_loc, locate_loc )
!     USE EIKONAL3D_TYPES, ONLY : localModelType
!     USE LOCATE_TYPES, ONLY : locateType
!     TYPE(localModelType), INTENT(IN) :: model_loc
!     TYPE(locateType), INTENT(INOUT) :: locate_loc
!     INTEGER i, nxyz_eik, nxyz_loc
!     ! verify the nodes are ordered
!     nxyz_eik = model_loc%nx*model_loc%ny*model_loc%nz
!     DO 1 i=2,nxyz_eik
!        IF (model_loc%l2g_node(i) <= model_loc%l2g_node(i-1)) THEN
!           lord = .FALSE.
!           EXIT 
!        ENDIF
!   1 CONTINUE
!     DO 2 i=2,nxyz_loc
!        IF (locate_loc%l2g_node(i) <= locate_loc%l2g_node(i-1)) THEN
!           lord = .FALSE.
!           EXIT
!        ENDIF
!   2 CONTINUE
!     RETURN
!     END

      SUBROUTINE LOCATE_SET_COMM(comm, nx, ny, nz, ndivx, ndivy, ndivz, &
                                 dx, dy, dz, x0, y0, z0, locate, ierr)
      USE MPIUTILS_MODULE, ONLY : MPIUTILS_GRD2IJK
      USE LOCATE_TYPES, ONLY : locateType
      DOUBLE PRECISION, INTENT(IN) :: dx, dy, dz, x0, y0, z0
      INTEGER, INTENT(IN) :: comm, nx, ny, nz, ndivx, ndivy, ndivz
      TYPE(locateType), INTENT(OUT) :: locate
      INTEGER, INTENT(OUT) :: ierr
      INTEGER i1, i2, imbx, imby, imbz, j1, j2, k1, k2
      CALL MPI_COMM_RANK(comm, myid, mpierr)
      CALL MPIUTILS_GRD2IJK(myid, ndivx, ndivy, ndivz, imbx, imby, imbz, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'Error finding process block'
         RETURN
      ENDIF
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
      nx_loc = i2 - i1 + 1
      ny_loc = i2 - i1 + 1
      nz_loc = i2 - i1 + 1
      locate%ngrd = nx_loc*ny_loc*nz_loc
      locate%nxLoc = nx_loc
      locate%nyLoc = ny_loc
      locate%nzLoc = nz_loc
      locate%nx = nx
      locate%ny = ny
      locate%nz = nz
      locate%ix0 = i1
      locate%iy0 = j1
      locate%iz0 = k1
      ALLOCATE(locate%l2g_node(locate%ngrd))
      ALLOCATE(locate%xlocs(locate%ngrd))
      ALLOCATE(locate%ylocs(locate%ngrd))
      ALLOCATE(locate%zlocs(locate%ngrd))     
      DO 1 iz=k1,k2
         DO 2 iy=j1,j2
            DO 3 ix=i1,i2
               i = (iz - k1)*nx_loc*ny_loc &
                 + (iy - j1)*nx_loc        &
                 + (ix - i1) + 1 
               locate%l2g_node(i) = (iz - 1)*nx*ny + (iy - 1)*nx + ix
               locate%xlocs(i) = x0 + FLOAT(ix - 1)*dx
               locate%ylocs(i) = y0 + FLOAT(iy - 1)*dy
               locate%zlocs(i) = z0 + FLOAT(iz - 1)*dz
    3       CONTINUE
    2    CONTINUE
    1 CONTINUE
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE LOCATE3D_INITIALIZE(comm, iverb, nx, ny, nz, &
                                     ndivx, ndivy, ndivz,     &
                                     ierr)                    &
                 BIND(C,NAME='locate3d_initialize')
      USE MPI
      USE MPIUTILS_MODULE, ONLY : MPIUTILS_INITIALIZE3D, linitComm
      USE LOCATE_MODULE, ONLY : parms
      USE ISO_C_BINDING
      INTEGER(C_INT), INTENT(IN) :: comm, iverb, nx, ny, nz, ndivx, ndivy, ndivz
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER, PARAMETER :: master = 0
      INTEGER, PARAMETER :: ireord = 1
      INTEGER, PARAMETER :: iwt = 0
      ierr = 0
      CALL MPI_COMM_RANK(comm, myid, mpierr)
      IF (myid == master) THEN
         parms%iverb = iverb
         parms%nx = nx
         parms%ny = ny
         parms%nz = nz
         parms%ndivx = ndivx
         parms%ndivy = ndivy
         parms%ndivz = ndivz
      ENDIF
      CALL MPI_BCAST(parms%iverb, 1, MPI_INTEGER, master, comm, mpierr)
      CALL MPI_BCAST(parms%nx,    1, MPI_INTEGER, master, comm, mpierr)
      CALL MPI_BCAST(parms%ny,    1, MPI_INTEGER, master, comm, mpierr)
      CALL MPI_BCAST(parms%nz,    1, MPI_INTEGER, master, comm, mpierr)
      CALL MPI_BCAST(parms%ndivx, 1, MPI_INTEGER, master, comm, mpierr)
      CALL MPI_BCAST(parms%ndivy, 1, MPI_INTEGER, master, comm, mpierr)
      CALL MPI_BCAST(parms%ndivz, 1, MPI_INTEGER, master, comm, mpierr)
      IF (.NOT.linitComm) THEN
         WRITE(*,*) 'locate_initialize3d: Splitting communicator...'
         CALL MPIUTILS_INITIALIZE3D(comm, ireord, iwt,                     &
                                    parms%ndivx, parms%ndivy, parms%ndivz, &
                                    ierr)
         IF (ierr /= 0) THEN
            WRITE(*,*) 'locate_initialize3d: Error splitting communicator'
            RETURN
         ENDIF
      ENDIF
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !

!     SUBROUTINE LOCATE_INITIALIZE(comm, iverb, nx, ny, nz, &
!                                  ndivx, ndivy, ndivz )
 
!     RETURN
!     END 
! 
!     SUBROUTINE LOCATE_UPDATE_PDF(nobs, )
!     DO 1 i=1,
!     RETURN
!     END

      subroutine delete_me()

      USE MPI
      USE MPIUTILS_MODULE, ONLY : global_comm, intra_table_comm, inter_table_comm 
      USE MPIUTILS_MODULE, ONLY : MPIUTILS_INITIALIZE3D
      CALL MPI_INIT(mpierr)

      ndivx = 2
      ndivy = 2
      ndivz = 2
      CALL MPIUTILS_INITIALIZE3D(MPI_COMM_WORLD, 1, 1,    &
                                 ndivx, ndivy, ndivz, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'Error initializing MPI communicators'
         GOTO 500
      ENDIF


  500 CONTINUE 
      CALL MPI_COMM_FREE(intra_table_comm, mpierr)
      CALL MPI_COMM_FREE(inter_table_comm, mpierr)
      CALL MPI_COMM_FREE(global_comm, mpierr)
      CALL MPI_FINALIZE(mpierr)
      STOP
      END

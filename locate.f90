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

      SUBROUTINE VERIF_MODEL(iphase, nx, ny, nz, ttimes4)
      USE ISO_C_BINDING
      INTEGER, INTENT(IN) :: iphase, nx, ny, nz
      REAL, INTENT(IN) :: ttimes4(nx*ny*nz)
      INTERFACE
         INTEGER(C_INT) FUNCTION computeHomogeneousTraveltimes_finter( &
                                    nx, ny, nz, &
                                    x0, y0, z0, &
                                    dx, dy, dz, &
                                    xs, ys, zs, &
                                    vel, ttimes) &
                        BIND(C, NAME='computeHomogeneousTraveltimes_finter')
         USE ISO_C_BINDING 
         INTEGER(C_INT), INTENT(IN) :: nx, ny, nz
         REAL(C_DOUBLE), INTENT(IN) :: x0, y0, z0, dx, dy, dz, xs, ys, zs, vel
         REAL(C_DOUBLE), INTENT(OUT) :: ttimes(nx*ny*nz)
         END
      END INTERFACE 
      REAL(C_DOUBLE), ALLOCATABLE :: ttimes(:)
      REAL(C_DOUBLE) dx, dy, dz, vel, x0, y0, z0, xs, ys, zs
      INTEGER(C_INT) ierr
      INTEGER i, igrd, ix, iy, iz
      dx = 1.d3
      dy = 1.d3
      dz = 1.d3
      x0 = 0.d0
      y0 = 0.d0
      z0 = 0.d0
      IF (iphase == 1) THEN
         vel = 2.d3
      ELSE
         vel = 2.d3/DSQRT(3.d0)
      ENDIF
      xs = 0.d0
      ys = 0.d0
      zs = 0.d0
      igrd = MINLOC(ttimes4, 1)
      DO iz=1,nz
         DO iy=1,ny
            DO ix=1,nx
               IF ((iz-1)*nx*ny + (iy-1)*nx + ix == igrd) THEN
                  xs = x0 + (ix - 1)*dx
                  ys = y0 + (iy - 1)*dy
                  zs = z0 + (iz - 1)*dz
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      ALLOCATE(ttimes(nx*ny*nz))
      ierr = computeHomogeneousTraveltimes_finter(nx, ny, nz, &
                                                  x0, y0, z0, &
                                                  dx, dy, dz, &
                                                  xs, ys, zs, &
                                                  vel, ttimes)
      DO i=1,nx*ny*nz
         IF (ABS(ttimes(i) - ttimes4(i)) > 1.d-5) THEN
            WRITE(*,*) 'model mismatch:', ttimes(i), ttimes4(i)
            ierr = 1
            RETURN
         ENDIF
      ENDDO
      WRITE(*,*) 'Model is okay', MAXVAL(ABS(ttimes - ttimes4))
      DEALLOCATE(ttimes) 
      RETURN
      END

      SUBROUTINE LOCATE3D_SET_TTABLES_FROM_FILE( )

      END

      SUBROUTINE LOCATE3D_SET_TTABLE_FROM_FILE(ngrd, iphase, istat, &
                                               model, test4, ierr)
      USE MPIUTILS_MODULE, ONLY : global_comm, inter_table_comm, intra_table_comm
      USE H5IO_MODULE, ONLY : H5IO_READ_TRAVELTIMESF
      USE LOCATE_MODULE, ONLY : locate, parms
      INTEGER, INTENT(IN) :: ngrd, iphase, istat, model
      INTEGER, INTENT(OUT) :: ierr
      REAL, INTENT(OUT) :: test4(ngrd)
      INTEGER mpierr, myid
      !----------------------------------------------------------------------------------!
      ierr = 0
      CALL MPI_COMM_RANK(intra_table_comm, myid, mpierr)
      IF (iphase < 1 .OR. iphase > 2) THEN
         WRITE(*,*) 'locate3d_set_ttable_from_file: Invalid phase', iphase, myid
         ierr = 1
      ENDIF
      IF (model < 1) THEN
         WRITE(*,*) 'locate3d_set_ttable_from_file: Invalid model', model, myid
         ierr = 1
      ENDIF
      IF (istat < 1) THEN
         WRITE(*,*) 'locate3d_set_ttable_from_file: Invalid station', istat, myid
         ierr = 1
      ENDIF
      IF (ierr /= 0) RETURN
      CALL H5IO_READ_TRAVELTIMESF(intra_table_comm, parms%tttFileID,        &
                                  istat, model, iphase,                     &
                                  locate%ix0, locate%iy0, locate%iz0,       &
                                  locate%nxLoc, locate%nyLoc, locate%nzLoc, &
                                  test4, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'locate3d_set_ttable_from_file: Error reading traveltimes'
         ierr = 1
      ENDIF
      RETURN
      END SUBROUTINE

      SUBROUTINE LOCATE3D_GRIDSEARCH_new( )
      USE MPI
      LOGICAL, ALLOCATABLE :: lhasTable(:)
      INTEGER stat(MPI_STATUS_SIZE), nevGroups
      INTEGER, PARAMETER :: master = 0
      ntables = 4
      nevGroups = 1
      myevGroup = 0
      IF (nevGroups == 1) THEN

      ELSE
         IF (myevGroup == master) THEN
            ALLOCATE(lhasTable(nevGroups*ntables))
            ALLOCATE(needTables(ntables))
         ELSE
            ALLOCATE(lhasTable(ntables))
            ALLOCATE(needTables(ntables))
         ENDIF
         lhasTable(:) = .FALSE.
         ! parallel loop on events
         IF (myevGroup == master) THEN
            DO 11 jev=1,nevents,nevGroups
               ! 
               IF (nevGroups > 1) THEN

               ! just do it all myself
               ELSE

               ENDIF 
   11       CONTINUE
         ELSE
            ! receive the event number
            mysrc = myTableID
            DO WHILE (.TRUE.)
               CALL MPI_RECV(iev, 1, MPI_INTEGER, mysrc, MPI_ANY_TAG, &
                             globalComm, stat, mpierr)
               IF (iev < 0) EXIT ! i'm done
               ! look through the catalog and see what tables i need
               needTable(:) = 0
               iobs1 = catptr(iev)
               iobs2 = catptr(iev+1) - 1
               nobs = iobs2 - iobs1 + 1
               nneed = 0
               DO 21 iobs=iobs1,iobs2
                  itable = tablePtr(iobs)
                  IF (lhasTable(itable)) THEN
                     nneed = nneed + 1
                     needTable(nneed) = itable
                  ENDIF
   21          CONTINUE
               nneed = SUM(needTable) 
               ! tell the master what tables i need
               CALL MPI_SEND(nneed, 1, MPI_INTEGER)
               IF (nneed > 0) THEN
                  CALL MPI_SEND(needTable, nneed, MPI_INTEGER)
               ENDIF
               ! now receive my tables
               DO 22 ineed=1,nneed
                  itable = needTable(ineed)
                  igrd1 = (itable - 1)*ngrd + 1
                  igrd2 = itable*ngrd
                  CALL MPI_RECV( )
   22          CONTINUE
            ENDDO
         ENDIF
         IF (ALLOCATED(lhasTable)) DEALLOCATE(lhasTable)
      ENDIF
      END SUBROUTINE
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Driver routine for earthquake location
!>
!>    @param[in] model        model number
!>    @param[in] job          if job == 0 then compute the locations only.
!>                            if job == 1 then compute the locations and origin times.
!>    @param[in] ngrd         number of grid points in domain
!>    @param[in] nobs         max number of observations for all the events
!>    @param[in] nevents      number of events to locate
!>    @param[in] luseObs      if 1 then use the iobs'th observation for the isrc'th
!>                            source [nobs*nevents]
!>    @param[in] statCor      station static correction for this observation type [nobs]
!>    @param[in] tori         origin time (epochal seconds)
!>    @param[in] varobs       variance in the iobs'th observation for the isrc'th
!>                            source [nobs*nevents] 
!>    @param[in] tobs         observed traveltime (pick time) for the iobs'th observation
!>                            for the isrc'th source [nobs*nevents]
!>
!>    @param[out] hypo        on successful output contains the optimal locations
!>                            (x,y,z,t0) for all sources [4*nevents]
!>    @param[out] ierr        0 indicates success
!>
!>    @author Ben Baker
!>
!>    @copyright Apache 2
!>
      SUBROUTINE LOCATE3D_GRIDSEARCH(model,                 &
                                     job, nobs, nevents, &
                                     luseObs, statPtr,  &
                                     pickType, statCor,      &
                                     tori, varobs,          &
                                     tobs, test,            &
                                     hypo, ierr)      &
                 BIND(C,NAME='locate3d_gridsearch')
      USE MPI
      USE H5IO_MODULE, ONLY : H5IO_READ_TRAVELTIMESF
      USE LOCATE_MODULE, ONLY : COMPUTE_LOCATION_ONLY,            &
                                COMPUTE_LOCATION_AND_ORIGIN_TIME, &
                                COMPUTE_LOCATION_AND_STATICS,     &
                                COMPUTE_LOCATION_ALL, zero, one, sqrt2i
      USE MPIUTILS_MODULE, ONLY : global_comm, inter_table_comm, intra_table_comm
      USE LOCATE_MODULE, ONLY : locate, parms
      USE ISO_C_BINDING
      INTERFACE
          SUBROUTINE DSCAL(n, alpha, x, incx)
          INTEGER, INTENT(IN) :: n, incx
          DOUBLE PRECISION, INTENT(IN) :: alpha 
          DOUBLE PRECISION, INTENT(INOUT) :: x(n)
          END SUBROUTINE DSCAL
      END INTERFACE
      INTEGER(C_INT), INTENT(IN) :: job, model, nobs, nevents
      INTEGER(C_INT), INTENT(IN) :: luseObs(nevents*nobs), pickType(nevents*nobs), &
                                    statPtr(nevents*nobs)
      REAL(C_DOUBLE), INTENT(IN) :: statCor(nobs), tori(nevents),       &
                                    varobs(nevents*nobs), tobs(nevents*nobs)
      REAL(C_DOUBLE), INTENT(OUT) :: hypo(4*nevents), test(nevents*nobs)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      DOUBLE PRECISION, ALLOCATABLE :: logPDF(:), logPDFbuf(:), objbuf(:), objloc(:), t0(:), t0buf(:)
      REAL, ALLOCATABLE :: test4(:)
      DOUBLE PRECISION hypo4(4), res, tobs_i, xnorm, wt_i, wt_i_sqrt2i
      INTEGER stat(MPI_STATUS_SIZE), igrd, isrc, globalID, myblockID, ntableGroups, tableID
      INTEGER, PARAMETER :: master = 0
      !----------------------------------------------------------------------------------!
      ierr = 0
!     hypo(1:4*nevents) = zero
      npdomain = 1
      mydomain_id = 0
      myobs_id = 0
      CALL MPI_COMM_SIZE(inter_table_comm, ntableGroups, mpierr)
      CALL MPI_COMM_RANK(inter_table_comm, tableID, mpierr)
      CALL MPI_COMM_SIZE(intra_table_comm, nblocks, mpierr)
      CALL MPI_COMM_RANK(intra_table_comm, myblockID, mpierr)
      !CALL MPI_COMM_RANK(obs_comm,    myobs_id,    mpierr) ! 
      !CALL MPI_COMM_RANK(domain_comm, mydomain_id, mpierr) ! block of domain 
      !CALL MPI_COMM_SIZE(obs_comm,    nobs_groups, mpierr) ! observations groups
      !CALL MPI_COMM_SIZE(domain_comm, npdomain,    mpierr) ! processors in domain
      ngrd = locate%ngrd
      ALLOCATE(logPDF(ngrd))
      ALLOCATE(logPDFbuf(ngrd))
      ALLOCATE(test4(ngrd))
      ALLOCATE(t0(ngrd))
      ALLOCATE(objloc(nblocks))
      ALLOCATE(objbuf(nblocks))
      logPDF(:) = zero
      logPDFbuf(:) = zero
      t0(:) = zero
      ! classify the job
      IF (job == COMPUTE_LOCATION_ONLY        .OR.    &
          job == COMPUTE_LOCATION_AND_ORIGIN_TIME) THEN
         DO 1 isrc=1,nevents
            ! compute the analytic origin time (e.g. moser 1992 eqn 19)
            IF (job == COMPUTE_LOCATION_AND_ORIGIN_TIME) THEN
               IF (ALLOCATED(t0buf)) deallocate(t0buf)
               ALLOCATE(t0buf(ngrd))
               t0buf(:) = zero 
               t0(:) = zero
               ! tabulate the common residual
               DO 2 iobs=1,nobs,ntableGroups
                  myobs = (isrc - 1)*nobs + iobs
                  IF (iobs + tableID > nobs) CYCLE
                  IF (luseObs(iobs) == 0) CYCLE
                  ! load the estimates from disk
                  istat = statPtr(myobs)
                  iphase = pickType(myobs)
                  CALL H5IO_READ_TRAVELTIMESF(intra_table_comm, parms%tttFileID,        &
                                              istat, model, iphase,                     &
                                              locate%ix0, locate%iy0, locate%iz0,       &
                                              locate%nxLoc, locate%nyLoc, locate%nzLoc, &
                                              test4, ierr)
                  IF (ierr /= 0) THEN
                     WRITE(*,*) 'Error reading observed traveltimes 1'
                     GOTO 500 
                  ENDIF
                  tobs_i =   tobs(myobs) - statCor(iobs)
                  wt_i   = varobs(myobs)
                  IF (wt_i > zero) THEN
                     !$OMP DO SIMD
                     DO 3 igrd=1,ngrd
                        t0buf(igrd) = t0buf(igrd) + wt_i*(tobs_i - DBLE(test4(igrd)))
    3                CONTINUE
                     !$OMP END DO SIMD
                  ENDIF
    2          CONTINUE
               ! normalize by the sum of the weights which is equivalent to scaling
               ! by the sum of the data variances
               xnorm = zero
               iobs1 = (isrc-1)*nobs+1
               iobs2 = isrc*nobs
               IF (SUM(luseObs(iobs1:iobs2)) > 0) THEN
                  xnorm = SUM(varobs(iobs1:iobs2))
                  CALL DSCAL(ngrd, xnorm, t0buf, 1) ! divide by the sum of the weights
               ENDIF
               ! add up the common residual (which is the origin time at each grid point)
               CALL MPI_ALLREDUCE(t0buf, t0, ngrd, MPI_DOUBLE_PRECISION, MPI_SUM, &
                                  inter_table_comm, mpierr)
            ELSE
               t0(:) = tori(isrc)
            ENDIF
            ! stack the residuals for this event into the PDF
            logPDFbuf(:) = zero
            DO 4 iobs=1,nobs,ntableGroups
               myobs = (isrc - 1)*nobs + iobs
               IF (iobs + tableID > nobs) CYCLE
               IF (luseObs(iobs) == 0) CYCLE
               ! load the estimates from disk
               istat = statPtr(myobs)
               iphase = pickType(myobs)
               CALL H5IO_READ_TRAVELTIMESF(intra_table_comm, parms%tttFileID,        &
                                           istat, model, iphase,                     &
                                           locate%ix0, locate%iy0, locate%iz0,       &
                                           locate%nxLoc, locate%nyLoc, locate%nzLoc, &
                                           test4, ierr)
               IF (ierr /= 0) THEN
                  WRITE(*,*) 'Error reading observed traveltimes 1'
                  GOTO 500 
               ENDIF
               ! set the observations and weights 
               tobs_i = tobs(myobs) - statCor(iobs)
               wt_i   = one/varobs(myobs)
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
            CALL MPI_REDUCE(logPDFbuf, logPDF, ngrd, MPI_DOUBLE_PRECISION, MPI_SUM, &
                            master, inter_table_comm, mpierr)
            ! now the head group can locate the event
            IF (tableID == master) THEN
               ioptLoc = MAXLOC(logPDF, 1)
               objbuf(:) =-HUGE(one)
               hypo(:) = zero
               objbuf(myblockID+1) = logPDF(ioptLoc)
               CALL MPI_ALLREDUCE(objbuf, objloc, nblocks, MPI_DOUBLE_PRECISION, &
                                 MPI_MAX, intra_table_comm, mpierr) 
               ioptBlock = MAXLOC(objloc, 1) - 1
               IF (ioptBlock == myblockID) THEN
                  print *, logPDF(ioptLoc), ioptloc
                  hypo4(1) = locate%xlocs(ioptLoc)
                  hypo4(2) = locate%ylocs(ioptLoc)
                  hypo4(3) = locate%zlocs(ioptLoc)
                  hypo4(4) = t0(ioptLoc)
               ENDIF
               ! get the optimal hypocenter back on the master
               IF (myblockID == master) THEN
                  IF (ioptBlock /= myblockID) THEN
print *, 'receiving from:', ioptBlock
                     CALL MPI_RECV(hypo4, 4, MPI_DOUBLE_PRECISION,   &
                                   ioptBlock, MPI_ANY_TAG, intra_node_comm, stat, &
                                   mpierr)
                  ENDIF
                  hypo(4*(isrc-1)+1:4*(isrc-1)+4) = hypo4(1:4)
print *, hypo4
               ELSE
                  IF (ioptBlock == myblockID) THEN
                     CALL MPI_SEND(hypo4, 4, MPI_DOUBLE_PRECISION, master,       &
                                   ioptBlock, intra_node_comm, mpierr)
                  ENDIF
               ENDIF
            ENDIF
    1    CONTINUE ! loop on sources
      ! locate many sources
      ELSEIF (job == COMPUTE_LOCATION_AND_STATICS) THEN
         WRITE(*,*) 'Not yet done'
         ierr = 1
!        DO 100 iobs=1,nobs
!           DO 200 isrc=1,nevents

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

!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Initializes the 3D location strutures
!>
!>    @param[in] comm       communicator to split in location grid
!>    @param[in] tttFileID  traveltime HDF5 file handles 
!>    @param[in] locFileID  location HDF5 file handles
!>    @param[in] ndivx      number of blocks to divide domain in x
!>    @param[in] ndivy      number of blocks to divide domain in y
!>    @param[in] ndivz      number of blocks to divide domain in z 
!>
!>    @param[out] ierr      0 indicates success
!>
!>    @author Ben Baker
!>
!>    @copyright Apache 2
!>
      SUBROUTINE LOCATE3D_INITIALIZE(comm, iverb,                 &
                                     tttFileID, locFileID,        &
                                     ndivx, ndivy, ndivz,         &
                                     ierr)  &
                 BIND(C,NAME='locate3d_initialize')
      USE MPI
      USE H5IO_MODULE, ONLY : H5IO_GET_MODEL_DIMENSIONSF, H5IO_READ_MODELF
      USE MPIUTILS_MODULE, ONLY : MPIUTILS_INITIALIZE3D, MPIUTILS_GRD2IJK, &
                                  global_comm, intra_table_comm, linitComm
      USE LOCATE_MODULE, ONLY : locate, parms, zero
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_LONG), INTENT(IN) :: locFileID, tttFileID
      INTEGER(C_INT), INTENT(IN) :: comm, iverb, ndivx, ndivy, ndivz
      INTEGER(C_INT), INTENT(OUT) :: ierr
      INTEGER(C_INT) nx, ny, nz
      INTEGER i1, i2, imbx, imby, imbz, j1, j2, k1, k2, &
              mpierr, myid, nall, ndx, ndy, ndz, ngrdAll
      INTEGER, PARAMETER :: master = 0
      INTEGER, PARAMETER :: ireord = 1
      INTEGER, PARAMETER :: iwt = 0
      ierr = 0
      CALL MPI_COMM_RANK(comm, myid, mpierr)
      IF (myid == master) THEN
         CALL H5IO_GET_MODEL_DIMENSIONSF(tttFileID,      &
                                         nx, ny, nz, ierr)
         IF (ierr /= 0) THEN
            WRITE(*,*) 'locate3d_initialize: Error getting dimensions'
         ENDIF
         parms%iverb = iverb
         parms%nx = nx
         parms%ny = ny
         parms%nz = nz
         parms%ndivx = ndivx
         parms%ndivy = ndivy
         parms%ndivz = ndivz
         parms%locFileID = locFileID
         parms%tttFileID = tttFileID
      ENDIF
      CALL MPI_BCAST(parms%iverb, 1, MPI_INTEGER, master, comm, mpierr)
      CALL MPI_BCAST(parms%nx,    1, MPI_INTEGER, master, comm, mpierr)
      CALL MPI_BCAST(parms%ny,    1, MPI_INTEGER, master, comm, mpierr)
      CALL MPI_BCAST(parms%nz,    1, MPI_INTEGER, master, comm, mpierr)
      CALL MPI_BCAST(parms%ndivx, 1, MPI_INTEGER, master, comm, mpierr)
      CALL MPI_BCAST(parms%ndivy, 1, MPI_INTEGER, master, comm, mpierr)
      CALL MPI_BCAST(parms%ndivz, 1, MPI_INTEGER, master, comm, mpierr)
      CALL MPI_BCAST(parms%locFileID, 1, MPI_LONG_LONG, master, comm, mpierr)
      CALL MPI_BCAST(parms%tttFileID, 1, MPI_LONG_LONG, master, comm, mpierr)
      IF (.NOT.linitComm) THEN
         WRITE(*,*) 'locate_initialize3d: Splitting communicator...'
         CALL MPIUTILS_INITIALIZE3D(comm, ireord, iwt,                     &
                                    parms%ndivx, parms%ndivy, parms%ndivz, &
                                    ierr)
         IF (ierr /= 0) THEN
            WRITE(*,*) 'locate_initialize3d: Error splitting communicator'
            ierr = 1
            RETURN
         ENDIF
      ENDIF
      CALL MPI_BARRIER(global_comm, mpierr)
      ! figure out my grid sizes
      CALL MPIUTILS_GRD2IJK(myid, parms%ndivx, parms%ndivy, parms%ndivz, &
                            imbx, imby, imbz, ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'Error finding process block'
         RETURN
      ENDIF
      ndx = MAX(parms%nx/ndivx, 1)
      ndy = MAX(parms%ny/ndivy, 1)
      ndz = MAX(parms%nz/ndivz, 1)
      i1 = ndx*imbx + 1 
      i2 = ndx*(imbx + 1)
      j1 = ndy*imby + 1 
      j2 = ndy*(imby + 1)  
      k1 = ndz*imbz + 1 
      k2 = ndz*(imbz + 1)
      IF (imbx + 1 == ndivx) i2 = parms%nx
      IF (imby + 1 == ndivy) j2 = parms%ny
      IF (imbz + 1 == ndivz) k2 = parms%nz
      locate%nxLoc = i2 - i1 + 1
      locate%nyLoc = j2 - j1 + 1 
      locate%nzLoc = k2 - k1 + 1 
      locate%nx = parms%nx
      locate%ny = parms%ny
      locate%nz = parms%nz
      locate%ix0 = i1
      locate%iy0 = j1
      locate%iz0 = k1
      ! verify i got 'em all
      locate%ngrd = locate%nxLoc*locate%nyLoc*locate%nzLoc
      CALL MPI_REDUCE(locate%ngrd, ngrdAll, 1, MPI_INTEGER, MPI_SUM, master, comm, mpierr)
      IF (myid == master) THEN
         IF (nx*ny*nz /= ngrdAll) THEN
            WRITE(*,*) 'locate3d_initialize: Failed to split grid', nall, ngrdAll 
            ierr = 1
         ENDIF
      ENDIF
      ! load the model grid - these are my locations
      ALLOCATE(locate%xlocs(MAX(1, locate%ngrd)))
      ALLOCATE(locate%ylocs(MAX(1, locate%ngrd)))
      ALLOCATE(locate%zlocs(MAX(1, locate%ngrd)))
      locate%xlocs(:) = REAL(zero)
      locate%ylocs(:) = REAL(zero)
      locate%zlocs(:) = REAL(zero)
      CALL H5IO_READ_MODELF(intra_table_comm, tttFileID,                &
                            locate%ix0, locate%iy0, locate%iz0,         &
                            locate%nxLoc, locate%nyLoc, locate%nzLoc,   &
                            locate%xlocs, locate%ylocs, locate%zlocs,   &
                            ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'locate3d_initialize: Error reading model'
         ierr = 1
      ENDIF
      parms%linit = .TRUE.
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE LOCATE3D_FINALIZE() &
                 BIND(C,NAME='locate3d_finalize')
      USE LOCATE_MODULE, ONLY : locate, parms
      IF (ALLOCATED(locate%xlocs)) DEALLOCATE(locate%xlocs)
      IF (ALLOCATED(locate%ylocs)) DEALLOCATE(locate%ylocs)
      IF (ALLOCATED(locate%zlocs)) DEALLOCATE(locate%zlocs)
      parms%linit = .FALSE.
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

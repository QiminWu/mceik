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
!>    @brief Driver routine for earthquake location
!>
!>    @param[in] tttFileID    HDF5 handle for traveltime tables
!>    @param[in] locFileID    HDF5 location file handle
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
      SUBROUTINE LOCATE3D_GRIDSEARCH(tttFileID, locFileID,  &
                                     model,                 &
                                     job, nobs, nevents, &
                                     luseObs, statCor,      &
                                     tori, varobs,          &
                                     tobs, test,            &
                                     hypo, ierr)      &
                 BIND(C,NAME='locate3d_gridSearch')
      USE MPI
      USE HDF5
      USE H5IO_MODULE, ONLY : H5IO_READ_TRAVELTIMESF
      USE LOCATE_MODULE, ONLY : COMPUTE_LOCATION_ONLY,            &
                                COMPUTE_LOCATION_AND_ORIGIN_TIME, &
                                COMPUTE_LOCATION_AND_STATICS,     &
                                COMPUTE_LOCATION_ALL, zero, one, sqrt2i
      USE LOCATE_TYPES, ONLY : locateType
      USE ISO_C_BINDING
      INTEGER(C_INT), INTENT(IN) :: job, locfileID, model, nobs, nevents, tttFileID
      INTEGER(C_INT), INTENT(IN) :: luseObs(nevents*nobs)
      REAL(C_DOUBLE), INTENT(IN) :: statCor(nobs), tori(nevents),       &
                                    varobs(nevents*nobs), tobs(nevents*nobs)
      REAL(C_DOUBLE), INTENT(OUT) :: hypo(4*nevents), test(nevents*nobs)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      DOUBLE PRECISION, ALLOCATABLE :: logPDF(:), logPDFbuf(:), t0(:), t0buf(:)
      REAL, ALLOCATABLE :: test4(:)
      DOUBLE PRECISION res, tobs_i, xnorm, wt_i, wt_i_sqrt2i
      INTEGER igrd, isrc
      INTEGER, PARAMETER :: master = 0
TYPE(locateType) :: locate
      !----------------------------------------------------------------------------------!
      ierr = 0
      hypo(1:4*nevents) = zero
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
         DO 1 isrc=1,nevents
            logPDF(:) = zero
            logPDFbuf(:) = zero
cycle
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
                  myobs = (isrc - 1)*nobs + iobs
                  IF (luseObs(iobs) == 0) CYCLE
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
               IF (SUM(luseObs((isrc-1)*nobs+1:isrc*nobs)) == 0) THEN
                  xnorm = SUM(varobs((isrc-1)*nobs+1:isrc*nobs))
                  CALL DSCAL(ngrd, xnorm, t0buff, 1) ! divide by the sum of the weights
               ENDIF
               ! add up the common residual (which is the origin time at each grid point)
               !CALL MPI_ALLREDUCE(t0buf, t0, MPI_DOUBLE_PRECISION, MPI_SUM, &
               !                   domain_comm, mpierr)
            ELSE
               t0(:) = tori(isrc)
            ENDIF
            ! stack the residuals for this event
            DO 4 iobs=1,nobs,nobs_groups
               myobs = (isrc - 1)*nobs + iobs
               IF (luseObs(iobs) == 0) CYCLE
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
            !CALL MPI_REDUCE(logPDFbuf, logPDF, ngrd, MPI_DOUBLE_PRECISION, MPI_SUM, master, &
            !                cross_comm, mpierr)
            ! now the head group can locate the event
            IF (mydomain_id == master) THEN
!              hypo(4*(isrc-1)+1) = xhypo
!              hypo(4*(isrc-1)+2) = yhypo
!              hypo(4*(isrc-1)+3) = zhypo
!              hypo(4*(isrc-1)+4) = thypo
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
         !parms%locFileID = locFileID
         !parms%tttFileID = tttFileID
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
      nall = locate%nxLoc*locate%nyLoc*locate%nzLoc
      CALL MPI_REDUCE(nall, ngrdAll, 1, MPI_INTEGER, MPI_SUM, master, comm, mpierr)
      IF (myid == master) THEN
         IF (nx*ny*nz /= ngrdAll) THEN
            WRITE(*,*) 'locate3d_initialize: Failed to split grid', nall, ngrdAll 
            ierr = 1
         ENDIF
      ENDIF
      ! load the model grid - these are my locations
      ALLOCATE(locate%xlocs(MAX(1, locate%nxLoc)))
      ALLOCATE(locate%ylocs(MAX(1, locate%nyLoc)))
      ALLOCATE(locate%zlocs(MAX(1, locate%nzLoc)))
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
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE LOCATE3D_FINALIZE() &
                 BIND(C,NAME='locate3d_finalize')

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

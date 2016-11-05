      INTERFACE
         SUBROUTINE LOCATE3D_GRIDSEARCH_DOUBLE64(ldgrd, ngrd, nobs, iwantOT, mask,  &
                                                 tobs, varobs, test, logPDF, ierr)  &
                    BIND(C, NAME='locate3d_gridsearch__double64')
         USE ISO_C_BINDING
         IMPLICIT NONE
         INTEGER(C_INT), INTENT(IN) :: iwantOT, ldgrd, ngrd, nobs
         REAL(C_DOUBLE), INTENT(IN) :: tobs(nobs), varobs(nobs), test(ngrd*nobs)
         INTEGER(C_INT), INTENT(IN) :: mask(nobs) 
         REAL(C_DOUBLE), INTENT(OUT) :: logPDF(ngrd)
         INTEGER(C_INT), INTENT(OUT) :: ierr
         !DIR$ ATTRIBUTES ALIGN: 64 :: logPDF, test, tobs, varobs
         END SUBROUTINE LOCATE3D_GRIDSEARCH_DOUBLE64

         SUBROUTINE LOCATE3D_GRIDSEARCH_FLOAT64(ldgrd, ngrd, nobs, iwantOT, mask,  &
                                                tobs, varobs, test, logPDF, ierr)  &
                    BIND(C, NAME='locate3d_gridsearch__float64')
         USE ISO_C_BINDING
         IMPLICIT NONE
         INTEGER(C_INT), INTENT(IN) :: iwantOT, ldgrd, ngrd, nobs
         REAL(C_FLOAT), INTENT(IN) :: tobs(nobs), varobs(nobs), test(ngrd*nobs)
         INTEGER(C_INT), INTENT(IN) :: mask(nobs)
         REAL(C_FLOAT), INTENT(OUT) :: logPDF(ngrd)
         INTEGER(C_INT), INTENT(OUT) :: ierr
         !DIR$ ATTRIBUTES ALIGN: 64 :: logPDF test, tobs, varobs
         END SUBROUTINE LOCATE3D_GRIDSEARCH_FLOAT64

      END INTERFACE
      double precision, allocatable :: logPDF(:), test(:), tobs(:), varobs(:)
      real, allocatable :: logPDF4(:), test4(:), tobs4(:), varobs4(:) 
      integer, allocatable :: mask(:)
      integer nobs, ngrd
      integer, parameter :: iwantOT = 1
      !DIR$ ATTRIBUTES ALIGN: 64 :: logPDF, test, tobs, varobs 
      !DIR$ ATTRIBUTES ALIGN: 64 :: logPDF4, test4, tobs4, varobs4
      double precision dx, dy, dz, xsrc, ysrc, zsrc
      double precision, allocatable :: xrec(:), yrec(:), zrec(:)
dx = 1.d3
dy = 1.d3
dz = 1.d3
nobs = 20
nx = 79
ny = 71
nz = 15
ixsrc = 31
iysrc = 55
izsrc = 4
xsrc = (ixsrc-1)*dx 
ysrc = (iysrc-1)*dy
zsrc = (izsrc-1)*dz
print *, nx*ny*nz, (izsrc-1)*nx*ny + (iysrc-1)*nx + ixsrc, ixsrc, iysrc, izsrc 
print *, 'true optimum', (izsrc-1)*nx*ny + (iysrc-1)*nx + ixsrc
nobs = 14
allocate(xrec(nobs))
allocate(yrec(nobs))
allocate(zrec(nobs))
allocate(mask(nobs))
mask(:) = 0
call random_number(xrec)
call random_number(yrec)
call random_number(zrec)
      ngrd = nx*ny*nz 
      ldgrd = ngrd + 64 - MOD(ngrd,64)
print *, ldgrd, ngrd, mod(ngrd,64)
      allocate(logPDF(ngrd))
      allocate(test(nobs*ldgrd))
      allocate(tobs(nobs))
      allocate(varobs(nobs)) 
      test(:) = 0.d0
      tobs(:) = 0.d0
do i=1,nobs
   xrec(i) = xrec(i)*dble(nx-1)*dx
   yrec(i) = yrec(i)*dble(ny-1)*dy
   zrec(i) = zrec(i)*dble(nz-1)*dz
   i1 = (i - 1)*ldgrd + 1
   i2 = i*ldgrd 
   call makeTest(nx, ny, nz, ngrd, dx, dy, dz, xrec(i), yrec(i), zrec(i), test(i1:i2))
   call makeObs(xrec(i), yrec(i), zrec(i), xsrc, ysrc, zsrc, tobs(i))
enddo
tobs(:) = tobs(:) + 4.d0
varobs(:) = 1.d0

print *, 'here'
!do i=1,nobs
      CALL LOCATE3D_GRIDSEARCH_DOUBLE64(ldgrd, ngrd, nobs, iwantOT, mask, &
                                        tobs, varobs, test, logPDF, ierr)
print *, 'optimum', minloc(logPDF), 'worst', maxloc(logPDF)
!enddo
      DEALLOCATE(logPDF)
!     DEALLOCATE(test)
!     DEALLOCATE(tobs)
!     DEALLOCATE(varobs)

      allocate(logPDF4(ngrd))
      allocate(test4(nobs*ldgrd))
      allocate(tobs4(nobs))
      allocate(varobs4(nobs)) 
      test4(:) = sngl(test(:))
      tobs4(:) = sngl(tobs(:)) 
      varobs4(:) = sngl(varobs(:))
print *, varobs4
  deallocate(test)
  deallocate(tobs)
  deallocate(varobs) 
print *, 'here2'
!do i=1,nobs
     CALL LOCATE3D_GRIDSEARCH_FLOAT64(ldgrd, ngrd, nobs, iwantOT, mask, &
                                      tobs4, varobs4, test4, logPDF4, ierr )
!enddo
print *, 'done'
      DEALLOCATE(logPDF4)
      DEALLOCATE(test4)
      DEALLOCATE(tobs4)
      DEALLOCATE(varobs4)
      stop
      end

      SUBROUTINE makeTest(nx, ny, nz, ngrd, dx, dy, dz, xsrc, ysrc, zsrc, test)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ngrd, nx, ny, nz
      DOUBLE PRECISION, INTENT(IN) :: dx, dy, dz, xsrc, ysrc, zsrc 
      DOUBLE PRECISION, INTENT(OUT) :: test(ngrd) 
      DOUBLE PRECISION d2, x, y, z
      INTEGER igrd, i, j, k
      DOUBLE PRECISION, PARAMETER :: slow = 1.d0/5.d3
      DO 1 k=1,nz
         DO 2 j=1,ny
            DO 3 i=1,nx
               igrd = (k-1)*nx*ny + (j-1)*nx + i
               x = DBLE(i - 1)*dx
               y = DBLE(j - 1)*dy
               z = DBLE(k - 1)*dz
               d2 = (x - xsrc)**2 + (y - ysrc)**2 + (z - zsrc)**2
               test(igrd) = SQRT(d2)*slow 
    3       CONTINUE
    2    CONTINUE
    1 CONTINUE

      RETURN
      END 

      SUBROUTINE makeObs(x, y, z, xsrc, ysrc, zsrc, tobs)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: x, y, z, xsrc, ysrc, zsrc
      DOUBLE PRECISION, INTENT(OUT) :: tobs
      DOUBLE PRECISION d2
      DOUBLE PRECISION, PARAMETER :: slow = 1.d0/5.d3
      d2 = (x - xsrc)**2 + (y - ysrc)**2 + (z - zsrc)**2
      tobs = SQRT(d2)*slow 
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Stacks the weighted residuals into the analytic origin time
!>           computation.  This is for analytic removal of the sourc time
!>           in an L2 optimization as described by Moser et. al. 1992 Eqn 19
!>           for diagonal weights.
!> 
!>    @param[in] ngrd     number of points in grid search
!>    @param[in] tobs_i   observed pick time (seconds)
!>    @param[in] xnorm    normalization factor which is the inverse of the
!>                        sum of the data variances for all observations.
!>    @param[in] var_i    variance in the pick time (seconds)
!>    @param[in] test     estimate arrival times at all points in grid
!>
!>    @param[in,out] t0   on input contains the current weighted residual
!>                        sum for previous observations
!>                        on output contains the updated weighted residual
!>                        sum which incorporates the current observation.
!>
!>    @author Ben Baker
!>
!>    @copyright Apache 2
!> 
      SUBROUTINE LOCATE3D_STACK_T0_DOUBLE64(ngrd, tobs_i, xnorm, var_i, test, t0)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ngrd
      DOUBLE PRECISION, INTENT(IN) :: test(ngrd), xnorm, tobs_i, var_i
      DOUBLE PRECISION, INTENT(INOUT) :: t0(ngrd)
      DOUBLE PRECISION tobs, wt
      INTEGER igrd
      DOUBLE PRECISION, PARAMETER :: one = 1.d0
      !DIR$ ATTRIBUTES ALIGN: 64 :: igrd, tobs, wt
      wt = one/(var_i*xnorm)
      tobs = tobs_i
      !$OMP PARALLEL DO SIMD FIRSTPRIVATE(tobs, wt), &
      !$OMP SHARED(ngrd, t0, test), DEFAULT(none)
      DO 1 igrd=1,ngrd
         t0(igrd) = t0(igrd) + wt*(tobs - test(igrd)) 
    1 CONTINUE
      !$OMP END PARALLEL DO SIMD
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Stacks the weighted residuals into the analytic origin time
!>           computation.  This is for analytic removal of the sourc time
!>           in an L2 optimization as described by Moser et. al. 1992 Eqn 19
!>           for diagonal weights. 
!> 
!>    @param[in] ngrd     number of points in grid search
!>    @param[in] tobs_i   observed pick time (seconds)
!>    @param[in] xnorm    normalization factor which is the inverse of the
!>                        sum of the data variances for all observations.
!>    @param[in] var_i    variance in the pick time (seconds)
!>    @param[in] test     estimate arrival times at all points in grid
!>
!>    @param[in,out] t0   on input contains the current weighted residual
!>                        sum for previous observations
!>                        on output contains the updated weighted residual
!>                        sum which incorporates the current observation.
!>
!>    @author Ben Baker
!>
!>    @copyright Apache 2
!>
      SUBROUTINE LOCATE3D_STACK_T0_FLOAT64(ngrd, tobs_i, xnorm, var_i, test, t0) 
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ngrd
      REAL, INTENT(IN) :: test(ngrd), xnorm, tobs_i, var_i
      REAL, INTENT(INOUT) :: t0(ngrd) 
      REAL tobs, wt
      INTEGER igrd
      REAL, PARAMETER :: one = 1.0
      !DIR$ ATTRIBUTES ALIGN: 64 :: igrd, tobs, wt
      wt = one/(var_i*xnorm)
      tobs = tobs_i
      !$OMP PARALLEL DO SIMD FIRSTPRIVATE(tobs, wt), &
      !$OMP SHARED(ngrd, t0, test), DEFAULT(none)
      DO 1 igrd=1,ngrd
         t0(igrd) = t0(igrd) + wt*(tobs - test(igrd))
    1 CONTINUE
      !$OMP END PARALLEL DO SIMD
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Stacks the squared residuals into the L2 penalty function.  This
!>           is for diagonal weighting.
!>
!>    @param[in] ngrd     number of points in grid search
!>    @param[in] tobs_i   i'th pick time (seconds)
!>    @param[in] var_i    variance (seconds) in i'th observation
!>    @param[in] test     theoretical traveltime to each point in grid
!>                        (seconds) [ngrd]
!>    @param[in] t0       analytic origin time (seconds) at all points
!>                        in grid [ngrd]
!>
!>    @param[in,out]      on input contains the weighted squared residuals
!>                        at all points in the grid.
!>                        on output contains the contribution of the current
!>                        squared residual to all points in the grid. [ngrd]
!>
!>    @author Ben Baker
!>
!>    @copyright Apache 2
!>
      SUBROUTINE LOCATE3D_STACK_LOGPDF_DOUBLE64(ngrd, tobs_i, var_i, test, &
                                                t0, logPDF)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ngrd
      DOUBLE PRECISION, INTENT(IN) :: t0(ngrd), test(ngrd), tobs_i, var_i
      DOUBLE PRECISION, INTENT(INOUT) :: logPDF(ngrd)
      DOUBLE PRECISION res, tobs, wt
      INTEGER igrd
      DOUBLE PRECISION, PARAMETER :: one = 1.d0
      DOUBLE PRECISION, PARAMETER :: two = 2.d0
      DOUBLE PRECISION, PARAMETER :: sqrt2i = one/SQRT(two)
      !DIR$ ATTRIBUTES ALIGN: 64 :: igrd, res, tobs, wt
      wt = sqrt2i/var_i ! accounts for 1/2 half scaling in L2 norm
      tobs = tobs_i
      !$OMP PARALLEL DO SIMD FIRSTPRIVATE(tobs, wt), PRIVATE(res), &
      !$OMP SHARED(logPDF, ngrd, t0, test), DEFAULT(none)
      DO 1 igrd=1,ngrd
         res = wt*(tobs - (test(igrd) + t0(igrd)))
         logPDF(igrd) = logPDF(igrd) + res*res
    1 CONTINUE
      !$OMP END PARALLEL DO SIMD
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Stacks the squared residuals into the L2 penalty function.  This
!>           is for diagonal weighting.
!>
!>    @param[in] ngrd     number of points in grid search
!>    @param[in] tobs_i   i'th pick time (seconds)
!>    @param[in] var_i    variance (seconds) in i'th observation
!>    @param[in] test     theoretical traveltime to each point in grid
!>                        (seconds) [ngrd]
!>    @param[in] t0       analytic origin time (seconds) at all points
!>                        in grid [ngrd]
!>
!>    @param[in,out]      on input contains the weighted squared residuals
!>                        at all points in the grid.
!>                        on output contains the contribution of the current
!>                        squared residual to all points in the grid. [ngrd]
!>
!>    @author Ben Baker
!>
!>    @copyright Apache 2
!>
      SUBROUTINE LOCATE3D_STACK_LOGPDF_FLOAT64(ngrd, tobs_i, var_i, test, &
                                               t0, logPDF)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ngrd
      REAL, INTENT(IN) :: t0(ngrd), test(ngrd), tobs_i, var_i
      REAL, INTENT(OUT) :: logPDF(ngrd)
      REAL res, tobs, wt
      INTEGER igrd
      REAL, PARAMETER :: one = 1.0
      REAL, PARAMETER :: two = 2.0
      REAL, PARAMETER :: sqrt2i = one/SQRT(two)
      !DIR$ ATTRIBUTES ALIGN: 64 :: igrd, res, tobs, wt
      wt = sqrt2i/var_i ! accounts for 1/2 half scaling in L2 norm
      tobs = tobs_i
      !$OMP PARALLEL DO SIMD FIRSTPRIVATE(tobs, wt), PRIVATE(res) &
      !$OMP SHARED(logPDF, ngrd, t0, test), DEFAULT(none)
      DO 1 igrd=1,ngrd
         res = wt*(tobs - (test(igrd) + t0(igrd)))
         logPDF(igrd) = logPDF(igrd) + res*res
    1 CONTINUE
      !$OMP END PARALLEL DO SIMD
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Zeros out a 64 bit aligned array x
!>
!>    @param[in] n    number of points in array x
!>
!>    @param[out] x   nulled out array [n]
!>
      SUBROUTINE LOCATE3D_NULL_DOUBLE64(n, x)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      DOUBLE PRECISION, INTENT(OUT) :: x(n)
      INTEGER i
      DOUBLE PRECISION, PARAMETER :: zero = 0.d0
      !DIR$ ATTRIBUTES ALIGN: 64 :: i, zero
      !$OMP PARALLEL DO SIMD SHARED(n, x), DEFAULT(none)
      DO 1 i=1,n
         x(i) = zero
    1 CONTINUE
      !$OMP END PARALLEL DO SIMD 
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Zeros out a 64 bit aligned array x
!>
!>    @param[in] n    number of points in array x
!>
!>    @param[out] x   nulled out array [n]
!>
      SUBROUTINE LOCATE3D_NULL_FLOAT64(n, x)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      REAL, INTENT(OUT) :: x(n)
      INTEGER i
      REAL, PARAMETER :: zero = 0.0
      !DIR$ ATTRIBUTES ALIGN: 64 :: i, zero
      !$OMP PARALLEL DO SIMD SHARED(n, x), DEFAULT(none)
      DO 1 i=1,n
         x(i) = zero
    1 CONTINUE
      !$OMP END PARALLEL DO SIMD 
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!>    @brief Performs the grid search over all estimates for aligned data
!>
!>    
      SUBROUTINE LOCATE3D_GRIDSEARCH_DOUBLE64(ldgrd, ngrd, nobs, iwantOT, mask,  &
                                              tobs, varobs, test, logPDF, ierr)  &
                 BIND(C, NAME='locate3d_gridsearch__double64')
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: iwantOT, ldgrd, ngrd, nobs
      REAL(C_DOUBLE), INTENT(IN) :: tobs(nobs), varobs(nobs), test(ngrd*nobs)
      INTEGER(C_INT), INTENT(IN) :: mask(nobs) 
      REAL(C_DOUBLE), INTENT(OUT) :: logPDF(ngrd)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      !DIR$ ATTRIBUTES ALIGN: 64 :: logPDF, test, tobs, varobs
      ! local variables
      DOUBLE PRECISION, ALLOCATABLE :: t0(:)
      DOUBLE PRECISION tobs_i, xnorm, var_i
      INTEGER ibeg, igrd, iobs, iopt
      DOUBLE PRECISION, PARAMETER :: zero = 0.d0
      DOUBLE PRECISION, PARAMETER :: one = 1.d0
      !DIR$ ATTRIBUTES ALIGN: 64 :: t0, xnorm, zero
      !----------------------------------------------------------------------------------!
      !
      ! check input variables
      ierr = 0
      IF (MOD(ldgrd, 64) /= 0) THEN
         WRITE(*,*) 'locate3d_gridsearch_double64: Require arrays be 64 byte aligned', &
                    ldgrd, MOD(8*ldgrd, 64)
         ierr = 1
         RETURN
      ENDIF
      IF (ngrd > ldgrd) THEN
         WRITE(*,*) 'locate3d_gridsearch_double64: ngrd cannot be greater than ldgrd'
         ierr = 1
         RETURN
      ENDIF
      IF (SUM(mask) == nobs) THEN
         WRITE(*,*) 'locate3d_gridsearch_double64: No observations'
         ierr = 1
         RETURN
      ENDIF
      IF (ABS(SUM(varobs) - zero) < EPSILON(one)) THEN
         WRITE(*,*) 'locate3d_gridsearch_double64: Will be division by zero'
         ierr = 1
         RETURN
      ENDIF
      ! initialize and set space
      ALLOCATE(t0(ngrd))
      CALL LOCATE3D_NULL_DOUBLE64(ngrd, t0)
      CALL LOCATE3D_NULL_DOUBLE64(ngrd, logPDF)
      ! compute the origin time in this grid
      IF (iwantOT == 1) THEN
         xnorm = zero
         DO 3 iobs=1,nobs
            IF (mask(iobs) /= 1) xnorm = xnorm + varobs(iobs) 
    3    CONTINUE
         DO 4 iobs=1,nobs
            IF (mask(iobs) == 1) CYCLE
            tobs_i = tobs(iobs)
            var_i = varobs(iobs)
            ibeg = ldgrd*(iobs - 1) + 1
            CALL LOCATE3D_STACK_T0_DOUBLE64(ngrd, tobs_i, xnorm, var_i, test(ibeg), t0)
    4    CONTINUE 
      ENDIF
      ! now compute the locations with the origin times at each grid point
      DO 11 iobs=1,nobs
         IF (mask(iobs) == 1) CYCLE
         tobs_i = tobs(iobs)
         var_i = varobs(iobs)
         ibeg = ldgrd*(iobs - 1) + 1
         CALL LOCATE3D_STACK_LOGPDF_DOUBLE64(ngrd, tobs_i, var_i, test(ibeg), &
                                             t0, logPDF)
   11 CONTINUE
      ! get the origin time
      IF (iwantOT == 1) THEN
         iopt = MINLOC(logPDF, 1)
print *, iopt, t0(iopt)
      ENDIF
      DEALLOCATE(t0)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE LOCATE3D_GRIDSEARCH_FLOAT64(ldgrd, ngrd, nobs, iwantOT, mask,  &
                                             tobs, varobs, test, logPDF, ierr)  &
                 BIND(C, NAME='locate3d_gridsearch__float64')
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), INTENT(IN) :: iwantOT, ldgrd, ngrd, nobs
      REAL(C_FLOAT), INTENT(IN) :: tobs(nobs), varobs(nobs), test(ngrd*nobs)
      INTEGER(C_INT), INTENT(IN) :: mask(nobs)
      REAL(C_FLOAT), INTENT(OUT) :: logPDF(ngrd)
      INTEGER(C_INT), INTENT(OUT) :: ierr
      !DIR$ ATTRIBUTES ALIGN: 64 :: logPDF, test, tobs, varobs 
      ! local variables
      REAL, ALLOCATABLE :: t0(:)
      REAL tobs_i, xnorm, var_i
      INTEGER ibeg, igrd, iobs, iopt
      REAL, PARAMETER :: zero = 0.0
      REAL, PARAMETER :: one = 1.0
      !DIR$ ATTRIBUTES ALIGN: 64 :: t0, xnorm, zero
      !----------------------------------------------------------------------------------!
      !   
      ! check input variables
      ierr = 0 
      IF (MOD(ldgrd, 64) /= 0) THEN
         WRITE(*,*) 'locate3d_gridsearch_float64: Require arrays be 64 byte aligned', &
                    ldgrd, MOD(ldgrd, 64)
         ierr = 1 
         RETURN
      ENDIF
      IF (ngrd > ldgrd) THEN
         WRITE(*,*) 'locate3d_gridsearch_float64: ngrd cannot be greater than ldgrd'
         ierr = 1
         RETURN
      ENDIF
      IF (SUM(mask) == nobs) THEN
         WRITE(*,*) 'locate3d_gridsearch_float64: No observations'
         ierr = 1
         RETURN
      ENDIF
      IF (ABS(SUM(varobs) - zero) < EPSILON(one)) THEN
         WRITE(*,*) 'locate3d_gridsearch_float64: Will be division by zero'
         ierr = 1
         RETURN
      ENDIF
      ! initialize and set space
      ALLOCATE(t0(ngrd))
      CALL LOCATE3D_NULL_FLOAT64(ngrd, t0)
      CALL LOCATE3D_NULL_FLOAT64(ngrd, logPDF)
      ! compute the origin time in this grid
      IF (iwantOT == 1) THEN
         xnorm = zero
         DO 3 iobs=1,nobs
            IF (mask(iobs) /= 1) xnorm = xnorm + varobs(iobs)
    3    CONTINUE
         DO 4 iobs=1,nobs
            IF (mask(iobs) == 1) CYCLE
            tobs_i = tobs(iobs)
            var_i = varobs(iobs)
            ibeg = ldgrd*(iobs - 1) + 1
            CALL LOCATE3D_STACK_T0_FLOAT64(ngrd, tobs_i, xnorm, var_i, test(ibeg), t0)
    4    CONTINUE
      ENDIF
      ! now compute the locations with the origin times at each grid point
      DO 11 iobs=1,nobs
         IF (mask(iobs) == 1) CYCLE
         tobs_i = tobs(iobs)
         var_i = varobs(iobs)
         ibeg = ldgrd*(iobs - 1) + 1
         CALL LOCATE3D_STACK_LOGPDF_FLOAT64(ngrd, tobs_i, var_i, test(ibeg), &
                                            t0, logPDF)
   11 CONTINUE
      ! get the origin time
      IF (iwantOT == 1) THEN
         iopt = MINLOC(logPDF, 1)
print *, iopt, t0(iopt)
      ENDIF
      DEALLOCATE(t0)
      RETURN
      END

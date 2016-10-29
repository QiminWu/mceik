      double precision, allocatable :: logPDF(:), test(:), tobs(:), varobs(:)
      real, allocatable :: logPDF4(:), test4(:), tobs4(:), varobs4(:) 
      integer nobs, ngrd
      integer, parameter :: iwantOT = 1
      !DIR$ ATTRIBUTES ALIGN: 64 :: logPDF, test, tobs, varobs 
      !DIR$ ATTRIBUTES ALIGN: 64 :: logPDF4, test4, tobs4, varobs4
      nobs = 20
      ngrd = 14*71*78
      ldgrd = ngrd + 64 - MOD(ngrd,64)
print *, ldgrd, ngrd, mod(ngrd,64)
      allocate(logPDF(ngrd))
      allocate(test(nobs*ldgrd))
      allocate(tobs(nobs))
      allocate(varobs(nobs)) 
      test(:) = 0.d0
      tobs(:) = 0.d0
      varobs(:) = 1.d0
print *, 'here'
do i=1,40
      CALL LOCATE3D_GRIDSEARCH_DOUBLE64(ldgrd, ngrd, nobs, iwantOT, &
                                        tobs, varobs, test, logPDF, ierr )
enddo
      DEALLOCATE(logPDF)
      DEALLOCATE(test)
      DEALLOCATE(tobs)
      DEALLOCATE(varobs)

      allocate(logPDF4(ngrd))
      allocate(test4(nobs*ldgrd))
      allocate(tobs4(nobs))
      allocate(varobs4(nobs)) 
      test4(:) = 0.d0
      tobs4(:) = 0.d0
      varobs4(:) = 1.d0
print *, 'here2'
do i=1,40
     CALL LOCATE3D_GRIDSEARCH_FLOAT64(ldgrd, ngrd, nobs, iwantOT, &
                                      tobs4, varobs4, test4, logPDF4, ierr )
enddo
print *, 'done'
      DEALLOCATE(logPDF4)
      DEALLOCATE(test4)
      DEALLOCATE(tobs4)
      DEALLOCATE(varobs4)
      stop
      end

      SUBROUTINE LOCATE3D_STACK_T0_DOUBLE64(ngrd, tobs_i, var_i, test, t0)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ngrd
      DOUBLE PRECISION, INTENT(IN) :: test(ngrd), tobs_i, var_i
      DOUBLE PRECISION, INTENT(INOUT) :: t0(ngrd) 
      DOUBLE PRECISION res, tobs, wt
      INTEGER igrd
      DOUBLE PRECISION, PARAMETER :: one = 1.d0
      !DIR$ ATTRIBUTES ALIGN: 64 :: igrd, res, tobs, wt
      wt = one/var_i
      tobs = tobs_i
      !$OMP PARALLEL DO SIMD
      DO 1 igrd=1,ngrd
         res = wt*(tobs - test(igrd))
         t0(igrd) = t0(igrd) + res*res
    1 CONTINUE
      !$OMP END PARALLEL DO SIMD
      RETURN
      END

      SUBROUTINE LOCATE3D_STACK_T0_FLOAT64(ngrd, tobs_i, var_i, test, t0) 
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ngrd
      REAL, INTENT(IN) :: test(ngrd), tobs_i, var_i
      REAL, INTENT(INOUT) :: t0(ngrd) 
      REAL res, tobs, wt
      INTEGER igrd
      REAL, PARAMETER :: one = 1.0
      !DIR$ ATTRIBUTES ALIGN: 64 :: igrd, res, tobs, wt
      wt = one/var_i
      tobs = tobs_i
      !$OMP PARALLEL DO SIMD
      DO 1 igrd=1,ngrd
         res = wt*(tobs - test(igrd))
         t0(igrd) = t0(igrd) + res*res
    1 CONTINUE
      !$OMP END PARALLEL DO SIMD
      RETURN
      END 

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
      !$OMP PARALLEL DO SIMD
      DO 1 igrd=1,ngrd
         res = wt*(tobs - (test(igrd) + t0(igrd)))
         logPDF(igrd) = logPDF(igrd) - res*res
    1 CONTINUE
      !$OMP END PARALLEL DO SIMD
      RETURN
      END

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
      !$OMP PARALLEL DO SIMD
      DO 1 igrd=1,ngrd
         res = wt*(tobs - (test(igrd) + t0(igrd)))
         logPDF(igrd) = logPDF(igrd) - res*res
    1 CONTINUE
      !$OMP END PARALLEL DO SIMD
      RETURN
      END


      SUBROUTINE LOCATE3D_NULL_DOUBLE64(n, x)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      DOUBLE PRECISION, INTENT(OUT) :: x(n)
      INTEGER i
      DOUBLE PRECISION, PARAMETER :: zero = 0.d0
      !DIR$ ATTRIBUTES ALIGN: 64 :: i, zero
      !$OMP PARALLEL DO SIMD
      DO 1 i=1,n
         x(i) = zero
    1 CONTINUE
      !$OMP END PARALLEL DO SIMD 
      RETURN
      END

      SUBROUTINE LOCATE3D_NULL_FLOAT64(n, x)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      REAL, INTENT(OUT) :: x(n)
      INTEGER i
      REAL, PARAMETER :: zero = 0.0
      !DIR$ ATTRIBUTES ALIGN: 64 :: i, zero
      !$OMP PARALLEL DO SIMD
      DO 1 i=1,n
         x(i) = zero
    1 CONTINUE
      !$OMP END PARALLEL DO SIMD 
      RETURN
      END


      SUBROUTINE LOCATE3D_GRIDSEARCH_DOUBLE64(ldgrd, ngrd, nobs, iwantOT, &
                                              tobs, varobs, test, logPDF, ierr )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iwantOT, ldgrd, ngrd, nobs
      DOUBLE PRECISION, INTENT(IN) :: tobs(nobs), varobs(nobs), test(ngrd*nobs)
      DOUBLE PRECISION, INTENT(OUT) :: logPDF(ngrd)
      INTEGER, INTENT(OUT) :: ierr
      ! local variables
      DOUBLE PRECISION, ALLOCATABLE :: t0(:)
      DOUBLE PRECISION tobs_i, xnorm, var_i
      INTEGER ibeg, igrd, iobs
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
      ! initialize and set space
      ALLOCATE(t0(ldgrd))
      CALL LOCATE3D_NULL_DOUBLE64(ngrd, t0)
      CALL LOCATE3D_NULL_DOUBLE64(ngrd, logPDF)
      ! compute the origin time in this grid
      IF (iwantOT == 1) THEN
         DO 3 iobs=1,nobs
            tobs_i = tobs(iobs)
            var_i = varobs(iobs)
            ibeg = ldgrd*(iobs - 1)
            CALL LOCATE3D_STACK_T0_DOUBLE64(ngrd, tobs_i, var_i, test(ibeg), t0)
    3    CONTINUE 
         xnorm = SUM(varobs)
         ! normalize by the sum of the weights which is equivalent to scaling by the
         ! sum of the variances 
         !$OMP DO SIMD ALIGNED(t0, xnorm: 64)
         DO 4 igrd=1,ngrd
            t0(igrd) = xnorm*t0(igrd)
    4    CONTINUE
         !$OMP END DO SIMD 
      ENDIF
      ! now compute the locations with the origin times at each grid point
      DO 11 iobs=1,nobs
         tobs_i = tobs(iobs)
         var_i = varobs(iobs)
         CALL LOCATE3D_STACK_LOGPDF_DOUBLE64(ngrd, tobs_i, var_i, test(ibeg), &
                                             t0, logPDF)
   11 CONTINUE
      DEALLOCATE(t0)
      RETURN
      END

      SUBROUTINE LOCATE3D_GRIDSEARCH_FLOAT64(ldgrd, ngrd, nobs, iwantOT, &
                                             tobs, varobs, test, logPDF, ierr )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iwantOT, ldgrd, ngrd, nobs
      REAL, INTENT(IN) :: tobs(nobs), varobs(nobs), test(ngrd*nobs)
      REAL, INTENT(OUT) :: logPDF(ngrd)
      INTEGER, INTENT(OUT) :: ierr
      ! local variables
      REAL, ALLOCATABLE :: t0(:)
      REAL tobs_i, xnorm, var_i
      INTEGER ibeg, igrd, iobs
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
      ! initialize and set space
      ALLOCATE(t0(ldgrd))
      CALL LOCATE3D_NULL_FLOAT64(ngrd, t0)
      CALL LOCATE3D_NULL_FLOAT64(ngrd, logPDF)
      ! compute the origin time in this grid
      IF (iwantOT == 1) THEN
         DO 3 iobs=1,nobs
            tobs_i = tobs(iobs)
            var_i = varobs(iobs)
            ibeg = ldgrd*(iobs - 1)
            CALL LOCATE3D_STACK_T0_FLOAT64(ngrd, tobs_i, var_i, test(ibeg), t0)
    3    CONTINUE
         xnorm = SUM(varobs)
         ! normalize by the sum of the weights which is equivalent to scaling by the
         ! sum of the variances 
         !$OMP DO SIMD ALIGNED(t0, xnorm: 64)
         DO 4 igrd=1,ngrd
            t0(igrd) = xnorm*t0(igrd)
    4    CONTINUE
         !$OMP END DO SIMD 
      ENDIF
      ! now compute the locations with the origin times at each grid point
      DO 11 iobs=1,nobs
         tobs_i = tobs(iobs)
         var_i = varobs(iobs)
         CALL LOCATE3D_STACK_LOGPDF_FLOAT64(ngrd, tobs_i, var_i, test(ibeg), &
                                            t0, logPDF)
   11 CONTINUE
      DEALLOCATE(t0)
      RETURN
      END

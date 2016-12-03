#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>
#include <omp.h>
#include <lapacke_utils.h>
#include <mpi.h>

#define memory_isAligned(POINTER, BYTE_COUNT) \
    (((uintptr_t)(const void *)(POINTER)) % (BYTE_COUNT) == 0)

int locate_l2_gridSearch__double64(const int ldgrd,
                                   const int ngrd,
                                   const int nobs,
                                   const int iwantOT,
                                   const double t0use,
                                   const int *__restrict__ mask,
                                   const double *__restrict__ tobs,
                                   const double *__restrict__ tcorr,
                                   const double *__restrict__ varobs,
                                   const double *__restrict__ test,
                                   double *__restrict__ t0, 
                                   double *__restrict__ objfn);
int locate_l2_gridSearch__float64(const int ldgrd,
                                  const int ngrd,
                                  const int nobs,
                                  const int iwantOT,
                                  const float t0use,
                                  const int *__restrict__ mask,
                                  const float *__restrict__ tobs,
                                  const float *__restrict__ tcorr,
                                  const float *__restrict__ varobs,
                                  const float *__restrict__ test,
                                  float *__restrict__ t0, 
                                  float *__restrict__ objfn);
int locate_l1_gridSearch__double64(const int ldgrd,
                                   const int ngrd,
                                   const int nobs,
                                   const int iwantOT,
                                   const double t0use,
                                   const int *__restrict__ mask,
                                   const double *__restrict__ tobs,
                                   const double *__restrict__ varobs,
                                   const double *__restrict__ test,
                                   double *__restrict__ t0,
                                   double *__restrict__ objfn);
static void locate_l2_stackT0__double64(const int ngrd,
                                        const double tobs_i,
                                        const double xnorm, 
                                        const double wt_i,
                                        const double *__restrict__ test,
                                        double *__restrict__ t0);
static void locate_l2_stackT0__float64(const int ngrd,
                                       const float tobs_i,
                                       const float xnorm,
                                       const float wt_i,
                                       const float *__restrict__ test,
                                       float *__restrict__ t0);
static void locate_l2_stackObjfn__double64(const int ngrd,
                                           const double tobs_i,
                                           const double wt_i,
                                           const double *__restrict__ test,
                                           const double *__restrict__ t0,
                                           double *__restrict__ objfn);
static void locate_l2_stackObjfn__float64(const int ngrd,
                                          const float tobs_i,
                                          const float wt_i,
                                          const float *__restrict__ test,
                                          const float *__restrict__ t0, 
                                          float *__restrict__ objfn);
double weightedMedian__double(const int n,
                              const double *__restrict__ x,
                              const double *__restrict__ w,
                              int *__restrict__ perm,
                              bool *lsort, int *ierr);
void locate_setDouble64(const int n, const double x0in,
                         double *__restrict__ x);
void locate_setFloat64(const int n, const float x0in,
                       float *__restrict__ x);
static void locate_nullDouble64(const int n, double *__restrict__ x);
void locate_nullFloat64(const int n, float *__restrict__ x);
double locate_sumDouble64(const int n, const double *__restrict__ x);
float locate_sumFloat64(const int n, const float *__restrict__ x);
int locate_sumInt(const int n, const int *__restrict__ x);
int locate_minLocDouble64(const int n, const double *__restrict__ x);
int locate_minLocFloat64(const int n, const float *__restrict__ x);
int memory_malloc__double(const int n, double **x);
int memory_malloc__float(const int n, float **x);
int memory_malloc__int(const int n, int **x);
double *memory_calloc__double(const int n);
float *memory_calloc__float(const int n);
int *memory_calloc__int(const int n);
void makeTest(const int nx, const int ny, const int nz,
              const double dx, const double dy, const double dz,
              const double xsrc, const double ysrc, const double zsrc,
              double *__restrict__ test);
double makeObs(const double x, const double y, const double z,
               const double xsrc, const double ysrc, const double zsrc);

int main2()
{

    return 0;
}

int main()
{
    double *objfn, *xrec, *yrec, *zrec, *t0, *tcorr, *test, *tobs, *varobs;
    float *objfn4, *t04, *tcorr4, *test4, *tobs4, *varobs4;
    int *mask;
    double dx, dy, dz, t0use, xsrc, ysrc, zsrc;
    float t0use4;
    int igrd, iobs, iopt, ixsrc, iysrc, izsrc, iwantOT,
        ldgrd, ngrd, nobs, nx, ny, nz;
    bool lsort;
    srand(4042);
    nobs = 20;
    iwantOT = 1;
    t0use = 4.0;
    dx = 1.e3;
    dy = 1.e3;
    dz = 1.e3;
    nx = 145;
    ny = 145;
    nz = 45;
    ngrd = nx*ny*nz;
    ldgrd = ngrd + 64 - ngrd%64;
    ixsrc = 12;
    iysrc = 15;
    izsrc = 5;
printf("true: %d\n", izsrc*nx*ny + iysrc*nx + ixsrc);
    xsrc = (double) ixsrc*dx;
    ysrc = (double) iysrc*dy;
    zsrc = (double) izsrc*dz;
    xrec = memory_calloc__double(nobs);
    yrec = memory_calloc__double(nobs);
    zrec = memory_calloc__double(nobs);

    tobs = memory_calloc__double(nobs);
    tcorr = memory_calloc__double(nobs);
    varobs = memory_calloc__double(nobs);
    test = memory_calloc__double(nobs*ldgrd);
    mask = (int *)calloc((size_t) nobs, sizeof(int));
    // Make the receiver locations
    for (iobs=0; iobs<nobs; iobs++)
    {
        xrec[iobs] = ((double) (rand()))/RAND_MAX;
        yrec[iobs] = ((double) (rand()))/RAND_MAX;
        zrec[iobs] = ((double) (rand()))/RAND_MAX;
        xrec[iobs] = xrec[iobs]*(double) (nx - 1)*dx;
        yrec[iobs] = yrec[iobs]*(double) (ny - 1)*dy;
        zrec[iobs] = zrec[iobs]*(double) (nz - 1)*dz;
        tobs[iobs] = makeObs(xrec[iobs], yrec[iobs], zrec[iobs],
                             xsrc, ysrc, zsrc);
        varobs[iobs] = ((double) (rand()))/RAND_MAX; //0.2;
        makeTest(nx, ny, nz, dx, dy, dz,
                 xrec[iobs], yrec[iobs], zrec[iobs],
                 &test[ldgrd*iobs]);
    }
    // Add the origin time
    if (iwantOT == 1)
    {
        for (iobs=0; iobs<nobs; iobs++)
        {
            tobs[iobs] = tobs[iobs] + t0use;
        }
    }
    t0 = memory_calloc__double(ngrd);
    objfn = memory_calloc__double(ngrd);
    locate_l1_gridSearch__double64(ldgrd, ngrd, nobs, iwantOT,
                                   t0use, mask, tobs, varobs,
                                   test, t0, objfn);
    iopt = locate_minLocDouble64(ngrd, objfn);
    printf("l1 first estimate: %d %f\n", iopt, t0[iopt]);

    locate_l2_gridSearch__double64(ldgrd, ngrd, nobs, iwantOT, 
                                   t0use, mask, tobs, tcorr, varobs,
                                   test, t0, objfn);
    iopt = locate_minLocDouble64(ngrd, objfn);
    printf("double estimate: %d %f\n", iopt, t0[iopt]);
    locate_l1_gridSearch__double64(ldgrd, ngrd, nobs, iwantOT,
                                   t0use, mask, tobs, varobs,
                                   test, t0, objfn);
    iopt = locate_minLocDouble64(ngrd, objfn);
    printf("double l1 testimate: %d %f\n", iopt, t0[iopt]);
    free(objfn);
    free(t0);
    // Set the float problem
    t0use4 = (float) t0use;
    tobs4 = memory_calloc__float(nobs);
    tcorr4 = memory_calloc__float(nobs);
    varobs4 = memory_calloc__float(nobs);
    test4 = memory_calloc__float(ldgrd*nobs);
    for (iobs=0; iobs<nobs; iobs++)
    {
        tobs4[iobs] = (float) tobs[iobs];
        varobs4[iobs] = (float) varobs[iobs];
        for (igrd=0; igrd<ngrd; igrd++)
        {
            test4[ldgrd*iobs+igrd] = (float) test[ldgrd*iobs+igrd];
        }
    }
    free(test);
    free(tobs);
    free(tcorr);
    free(varobs);
    objfn4 = memory_calloc__float(ngrd);
    t04 = memory_calloc__float(ngrd);
    locate_l2_gridSearch__float64(ldgrd, ngrd, nobs, iwantOT,
                                  t0use4, mask, tobs4, tcorr4, varobs4,  
                                  test4, t04, objfn4);
    iopt = locate_minLocFloat64(ngrd, objfn4);
    printf("float estimate: %d\n", iopt);

    free(test4);
    free(objfn4);
    free(t04);
    free(tobs4);
    free(tcorr4);
    free(varobs4);

    free(mask);
    free(xrec);
    free(yrec);
    free(zrec);
double xs[7] = {0.1, 0.35, 0.05, 0.1, 0.15, 0.05, 0.2};
double ws[7] = {0.1, 0.35, 0.05, 0.1, 0.15, 0.05, 0.2};
double w1[7] = {1.0, 1.0,  1.0,  1.0,  1.0, 1.0,  1.0};
int perm[7];
double wmed;
int ierr;
wmed = weightedMedian__double(7, xs, ws, perm, &lsort, &ierr);
printf("%f\n", wmed);
wmed = weightedMedian__double(7, xs, w1, perm, &lsort, &ierr);
printf("%f\n", wmed);
    return 0;
}

void makeTest(const int nx, const int ny, const int nz,
              const double dx, const double dy, const double dz,
              const double xsrc, const double ysrc, const double zsrc,
              double *__restrict__ test)
{
    double d2, x, y, z;
    int i, igrd, j, k;
    const double slow = 1.0/5000.0;
    for (k=0; k<nz; k++)
    {
        for (j=0; j<ny; j++)
        {
            for (i=0; i<nx; i++)
            {
                igrd = k*nx*ny + j*nx + i;
                x = (double) i*dx;
                y = (double) j*dy;
                z = (double) k*dz;
                d2 = pow(x-xsrc, 2) + pow(y-ysrc, 2) + pow(z-zsrc, 2);
                test[igrd] = sqrt(d2)*slow;
            }
        }
    }
    return;
}

double makeObs(const double x, const double y, const double z,
               const double xsrc, const double ysrc, const double zsrc)
{
    const double slow = 1.0/5000.0;
    double d2, tobs;
    d2 = pow(x-xsrc, 2) + pow(y-ysrc, 2) + pow(z-zsrc, 2);
    tobs = sqrt(d2)*slow;
    return tobs; 
}

int computeAngles(const int nx, const int ny, const int nz,
                  const bool lbackProj,
                  const double *__restrict__ ttimes,
                  double *__restrict__ az,
                  double *__restrict__ aoi)
{
    int i, i1, i2, ic, igrd, ix, iy, iz, j, k, ngrdx, ngrdx3, nxy, y3, z3;
    double *work, *workxm, *workxp, *workym, *workyp, *workzm, *workzp,
           dxh, dxh2, dyh, dyh2, dzh, dzh2,
           fxmh, fxph, fymh, fyph, fzmh, fzph, phi, rho, theta;
    const int chunkx = 16;
    const double half = 0.5; 
    const double pi180i = 180.0/M_PI;
    ngrdx = 16; 
    ngrdx3 = ngrdx*3;
    nxy = nx*ny;
    work = memory_calloc__double((chunkx+1)*9);
    workxp = memory_calloc__double(2*(chunkx+1));
    workyp = memory_calloc__double(2*(chunkx+1));
    workzp = memory_calloc__double(2*(chunkx+1));
    workxm = memory_calloc__double(chunkx+1);
    workym = memory_calloc__double(chunkx+1);
    workzm = memory_calloc__double(chunkx+1);
    for (k=1; k<nz-1; k++)
    {
        for (j=1; j<ny-1; j++)
        {
            for (i=1; i<nx-1; i=i+chunkx)
            {
                // prefetch the values in the grid
                i1 = i-1;
                i2 = MIN(nx-2, i+chunkx);

                // differentiate pre-fetched values
                for (z3=-1; z3<=1; z3++)
                {
                    for (y3=-1; y3<=1; y3++)
                    {
                        iz = k + z3; 
                        iy = j + y3;
                        for (ix=i1; ix<i2; ix++)
                        {
                            igrd = iz*nxy + iy*nx + ix;
                            // central difference
                            ic = (z3 + 1)*ngrdx3
                               + (y3 + 1)*ngrdx
                               + (ix - i1 + 1);
                            ic = ix;
                            fxph = half*(workxp[ic] + workxm[ic]); //half*(work[ic+1] + work[ic]);
                            //fxmh = half*(workx[ic] + workx[ic-4]); //half*(work[ic] + work[ic-1]);
                            fyph = half*(workyp[ic] + workym[ic]); //half*(work[ic+ngrdx] + work[ic]);
                            //fymh = half*(worky[ic] + worky[ic-4]); //half*(work[ic] + work[ic-ngrdx]);
                            fzph = half*(workzp[ic] + workzm[ic]); //half*(work[ic+ngrdx3] + work[ic]);
                            //fzph = half*(workz[ic] + workz[ic-4]); //half*(work[ic] + work[ic-ngrdx3]);
                            dxh = fxph;// - fxmh;
                            dyh = fyph;// - fymh;
                            dzh = fzph;// - fzmh; 
                            dxh2 = dxh*dxh;
                            dyh2 = dyh*dyh;
                            dzh2 = dzh*dzh;
                            rho = sqrt(dxh2 + dyh2 + dzh2);
                            phi = acos(dzh/rho)*pi180i;
                            theta = atan2(dyh, dxh)*pi180i;
                             
                            // the modeling coordinate system is right handed up
                            // while the observation coordinate system is + from
                            // north.  hence, make +x 'north' and -y 'east'.
                            az[igrd] = 90.0 - theta;
                            if (az[igrd] < 0.0){az[igrd] = az[igrd] + 360.0;}
                            aoi[igrd] = theta - 180.0;
                        }
                    }
                } 
            }
        }
    }
    free(work);
    free(workxm);
    free(workym);
    free(workzm);
    free(workxp);
    free(workyp);
    free(workzp);
    return 0;
}

/*!
 * @brief Stacks the weighted residuals into the analytic origin time
 *        computation.  This is for the analytic removal of the source
 *        time in a least squares optimization as descxribed by Moser
 *        et. al. 1992 Eqn 19 for diagonally weighted matrices.
 *
 * @param[in] ngrd    number of grid points in the grid
 * @param[in] tobs_i  i'th observed pick time (seconds)
 * @param[in] xnorm   normalization factor which is the inverse of the
 *                    sum of the data weights for all observations
 *                    (1/seconds)
 * @param[in] wt_i    i'th observed pick time weight (1/seconds)
 * @param[in] test    estimate arrival times (seconds) at all points in the
 *                    grid. [ngrd]
 *
 * @param[in,out] t0  on input contains the currented weighted residual
 *                    sum for previous observations.
 *                    on output contains the updated weighted residual
 *                    sum which incorporates the current observation. [ngrd]
 *
 * @author Ben Baker
 *
 * @copyright Apache 2
 *
 */
static void locate_l2_stackT0__double64(const int ngrd,
                                        const double tobs_i,
                                        const double xnorm, 
                                        const double wt_i,
                                        const double *__restrict__ test,
                                        double *__restrict__ t0)
{
    double wt __attribute__ ((aligned(64))) = 0.0;
    double tobs __attribute__ ((aligned(64))) = 0.0;
    int igrd __attribute__ ((aligned(64))) = 0;
    //------------------------------------------------------------------------//
    wt = wt_i/xnorm;
    tobs = tobs_i;
#ifdef __INTEL_COMPILER
    __assume_aligned(test, 64); 
    __assume_aligned(t0, 64);
#else
    #pragma omp simd aligned(t0, test:64)
#endif
    for (igrd=0; igrd<ngrd; igrd++)
    {
        t0[igrd] = t0[igrd] + wt*(tobs - test[igrd]);
    }
    return;
}
//============================================================================//
/*!
 * @brief Stacks the weighted residuals into the analytic origin time
 *        computation.  This is for the analytic removal of the source
 *        time in a least squares optimization as descxribed by Moser
 *        et. al. 1992 Eqn 19 for diagonally weighted matrices.
 *
 * @param[in] ngrd    number of grid points in the grid
 * @param[in] tobs_i  i'th observed pick time (seconds)
 * @param[in] xnorm   normalization factor which is the inverse of the
 *                    sum of the data weights for all observations
 *                    (1/seconds)
 * @param[in] wt_i    i'th observed pick time weight (1/seconds)
 * @param[in] test    estimate arrival times (seconds) at all points in the
 *                    grid. [ngrd]
 *
 * @param[in,out] t0  on input contains the currented weighted residual
 *                    sum for previous observations.
 *                    on output contains the updated weighted residual
 *                    sum which incorporates the current observation. [ngrd]
 *
 * @author Ben Baker
 *
 * @copyright Apache 2
 *
 */
static void locate_l2_stackT0__float64(const int ngrd,
                                       const float tobs_i,
                                       const float xnorm, 
                                       const float wt_i,
                                       const float *__restrict__ test,
                                       float *__restrict__ t0) 
{
    float wt __attribute__ ((aligned(64))) = 0.0f;
    float tobs __attribute__ ((aligned(64))) = 0.0f;
    int igrd __attribute__ ((aligned(64))) = 0;
    //------------------------------------------------------------------------//
    wt = wt_i/xnorm;
    tobs = tobs_i;
#ifdef __INTEL_COMPILER
    __assume_aligned(test, 64); 
    __assume_aligned(t0, 64);
#else
    #pragma omp simd aligned(t0, test:64)
#endif
    for (igrd=0; igrd<ngrd; igrd++)
    {   
        t0[igrd] = t0[igrd] + wt*(tobs - test[igrd]);
    }   
    return;
}
//============================================================================//
/*!
 * @brief Stacks the squared residuals into the least squares penalty
 *        function.  This assumes diagonal weighting.
 *
 * @param[in] ngrd        number of poitns in grid search
 * @param[in] tobs_i      i'th observed pick time (seconds)
 * @param[in] wt_i        i'th observed pick time weight (1/seconds)
 * @param[in] test        estimate arrival times (seconds) at all points in the
 *                        grid. [ngrd]
 * @param[in] t0          origin time (seconds) at all points in grid [ngrd]
 *
 * @param[in,out] objfn   on input contains the weighted squared residuals
 *                        at all poitns in the grid.
 *                        on output contains the contribution of this
 *                        observation to all the squared residuals at all
 *                        points in the grid. [ngrd] 
 *
 * @author Ben Baker
 *
 * @copyright Apache 2
 *
 */
static void locate_l2_stackObjfn__double64(const int ngrd,
                                           const double tobs_i,
                                           const double wt_i,
                                           const double *__restrict__ test,
                                           const double *__restrict__ t0,
                                           double *__restrict__ objfn)
{
    double wt __attribute__ ((aligned(64))) = 0.0;
    double res __attribute__ ((aligned(64))) = 0.0;
    double tobs __attribute__ ((aligned(64))) = 0.0; 
    int igrd __attribute__ ((aligned(64))) = 0; 
    const double sqrt2i = 0.7071067811865475; //one/sqrt(two);
    //------------------------------------------------------------------------//
    wt = wt_i*sqrt2i;
    tobs = tobs_i;
#ifdef __INTEL_COMPILER
    __assume_aligned(test, 64);
    __assume_aligned(t0, 64);
    __assume_aligned(objfn, 64);
#else
    #pragma omp simd aligned(t0, test, objfn: 64)
#endif
    for (igrd=0; igrd<ngrd; igrd++)
    {
        res = wt*(tobs - (test[igrd] + t0[igrd]));
        objfn[igrd] = objfn[igrd] + res*res;
    }
    return;
}
//============================================================================//
/*!
 * @brief Stacks the squared residuals into the least squares penalty
 *        function.  This assumes diagonal weighting.
 *
 * @param[in] ngrd        number of poitns in grid search
 * @param[in] tobs_i      i'th observed pick time (seconds)
 * @param[in] wt_i        i'th observed pick time weight (1/seconds)
 * @param[in] test        estimate arrival times (seconds) at all points in the
 *                        grid. [ngrd]
 * @param[in] t0          origin time (seconds) at all points in grid [ngrd]
 *
 * @param[in,out] objfn   on input contains the weighted squared residuals
 *                        at all poitns in the grid.
 *                        on output contains the contribution of this
 *                        observation to all the squared residuals at all
 *                        points in the grid. [ngrd] 
 *
 * @author Ben Baker
 *
 * @copyright Apache 2
 *
 */
static void locate_l2_stackObjfn__float64(const int ngrd,
                                          const float tobs_i,
                                          const float wt_i,
                                          const float *__restrict__ test,
                                          const float *__restrict__ t0, 
                                          float *__restrict__ objfn)
{
    float wt __attribute__ ((aligned(64))) = 0.0f;
    float res __attribute__ ((aligned(64))) = 0.0f;
    float tobs __attribute__ ((aligned(64))) = 0.0f;
    int igrd __attribute__ ((aligned(64))) = 0;  
    const float sqrt2i = 0.7071067811865475f; //one/sqrt(two);
    //------------------------------------------------------------------------//
    wt = wt_i*sqrt2i;
    tobs = tobs_i;
#ifdef __INTEL_COMPILER
    __assume_aligned(test, 64);
    __assume_aligned(t0, 64);
    __assume_aligned(objfn, 64);
#else
    #pragma omp simd aligned(t0, test, objfn:64)
#endif
    for (igrd=0; igrd<ngrd; igrd++)
    {
        res = wt*(tobs - (test[igrd] + t0[igrd]));
        objfn[igrd] = objfn[igrd] + res*res;
    }   
    return;
}
//============================================================================//
double *memory_calloc__double(const int n)
{
    double *x = NULL;
    int i, ierr;
    const double zero __attribute__ ((aligned(64))) = 0.0; 
    ierr = memory_malloc__double(n, &x);
#ifdef __INTEL_COMPILER
    __assume_aligned(x, 64);
#endif
    //memset(x, 0.0, (size_t) n*sizeof(double));
    for (i=0; i<n; i++)
    {
        x[i] = zero;
    }
    return x; 
}

float *memory_calloc__float(const int n)
{
    float *x = NULL;
    int i, ierr;
    const float zero __attribute__ ((aligned(64))) = 0.0f;
    ierr = memory_malloc__float(n, &x);
#ifdef __INTEL_COMPILER
    __assume_aligned(x, 64); 
#endif
    //memset(x, 0.0f, (size_t) n*sizeof(float));
    for (i=0; i<n; i++)
    {
        x[i] = zero;
    }
    return x;
}

int *memory_calloc__int(const int n)
{
    int *x = NULL;
    int i, ierr;
    const int zero __attribute__ ((aligned(64))) = 0;
    ierr = memory_malloc__int(n, &x);
#ifdef __INTEL_COMPILER 
    __assume_aligned(x, 64);
#endif
    //memset(x, 0.0f, (size_t) n*sizeof(float));
    for (i=0; i<n; i++)
    {   
        x[i] = zero;
    }
    return x;
}

int memory_malloc__double(const int n, double **x)
{
    const char *fcnm = "memory_malloc__double\0";
    size_t nbytes;
    int ierr;
    ierr = 0;
    nbytes = (size_t) (n + 64 - n%64)*sizeof(double);
    *x = (double *) aligned_alloc(64, nbytes);
    if (*x == NULL)
    {   
        printf("%s: Error allocating array\n", fcnm);
        ierr = 1;
    }   
    return ierr;
}

int memory_malloc__float(const int n, float **x)
{
    const char *fcnm = "memory_malloc__float\0";
    size_t nbytes;
    int ierr;
    ierr = 0;
    nbytes = (size_t) (n + 64 - n%64)*sizeof(float);
    *x = (float *) aligned_alloc(64, nbytes);
    if (*x == NULL)
    {   
        printf("%s: Error allocating array\n", fcnm);
        ierr = 1;
    }   
    return ierr;
}

int memory_malloc__int(const int n, int **x)
{
    const char *fcnm = "memory_malloc__int\0";
    size_t nbytes;
    int ierr;
    ierr = 0;
    nbytes = (size_t) (n + 64 - n%64)*sizeof(int);
    *x = (int *) aligned_alloc(64, nbytes);
    if (*x == NULL)
    {
        printf("%s: Error allocating array\n", fcnm);
        ierr = 1;
    }
    return ierr;
}
//============================================================================//
/*!
 * @brief Sets all elements of a 64 bit aligned array x to x0in
 *
 * @param[in] n      number of points in array x
 * @param[in] x0in   value to set
 *
 * @param[out] x     all values of array set to x0 [n]
 *
 * @author Ben Baker
 *
 * @copyright Apache 2
 *
 */
void locate_setDouble64(const int n, const double x0in,
                        double *__restrict__ x)
{
    double x0 __attribute__ ((aligned(64))) = x0in;
    int i __attribute__ ((aligned(64))) = 0;
#ifdef __INTEL_COMPILER
    __assume_aligned(x, 64);
#else
    #pragma omp simd aligned(x : 64)
#endif
    for (i=0; i<n; i++)
    {
        x[i] = x0;
    }
}
//============================================================================//
/*!
 * @brief Sets all elements of a 64 bit aligned array x to x0in
 *
 * @param[in] n      number of points in array x
 * @param[in] x0in   value to set
 *
 * @param[out] x     all values of array set to x0 [n]
 *
 * @author Ben Baker
 *
 * @copyright Apache 2
 *
 */
void locate_setFloat64(const int n, const float x0in,
                       float *__restrict__ x)
{
    float x0 __attribute__ ((aligned(64))) = x0in;
    int i __attribute__ ((aligned(64))) = 0;
#ifdef __INTEL_COMPILER
    __assume_aligned(x, 64);
#else
    #pragma omp simd aligned(x : 64)
#endif
    for (i=0; i<n; i++)
    {
        x[i] = x0;
    }
}
//============================================================================//
/*!
 * @brief Zeros out a 64 bit aligned array x
 *
 * @param[in] n    number of points in arra yx
 *
 * @param[out] x   nulled out array [n]
 *
 * @author Ben Baker
 *
 * @copyright Apache 2
 *
 */
static void locate_nullDouble64(const int n, double *__restrict__ x)
{
    double zero __attribute__ ((aligned(64))) = 0.0;
    int i __attribute__ ((aligned(64))) = 0;
#ifdef __INTEL_COMPILER 
    __assume_aligned(x, 64);
#else
    #pragma omp simd aligned(x : 64)
#endif
    for (i=0; i<n; i++)
    {   
        x[i] = zero;
    }
    return;
}
//============================================================================//
/*!
 * @brief Zeros out a 64 bit aligned array x
 *
 * @param[in] n    number of points in arra yx
 *
 * @param[out] x   nulled out array [n]
 *
 * @author Ben Baker
 *
 * @copyright Apache 2
 *
 */
void locate_nullFloat64(const int n, float *__restrict__ x)
{
    float zero __attribute__ ((aligned(64))) = 0.0f;
    int i __attribute__ ((aligned(64))) = 0;
#ifdef __INTEL_COMPILER
    __assume_aligned(x, 64);
#else
    #pragma omp simd aligned(x : 64)
#endif
    for (i=0; i<n; i++)
    {   
        x[i] = zero; 
    }
    return;
}
//============================================================================//
/*!
 * @brief Sums all elements in an array
 *
 * @param[in] n    number of points in array to sum
 * @param[in] x    array to sum [n]
 *
 * @result sum of all elements in array x
 *
 * @author Ben Baker
 *
 * @copyright Apache 2
 *
 */
double locate_sumDouble64(const int n, const double *__restrict__ x)
{
    double xsum __attribute__ ((aligned(64))) = 0.0;
    int i __attribute__ ((aligned(64))) = 0;
#ifdef __INTEL_COMPILER
    __assume_aligned(x, 64);
#else
    #pragma omp simd aligned(x : 64) reduction(+:xsum)
#endif
    for (i=0; i<n; i++)
    {
        xsum = xsum + x[i];
    }
    return xsum;
}
//============================================================================//
int locate_minLocDouble64(const int n, const double *__restrict__ x)
{
    double xmin __attribute__ ((aligned(64))) = 0.0;
    int imin __attribute__ ((aligned(64))) = 0;
    int i __attribute__ ((aligned(64))) = 0;
    xmin = x[0];
    imin = 0;
#ifdef __INTEL_COMPILER
    __assume_aligned(x, 64);
#endif
    for (i=1; i<n; i++)
    {
        if (x[i] < xmin)
        {
            imin = i;
            xmin = x[i];
        }
    }
    return imin;
}
//============================================================================//
int locate_minLocFloat64(const int n, const float *__restrict__ x)
{
    float xmin __attribute__ ((aligned(64))) = 0.0f;
    int imin __attribute__ ((aligned(64))) = 0;
    int i __attribute__ ((aligned(64))) = 0;
    xmin = x[0];
    imin = 0;
#ifdef __INTEL_COMPILER
    __assume_aligned(x, 64);
#endif
    for (i=1; i<n; i++)
    {   
        if (x[i] < xmin)
        {
            imin = i;
            xmin = x[i];
        }
    }   
    return imin;
}
//============================================================================//
/*!
 * @brief Sums all elements in an array
 *
 * @param[in] n    number of points in array to sum
 * @param[in] x    array to sum [n]
 *
 * @result sum of all elements in array x
 *
 * @author Ben Baker
 *
 * @copyright Apache 2
 *
 */
float locate_sumFloat64(const int n, const float *__restrict__ x)
{
    float xsum __attribute__ ((aligned(64))) = 0.0f;
    int i __attribute__ ((aligned(64))) = 0;
#ifdef __INTEL_COMPILER
    __assume_aligned(x, 64);
#else
    #pragma omp simd aligned(x : 64) reduction(+:xsum)
#endif
    for (i=0; i<n; i++)
    {   
        xsum = xsum + x[i];
    }   
    return xsum;
}
//============================================================================//
int locate_sumInt(const int n, const int *__restrict__ x)
{
    int i, xsum;
    xsum = 0;
    for (i=0; i<n; i++)
    {
        xsum = xsum + x[i];
    }
    return xsum;
}
//============================================================================//
/*!
 * @brief Performs the least squares gridsearch.  If the origin time is 
 *        desired then the locations correspond to an analytic integration
 *        of the origin time (e.g. Moser et al. 1992).
 *
 * @param[in] ldgrd    leading dimension of traveltime tables (>= ngrd)
 * @param[in] ngrd     number of grid points in traveltime tables
 * @param[in] nobs     number of observations
 * @param[in] iwantOT  if 1 then compute the origint time.
 *                     otherwise the origin time will be set by t0use.
 * @param[in] t0use    if iwantOT is not when this this is the source
 *                     origin time (seconds)
 * @param[in] mask     if the i'th observation is 1 then it is masked from
 *                     from the location [nobs]
 * @param[in] tobs     observed pick times (seconds) [nobs]
 * @param[in] tcorr    static corrections (seconds) for stations [nobs].
 *                     if NULL then it is ignored and assumed all zero.
 * @param[in] varobs   variance in pick times (seconds) [nobs]
 * @param[in] test     traveltime tables for each observation [ldgrd x nobs]
 *                     with leading dimension ldgrd.
 *
 * @param[out] t0      origin time (seconds) at each point in grid serach [ngrd]
 * @param[out] objfn   residual squared objective function at each point in
 *                     the grid [ngrd]. 
 *
 * @author Ben Baker
 *
 * @copyright Apache 2
 *
 */
int locate_l2_gridSearch__double64(const int ldgrd,
                                   const int ngrd,
                                   const int nobs,
                                   const int iwantOT,
                                   const double t0use,
                                   const int *__restrict__ mask,
                                   const double *__restrict__ tobs,
                                   const double *__restrict__ tcorr,
                                   const double *__restrict__ varobs,
                                   const double *__restrict__ test,
                                   double *__restrict__ t0,
                                   double *__restrict__ objfn)
{
    const char *fcnm = "locate_l2_gridSearch__double64\0";
    double *tobsCor, *wtUse;
    double xnorm __attribute__ ((aligned(64))) = 0.0; 
    double tobs_i __attribute__ ((aligned(64))) = 0.0;
    double wt_i __attribute__ ((aligned(64))) = 0.0;
    const double zero = 0.0;
    const double one = 1.0;
    int *obsPtr, ibeg, ierr, iobs, jobs, nobsUse;
    //------------------------------------------------------------------------//
    //
    // Error checking
    ierr = 0;
    if ((sizeof(double)*(size_t) ldgrd)%64 != 0 || ldgrd < ngrd || nobs < 1 ||
        mask == NULL || tobs == NULL || varobs == NULL ||
        test == NULL || t0 == NULL || objfn == NULL)
    {
        if ((sizeof(double)*(size_t) ldgrd)%64 != 0)
        {
            printf("%s: Error ldgrd must be divisible by 64\n", fcnm);
        }
        if (ldgrd < ngrd){printf("%s: Error ldgrd < ngrd\n", fcnm);}
        if (mask == NULL){printf("%s: mask is null\n", fcnm);}
        if (tobs == NULL){printf("%s: tobs is null\n", fcnm);}
        if (varobs == NULL){printf("%s: varobs is null\n", fcnm);}
        if (test == NULL){printf("%s: test is null\n", fcnm);}
        if (t0 == NULL){printf("%s: t0 is null\n", fcnm);}
        if (objfn == NULL){printf("%s: objfn is null\n", fcnm);}
        ierr = 1;
        return ierr;
    }
    // Require the arrays be 64 bit aligned
    if (memory_isAligned(t0, 64) != 1 || 
        memory_isAligned(test, 64) != 1 ||
        memory_isAligned(objfn, 64) != 1)
    {
        printf("%s: Input arrays are not 64 bit aligned\n", fcnm);
        ierr = 1;
        return ierr;
    }
    // zero out the result
    locate_nullDouble64(ngrd, objfn);
    // Set the static corrections.  While it would seem more sensible
    // to add the correction to the estimate recall that we are ultimately
    // interested in residuals so we instead remove it from the observation
    // because t_obs - (t_est + t_stat) = t_obs - t_stat - t_est = t_cor - t_et
    obsPtr = memory_calloc__int(nobs);
    tobsCor = memory_calloc__double(nobs);
    wtUse = memory_calloc__double(nobs);
    nobsUse = 0;
    xnorm = zero;
    if (tcorr == NULL)
    {
        for (iobs=0; iobs<nobs; iobs++)
        {
            if (mask[iobs] == 0)
            {
                tobsCor[nobsUse] = tobs[iobs];
                obsPtr[nobsUse] = iobs;
                wtUse[nobsUse] = one/varobs[iobs];
                xnorm = xnorm + wtUse[nobsUse]; //one/varobs[iobs];
                nobsUse = nobsUse + 1;
            }
        }
    }
    else
    {
        for (iobs=0; iobs<nobs; iobs++)
        {
            if (mask[iobs] == 0)
            {
                tobsCor[nobsUse] = tobs[iobs] - tcorr[iobs];
                obsPtr[nobsUse] = iobs;
                wtUse[nobsUse] = one/varobs[iobs];
                xnorm = xnorm + wtUse[nobsUse]; //one/varobs[iobs];
                nobsUse = nobsUse + 1;
            }
        }
    }
    // Compute the least-squares origin time which is the average reisdual
    if (iwantOT == 1)
    {
        locate_nullDouble64(ngrd, t0);
        for (jobs=0; jobs<nobsUse; jobs++)
        {
            iobs = obsPtr[jobs];
            tobs_i = tobsCor[jobs];
            wt_i = wtUse[jobs]; //one/varUse[jobs];
            ibeg = ldgrd*iobs;
            locate_l2_stackT0__double64(ngrd, tobs_i, xnorm,  wt_i,
                                        &test[ibeg], t0);
        }
    }
    // Set the desired residual
    else
    {
        locate_setDouble64(ngrd, t0use, t0);
    }
    // Compute the locations with the origin times at each grid point
    for (jobs=0; jobs<nobsUse; jobs++)
    {
        iobs = obsPtr[jobs];
        tobs_i = tobsCor[jobs];
        wt_i = wtUse[jobs]; //one/varUse[jobs];
        ibeg = ldgrd*iobs;
        locate_l2_stackObjfn__double64(ngrd, tobs_i, wt_i,
                                       &test[ibeg], t0, objfn);
    }
    free(obsPtr);
    free(tobsCor);
    free(wtUse);
    return ierr;     
}
//============================================================================//
/*!
 * @brief Performs the least squares gridsearch.  If the origin time is 
 *        desired then the locations correspond to an analytic integration
 *        of the origin time (e.g. Moser et al. 1992).
 *
 * @param[in] ldgrd    leading dimension of traveltime tables (>= ngrd)
 * @param[in] ngrd     number of grid points in traveltime tables
 * @param[in] nobs     number of observations
 * @param[in] iwantOT  if 1 then compute the origint time.
 *                     otherwise the origin time will be set by t0use.
 * @param[in] t0use    if iwantOT is not when this this is the source
 *                     origin time (seconds)
 * @param[in] mask     if the i'th observation is 1 then it is masked from
 *                     from the location [nobs]
 * @param[in] tobs     observed pick times (seconds) [nobs]
 * @param[in] tcorr    static corrections (seconds) for stations [nobs].
 *                     if NULL then it is ignored and assumed all zero.
 * @param[in] varobs   variance in pick times (seconds) [nobs]
 * @param[in] test     traveltime tables for each observation [ldgrd x nobs]
 *                     with leading dimension ldgrd.
 *
 * @param[out] t0      origin time (seconds) at each point in grid serach [ngrd]
 * @param[out] objfn   residual squared objective function at each point in
 *                     the grid [ngrd]. 
 *
 * @author Ben Baker
 *
 * @copyright Apache 2
 *
 */
int locate_l2_gridSearch__float64(const int ldgrd,
                                  const int ngrd,
                                  const int nobs,
                                  const int iwantOT,
                                  const float t0use,
                                  const int *__restrict__ mask,
                                  const float *__restrict__ tobs,
                                  const float *__restrict__ tcorr, 
                                  const float *__restrict__ varobs,
                                  const float *__restrict__ test,
                                  float *__restrict__ t0,
                                  float *__restrict__ objfn)
{
    const char *fcnm = "locate_l2_gridSearch__float64\0";
    float *tobsCor, *wtUse;
    float xnorm __attribute__ ((aligned(64))) = 0.0f;
    float tobs_i __attribute__ ((aligned(64))) = 0.0f;
    float wt_i __attribute__ ((aligned(64))) = 0.0f;
    const float zero = 0.0f;
    const float one = 1.0f;
    int ibeg, ierr, iobs, jobs, nobsUse, *obsPtr;
    //------------------------------------------------------------------------//
    //   
    // Error checking
    ierr = 0;
    if ((sizeof(float)*(size_t) ldgrd)%64 != 0 || ldgrd < ngrd || nobs < 1 ||
        mask == NULL || tobs == NULL || varobs == NULL ||
        test == NULL || t0 == NULL || objfn == NULL)
    {
        if ((sizeof(float)*(size_t) ldgrd)%64 != 0)
        {
            printf("%s: Error ldgrd must be divisible by 64\n", fcnm);
        }
        if (ldgrd < ngrd){printf("%s: Error ldgrd < ngrd\n", fcnm);}
        if (mask == NULL){printf("%s: mask is null\n", fcnm);}
        if (tobs == NULL){printf("%s: tobs is null\n", fcnm);}
        if (varobs == NULL){printf("%s: varobs is null\n", fcnm);}
        if (test == NULL){printf("%s: test is null\n", fcnm);}
        if (t0 == NULL){printf("%s: t0 is null\n", fcnm);}
        if (objfn == NULL){printf("%s: objfn is null\n", fcnm);}
        ierr = 1;
        return ierr;
    }
    // Require the arrays be 64 bit aligned
    if (memory_isAligned(t0, 64) != 1 ||
        memory_isAligned(test, 64) != 1 ||
        memory_isAligned(objfn, 64) != 1)
    {
        printf("%s: Input arrays are not 64 bit aligned\n", fcnm);
        ierr = 1;
        return ierr;
    }
    // zero out the result
    locate_nullFloat64(ngrd, objfn);
    // Set the static corrections.  While it would seem more sensible
    // to add the correction to the estimate recall that we are ultimately
    // interested in residuals so we instead remove it from the observation
    // because t_obs - (t_est + t_stat) = t_obs - t_stat - t_est = t_cor - t_est
    obsPtr = memory_calloc__int(nobs);
    tobsCor = memory_calloc__float(nobs);
    wtUse = memory_calloc__float(nobs);
    nobsUse = 0;
    xnorm = zero;
    if (tcorr == NULL)
    {
        for (iobs=0; iobs<nobs; iobs++)
        {
            if (mask[iobs] == 0)
            {
                tobsCor[nobsUse] = tobs[iobs];
                obsPtr[nobsUse] = iobs;
                wtUse[nobsUse] = one/varobs[iobs];
                xnorm = xnorm + wtUse[nobsUse]; //one/varobs[iobs];
                nobsUse = nobsUse + 1;
            }
        }
    }
    else
    {
        for (iobs=0; iobs<nobs; iobs++)
        {
            if (mask[iobs] == 0)
            {
                tobsCor[nobsUse] = tobs[iobs] - tcorr[iobs];
                obsPtr[nobsUse] = iobs;
                wtUse[nobsUse] = one/varobs[iobs];
                xnorm = xnorm + wtUse[nobsUse]; //one/varobs[iobs];
                nobsUse = nobsUse + 1;
            }
        }
    }
    // Compute the least-squares origin time which is the average reisdual
    if (iwantOT == 1)
    {
        locate_nullFloat64(ngrd, t0);
        for (jobs=0; jobs<nobsUse; jobs++)
        {
            iobs = obsPtr[jobs];
            tobs_i = tobsCor[jobs];
            wt_i = wtUse[jobs]; //one/varUse[jobs];
            ibeg = ldgrd*iobs;
            locate_l2_stackT0__float64(ngrd, tobs_i, xnorm,  wt_i,
                                       &test[ibeg], t0);
        }
    }
    // Set the desired residual
    else
    {
        locate_setFloat64(ngrd, t0use, t0);
    }
    // Compute the locations with the origin times at each grid point
    for (jobs=0; jobs<nobsUse; jobs++)
    {
        iobs = obsPtr[jobs];
        tobs_i = tobsCor[jobs];
        wt_i = wtUse[jobs]; //one/varUse[jobs];
        ibeg = ldgrd*iobs;
        locate_l2_stackObjfn__float64(ngrd, tobs_i, wt_i,
                                      &test[ibeg], t0, objfn);
    }
    free(obsPtr);
    free(tobsCor);
    free(wtUse);
    return ierr;
}
//============================================================================//
int locate_l1_gridSearch__double64(const int ldgrd,
                                   const int ngrd,
                                   const int nobs,
                                   const int iwantOT,
                                   const double t0use,
                                   const int *__restrict__ mask,
                                   const double *__restrict__ tobs,
                                   const double *__restrict__ varobs,
                                   const double *__restrict__ test,
                                   double *__restrict__ t0,
                                   double *__restrict__ objfn)
{
    double *testPerm, *tobsUse, *wt, *wtSort, *res;//, *resSort;
    double rsum __attribute__ ((aligned(64))) = 0.0; 
    double t0opt __attribute__ ((aligned(64))) = 0.0;
    double tobs_i __attribute__ ((aligned(64))) = 0.0;
    double wt_i __attribute__ ((aligned(64))) = 0.0;
    const double one __attribute__ ((aligned(64))) = 1.0;
    const double zero __attribute__ ((aligned(64))) = 0.0;
    double wtsum, wtsumi;
    int *obsPtr, *perm, ierr, igrd, igrd1, igrd2, iobs, indx, jgrd,
        jndx, jobs, ldobs, nobsUse, nsort;
    bool lsort;
    const int nchunk = 256;

    ierr = 0;
    perm = memory_calloc__int(nobs);
    obsPtr = memory_calloc__int(nobs);
    tobsUse = memory_calloc__double(nobs);
    wt = memory_calloc__double(nobs);
    wtsum = 0.0;
    nobsUse = 0;
    for (iobs=0; iobs<nobs; iobs++)
    {
        if (mask[iobs] == 0)
        {
            wt[nobsUse] = one/varobs[iobs];
            tobsUse[nobsUse] = tobs[iobs];
            obsPtr[iobs] = iobs;
            nobsUse = nobsUse + 1;
        }
        wtsum = wtsum + wt[iobs]; 
    }
    ldobs = nobsUse + 8 - nobsUse%8;
    locate_nullDouble64(ngrd, objfn);
    if (iwantOT == 1)
    {
        lsort = false;
        nsort = 0;
        // Zero out the origin time
        locate_nullDouble64(ngrd, t0);
        // Set space 
        res = memory_calloc__double(nobsUse);
        //resSort = memory_calloc__double(nobsUse);
        wtSort = memory_calloc__double(nobsUse);
        // Normalize the weights?
        wtsumi = one/wtsum;
        if (fabs(wtsum - one) > 1.e-14)
        {
            for (iobs=0; iobs<nobsUse; iobs++)
            {
                wt[iobs] = wt[iobs]*wtsumi;
            }    
        }
        for (iobs=0; iobs<nobsUse; iobs++){perm[iobs] = iobs;}
/*
        #pragma omp parallel for \
         firstprivate (obsPtr, perm, res, tobsUse, wt, wtSort) \
         private (ierr, igrd, igrd1, igrd2, iobs, jgrd, jndx, jobs, lsort, rsum, t0opt) \
         shared (nchunk, ngrd, nobsUse, objfn, t0, test, zero) \
         default (none) reduction (+:nsort)
*/
        // Loop on chunks
        for (jgrd=0; jgrd<ngrd; jgrd=jgrd+nchunk)
        {
            igrd1 = jgrd; 
            igrd2 = MIN(ngrd, jgrd+nchunk);
            // Loop on the subgrid 
            for (igrd=igrd1; igrd<igrd2; igrd++)
            {
                // Extract the estimates
                for (iobs=0; iobs<nobsUse; iobs++)
                {
                    jobs = obsPtr[iobs];
                    jndx = jobs*ldgrd + igrd;
                    res[iobs] = tobsUse[iobs] - test[jndx];
                    wtSort[iobs] = wt[iobs];
                }
                t0opt = weightedMedian__double(nobsUse, res,
                                               wtSort, perm, &lsort, &ierr);
                if (lsort){nsort = nsort + 1;}
                t0[igrd] = t0opt;
            }
        }
        free(res);
printf("%d %d\n", nsort, ngrd);
    }
    // Simply set the origin time to t0
    else
    {
        locate_setDouble64(ngrd, t0use, t0);
    }
    // Compute the L1 misfit function at all points - loop on chunks
    #pragma omp parallel for \
     private (igrd, igrd1, igrd2, iobs, jgrd, jobs, jndx, tobs_i, wt_i) \
     shared (ldgrd, nchunk, nobsUse, objfn, obsPtr, test, t0, tobsUse, wt) \
     default (none)
    for (jgrd=0; jgrd<ngrd; jgrd=jgrd+nchunk)
    {
        igrd1 = jgrd;
        igrd2 = MIN(ngrd, jgrd+nchunk);
        // Compute the L1 norm for all observations
        for (iobs=0; iobs<nobsUse; iobs++)
        {
            jobs = obsPtr[iobs];
            tobs_i = tobsUse[iobs]; 
            wt_i = wt[iobs];
            // Extract the estimate at each grid point in chunk
            for (igrd=igrd1; igrd<igrd2; igrd++)
            {
                jndx = jobs*ldgrd + igrd;
                objfn[igrd] = objfn[igrd]
                            + wt_i*fabs(tobs_i - test[jndx] - t0[igrd]);
            }
        }
    }
    //free(testPerm);
    free(wt);
    free(perm);
    return ierr;
}

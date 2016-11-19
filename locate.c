#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>

#define memory_isAligned(POINTER, BYTE_COUNT) \
    (((uintptr_t)(const void *)(POINTER)) % (BYTE_COUNT) == 0)

int locate_l2_gridSearch__double64(const int ldgrd,
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
int locate_l2_gridSearch__float64(const int ldgrd,
                                  const int ngrd,
                                  const int nobs,
                                  const int iwantOT,
                                  const float t0use,
                                  const int *__restrict__ mask,
                                  const float *__restrict__ tobs,
                                  const float *__restrict__ varobs,
                                  const float *__restrict__ test,
                                  float *__restrict__ t0, 
                                  float *__restrict__ objfn);
int locate_l1_gridSearch__double64(const int ldobs,
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
double locate_sumDouble64(const int n, double *__restrict__ x);
float locate_sumFloat64(const int n, float *__restrict__ x);
int locate_minLocDouble64(const int n, double *__restrict__ x);
int locate_minLocFloat64(const int n, float *__restrict__ x);
int memory_malloc__double(const int n, double **x);
int memory_malloc__float(const int n, float **x);
int memory_calloc__double(const int n, double **x);
int memory_calloc__float(const int n, float **x);
void makeTest(const int nx, const int ny, const int nz,
              const double dx, const double dy, const double dz,
              const double xsrc, const double ysrc, const double zsrc,
              double *__restrict__ test);
double makeObs(const double x, const double y, const double z,
               const double xsrc, const double ysrc, const double zsrc);

int main()
{
    double *objfn, *xrec, *yrec, *zrec, *t0, *test, *tobs, *varobs;
    float *objfn4, *t04, *test4, *tobs4, *varobs4;
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
    nz = 55;
    ngrd = nx*ny*nz;
    ldgrd = ngrd + 64 - ngrd%64;
    ixsrc = 12;
    iysrc = 15;
    izsrc = 5;
printf("true: %d\n", izsrc*nx*ny + iysrc*nx + ixsrc);
    xsrc = (double) ixsrc*dx;
    ysrc = (double) iysrc*dy;
    zsrc = (double) izsrc*dz;
    memory_calloc__double(nobs, &xrec);
    memory_calloc__double(nobs, &yrec);
    memory_calloc__double(nobs, &zrec);

    memory_calloc__double(nobs, &tobs);
    memory_calloc__double(nobs, &varobs);
    memory_calloc__double(nobs*ldgrd, &test);
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
    memory_calloc__double(ngrd, &t0);
    memory_calloc__double(ngrd, &objfn);
    locate_l2_gridSearch__double64(ldgrd, ngrd, nobs, iwantOT, 
                                   t0use, mask, tobs, varobs,
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
    memory_calloc__float(nobs, &tobs4);
    memory_calloc__float(nobs, &varobs4);
    memory_calloc__float(ldgrd*nobs, &test4);
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
    free(varobs);
    memory_calloc__float(ngrd, &objfn4);
    memory_calloc__float(ngrd, &t04);
    locate_l2_gridSearch__float64(ldgrd, ngrd, nobs, iwantOT,
                                  t0use4, mask, tobs4, varobs4,  
                                  test4, t04, objfn4);
    iopt = locate_minLocFloat64(ngrd, objfn4);
    printf("float estimate: %d\n", iopt);

    free(test4);
    free(objfn4);
    free(t04);
    free(tobs4);
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

struct double2d_struct
{
    double val;
    int indx;
    char pad[4];
};

struct float2d_struct
{
    float val;
    int indx;
};

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wcast-qual"
#endif
static int cmp_double_array(const void *x, const void *y)
{
    struct double2d_struct xx = *(struct double2d_struct *) x;
    struct double2d_struct yy = *(struct double2d_struct *) y;  
    if (xx.val < yy.val) return -1; 
    if (xx.val > yy.val) return  1;  
    return 0;
} 

static int cmp_float_array(const void *x, const void *y) 
{
    struct float2d_struct xx = *(struct float2d_struct *) x;
    struct float2d_struct yy = *(struct float2d_struct *) y;  
    if (xx.val < yy.val) return -1; 
    if (xx.val > yy.val) return  1;  
    return 0;
}
#ifdef __clang__
#pragma clang diagnostic pop
#endif

/*!
 * @brief Weighted median algorithm for L1 approximation - Gurwitz 1990.
 *        More implementation details can be found at: 
 *        https://www.mathworks.com/matlabcentral/fileexchange/23077-weighted-median
 *
 * @param[in] n       number of elements (must be less than 256)
 * @param[in] x       array of values [n]
 * @param[in] w       weights [n].
 *
 * @param[out] perm   permutation which sorts x in ascending order [n]
 * @param[out] lsort  if true then this function had to perform a sort.
 * @param[out] ierr   0 indicates success
 *
 * @result the weighted median of a list
 *
 * @copyright Apache 2
 *
 */
double weightedMedian__double(const int n,
                              const double *__restrict__ x,
                              const double *__restrict__ w,
                              int *__restrict__ perm,
                              bool *lsort, int *ierr)
{
    const char *fcnm = "weightedMedian__double\0";
    struct double2d_struct vals[256] __attribute__ ((aligned(64)));
    double wt[256] __attribute__ ((aligned(64)));
    double xnorm, wt0, wtheap, wtsum1, wtsum3;
    int i, j, n2;
    bool lmed, lnorm;
    *ierr = 0;
    wtheap = 0.0;
    if (n > 256)
    {
        printf("%s: insufficient space\n", fcnm);
        *ierr = 1;
        return wtheap;
    }
    if (n == 1){return x[0];}
    lmed = true;
    *lsort = false;
    wtsum3 = 0.0;
    wt0 = w[0];
    // put x_i in heap
    for (i=0; i<n; i++)
    {
        wt[i] = w[i];
        wtsum3 = wtsum3 + wt[i];
        if (fabs(wt[i] - wt0) > 1.e-14){lmed = false;}
        if (i > 0)
        {
            if (x[i] < x[i-1]){*lsort = true;}
        }
        vals[i].val = x[i];
        vals[i].indx = i;
    }
    // avoid a division by zero
    if (wtsum3 == 0.0)
    {
        printf("%s: Sum of weights is zero!\n", fcnm);
        *ierr = 1;
        return wtheap;
    }
    // do i have to normalize?
    xnorm = 1.0;
    lnorm = false;
    if (fabs(wtsum3 - 1.0) > 1.e-14)
    {
        xnorm = 1.0/wtsum3;
        lnorm = true;
    }
    if (*lsort)
    {
        // Sort x in ascending order
        qsort((void *) vals, (size_t) n,
              sizeof(struct double2d_struct), cmp_double_array);
        // Get the permutation
        for (i=0; i<n; i++)
        {
            perm[i] = vals[i].indx;
        }
    }
    else
    {
        for (i=0; i<n; i++)
        {
            perm[i] = i;
        }
    }
    // Might be a straight median
    if (lmed)
    {
        n2 = n/2;
        if (n%2 == 0)
        {
            wtheap = 0.5*(x[perm[n2-1]] + x[perm[n2]]); 
        }
        else
        {
            wtheap = x[perm[n2]];
        }
    }
    // It's a weighted median
    else
    {
        i = 0;
        wtsum1 = 0.0;
        // Handle the normalization on the fly
        if (lnorm)
        {
            while (wtsum1 < 0.5)
            {
                j = perm[i]; 
                wtheap = x[j];
                wtsum1 = wtsum1 + wt[j]*xnorm;
                i = i + 1;
            }
        }
        // No normalization required
        else
        {
            while (wtsum1 < 0.5)
            {
                j = perm[i];
                wtheap = x[j];
                wtsum1 = wtsum1 + wt[j];
                i = i + 1;
            }
        }
    }
    return wtheap;
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
int memory_calloc__double(const int n, double **x)
{
    int ierr;
    ierr = memory_malloc__double(n, x);
    memset(*x, 0.0, (size_t) n*sizeof(double));
    return ierr;
}

int memory_calloc__float(const int n, float **x)
{
    int ierr;
    ierr = memory_malloc__float(n, x);
    memset(*x, 0.0f, (size_t) n*sizeof(float));
    return ierr;
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
//============================================================================//
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
    //__assume_aligned(x, 64);
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
double locate_sumDouble64(const int n, double *__restrict__ x)
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
int locate_minLocDouble64(const int n, double *__restrict__ x)
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
int locate_minLocFloat64(const int n, float *__restrict__ x)
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
float locate_sumFloat64(const int n, float *__restrict__ x)
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
                                   const double *__restrict__ varobs,
                                   const double *__restrict__ test,
                                   double *__restrict__ t0,
                                   double *__restrict__ objfn)
{
    const char *fcnm = "locate_l2_gridSearch__double64\0";
    double xnorm __attribute__ ((aligned(64))) = 0.0; 
    double tobs_i __attribute__ ((aligned(64))) = 0.0;
    double wt_i __attribute__ ((aligned(64))) = 0.0;
    const double zero = 0.0;
    const double one = 1.0;
    int ibeg, ierr, iobs;
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
    locate_nullDouble64(ngrd, objfn);
    // Compute the least-squares origin time which is the average reisdual
    if (iwantOT == 1)
    {
        locate_nullDouble64(ngrd, t0);
        xnorm = zero;
        for (iobs=0; iobs<nobs; iobs++)
        {
            if (mask[iobs] != 1){xnorm = xnorm + one/varobs[iobs];}
        }
        for (iobs=0; iobs<nobs; iobs++)
        {
            if (mask[iobs] == 1){continue;}
            tobs_i = tobs[iobs];
            wt_i = one/varobs[iobs];
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
    for (iobs=0; iobs<nobs; iobs++)
    {
        if (mask[iobs] == 1){continue;}
        tobs_i = tobs[iobs];
        wt_i = one/varobs[iobs];
        ibeg = ldgrd*iobs;
        locate_l2_stackObjfn__double64(ngrd, tobs_i, wt_i,
                                       &test[ibeg], t0, objfn);
    }
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
                                  const float *__restrict__ varobs,
                                  const float *__restrict__ test,
                                  float *__restrict__ t0,
                                  float *__restrict__ objfn)
{
    const char *fcnm = "locate_l2_gridSearch__float64\0";
    float xnorm __attribute__ ((aligned(64))) = 0.0f;
    float tobs_i __attribute__ ((aligned(64))) = 0.0f;
    float wt_i __attribute__ ((aligned(64))) = 0.0f;
    const float zero = 0.0f;
    const float one = 1.0f;
    int ibeg, ierr, iobs;
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
    locate_nullFloat64(ngrd, objfn);
    // Compute the least-squares origin time which is the average reisdual
    if (iwantOT == 1)
    {
        locate_nullFloat64(ngrd, t0);
        xnorm = zero;
        for (iobs=0; iobs<nobs; iobs++)
        {
            if (mask[iobs] != 1){xnorm = xnorm + one/varobs[iobs];}
        }
        for (iobs=0; iobs<nobs; iobs++)
        {
            if (mask[iobs] == 1){continue;}
            tobs_i = tobs[iobs];
            wt_i = one/varobs[iobs];
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
    for (iobs=0; iobs<nobs; iobs++)
    {
        if (mask[iobs] == 1){continue;}
        tobs_i = tobs[iobs];
        wt_i = one/varobs[iobs];
        ibeg = ldgrd*iobs;
        locate_l2_stackObjfn__float64(ngrd, tobs_i, wt_i,
                                      &test[ibeg], t0, objfn);
    }
    return ierr;
}
//============================================================================//
int locate_l1_gridSearch__double64(const int ldgrd, //const int ldobs,
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
    double *wt __attribute__ ((aligned(64)));
    double *wtSort __attribute__ ((aligned(64)));
    double *res __attribute__ ((aligned(64)));
    double *resSort __attribute__ ((aligned(64)));
    bool *lmask __attribute__ ((aligned(64)));
    double rsum __attribute__ ((aligned(64))) = 0.0; 
    double wtres __attribute__ ((aligned(64))) = 0.0;
    const double one __attribute__ ((aligned(64))) = 1.0;
    const double zero __attribute__ ((aligned(64))) = 0.0;
    double wtsum, wtsumi;
    int *perm, ierr, igrd, iobs, indx, nobsUse, nsort;
    bool lsort;

    ierr = 0;
    lmask = (bool *) calloc(nobs, sizeof(bool));
    memory_calloc__double(nobs, &wt);
    wtsum = 0.0;
    nobsUse = 0;
    for (iobs=0; iobs<nobs; iobs++)
    {
        if (mask[iobs] == 1)
        {
            lmask[iobs] = true;
        }
        else
        {
            wt[nobsUse] = one/varobs[iobs];
            nobsUse = nobsUse + 1;
        }
        wtsum = wtsum + wt[iobs]; 
    }
    locate_nullDouble64(ngrd, objfn);
    if (iwantOT == 1)
    {
        nsort = 0;
        memory_calloc__double(nobsUse, &res);
        memory_calloc__double(nobsUse, &resSort);
        memory_calloc__double(nobsUse, &wtSort);
        perm = (int *)calloc(nobsUse, sizeof(int));
        // initialize the permutation
        for (iobs=0; iobs<nobs; iobs++)
        {
            perm[iobs] = iobs;
        }
        // Normalize the weights? 
        wtsumi = one/wtsum;
        if (fabs(wtsum - one) > 1.e-14)
        {
            for (iobs=0; iobs<nobsUse; iobs++)
            {
                wt[iobs] = wt[iobs]*wtsumi;
            } 
        }
        // Zero out the origin time
        locate_nullDouble64(ngrd, t0);
        // Loop on grid
        for (igrd=0; igrd<ngrd; igrd++)
        {
            // Compute the residuals in ascending order
            for (iobs=0; iobs<nobsUse; iobs++)
            {
                indx = iobs*ldgrd + igrd; //igrd*ldobs + iobs;
                res[iobs] = tobs[iobs] - test[indx]; 
            }
            // Sort?
            for (iobs=0; iobs<nobs; iobs++)
            {
                resSort[iobs] = res[iobs]; //res[perm[iobs]];
                wtSort[iobs] = wt[iobs]; //wt[perm[iobs]];
            }
            // Optimize weighted median for origin time
            t0[igrd] = weightedMedian__double(nobsUse, resSort,
                                              wtSort, perm, &lsort, &ierr);
//for (iobs=0; iobs<nobs; iobs++){printf("%d\n", perm[iobs]);}
//getchar();
            if (lsort){nsort = nsort + 1;}
            // Compute L1 norm
            for (iobs=0; iobs<nobsUse; iobs++)
            {
                indx = iobs*ldgrd + igrd; //igrd*ldobs + iobs;
                objfn[igrd] = objfn[igrd]
                            + wt[iobs]*fabs(tobs[iobs] - test[indx] - t0[igrd]);
            } 
        }
        free(res);
        free(resSort);
        free(wtSort);
        free(perm);
    }
    else
    {
        // Set the origin time
        locate_setDouble64(ngrd, t0use, t0);
        // Loop on grid points
        //__assume_aligned(wt , 64);
        for (igrd=0; igrd<ngrd; igrd++)
        {
            // Optimize each point
            rsum = zero;
            for (iobs=0; iobs<nobsUse; iobs++)
            {
                indx = iobs*ldgrd + igrd; //igrd*ldobs + iobs;
                rsum = rsum + wt[iobs]*fabs(tobs[iobs] - test[indx] - t0use);
            }
            objfn[igrd] = rsum;
        }
    }
    free(wt);
    free(lmask);
    return ierr;
}

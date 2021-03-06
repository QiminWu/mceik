#ifndef _mceik_struct_h__
#define _mceik_struct_h__ 1

enum pick_type_enum
{
    P_PRIMARY_PICK = 1,  /*!< P first arrival phase pick */
    S_PRIMARY_PICK = 2,  /*!< S first arrival phase pick */
};

struct mceik_catalog_struct
{
    double *xsrc;  /*!< x source ordinate (m) or longitude (deg) [nevents] */
    double *ysrc;  /*!< y source ordinate (m) or latitude (deg) [nevents] */
    double *zsrc;  /*!< station elevation (m) above base of model [nevents] */
    double *tori;  /*!< Epochal origin time (s) [nevents] */
    double *tobs;  /*!< Observed epochal pick times (s) [obsptr[nevents]-1] */
    double *test;  /*!< Estimate epochal pick times (s) corresponding to the
                        source at (xsrc,ysrc,zsrc,t0) [obsptr[nevents]-1] */
    double *varObs;   /*!< Variance (s) in the i'th observations
                           [obsptr[nevents-1] */
    int *luseObs;  /*!< If 0 then the i'th observation is not used 
                        [obsptr[nevents-1] */
    int *pickType; /*!< Phase identifier (P=1) or (S=2) for i'th observation
                        [obsptr[nevents-1] */
    int *statPtr;  /*!< Maps from the i'th observation to the station 
                        information [obsptr[nevents-1]] */
    int *obsPtr;   /*!< Maps from i'th event to start of observations
                        [nevents + 1]*/
    int nevents;   /*!< Number of events */
};

struct mceik_stations_struct
{
    char **netw;     /*!< Network code [nstat] */
    char **stnm;     /*!< Station name [nstat] */
    char **chan;     /*!< Channel (e.g. BH?) [nstat] */
    char **loc;      /*!< Location code [nstat] */
    double *xrec;    /*!< x station ordinate (m) or longitude (deg) [nstat] */
    double *yrec;    /*!< y station ordinate (m) or latitude (deg) [nstat] */
    double *zrec;    /*!< z station elevation (m) above base of model [nstat] */
    double *pcorr;   /*!< P static correction (s) [nstat] */
    double *scorr;   /*!< S static correction (s) [nstat] */
    int *lhasP;      /*!< If 1 then the i'th station has P-picks [nstat] */
    int *lhasS;      /*!< If 1 then the i'th station has S-picks [nstat] */
    int nstat;       /*!< Number of sites */
    int lcartesian;  /*!< If 1 this is a cartesian system */
};

struct catalog_struct
{
    int nevents;
};

struct mcmc_parms_struct
{
    char resdir[512]; /*!< Results directory */
    int nburnIn;      /*!< Number of burn in steps */
    int niter;        /*!< Number of iterations (total forward problems) */
    int keepK;        /*!< Keep the k'th iteration */
};

struct eik_parms_struct
{
    double tol;   /*!< Tolerance (seconds) for declaring convergence */
    int maxit;    /*!< Max number of sweeps */
};

struct mceik_parms_struct
{
    struct mcmc_parms_struct mcparms;   /*!< Markov-Chain Monte-Carlo
                                             parameters */
    struct eik_parms_struct eikparms;   /*!< Eikonal solver parameters */
    char projnm[128];                   /*!< Project name */
    char scratch_dir[512];              /*!< Scratch directory */
    double x0;                          /*!< x origin (m) */
    double y0;                          /*!< y origin (m) */
    double z0;                          /*!< z origin (m) */
    double dx;                          /*!< Grid spacing in x (m) */
    double dy;                          /*!< Grid spacing in y (m) */
    double dz;                          /*!< Grid spacing in z (m) */
    int ndivx;                          /*!< Number of x domain divisions */
    int ndivy;                          /*!< Number of y domain divisions */
    int ndivz;                          /*!< Number of z domain divisions */
    int nrefx;                /*!< Refinement factor in x for translating
                                   from inversion grid to eikonal grid */
    int nrefy;                /*!< Refinement factor in y for translating
                                   from inversion grid to eikonal grid */
    int nrefz;                /*!< Refinement factor in z for translating
                                   from inversion grid to eikonal grid */
};

#endif /* __MCEIK_STRUCT_H__ */

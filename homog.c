#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <mpi.h>
#include <hdf5.h>

#include "mceik_struct.h"
#include "mceik_broadcast.h"
#include "mpiutils.h"
#include "locate.h"
#include "h5io.h"


void freeStations(struct mceik_stations_struct *stations);
void freeCatalog(struct mceik_catalog_struct *catalog);
static void createHomogeneousModel(const int nx, const int ny, const int nz, 
                                   const double vel_ms, int *__restrict__ vmod);
int computeHomogeneousTraveltimes(
    const int nx, const int ny, const int nz, 
    double x0, double  y0, double z0, 
    const double dx, double dy, const double dz, 
    const double xs, const double ys, double zs, 
    const double vel, double *__restrict__ ttimes);
float *double2FloatArray(const int n, double *__restrict__ x);

/*!
 * @brief Homgeneous test case for earthquake location
 */
int main(int argc, char **argv)
{
    const char *fcnm = "xhomog\0";
    const char *projnm = "homog\0";
    char ttimeScratchFile[PATH_MAX];
    MPI_Comm globalComm, intraTableComm, interTableComm;
    struct mceik_catalog_struct catalog;
    struct mceik_stations_struct stations;
    double *ttimes;
    float *ttimes4;
    int *vpmod, *vsmod;
    double dist, dx, dy, dz, velUse, x0, x0Loc, x1, y0, y0Loc, y1, z0, z0Loc, z1;
    int myid, nprocs, nx, ny, nz;
    int imbx, imby, imbz;
    int *tableToStation, *tablePhase;
    int i, ierr, iphase, itable, ix, ix0, ix1, iy, iy0, iy1, iz0, iz1, k,
        nevents, nmodels, nxrec, nyrec,
        ndivx, ndivy, ndivz, ndx, ndy, ndz, nkeep, ntables,
        nwork, nxLoc, nyLoc, nzLoc;
    int globalCommInt, intraTableCommInt, interTableCommInt;
    bool lsaveScratch;
    hid_t tttFileID;
    const double const_vp = 2000.0; // Slower is harder for the solver
    const double const_vs = const_vp/sqrt(3.0);
    const double varVp = 0.25;
    const double varVs = 0.25;
    const int master = 0;
    const int model = 1;
    const int ireord = 1; // Reorder the communicator
    const int iwt = 1;    // Use heuristic waiting 
    //------------------------------------------------------------------------//
    //
    // Initialize mpi
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    // initializations
    memset(&catalog,  0, sizeof(struct mceik_catalog_struct));
    memset(&stations, 0, sizeof(struct mceik_stations_struct));
    vpmod = NULL;
    vsmod = NULL;
    x0 = 0.0;
    y0 = 0.0;
    z0 = 0.0;
    x1 = 31.e3;
    y1 = 28.e3;
    z1 = 25.e3;
    dx = 1000.0;
    dy = 1000.0;
    dz = 1000.0;
    ndivx = 2;
    ndivy = 1;
    ndivz = 1;
    nmodels = 1;
    lsaveScratch = false;
    nx = (int) ((x1 - x0)/dx + 0.5) + 1;
    ny = (int) ((y1 - y0)/dy + 0.5) + 1;
    nz = (int) ((z1 - z0)/dz + 0.5) + 1;
    if (myid == master){printf("%s: Splitting the commuicator...\n", fcnm);}
    globalCommInt = (int) (MPI_Comm_c2f(MPI_COMM_WORLD));
    mpiutils_initialize3d(&globalCommInt, &ireord, &iwt,
                          &ndivx, &ndivy, &ndivz, &ierr);
    if (ierr != 0)
    {
        printf("%s: Error splitting the communicators!\n", fcnm);
        MPI_Abort(MPI_COMM_WORLD, 30);
    }
    // get the communicator IDs
    mpiutils_getCommunicators(&globalCommInt, &intraTableCommInt,
                              &interTableCommInt, &ierr);
    if (ierr != 0)
    {
        printf("%s: Error getting communicators\n", fcnm);
        MPI_Abort(MPI_COMM_WORLD, 30);
    }
    globalComm = MPI_Comm_f2c((MPI_Fint) globalCommInt);
    intraTableComm = MPI_Comm_f2c((MPI_Fint) intraTableCommInt);
    interTableComm = MPI_Comm_f2c((MPI_Fint) interTableCommInt);
    MPI_Barrier(globalComm);
    if (myid == master)
    {
        srand(2016);
        //srand(time(0)); // will create random number
        nevents = 4;
        // Scatter the receivers across the free surface
        printf("%s: Creating station locations...\n", fcnm);
        nxrec = 2;
        nyrec = 3;
        stations.nstat = nxrec*nyrec;
        stations.lcartesian = 1; // Cartesian system
        stations.netw  = (char **)calloc(stations.nstat, sizeof(char *));
        stations.stnm  = (char **)calloc(stations.nstat, sizeof(char *));
        stations.chan  = (char **)calloc(stations.nstat, sizeof(char *));
        stations.loc   = (char **)calloc(stations.nstat, sizeof(char *));
        for (i=0; i<stations.nstat; i++) 
        {
            stations.netw[i] = (char *)calloc(64, sizeof(char));
            stations.stnm[i] = (char *)calloc(64, sizeof(char));
            stations.chan[i] = (char *)calloc(64, sizeof(char));
            stations.loc[i]  = (char *)calloc(64, sizeof(char));
            strcpy(stations.netw[i],  "NA\0");
            sprintf(stations.stnm[i], "RC%d", i+1);
            strcpy(stations.chan[i],  "HH?\0");
            strcpy(stations.loc[i],   "00\0");
        }
        // Randomly scatter stations
        stations.xrec  = (double *)calloc(stations.nstat, sizeof(double));
        stations.yrec  = (double *)calloc(stations.nstat, sizeof(double));
        stations.zrec  = (double *)calloc(stations.nstat, sizeof(double));
        for (iy=0; iy<nyrec; iy++)
        {
            for (ix=0; ix<nxrec; ix++)
            {
                stations.xrec[iy*nxrec+ix]
                     = x0 + ((int) ((double) rand()/RAND_MAX*(nx - 1)))*dx;
                stations.yrec[iy*nxrec+ix]
                     = y0 + ((int) ((double) rand()/RAND_MAX*(ny - 1)))*dy;
                //stations.xrec[iy*nxrec+ix] = x0 + (nx/4 + ix)*dx;
                //stations.yrec[iy*nxrec+ix] = y0 + (ny/4 + iy)*dy;
                stations.zrec[iy*nxrec+ix] = z1;
                if (stations.xrec[iy*nxrec+ix] < x0 ||
                    stations.xrec[iy*nxrec+ix] > x1 ||
                    stations.yrec[iy*nxrec+ix] < y0 ||
                    stations.yrec[iy*nxrec+ix] > y1 ||
                    stations.zrec[iy*nxrec+ix] < z0 || 
                    stations.zrec[iy*nxrec+ix] > z1)
                {
                    printf("%s: Station out of bounds!\n", fcnm);
                    MPI_Abort(MPI_COMM_WORLD, 20);
                }
                printf("%f %f %f\n", stations.xrec[iy*nxrec+ix],
                                     stations.yrec[iy*nxrec+ix],
                                     stations.zrec[iy*nxrec+ix]);
            }
        }
        stations.pcorr = (double *)
                         calloc((size_t) stations.nstat, sizeof(double));
        stations.scorr = (double *)
                         calloc((size_t) stations.nstat, sizeof(double));
        stations.lhasP = (int *)
                         calloc((size_t) stations.nstat, sizeof(int));
        stations.lhasS = (int *)
                         calloc((size_t) stations.nstat, sizeof(int));
        // Make some events
        printf("%s: Creating events...\n", fcnm);
        catalog.nevents = nevents;
        nwork = 2*catalog.nevents*stations.nstat;
        catalog.xsrc = (double *)
                       calloc((size_t) catalog.nevents, sizeof(double));
        catalog.ysrc = (double *)
                       calloc((size_t) catalog.nevents, sizeof(double));
        catalog.zsrc = (double *)
                       calloc((size_t) catalog.nevents, sizeof(double));
        catalog.t0 = (double *)
                     calloc((size_t) catalog.nevents, sizeof(double));
        catalog.tobs = (double *) calloc((size_t) nwork, sizeof(double));
        catalog.test = (double *) calloc((size_t) nwork, sizeof(double));
        catalog.var =  (double *) calloc((size_t) nwork, sizeof(double));
        catalog.luseObs  = (int *) calloc((size_t) nwork, sizeof(int));
        catalog.pickType = (int *) calloc((size_t) nwork, sizeof(int));
        catalog.statPtr  = (int *) calloc((size_t) nwork, sizeof(int));
        catalog.obsPtr   = (int *)
                           calloc((size_t) catalog.nevents + 1, sizeof(int));
        
        nkeep = 0;
        for (i=0; i<catalog.nevents; i++)
        {
            catalog.xsrc[i] = x0 + (x1 - x0)*(double) rand()/RAND_MAX;
            catalog.ysrc[i] = y0 + (y1 - y0)*(double) rand()/RAND_MAX;
            catalog.zsrc[i] = z0 + (z1 - z0)*(double) rand()/RAND_MAX;
            // Now attach some theoreticals
            for (k=0; k<stations.nstat; k++)
            {
                dist = sqrt( pow(stations.xrec[k] - catalog.xsrc[i], 2)
                           + pow(stations.yrec[k] - catalog.ysrc[i], 2)
                           + pow(stations.zrec[k] - catalog.zsrc[i], 2) );
                for (iphase=1; iphase<=2; iphase++)
                {
                    if (iphase == P_PRIMARY_PICK)
                    {
                        catalog.tobs[nkeep] = dist/const_vp;
                        catalog.var[nkeep] = varVp;
                    }
                    else if (iphase == S_PRIMARY_PICK)
                    {
                        catalog.tobs[nkeep] = dist/const_vs;
                        catalog.var[nkeep] = varVs;
                    }
                    else
                    {
                        printf("%s: Invalid phase %d\n", fcnm, iphase);
                        MPI_Abort(MPI_COMM_WORLD, 30);
                    }
                    catalog.luseObs[nkeep] = 1;
                    catalog.pickType[nkeep] = iphase;
                    catalog.statPtr[nkeep] = k;
                    nkeep = nkeep + 1;
                }
                // Do i have a P and S phase?
                if (i == 0)
                {
                    stations.lhasP[k] = 1;
                    stations.lhasS[k] = 1;
                }
            }
            catalog.obsPtr[i+1] = nkeep;
        }
    }

    // Cook up a list of bogus earthquakes
    // Initialize the eikonal solver
    if (myid == master)
    {
        printf("%s: Generating constant velocity model...\n", fcnm);
        vpmod = (int *)calloc((size_t) (nx*ny*nz), sizeof(int));
        vsmod = (int *)calloc((size_t) (nx*ny*nz), sizeof(int));
        createHomogeneousModel(nx, ny, nz, const_vp, vpmod);
        createHomogeneousModel(nx, ny, nz, const_vs, vsmod);
    }
    else
    {
        vpmod = (int *)calloc((size_t) (nx*ny*nz), sizeof(int));
        vsmod = (int *)calloc((size_t) (nx*ny*nz), sizeof(int));
    }
    // Distribute the inversion model to all
    MPI_Bcast(vpmod, nx*ny*nz, MPI_INTEGER, master, MPI_COMM_WORLD);
    MPI_Bcast(vsmod, nx*ny*nz, MPI_INTEGER, master, MPI_COMM_WORLD);
    // Distribute the station information
    broadcast_stations(MPI_COMM_WORLD, master, &stations);
    // Distribute the catalog
    broadcast_catalog(MPI_COMM_WORLD, master, &catalog);
    // Make the local model
   ix0 = 0;
   iy0 = 0;
   iz0 = 0;
    mpiutils_grd2ijk(&myid, &ndivx, &ndivy, &ndivz,
                     &imbx, &imby, &imbz, &ierr);
    if (ierr != 0){printf("%s: Failed to map rank to block\n", fcnm);}
    ndx = fmax(nx/ndivx, 1);
    ndy = fmax(ny/ndivy, 1);
    ndz = fmax(nz/ndivz, 1);
    ix0 = imbx*ndx;
    iy0 = imby*ndy;
    iz0 = imbz*ndz;
    ix1 = (imbx + 1)*ndx;
    iy1 = (imby + 1)*ndy;
    iz1 = (imbz + 1)*ndz;
    if (imbx + 1 == ndivx){ix1 = nx;}
    if (imby + 1 == ndivy){iy1 = ny;}
    if (imbz + 1 == ndivz){iz1 = nz;}
    nxLoc = ix1 - ix0;
    nyLoc = iy1 - iy0;
    nzLoc = iz1 - iz0;
    x0Loc = x0 + ix0*dx;
    y0Loc = y0 + iy0*dy;
    z0Loc = z0 + iz0*dz;

    // Initialize the HDF5 traveltime file
    int myTableID;
    MPI_Comm_rank(interTableComm, &myTableID);
    memset(ttimeScratchFile, 0, sizeof(ttimeScratchFile));
    sprintf(ttimeScratchFile, "%s_%d", projnm, myTableID+1); 
    ierr = eikonal_h5io_initialize(intraTableComm, //MPI_COMM_WORLD,
                                   "./\0",
                                   ttimeScratchFile, //"test\0",
                                   ix0, iy0, iz0,
                                   nx, ny, nz,
                                   nxLoc, nyLoc, nzLoc,
                                   nmodels,
                                   stations.nstat,
                                   lsaveScratch,
                                   x0, y0, z0,
                                   dx, dy, dz,
                                   &tttFileID);
    if (ierr != 0)
    {
        printf("%s: Error initializing H5 file\n", fcnm);
        return EXIT_FAILURE;
    }
    // Determine how many traveltime tables to make 
    ntables = 0;
    for (i=0; i<stations.nstat; i++)
    {
        if (stations.lhasP[i] == 1){ntables = ntables + 1;}
        if (stations.lhasS[i] == 1){ntables = ntables + 1;}
    }
    tableToStation = (int *) calloc((size_t) ntables, sizeof(int));
    tablePhase     = (int *) calloc((size_t) ntables, sizeof(int));
    itable = 0;
    for (i=0; i<stations.nstat; i++)
    {
        if (stations.lhasP[i] == 1)
        {
            tableToStation[itable] = i + 1;
            tablePhase[itable] = 1;
            itable = itable + 1;
        }
        if (stations.lhasS[i] == 1)
        {
            tableToStation[itable] = i + 1;
            tablePhase[itable] = 2;
            itable = itable + 1;
        }
    }
    if (myid == master)
    {
        printf("%s: Will compute %d travel time tables\n", fcnm, ntables);
    }
    // This is the parallel loop on tables
    ttimes = (double *) calloc((size_t) nxLoc*nyLoc*nzLoc, sizeof(double)); 
    MPI_Barrier(MPI_COMM_WORLD);
    for (itable=0; itable<ntables; itable++)
    {
        iphase = tablePhase[itable];
        k = tableToStation[itable];
        velUse = const_vp;
        if (iphase == S_PRIMARY_PICK)
        {
            velUse = const_vs;
        }
        else
        {
            if (iphase != P_PRIMARY_PICK)
            {
                printf("%s: Invalid phase\n", fcnm);
                MPI_Abort(MPI_COMM_WORLD, 30);
            }
        }
        // Compute traveltimes from station to all points in medium
        ierr = computeHomogeneousTraveltimes(nxLoc, nyLoc, nzLoc,
                                             x0Loc, y0Loc, z0Loc,
                                             dx, dy, dz,
                                             stations.xrec[k],
                                             stations.yrec[k],
                                             stations.zrec[k],
                                             velUse, ttimes);
        if (ierr != 0)
        {
            printf("%s: Error computing homogeneous traveltimes\n", fcnm);
            MPI_Abort(MPI_COMM_WORLD, 30);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        // Convert to float
        ttimes4 = double2FloatArray(nxLoc*nyLoc*nzLoc, ttimes);
        // Save the data
        ierr = eikonal_h5io_writeTravelTimes(intraTableComm, //MPI_COMM_WORLD,
                                             tttFileID,
                                             k, model,   
                                             iphase,
                                             ix0, iy0, iz0,
                                             nxLoc, nyLoc, nzLoc,
                                             ttimes4);
        if (ierr != 0)
        {
            printf("%s: Failed writing traveltimes %d\n", fcnm, myid);
            MPI_Abort(MPI_COMM_WORLD, 30);
        }
        // Verify
        memset(ttimes4, 0, (size_t) (nxLoc*nyLoc*nzLoc)*sizeof(float));
        ierr = eikonal_h5io_readTravelTimes(intraTableComm, //MPI_COMM_WORLD,
                                            tttFileID,
                                            k, model,
                                            iphase,
                                            ix0, iy0, iz0,
                                            nxLoc, nyLoc, nzLoc,
                                            ttimes4);
        if (ierr != 0)
        {
            printf("%s: Error loading traveltimes\n", fcnm);
            MPI_Abort(MPI_COMM_WORLD, 30);
        }
        double difMax = 0.0;
        int in;
        for (in=0; in<nxLoc*nyLoc*nzLoc; in++)
        {
            difMax = fmax(difMax, fabs(ttimes4[in] - (float) ttimes[in]));
        }
        if (difMax > 1.e-5)
        {
            printf("%s: Failed to read/write traveltime verification\n", fcnm);
            MPI_Abort(MPI_COMM_WORLD, 30);
        }
        free(ttimes4);
    }
    free(ttimes);
    // I am now ready to locate some earthquakes
   int iverb = 0;
    locate3d_initialize(&globalCommInt, &iverb, &nx, &ny, &nz,
                        &ndivx, &ndivy, &ndivz, &ierr);
    // Finalize
    eikonal_h5io_finalize(MPI_COMM_WORLD, &tttFileID);
    freeStations(&stations);
    freeCatalog(&catalog);
    if (vpmod != NULL){free(vpmod);}
    if (vsmod != NULL){free(vsmod);}
    mpiutils_finalize();
    MPI_Finalize();
    return EXIT_SUCCESS;
}
//============================================================================//
/*!
 * @brief Frees memory on the catalog structure
 *
 * @param[in,out] catalog    on exit all memory on catalog has been freed and
 *                           any scalars have been nulled out
 *
 * @author Ben Baker
 *
 */
void freeCatalog(struct mceik_catalog_struct *catalog)
{
    if (catalog->xsrc     != NULL){free(catalog->xsrc);}
    if (catalog->ysrc     != NULL){free(catalog->ysrc);}
    if (catalog->zsrc     != NULL){free(catalog->zsrc);}
    if (catalog->t0       != NULL){free(catalog->t0);}
    if (catalog->tobs     != NULL){free(catalog->tobs);}
    if (catalog->test     != NULL){free(catalog->test);}
    if (catalog->var      != NULL){free(catalog->var);}
    if (catalog->luseObs  != NULL){free(catalog->luseObs);}
    if (catalog->pickType != NULL){free(catalog->pickType);}
    if (catalog->statPtr  != NULL){free(catalog->statPtr);}
    if (catalog->obsPtr   != NULL){free(catalog->obsPtr);}
    memset(catalog, 0, sizeof(struct catalog_struct));
    return;
}
//============================================================================//
/*!
 * @brief Frees the station structure
 *
 * @param[in,out] station     on input contains the station list.
 *                            on output all memory has been released from
 *                            the station list and it has been reset.
 *
 * @author Ben Baker
 *
 */
void freeStations(struct mceik_stations_struct *stations)
{
    int i;
    for (i=0; i<stations->nstat; i++) 
    {
        if (stations->netw != NULL)
        {
            if (stations->netw[i] != NULL){free(stations->netw[i]);}
        }
        if (stations->stnm != NULL)
        {
            if (stations->stnm[i] != NULL){free(stations->stnm[i]);}
        }
        if (stations->chan != NULL)
        {
            if (stations->chan[i] != NULL){free(stations->chan[i]);}
        }
        if (stations->loc != NULL)
        {
            if (stations->loc[i]  != NULL){free(stations->loc[i]);}
        }
    }
    if (stations->netw  != NULL){free(stations->netw);}
    if (stations->stnm  != NULL){free(stations->stnm);}
    if (stations->chan  != NULL){free(stations->chan);}
    if (stations->loc   != NULL){free(stations->loc);}
    if (stations->xrec  != NULL){free(stations->xrec);}
    if (stations->yrec  != NULL){free(stations->yrec);}
    if (stations->zrec  != NULL){free(stations->zrec);}
    if (stations->pcorr != NULL){free(stations->pcorr);}
    if (stations->scorr != NULL){free(stations->scorr);}
    if (stations->lhasP != NULL){free(stations->lhasP);}
    if (stations->lhasS != NULL){free(stations->lhasS);}
    memset(stations, 0, sizeof(struct mceik_stations_struct));
    return;
}
//============================================================================//
/*!
 * @brief Sets a homogeneous model 
 *
 * @param[in] nx      number of x grid points in grid
 * @param[in] ny      number of y grid points in grid
 * @param[in] nz      number of z grid points in grid
 * @param[in] vel_ms  velocity (m/s)
 *
 * @param[in] vmod    const velocity model [nz x ny x nx]
 *
 * @author Ben Baker 
 *
 */
static void createHomogeneousModel(const int nx, const int ny, const int nz,
                                   const double vel_ms, int *__restrict__ vmod)
{
    int indx, ix, iy, iz, nxy;
    nxy = nx*ny;
    for (iz=0; iz<nz; iz++)
    {
        for (iy=0; iy<ny; iy++)
        {
            for (ix=0; ix<nx; ix++)
            {
                indx = iz*nxy + iy*nx + ix;
                vmod[indx] = (int) vel_ms;
            }
        }
    }
    return;
}
//============================================================================//
/*!
 * @brief Computes the traveltimes to all points in a constant gridded
 *        velocity model
 *
 * @param[in] nx        number of x grid points in model
 * @param[in] ny        number of y grid points in model
 * @param[in] nz        number of z grid points in model
 * @param[in] x0        x origin (m)
 * @param[in] y0        y origin (m)
 * @param[in] z0        z origin (m)
 * @param[in] dx        grid spacing in x (m)
 * @param[in] dy        grid spacing in y (m)
 * @param[in] dz        grid spacing in z (m)
 * @param[in] xs        x source position (m)
 * @param[in] ys        y source position (m)
 * @param[in] zs        z source position (m)
 * @param[in] vel       constant medium velocity (m/s)
 *
 * @param[out] ttimes   traveltimes (s) at each point in medium [nx*ny*nz].
 *                      the k'th index for the (ix,iy,iz)'th grid point
 *                      is accessed by k = iz*nx*ny + iy*nx + ix for 
 *                      ix=0,1,...,nx-1, iy=0,1,...,ny-1, iz=0,1,...,nz-1.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker
 *
 */ 
int computeHomogeneousTraveltimes(
    const int nx, const int ny, const int nz,
    double x0, double  y0, double z0,
    const double dx, double dy, const double dz,
    const double xs, const double ys, double zs,
    const double vel, double *__restrict__ ttimes)
{
    double dist, slow, x, y, z;
    int indx, ix, iy, iz, nxy;
    nxy = nx*ny;
    slow = 1.0/vel;
    for (iz=0; iz<nz; iz++)
    {
        for (iy=0; iy<ny; iy++)
        {
            for (ix=0; ix<nx; ix++)
            {
                x = x0 + (double) ix*dx;
                y = y0 + (double) iy*dy;
                z = z0 + (double) iz*dz;
                dist = sqrt(pow(xs-x, 2) + pow(ys-y, 2) + pow(zs-z, 2));
                indx = iz*nxy + iy*nx + ix;
                ttimes[indx] = dist*slow; 
            }
        }
    }
    return 0;
}
//============================================================================//

float *double2FloatArray(const int n, double *__restrict__ x)
{
    float *x4;
    int i;
    x4 = (float *) calloc((size_t) n, sizeof(float));
    #pragma omp simd
    for (i=0; i<n; i++)
    {
        x4[i] = (float) x[i];
    }
    return x4;
}

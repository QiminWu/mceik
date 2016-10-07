#include <stdio.h>
#include <stdlib.h>
#include "mceik_broadcast.h"

/*!
 * @brief Broadcasts the station list from root to all processes on the
 *        communicator
 *
 * @author Ben Baker
 *
 * @copyright Apache 2
 *
 */
void broadcast_stations(MPI_Comm comm, const int root,
                        struct mceik_stations_struct *stations)
{
    int k, myid;
    //------------------------------------------------------------------------//
    //
    // Get communicator info
    MPI_Comm_rank(comm, &myid);
    // Integers
    MPI_Bcast(&stations->nstat,      1, MPI_INT, root, comm);
    MPI_Bcast(&stations->lcartesian, 1, MPI_INT, root, comm); 
    if (stations->nstat < 1){return;}
    // Arrays
    if (myid != root)
    {
        stations->lhasP = (int *)calloc((size_t) stations->nstat, sizeof(int));
        stations->lhasS = (int *)calloc((size_t) stations->nstat, sizeof(int));
    }
    MPI_Bcast(stations->lhasP, stations->nstat, MPI_INT, root, comm);
    MPI_Bcast(stations->lhasS, stations->nstat, MPI_INT, root, comm);
    // Doubles
    if (myid != root)
    {
        stations->xrec = (double *)
                         calloc((size_t) stations->nstat, sizeof(double));
        stations->yrec = (double *)
                         calloc((size_t) stations->nstat, sizeof(double));
        stations->zrec = (double *)
                         calloc((size_t) stations->nstat, sizeof(double));
        stations->pcorr = (double *)
                          calloc((size_t) stations->nstat, sizeof(double));
        stations->scorr = (double *)
                          calloc((size_t) stations->nstat, sizeof(double));
    }
    MPI_Bcast(stations->xrec,  stations->nstat, MPI_DOUBLE, root, comm);
    MPI_Bcast(stations->yrec,  stations->nstat, MPI_DOUBLE, root, comm);
    MPI_Bcast(stations->zrec,  stations->nstat, MPI_DOUBLE, root, comm);
    MPI_Bcast(stations->pcorr, stations->nstat, MPI_DOUBLE, root, comm);
    MPI_Bcast(stations->scorr, stations->nstat, MPI_DOUBLE, root, comm);
    // Characters
    if (myid != root)
    {
        stations->netw = (char **)
                         calloc((size_t) stations->nstat, sizeof(char *));
        stations->stnm = (char **)
                         calloc((size_t) stations->nstat, sizeof(char *));
        stations->chan = (char **)
                         calloc((size_t) stations->nstat, sizeof(char *));
        stations->loc  = (char **)
                         calloc((size_t) stations->nstat, sizeof(char *));
        for (k=0; k<stations->nstat; k++)
        {
            stations->netw[k] = (char *)calloc(64, sizeof(char));
            stations->stnm[k] = (char *)calloc(64, sizeof(char));
            stations->chan[k] = (char *)calloc(64, sizeof(char));
            stations->loc[k]  = (char *)calloc(64, sizeof(char));
        }
    }
    for (k=0; k<stations->nstat; k++)
    {
        MPI_Bcast(stations->netw[k], 64, MPI_CHAR, root, comm);
        MPI_Bcast(stations->stnm[k], 64, MPI_CHAR, root, comm);
        MPI_Bcast(stations->chan[k], 64, MPI_CHAR, root, comm);
        MPI_Bcast(stations->loc[k],  64, MPI_CHAR, root, comm);
    }
    return;
} 
//============================================================================//
/*!
 * @brief Distributes the catalog to all processes on communicator from
 *        root process
 *
 * @param[in] comm         MPI communicator
 * @param[in] root         root process ID on communicator with catalog
 *
 * @param[in,out] catalog  on input this is the catalog on the root process.
 *                         on output this is the catalog on all the processes
 *                         on the communicator.
 *
 * @author Ben Baker
 *
 * @copyright Apache 2
 *
 */
void broadcast_catalog(MPI_Comm comm, const int root,
                       struct mceik_catalog_struct *catalog)
{
    int myid, nevents, nwork;
    //------------------------------------------------------------------------//
    //
    // Get communicator info
    MPI_Comm_rank(comm, &myid);
    // Sizes
    MPI_Bcast(&catalog->nevents, 1, MPI_INT, root, comm);
    nevents = catalog->nevents;
    if (nevents < 1){return;}
    if (myid == root){nwork = catalog->obsPtr[nevents];}
    MPI_Bcast(&nwork, 1, MPI_INT, root, comm);
    // Integers
    if (myid != root)
    {
        catalog->luseObs  = (int *) calloc((size_t) nwork,     sizeof(int));
        catalog->pickType = (int *) calloc((size_t) nwork,     sizeof(int));
        catalog->statPtr  = (int *) calloc((size_t) nwork,     sizeof(int));
        catalog->obsPtr   = (int *) calloc((size_t) nevents+1, sizeof(int));
    }
    MPI_Bcast(catalog->luseObs,  nwork,     MPI_INT, root, comm);
    MPI_Bcast(catalog->pickType, nwork,     MPI_INT, root, comm);
    MPI_Bcast(catalog->statPtr,  nwork,     MPI_INT, root, comm);
    MPI_Bcast(catalog->obsPtr,   nevents+1, MPI_INT, root, comm);
    // Doubles
    if (myid != root)
    {
        catalog->xsrc = (double *) calloc((size_t) nevents, sizeof(double));
        catalog->ysrc = (double *) calloc((size_t) nevents, sizeof(double));
        catalog->zsrc = (double *) calloc((size_t) nevents, sizeof(double));
        catalog->t0   = (double *) calloc((size_t) nevents, sizeof(double));
        catalog->tobs = (double *) calloc((size_t) nwork,   sizeof(double));
        catalog->test = (double *) calloc((size_t) nwork,   sizeof(double));
        catalog->var  = (double *) calloc((size_t) nwork,   sizeof(double));
    }
    MPI_Bcast(catalog->xsrc, nevents, MPI_DOUBLE, root, comm);
    MPI_Bcast(catalog->ysrc, nevents, MPI_DOUBLE, root, comm);
    MPI_Bcast(catalog->zsrc, nevents, MPI_DOUBLE, root, comm);
    MPI_Bcast(catalog->t0,   nevents, MPI_DOUBLE, root, comm);
    MPI_Bcast(catalog->tobs, nwork,   MPI_DOUBLE, root, comm);
    MPI_Bcast(catalog->test, nwork,   MPI_DOUBLE, root, comm);
    MPI_Bcast(catalog->var,  nwork,   MPI_DOUBLE, root, comm);
    return;
}

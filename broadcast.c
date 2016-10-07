#include <stdio.h>
#include <stdlib.h>
#include "mceik_broadcast.h"

/*!
 * @brief Broadcasts the station list from root to all processes on the
 *        communicator
 *
 * @author Ben Baker
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

#ifndef _mceik_broadcast_h__
#define _mceik_broadcast_h__ 1
#include <mpi.h>
#include "mceik_struct.h"

#ifdef __cplusplus
extern "C" 
{
#endif

/* Broadcasts the catalog */
void broadcast_catalog(MPI_Comm comm, const int root,
                       struct mceik_catalog_struct *catalog);
/* Broadcasts the stations list */
void broadcast_stations(MPI_Comm comm, const int root,
                        struct mceik_stations_struct *stations);

#ifdef __cplusplus
}
#endif
#endif /* _mceik_broadcast_h__ */


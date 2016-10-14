#ifndef _locate_h__
#define _locate_h__ 1
#include <mpi.h>

#ifdef __cplusplus
extern "C" 
{
#endif

void locate3d_initialize(const int *comm, const int *iverb,
                         const int *nx, const int *ny, const int *nz,
                         const int *ndivx, const int *ndivy, const int *ndivz,
                         int *ierr);
#ifdef __cpluplus
}
#endif
#endif /* _locate_h__ */


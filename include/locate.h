#ifndef _locate_h__
#define _locate_h__ 1
#include <mpi.h>

#ifdef __cplusplus
extern "C" 
{
#endif

void locate3d_initialize(const int *comm, const int *iverb,
                         const long *tttFileID, const long *locFileID,
                         const int *ndivx, const int *ndivy, const int *ndivz,
                         int *ierr);

void locate3d_gridsearch(const int *tttFileID, const int *locFileID,
                         const int *model,
                         const int *job, const int *ngrd, const int *nsrc,
                         const int *luseObs, const double *statCor,
                         const double *tori, const double *varobs,
                         const double *tobs, double *hypo, int *ierr);

void locate3d_finalize(void);

#ifdef __cpluplus
}
#endif
#endif /* _locate_h__ */


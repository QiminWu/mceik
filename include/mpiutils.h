#ifndef __MPIUTILS_H__
#define __MPIUTILS_H__

#ifdef __cplusplus
extern "C"
{
#endif

void mpiutils_getCommunicators(int *globalComm, int *intraTableComm,
                               int *interTableComm, int *ierr);
void mpiutils_grd2ijk(const int *igrd,
                      const int *nx, const int *ny, const int *nz,
                      int *i, int *j, int *k, int *ierr);
void mpiutils_initialize3d(const int *comm, const int *ireord, const int *iwt,
                           const int *ndivx, const int *ndivy, const int *ndivz,
                           int *ierr); 
void mpiutils_finalize(void); 
                            

#ifdef __cplusplus
}
#endif
#endif /*! __MPIUTILS_H__ */

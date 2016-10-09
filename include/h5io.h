#include <stdbool.h>
#include <limits.h>
#include <hdf5.h>
#include <mpi.h>
#ifndef _h5io_h__
#define _h5io_h__ 1
#ifdef __cplusplus
extern "C" 
{
#endif
/* Finialize the hdf5 io */
int eikonal_h5io_finalize(const MPI_Comm comm, hid_t *tttFileID);
/* Sets the traveltime table filename */
int eikonal_h5io_setFileName(const char *dirnm, const char *projnm,
                             char fileName[PATH_MAX]);
/* Sets the traveltime table dataset name for the given station and model */
void eikonal_h5io_setTravelTimeName(const int model, const int station,
                                    const bool isP, char dataSetName[512]);
/* Initialize the HDF5 traveltime table */
int eikonal_h5io_initialize(const MPI_Comm comm,
                            const char *dirnm,
                            const char *projnm,
                            const int ix0, const int iy0, const int iz0,
                            const int nx, const int ny, const int nz, 
                            const int nxLoc, const int nyLoc, const int nzLoc,
                            const int nmodels,
                            const int nstations,
                            const bool lsaveScratch,
                            const double x0, const double y0, const double z0, 
                            const double dx, const double dy, const double dz, 
                            hid_t *tttFileID);
/* Reads a traveltime grid */
int eikonal_h5io_readTravelTimes(const MPI_Comm comm,
                                 const hid_t tttFileID,
                                 const int station, const int model,
                                 const int iphase,
                                 const int ix0, const int iy0, const int iz0,
                                 const int nxLoc, const int nyLoc,
                                 const int nzLoc,
                                 float *__restrict__ ttimes);
/* Writes a traveltime grid */
int eikonal_h5io_writeTravelTimes(const MPI_Comm comm,
                                  const hid_t tttFileID,
                                  const int station, const int model,
                                  const int iphase,
                                  const int ix0, const int iy0, const int iz0,
                                  const int nxLoc, const int nyLoc,
                                  const int nzLoc,
                                  const float *__restrict__ ttimes);

#ifdef __cplusplus
}
#endif
#endif /* __H5IO_H__ */

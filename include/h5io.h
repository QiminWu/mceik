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

enum fileName_enum
{
    TRAVELTIME_FILE = 1, 
    LOCATION_FILE = 2,
};
/* Finialize the hdf5 io */
int eikonal_h5io_finalize(const MPI_Comm comm, hid_t *tttFileID);
/* Sets the traveltime table filename */
int eikonal_h5io_setFileName(enum fileName_enum job,
                             const char *dirnm, const char *projnm,
                             char fileName[PATH_MAX]);
/* Sets the traveltime table dataset name for the given station and model */
void eikonal_h5io_setTravelTimeName(const int model, const int station,
                                    const bool isP, char dataSetName[512]);
/* Get the model dimensions */
void eikonal_h5io_getModelDimensionsF(const long *inFileID,
                                      int *nx, int *ny, int *nz, 
                                      int *ierr);
int eikonal_h5io_getModelDimensions(const hid_t fileID,
                                    int *nx, int *ny, int *nz);
/* Read model */
void eikonal_h5io_readModelF(const int *comm,
                             const long *inFileID,
                             const int *ix0, const int *iy0, const int *iz0, 
                             const int *nxLoc, const int *nyLoc,
                             const int *nzLoc,
                             float *__restrict__ xlocs,
                             float *__restrict__ ylocs,
                             float *__restrict__ zlocs,
                             int *ierr);
int eikonal_h5io_readModel(const MPI_Comm comm,
                           const hid_t fileID,
                           const int ix0, const int iy0, const int iz0,
                           const int nxLoc, const int nyLoc, const int nzLoc,
                           float *__restrict__ xlocs,
                           float *__restrict__ ylocs,
                           float *__restrict__ zlocs);
/* Initialize the HDF5 locations */
int eikonal_h5io_initLocations(
    const MPI_Comm comm,
    const char *dirnm, const char *projnm,
    const int ix0, const int iy0, const int iz0,
    const int nx, const int ny, const int nz, 
    const int nxLoc, const int nyLoc, const int nzLoc,
    const int nmodels, const int nevents,
    const double x0, const double y0, const double z0, 
    const double dx, const double dy, const double dz, 
    hid_t *locFileID);
/* Initialize the HDF5 traveltime table */
int eikonal_h5io_initTTables(const MPI_Comm comm,
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
/* Makes the model group during initialization */
int eikonal_h5io_makeModelGroup(
    const MPI_Comm comm, const hid_t fileID,
    const int ix0, const int iy0, const int iz0,
    const int nxGlob, const int nyGlob, const int nzGlob,
    const int nxLoc, const int nyLoc, const int nzLoc,
    const int nxMax, const int nyMax, const int nzMax,
    const double dx, const double dy, const double dz, 
    const double x0, const double y0, const double z0);
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
/* Writes the logJPDF locations */
int eikonal_h5io_writeLocationLogJPDF(
    const MPI_Comm comm, const hid_t locFileID,
    const int model, const int event,
    const int ix0, const int iy0, const int iz0,
    const int nxLoc, const int nyLoc, const int nzLoc,
    const float *__restrict__ logJPDF);

#ifdef __cplusplus
}
#endif
#endif /* __H5IO_H__ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include <mpi.h>
#include "h5io.h"
#include "os.h"

int eikonal_h5io_setFileName(enum fileName_enum job,
                             const char *dirnm, const char *projnm,
                             char fileName[PATH_MAX])
{
    const char *fcnm = "eikonal_h5io_setFileName\0";
    int lendir;
    memset(fileName, 0, PATH_MAX*sizeof(char));
    if (dirnm == NULL)
    {
        strcpy(fileName, "./\0");
    }
    else
    {
        lendir = strlen(dirnm);
        if (lendir > 0)
        {
            strcpy(fileName, dirnm);
            // Make sure it ends with a "/"
            if (fileName[lendir-1] != '/'){strcat(fileName, "/\0");}
            // Require the directory exists
            if (!os_path_isdir(fileName))
            {

            }
        }
        else
        {
            strcpy(fileName, "./\0");
        }
    }
    if (projnm == NULL)
    {
        printf("%s: Project name must be defined\n", fcnm);
        return -1;
    }
    if (strlen(projnm) == 0)
    {
        printf("%s: Project can't be empty\n", fcnm);
        return -1;
    }
    strcat(fileName, projnm);
    if (job == TRAVELTIME_FILE)
    {
        strcat(fileName, "_ttimes.h5\0");
    }
    else if (job == 2)
    {
        strcat(fileName, "_locations.h5\0");
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Fortran interface for getting model dimensions
 *
 * @param[in] fileID  HDF5 file ID handle
 * 
 * @param[out] nx     number of x grid points
 * @param[out] ny     number of y grid points
 * @param[out] nz     number of z grid points
 * @param[out] ierr   0 indicates success
 *
 * @author Ben Baker
 *
 * @copyright Apache 2
 *
 */
void eikonal_h5io_getModelDimensionsF(const long *inFileID,
                                      int *nx, int *ny, int *nz,
                                      int *ierr)
{
    const char *fcnm = "eikonal_h5io_getModelDimensionsF\0";
    hid_t fileID = (hid_t) *inFileID;
    *ierr = eikonal_h5io_getModelDimensions(fileID, nx, ny, nz);
    if (*ierr != 0)
    {
        printf("%s: Error getting dimensions\n", fcnm);
        *ierr = 1;
    }
    return;
}
//============================================================================//
/*!
 * @brief Gets the model dimensions
 *
 * @param[in] fileID  HDF5 file ID handle
 * 
 * @param[out] nx     number of x grid points
 * @param[out] ny     number of y grid points
 * @param[out] nz     number of z grid points
 *
 * @result 0 indicates success
 *
 * @author Ben Baker
 *
 * @copyright Apache 2
 *
 */
int eikonal_h5io_getModelDimensions(const hid_t fileID,
                                    int *nx, int *ny, int *nz)
{
    const char *fcnm = "eikonal_h5io_getModelDimensions\0";
    char dataSetName[512];
    hid_t dataSetID, dataSpace;
    hsize_t dims[3];
    herr_t status;
    int ierr, rank;
    //------------------------------------------------------------------------//
    *nx = 0;
    *ny = 0;
    *nz = 0;
    memset(dataSetName, 0, sizeof(dataSetName));
    strcpy(dataSetName, "/Model/xlocs\0");
    if (H5Lexists(fileID, dataSetName, H5P_DEFAULT) != 1)
    {
        printf("%s: Error dataset %s doesn't exist\n", fcnm, dataSetName);
        return -1;
    }
    dataSetID = H5Dopen2(fileID, dataSetName, H5P_DEFAULT);
    dataSpace = H5Dget_space(dataSetID);
    rank = H5Sget_simple_extent_ndims(dataSpace);
    if (rank > 3)
    {
        printf("%s: Invalid rank %d\n", fcnm, rank);
        return -1;
    } 
    status = H5Sget_simple_extent_dims(dataSpace, dims, NULL);
    if (status < 0)
    {
        printf("%s: Error getting dimensions\n", fcnm);
        ierr = 1;
    }
    else
    {
        // Unstructured mesh 
        if (rank == 1)
        {
            *nx = dims[0];
            *ny = dims[0];
            *nz = dims[0];
        }
        // Structured mesh
        else
        {
            *nx = dims[0];
            *ny = dims[1];
            *nz = dims[2];
        }
        ierr = 0;
    }
    status = H5Sclose(dataSpace);
    status = H5Dclose(dataSetID);
    return ierr;
}
//============================================================================//
void eikonal_h5io_setTravelTimeName(const int model, const int station,
                                    const bool isP, char dataSetName[512])
{
    memset(dataSetName, 0, 512*sizeof(char));
    if (isP)
    {
        sprintf(dataSetName,
                "/TravelTimeTables/Model_%d/Station_%d/PTravelTimes",
                model, station);
    }
    else
    {
        sprintf(dataSetName,
                "/TravelTimeTables/Model_%d/Station_%d/STravelTimes",
                model, station);
    }
    return;
}
//============================================================================//
void eikonal_h5io_setLocationName(const int model, const int event,
                                  char dataSetName[512])
{
    memset(dataSetName, 0, 512*sizeof(char));
    sprintf(dataSetName,
            "/logJPDFs/Event_%d/Model_%d/logJPDF", event, model);
    return;
}
//============================================================================//
int eikonal_h5io_finalize(const MPI_Comm comm, hid_t *tttFileID)
{
    const char *fcnm = "eiknoal_h5io_finalize\0";
    herr_t status;
    status = 0;
    status = H5Fclose(*tttFileID);
    if (status != 0)
    {
        printf("%s: Failed closing travel time table file\n", fcnm);
        return -1;
    }
    return 0;
} 
//============================================================================//
/*!
 * @brief Initialize the HDF5 location file with structure
 *
 *        /ModelGeometry
 *          /xlocs
 *          /ylocs
 *          /zlocs
 *        /logJPDFs
 *          /Event_1
 *            /priorLocationModel
 *            /Model_1
 *            .
 *            .
 *            .
 *            /Model_nmodels
 *          /Event_2
 *          .
 *          .
 *          .
 *          /Event_nevents
 *
 * @author Ben Baker
 *
 * @copyright Apache 2
 *
 */
int eikonal_h5io_initLocations(
    const MPI_Comm comm,
    const char *dirnm, const char *projnm,
    const int ix0, const int iy0, const int iz0,
    const int nx, const int ny, const int nz,
    const int nxLoc, const int nyLoc, const int nzLoc,
    const int nmodels, const int nevents,
    const double x0, const double y0, const double z0,
    const double dx, const double dy, const double dz,
    hid_t *locFileID)
{
    const char *fcnm = "eikonal_h5io_initLocations\0";
    MPI_Info info = MPI_INFO_NULL;
    char h5name[PATH_MAX], dataSetName[512];
    hid_t dataSetID, dataSpace, groupID, memSpace, plistID;
    herr_t status;
    size_t blockSize;
    float *data;
    int event, i, ierr, k, model, myid, nxMax, nyMax, nzMax;
    const int rank = 3;
    hsize_t block[3], dimsLocal[3];
    hsize_t count[3] = {1, 1, 1};
    hsize_t dimsGlobal[3] = {nx, ny, nz};
    hsize_t offset[3] = {ix0, iy0, iz0};
    hsize_t stride[3] = {1, 1, 1};
    const bool lsaveRAM = false; // Not supported
    //------------------------------------------------------------------------//
    //  
    // Set the scratch file name
    ierr = 0;
    MPI_Comm_rank(comm, &myid);
    status = 0;
    ierr = eikonal_h5io_setFileName(LOCATION_FILE, dirnm, projnm, h5name); 
    if (ierr != 0)
    {   
        printf("%s: Error setting filename\n", fcnm);
        return -1; 
    }   
    // Get the chunk sizes 
    MPI_Allreduce(&nzLoc, &nzMax, 1, MPI_INTEGER, MPI_MAX, comm);
    MPI_Allreduce(&nyLoc, &nyMax, 1, MPI_INTEGER, MPI_MAX, comm);
    MPI_Allreduce(&nxLoc, &nxMax, 1, MPI_INTEGER, MPI_MAX, comm);
    dimsLocal[0] = nxMax;
    dimsLocal[1] = nyMax;
    dimsLocal[2] = nzMax;
    block[0] = nxMax;
    block[1] = nyMax;
    block[2] = nzMax;
    // Set the properties (RAM or disk + MPI)
    plistID = H5Pcreate(H5P_FILE_ACCESS);
    status = H5Pset_fapl_mpio(plistID, comm, info);
    // If `saving' in RAM decide if we want to write the disk when done 
    if (lsaveRAM)
    {   
        blockSize = (size_t) (nxMax*nyMax*nzMax*4*2 + 4*nxMax*nyMax*nzMax);
        status = H5Pset_fapl_core(plistID, blockSize, true);
        if (status < 0)
        {
            printf("%s: Failed to set RAM memory file open\n", fcnm);
        }
    }   
    *locFileID = H5Fcreate(h5name, H5F_ACC_TRUNC, H5P_DEFAULT, plistID);
    status = H5Pclose(plistID);
    // Make the model group
    ierr = eikonal_h5io_makeModelGroup(comm, *locFileID,
                                       ix0, iy0, iz0,
                                       nx, ny, nz,
                                       nxLoc, nyLoc, nzLoc,
                                       nxMax, nyMax, nzMax,
                                       dx, dy, dz,
                                       x0, y0, z0);
    if (ierr != 0)
    {
        printf("%s: Error making model group\n", fcnm);
        return -1;
    }
    // Make the uniform prior location model
    data = (float *)calloc((size_t) (nzMax*nyMax*nxMax), sizeof(float));
    for (k=0; k<nxMax*nyMax*nzMax; k++)
    {
        data[k] = 1.0;
    }
    memset(dataSetName, 0, sizeof(dataSetName));
    strcpy(dataSetName, "/Model/priorLocationModel");
    memSpace  = H5Screate_simple(rank, dimsLocal,  NULL);
    dataSpace = H5Screate_simple(rank, dimsGlobal, NULL);
    dataSetID = H5Dcreate(*locFileID, dataSetName,
                          H5T_NATIVE_FLOAT,
                          dataSpace, H5P_DEFAULT,
                          H5P_DEFAULT, H5P_DEFAULT);
    status = H5Sclose(dataSpace);
    // Open the chunked dataset
    dataSpace = H5Dget_space(dataSetID);
    // Select hyperslab in file
    status = H5Sselect_hyperslab(dataSpace, H5S_SELECT_SET, offset,
                                 stride, count, block);
    if (status < 0)
    {
        printf("%s: Error selecting hyperslab!\n", fcnm);
        return -1;
    }
    // Create a property list for collective dataset write
    plistID = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plistID, H5FD_MPIO_COLLECTIVE);
    // Write the data
    status = H5Dwrite(dataSetID, H5T_NATIVE_FLOAT, memSpace, dataSpace,
                      plistID, data);
    if (status < 0)
    {
        printf("%s: Error writing prior!\n", fcnm);
        return -1;
    }
    // Free space
    H5Pclose(plistID);
    H5Sclose(dataSpace);
    H5Dclose(dataSetID);
    H5Sclose(memSpace);
    // Make the location groups
    groupID = H5Gcreate2(*locFileID, "/logJPDFs\0", H5P_DEFAULT,
                         H5P_DEFAULT, H5P_DEFAULT);
    status = H5Gclose(groupID);
    MPI_Barrier(comm);
    for (k=0; k<nevents; k++)
    {
        event = k + 1;
        memset(dataSetName, 0, sizeof(dataSetName));
        sprintf(dataSetName, "/logJPDFs/Event_%d", event);
        groupID = H5Gcreate2(*locFileID, dataSetName, H5P_DEFAULT,
                             H5P_DEFAULT, H5P_DEFAULT);
        status = H5Gclose(groupID);
        if (status < 0)
        {
            printf("%s: Error making model group\n", fcnm);
            return -1;
        }
        // Loop on distinct models
        for (i=0; i<nmodels; i++)
        {
            model = i + 1;
            memset(dataSetName, 0, sizeof(dataSetName));
            sprintf(dataSetName, "/logJPDFs/Event_%d/Model_%d", event, model);
            groupID = H5Gcreate2(*locFileID, dataSetName, H5P_DEFAULT,
                                 H5P_DEFAULT, H5P_DEFAULT);
            status = H5Gclose(groupID);
            if (status < 0)
            {
                printf("%s: Error making model group\n", fcnm);
                return -1; 
            }
            // Create the dataset
            dataSpace = H5Screate_simple(rank, dimsGlobal, NULL);
            eikonal_h5io_setLocationName(model, event, dataSetName);
            dataSetID = H5Dcreate(*locFileID, dataSetName,
                                  H5T_NATIVE_FLOAT,
                                  dataSpace, H5P_DEFAULT,
                                  H5P_DEFAULT, H5P_DEFAULT);
            status = H5Dclose(dataSetID); 
            status = H5Sclose(dataSpace);
        }
    }
    MPI_Barrier(comm);
    // Write the bogus jPDFs
    memset(data, 0, (size_t) (nxMax*nyMax*nzMax)*sizeof(float));
    for (k=0; k<nevents; k++)
    {
        event = k + 1;
        for (i=0; i<nmodels; i++)
        {
            model = i + 1;
            ierr = eikonal_h5io_writeLocationLogJPDF(comm, *locFileID,
                                                     model, event,
                                                     ix0, iy0, iz0,
                                                     nxLoc, nyLoc, nzLoc,
                                                     data);
            if (ierr != 0)
            {
                printf("%s: Error writing null jpdfs\n", fcnm);
                return -1; 
            }
        } // Loop on stations
    } // Loop on models
    free(data);
    MPI_Barrier(comm);
    return 0;
}
//============================================================================//
int eikonal_h5io_makeModelGroup(
    const MPI_Comm comm, const hid_t fileID,
    const int ix0, const int iy0, const int iz0,
    const int nxGlob, const int nyGlob, const int nzGlob,
    const int nxLoc, const int nyLoc, const int nzLoc,
    const int nxMax, const int nyMax, const int nzMax,
    const double dx, const double dy, const double dz,
    const double x0, const double y0, const double z0)
{
    const char *fcnm = "eikonal_h5io_makeModelGroup\0";
    char dataSetName[512];
    hsize_t stride[3] = {1, 1, 1};
    hsize_t block[3] = {nxMax, nyMax, nzMax}; //{nzMax, nyMax, nxMax};
    hsize_t count[3] = {1, 1, 1};
    hsize_t dimsLocal[3] = {nxMax, nyMax, nzMax}; //{nzMax, nyMax, nxMax};
    hsize_t dimsGlobal[3] = {nxGlob, nyGlob, nzGlob}; //{nzGlob, nyGlob, nxGlob};
    hsize_t offset[3] = {ix0, iy0, iz0}; //{iz0, iy0, ix0};
    hid_t dataSetID, dataSpace, groupID, memSpace, plistID;
    herr_t status;
    int i, indx, ivar, j, k;
    float *data;
    const int rank = 3;
    //------------------------------------------------------------------------//
    // Make the model group
    groupID = H5Gcreate2(fileID, "/Model\0", H5P_DEFAULT,
                         H5P_DEFAULT, H5P_DEFAULT);
    status = H5Gclose(groupID);
    MPI_Barrier(comm);
    for (ivar=0; ivar<3; ivar++)
    {
        // Create dataspace for dataset
        dataSpace = H5Screate_simple(rank, dimsGlobal, NULL);
        memSpace  = H5Screate_simple(rank, dimsLocal,  NULL);
        // Create the chunked dataset
        //H5Pset_chunk(plistID, rank, dimsLocal);
        data = (float *)calloc((size_t) (nzMax*nyMax*nxMax), sizeof(float));
        memset(dataSetName, 0, sizeof(dataSetName));
        if (ivar == 0)
        {
            strcpy(dataSetName, "/Model/xlocs\0");
            for (k=0; k<nzLoc; k++)
            {
                for (j=0; j<nyLoc; j++)
                {
                    for (i=0; i<nxLoc; i++)
                    {
                        indx = k*nxMax*nyMax + j*nxMax + i;
                        data[indx] = (float) (x0 + (double) (i + ix0)*dx);
                    }
                }
            }
        }
        else if (ivar == 1)
        {
            strcpy(dataSetName, "/Model/ylocs\0");
            for (k=0; k<nzLoc; k++)
            {
                for (j=0; j<nyLoc; j++)
                {
                    for (i=0; i<nxLoc; i++)
                    {
                        indx = k*nxMax*nyMax + j*nxMax + i;
                        data[indx] = (float) (y0 + (double) (j + iy0)*dy);
                    }
                }
            }
        }
        else
        {
            strcpy(dataSetName, "/Model/zlocs\0");
            for (k=0; k<nzLoc; k++)
            {
                for (j=0; j<nyLoc; j++)
                {
                    for (i=0; i<nxLoc; i++)
                    {
                        indx = k*nxMax*nyMax + j*nxMax + i;
                        data[indx] = (float) (z0 + (double) (k + iz0)*dz);
                    }
                }
            }
        }
        dataSetID = H5Dcreate(fileID, dataSetName,
                              H5T_NATIVE_FLOAT,
                              dataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Sclose(dataSpace);
        // Open the chunked dataset
        dataSpace = H5Dget_space(dataSetID);
        // Select hyperslab in file
        status = H5Sselect_hyperslab(dataSpace, H5S_SELECT_SET, offset,
                                     stride, count, block);
        if (status < 0)
        {
            printf("%s: Error selecting hyperslab!\n", fcnm);
            return -1;
        }
        // Create a property list for collective dataset write
        plistID = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plistID, H5FD_MPIO_COLLECTIVE);
        // Write the data
        status = H5Dwrite(dataSetID, H5T_NATIVE_FLOAT, memSpace, dataSpace,
                          plistID, data);
        if (status < 0)
        {
            printf("%s: Error writing data!\n", fcnm);
            return -1;
        }
        // Free space
        H5Pclose(plistID);
        H5Sclose(dataSpace);
        H5Dclose(dataSetID);
        H5Sclose(memSpace);
        free(data);
    }
    MPI_Barrier(comm); 
    return 0;
}
//============================================================================//
/*!
 * @brief Initializes the HDF5 file file directory structure such that:
 *
 *        /ModelGeometry
 *          /xlocs
 *          /ylocs
 *          /zlocs
 *        /TravelTimeTables
 *          /Model_1
 *            /Station_1
 *              /PTravelTimes
 *              /STravelTimes
 *            /Station_2
 *            .
 *            .
 *            .
 *            /Station_nstations
 *          /Model_2
 *          .
 *          .
 *          /Model_nModels
 *
 */
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
                             hid_t *tttFileID)
{
    const char *fcnm = "eikonal_h5io_initTTables\0";
    MPI_Info info = MPI_INFO_NULL;
    char h5name[PATH_MAX], dataSetName[512], groupName[512];
    hid_t dataSetID, dataSpace, groupID, plistID;
    herr_t status;
    size_t blockSize;
    float *data;
    int i, ierr, iphase, k, model, myid, nxMax, nyMax, nzMax, station;
    const int rank = 3;
    hsize_t dimsGlobal[3] = {nx, ny, nz}; //{nz, ny, nx};
    const bool lsaveRAM = false; // Not supported
    //------------------------------------------------------------------------//
    //
    // Set the scratch file name
    MPI_Comm_rank(comm, &myid);
    status = 0;
    ierr = eikonal_h5io_setFileName(TRAVELTIME_FILE, dirnm, projnm, h5name); 
    if (ierr != 0)
    {
        printf("%s: Error setting filename\n", fcnm);
        return -1;
    }
    // Get the chunk sizes 
    MPI_Allreduce(&nzLoc, &nzMax, 1, MPI_INTEGER, MPI_MAX, comm);
    MPI_Allreduce(&nyLoc, &nyMax, 1, MPI_INTEGER, MPI_MAX, comm);
    MPI_Allreduce(&nxLoc, &nxMax, 1, MPI_INTEGER, MPI_MAX, comm);
    // Set the properties (RAM or disk + MPI)
    plistID = H5Pcreate(H5P_FILE_ACCESS);
    status = H5Pset_fapl_mpio(plistID, comm, info);
    // If `saving' in RAM decide if we want to write the disk when done 
    if (lsaveRAM)
    {
        blockSize = (size_t) (nxMax*nyMax*nzMax*4*2 + 4*nxMax*nyMax*nzMax);
        status = H5Pset_fapl_core(plistID, blockSize, lsaveScratch);
        if (status < 0)
        {
            printf("%s: Failed to set RAM memory file open\n", fcnm);
        }
    }
    *tttFileID = H5Fcreate(h5name, H5F_ACC_TRUNC, H5P_DEFAULT, plistID);
    status = H5Pclose(plistID);
    // Make the model group
    ierr = eikonal_h5io_makeModelGroup(comm, *tttFileID,
                                       ix0, iy0, iz0,
                                       nx, ny, nz,
                                       nxLoc, nyLoc, nzLoc,
                                       nxMax, nyMax, nzMax,
                                       dx, dy, dz,
                                       x0, y0, z0);
    if (ierr != 0)
    {
        printf("%s: Error making model group\n", fcnm);
        return -1;
    }
    // Make the traveltime table groups
    groupID = H5Gcreate2(*tttFileID, "/TravelTimeTables\0", H5P_DEFAULT,
                         H5P_DEFAULT, H5P_DEFAULT);
    status = H5Gclose(groupID);
    MPI_Barrier(comm);
    // Make the traveltime table for the model and then the stations.  
    // Since this is a pure initialization phase we will set the space
    // and update throughout the program execution
    for (i=0; i<nmodels; i++)
    {
        model = i + 1;
        memset(groupName, 0, sizeof(groupName));
        sprintf(groupName, "/TravelTimeTables/Model_%d", model);
        groupID = H5Gcreate2(*tttFileID, groupName, H5P_DEFAULT,
                             H5P_DEFAULT, H5P_DEFAULT);
        status = H5Gclose(groupID);
        if (status < 0)
        {
            printf("%s: Failed to create group %s\n", fcnm, groupName);
            return -1; 
        }
        for (k=0; k<nstations; k++)
        {
            station = k + 1;
            memset(groupName, 0, sizeof(groupName));
            sprintf(groupName, "/TravelTimeTables/Model_%d/Station_%d",
                    model, station);
            groupID = H5Gcreate2(*tttFileID, groupName, H5P_DEFAULT,
                                 H5P_DEFAULT, H5P_DEFAULT);
            status = H5Gclose(groupID);
            if (status < 0)
            {
                printf("%s: Failed to create group %s\n", fcnm, groupName);
                return -1;
            }
            // Create dataspace for dataset
            for (iphase=1; iphase<=2; iphase++)
            {
                dataSpace = H5Screate_simple(rank, dimsGlobal, NULL);
                if (iphase == 1)
                {
                    eikonal_h5io_setTravelTimeName(model, station,
                                                   true, dataSetName);
                }
                else
                {
                    eikonal_h5io_setTravelTimeName(model, station,
                                                   false, dataSetName);
                }
                dataSetID = H5Dcreate(*tttFileID, dataSetName,
                                      H5T_NATIVE_FLOAT,
                                      dataSpace, H5P_DEFAULT,
                                      H5P_DEFAULT, H5P_DEFAULT);
                status = H5Dclose(dataSetID); 
                status = H5Sclose(dataSpace);
            } // Loop on phases
        } // loop on stations
    } // loop on models
    MPI_Barrier(comm);
    // Write the null traveltimes
    data = (float *)calloc((size_t) (nzMax*nyMax*nxMax), sizeof(float));
    for (i=0; i<nmodels; i++)
    {
        model = i + 1;
        for (k=0; k<nstations; k++)
        {
            station = k + 1;
            for (iphase=1; iphase<=2; iphase++)
            {
                ierr = eikonal_h5io_writeTravelTimes(comm, *tttFileID,
                                                     station, model,
                                                     iphase,
                                                     ix0, iy0, iz0,
                                                     nxLoc, nyLoc, nzLoc,
                                                     data);
                if (ierr != 0)
                {
                    printf("%s: Error writing null ttimes\n", fcnm);
                    return -1;
                }
            } // Loop on phases
        } // Loop on stations
    } // Loop on models
    free(data);
    MPI_Barrier(comm);
    return 0;
}
//============================================================================//
int eikonal_h5io_writeLocationLogJPDF(
    const MPI_Comm comm, const hid_t locFileID,
    const int model, const int event,
    const int ix0, const int iy0, const int iz0,
    const int nxLoc, const int nyLoc, const int nzLoc,
    const float *__restrict__ logJPDF)
{
    const char *fcnm = "eikonal_h5io_writeLocationLogJPDF\0";
    char dataSetName[512];
    hid_t dataSetID, dataSpace, memSpace, plistID;
    hsize_t dimsLocal[3], block[3];
    hsize_t count[3] = {1, 1, 1};
    hsize_t stride[3] = {1, 1, 1};
    hsize_t offset[3] = {ix0, iy0, iz0}; //{iz0, iy0, ix0};
    herr_t status;
    int i, indx, j, jndx, k, nxMax, nyMax, nzMax;
    const int rank = 3;
    float *data;
    //------------------------------------------------------------------------//
    //
    // Make the dataset name and check that it exists
    status = 0;
    eikonal_h5io_setLocationName(model, event, dataSetName);
    if (H5Lexists(locFileID, dataSetName, H5P_DEFAULT) != 1)
    {
        printf("%s: Error dataset %s doesn't exist\n", fcnm, dataSetName);
        return -1;
    }
    // Set the chunk sizes 
    MPI_Allreduce(&nzLoc, &nzMax, 1, MPI_INTEGER, MPI_MAX, comm);
    MPI_Allreduce(&nyLoc, &nyMax, 1, MPI_INTEGER, MPI_MAX, comm);
    MPI_Allreduce(&nxLoc, &nxMax, 1, MPI_INTEGER, MPI_MAX, comm);
    dimsLocal[0] = nxMax; //nzMax;
    dimsLocal[1] = nyMax; //nyMax;
    dimsLocal[2] = nzMax; //nxMax;
    block[0] = nxMax; //nzMax;
    block[1] = nyMax; //nyMax;
    block[2] = nzMax; //nxMax;
    // Set the data
    data = (float *)calloc((size_t) (nxMax*nyMax*nzMax), sizeof(float));
    for (k=0; k<nzLoc; k++)
    {
        for (j=0; j<nyLoc; j++)
        {
            for (i=0; i<nxLoc; i++)
            {
                indx = k*nxMax*nyMax + j*nxMax + i;
                jndx = k*nxLoc*nyLoc + j*nxLoc + i;
                data[indx] = logJPDF[jndx];
            }
        }
    }
    // Open the memory spaces
    memSpace  = H5Screate_simple(rank, dimsLocal,  NULL);
    // Open the dataset
    dataSetID = H5Dopen2(locFileID, dataSetName, H5P_DEFAULT);
    // Open the chunked dataset
    dataSpace = H5Dget_space(dataSetID);
    // Select hyperslab in file
    status = H5Sselect_hyperslab(dataSpace, H5S_SELECT_SET, offset,
                                 stride, count, block);
    if (status < 0)
    {
        printf("%s: Error selecting hyperslab!\n", fcnm);
        return -1;
    }
    // Create a property list for collective dataset write
    plistID = H5Pcreate(H5P_DATASET_XFER);
    status = H5Pset_dxpl_mpio(plistID, H5FD_MPIO_COLLECTIVE);
    // Write the data
    status = H5Dwrite(dataSetID, H5T_NATIVE_FLOAT, memSpace,
                      dataSpace, plistID, data);
    if (status < 0)
    {
        printf("%s: Error writing dataset: %s!\n", fcnm, dataSetName);
        return -1;
    }
    // Free space
    free(data);
    status = H5Pclose(plistID);
    if (status < 0)
    {
        printf("%s: Error closing plistID\n", fcnm);
        return -1;
    }
    status = H5Sclose(dataSpace);
    if (status < 0)
    {
        printf("%s: Error closing dataSpace\n", fcnm);
        return -1;
    }
    status = H5Dclose(dataSetID);
    if (status < 0)
    {
        printf("%s: Error closing dataSetID\n", fcnm);
        return -1;
    }
    status = H5Sclose(memSpace);
    if (status < 0)
    {
        printf("%s: Error closing memSpace\n", fcnm);
        return -1;
    }
    MPI_Barrier(comm);
    return 0;
}
                             
//============================================================================//
/*!
 * @brief Writes the traveltimes computed on this communicator to the
 *        tttfileID HDF5 file.  If a dataset exists then it will be
 *        overwritten.
 *
 * @param[in] comm       MPI communicator
 * @param[in] tttfileID  HDF5 traveltime file handle.  The file must be
 *                       opened in read/write mode. 
 * @param[in] station    staion number
 * @param[in] model      earth model number
 * @param[in] iphase     if 1 then this is a P phase.
 *                       if 2 then this is an S phase.
 * @param[in] ix0        processes first global x index in grid
 * @param[in] iy0        processes first global y index in grid
 * @param[in] iz0        processes first global z index in grid
 * @param[in] nxLoc      number of local x grid points.
 * @param[in] nyLoc      number of local y grid points.
 * @param[in] nzLoc      number of local z grid points.
 * @param[in] ttimes     traveltimes to write.  the (ix,iy,iz)'th grid point
 *                       is given by (iz-1)*nxLoc*nyLoc + (iy-1)*nxLoc + ix
 *                       [nzLoc*nyLoc*nxLoc]
 *
 * @result 0 indicates success
 *
 * @author Ben Baker
 *
 * @copyright Apache 2
 *
 */
int eikonal_h5io_writeTravelTimes(
    const MPI_Comm comm, const hid_t tttFileID,
    const int station, const int model, const int iphase,
    const int ix0, const int iy0, const int iz0,
    const int nxLoc, const int nyLoc, const int nzLoc,
    const float *__restrict__ ttimes)
{
    const char *fcnm = "eikonal_h5io_writeTravelTimes\0";
    char dataSetName[512];
    hid_t dataSetID, dataSpace, memSpace, plistID;
    hsize_t dimsLocal[3], block[3];
    hsize_t count[3] = {1, 1, 1};
    hsize_t stride[3] = {1, 1, 1};
    hsize_t offset[3] = {ix0, iy0, iz0}; //{iz0, iy0, ix0};
    herr_t status;
    int i, indx, j, jndx, k, nxMax, nyMax, nzMax;
    bool isP;
    const int rank = 3;
    float *data;
    //------------------------------------------------------------------------//
    //
    // Make the dataset name and check that it exists
    status = 0;
    isP = true;
    if (iphase == 2){isP = false;}
    eikonal_h5io_setTravelTimeName(model, station, isP, dataSetName);
    if (H5Lexists(tttFileID, dataSetName, H5P_DEFAULT) != 1)
    {
        printf("%s: Error dataset %s doesn't exist\n", fcnm, dataSetName);
        return -1;
    }
    // Set the chunk sizes 
    MPI_Allreduce(&nzLoc, &nzMax, 1, MPI_INTEGER, MPI_MAX, comm);
    MPI_Allreduce(&nyLoc, &nyMax, 1, MPI_INTEGER, MPI_MAX, comm);
    MPI_Allreduce(&nxLoc, &nxMax, 1, MPI_INTEGER, MPI_MAX, comm);
    dimsLocal[0] = nxMax; //nzMax;
    dimsLocal[1] = nyMax; //nyMax;
    dimsLocal[2] = nzMax; //nxMax;
    block[0] = nxMax; //nzMax;
    block[1] = nyMax; //nyMax;
    block[2] = nzMax; //nxMax;
    // Set the data
    data = (float *)calloc((size_t) (nxMax*nyMax*nzMax), sizeof(float));
    for (k=0; k<nzLoc; k++)
    {   
        for (j=0; j<nyLoc; j++)
        {   
            for (i=0; i<nxLoc; i++)
            {   
                indx = k*nxMax*nyMax + j*nxMax + i;
                jndx = k*nxLoc*nyLoc + j*nxLoc + i;
                data[indx] = ttimes[jndx];
            }   
        }   
    }
    // Open the memory spaces
    memSpace  = H5Screate_simple(rank, dimsLocal,  NULL);
    // Open the dataset
    dataSetID = H5Dopen2(tttFileID, dataSetName, H5P_DEFAULT);
    // Open the chunked dataset
    dataSpace = H5Dget_space(dataSetID);
    // Select hyperslab in file
    status = H5Sselect_hyperslab(dataSpace, H5S_SELECT_SET, offset,
                                 stride, count, block);
    if (status < 0)
    {
        printf("%s: Error selecting hyperslab!\n", fcnm);
        return -1;
    }
    // Create a property list for collective dataset write
    plistID = H5Pcreate(H5P_DATASET_XFER);
    status = H5Pset_dxpl_mpio(plistID, H5FD_MPIO_COLLECTIVE);
    // Write the data
    status = H5Dwrite(dataSetID, H5T_NATIVE_FLOAT, memSpace,
                      dataSpace, plistID, data);
    if (status < 0)
    {
        printf("%s: Error writing dataset: %s!\n", fcnm, dataSetName);
        return -1;
    }
    // Free space
    free(data);
    status = H5Pclose(plistID);
    if (status < 0)
    {
        printf("%s: Error closing plistID\n", fcnm);
        return -1;
    }
    status = H5Sclose(dataSpace);
    if (status < 0)
    {
        printf("%s: Error closing dataSpace\n", fcnm);
        return -1;
    } 
    status = H5Dclose(dataSetID);
    if (status < 0)
    {
        printf("%s: Error closing dataSetID\n", fcnm);
        return -1;
    }
    status = H5Sclose(memSpace);
    if (status < 0)
    {
        printf("%s: Error closing memSpace\n", fcnm);
        return -1;
    }
    return 0;
}
//============================================================================//
void eikonal_h5io_readModelF(const int *comm, const long *inFileID,
                             const int *ix0, const int *iy0, const int *iz0, 
                             const int *nxLoc, const int *nyLoc,
                             const int *nzLoc,
                             float *__restrict__ xlocs,
                             float *__restrict__ ylocs,
                             float *__restrict__ zlocs,
                             int *ierr)
{
    const char *fcnm = "eikonal_h5io_readModelF\0";
    MPI_Comm mpiComm = MPI_Comm_f2c(*comm);
    hid_t fileID = (hid_t) *inFileID;
    int ix, iy, iz;
    *ierr = 0;
    ix = *ix0 - 1;
    iy = *iy0 - 1;
    iz = *iz0 - 1;
    *ierr = eikonal_h5io_readModel(mpiComm, fileID,
                                   ix, iy, iz,
                                   *nxLoc, *nyLoc, *nzLoc,
                                   xlocs, ylocs, zlocs);
    if (*ierr != 0)
    {
        printf("%s: Error reading model\n", fcnm);
        *ierr = 1;
    }
    return;
}
//============================================================================//
/*!
 * @brief Processes each read a chunk of a traveltime table for the given
 *        station and model 
 *
 * @param[in] comm       MPI communicator
 * @param[in] tttfileID  HDF5 traveltime file handle.  The file must be
 *                       opened in read/write mode. 
 * @param[in] station    staion number
 * @param[in] model      earth model number
 * @param[in] iphase     if 1 then this is a P phase.
 *                       if 2 then this is an S phase.
 * @param[in] ix0        processes first global x index in grid
 * @param[in] iy0        processes first global y index in grid
 * @param[in] iz0        processes first global z index in grid
 * @param[in] nxLoc      number of local x grid points.
 * @param[in] nyLoc      number of local y grid points.
 * @param[in] nzLoc      number of local z grid points.
 *
 * @param[out] xlocs     x locations defining grid [nxLoc*nyLoc*nzLoc]
 * @param[out] ylocs     y locations defining grid [nxLoc*nyLoc*nzLoc]
 * @param[out] zlocs     z locations defining grid [nxLoc*nyLoc*nzLoc]
 *
 * @result 0 indicates success
 * 
 * @author Ben Baker
 *
 * @copyright Apache 2
 *
 */
int eikonal_h5io_readModel(const MPI_Comm comm,
                           const hid_t fileID,
                           const int ix0, const int iy0, const int iz0,
                           const int nxLoc, const int nyLoc, const int nzLoc,
                           float *__restrict__ xlocs,
                           float *__restrict__ ylocs,
                           float *__restrict__ zlocs)
{
    const char *fcnm = "eikonal_h5io_readModel\0";
    char dataSetName[512];
    hid_t dataSetID, dataSpace, memSpace, plistID;
    hsize_t dimsLocal[3], block[3];
    hsize_t count[3] = {1, 1, 1};
    hsize_t stride[3] = {1, 1, 1};
    hsize_t offset[3] = {ix0, iy0, iz0};
    herr_t status;
    int i, indx, ivar, j, jndx, k, nxMax, nyMax, nzMax;
    const int rank = 3;
    float *data;
    //------------------------------------------------------------------------//
    //  
    // Set the chunk sizes 
    status = 0;
    MPI_Allreduce(&nzLoc, &nzMax, 1, MPI_INTEGER, MPI_MAX, comm);
    MPI_Allreduce(&nyLoc, &nyMax, 1, MPI_INTEGER, MPI_MAX, comm);
    MPI_Allreduce(&nxLoc, &nxMax, 1, MPI_INTEGER, MPI_MAX, comm);
    dimsLocal[0] = nxMax;
    dimsLocal[1] = nyMax;
    dimsLocal[2] = nzMax;
    block[0] = nxMax;
    block[1] = nyMax;
    block[2] = nzMax;
    // Set space
    data = (float *)calloc((size_t) (nxMax*nyMax*nzMax), sizeof(float));
    for (ivar=0; ivar<3; ivar++)
    {
        memset(dataSetName, 0, sizeof(dataSetName));
        if (ivar == 0)
        {
            strcpy(dataSetName, "/Model/xlocs\0");
        }
        else if (ivar == 1)
        {
            strcpy(dataSetName, "/Model/ylocs\0");
        }
        else if (ivar == 2)
        {
            strcpy(dataSetName, "/Model/zlocs\0");
        }
        if (H5Lexists(fileID, dataSetName, H5P_DEFAULT) != 1)
        {
            printf("%s: Error dataset %s doesn't exist\n", fcnm, dataSetName);
            return -1;
        }
        // Open the memory spaces
        memSpace  = H5Screate_simple(rank, dimsLocal,  NULL);
        // Open the dataset
        dataSetID = H5Dopen2(fileID, dataSetName, H5P_DEFAULT);
        // Open the chunked dataset
        dataSpace = H5Dget_space(dataSetID);
        // Define the block to read
        status = H5Sselect_hyperslab(dataSpace, H5S_SELECT_SET, offset,
                                     stride, count, block);
        if (status < 0)
        {
            printf("%s: Error setting slab dimensions\n", fcnm);
            return -1;
        }
        // Create a property list for collective dataset write
        plistID = H5Pcreate(H5P_DATASET_XFER);
        status = H5Pset_dxpl_mpio(plistID, H5FD_MPIO_COLLECTIVE);
        if (status < 0)
        {
            printf("%s: Error modifying plistID\n", fcnm);
            return -1;
        }
        status = H5Dread(dataSetID, H5T_NATIVE_FLOAT, memSpace,
                         dataSpace, plistID, data);
        if (status < 0)
        {
            printf("%s: Error reading dataset %s\n", fcnm, dataSetName);
            return -1;
        }
        // Copy result
        if (ivar == 0)
        {
            for (k=0; k<nzLoc; k++)
            {
                for (j=0; j<nyLoc; j++)
                {
                    for (i=0; i<nxLoc; i++)
                    {
                        indx = k*nxMax*nyMax + j*nxMax + i;
                        jndx = k*nxLoc*nyLoc + j*nxLoc + i;
                        xlocs[jndx] = data[indx];
                    }
                }
            }
        }
        else if (ivar == 1)
        {
            for (k=0; k<nzLoc; k++)
            {
                for (j=0; j<nyLoc; j++)
                {   
                    for (i=0; i<nxLoc; i++)
                    {   
                        indx = k*nxMax*nyMax + j*nxMax + i;
                        jndx = k*nxLoc*nyLoc + j*nxLoc + i;
                        ylocs[jndx] = data[indx];
                    }
                }
            }
        }
        else if (ivar == 2)
        {
            for (k=0; k<nzLoc; k++)
            {
                for (j=0; j<nyLoc; j++)
                {   
                    for (i=0; i<nxLoc; i++)
                    {   
                        indx = k*nxMax*nyMax + j*nxMax + i;
                        jndx = k*nxLoc*nyLoc + j*nxLoc + i;
                        zlocs[jndx] = data[indx];
                    }
                }
            }
        }
        // Free space
        status = H5Pclose(plistID);
        status = H5Sclose(dataSpace);
        status = H5Dclose(dataSetID);
        status = H5Sclose(memSpace);
        MPI_Barrier(comm);
    }
    free(data);
    return 0;
}

//============================================================================//
/*!
 * @brief Fortran interface to H5 reader
 */
void eikonal_h5io_readTraveltimesF(const int *comm,
                                   const long *tttFileID,
                                   const int *station, const int *model,
                                   const int *iphase,
                                   const int *ix0f, const int *iy0f,
                                   const int *iz0f,
                                   const int *nxLoc, const int *nyLoc,
                                   const int *nzLoc,
                                   float *ttimes, int *ierr)
{
    const char *fcnm = "eikonal_h5io_readTraveltimesF\0";
    MPI_Comm mpiComm = MPI_Comm_f2c(*comm);
    hid_t fileID = (hid_t) *tttFileID;
    int ix0, iy0, iz0;
    ix0 = *ix0f - 1;
    iy0 = *iy0f - 1;
    iz0 = *iz0f - 1;
    *ierr = 0;
    *ierr = eikonal_h5io_readTravelTimes(mpiComm,
                                         fileID,
                                         *station, *model, *iphase,
                                         ix0, iy0, iz0,
                                         *nxLoc, *nyLoc, *nzLoc,
                                         ttimes);
    if (*ierr != 0)
    {
        printf("%s: Error calling eikonal_h5io_readTravelTimes\n", fcnm);
        *ierr = 1;
        memset(ttimes, 0.0, (size_t) (*nxLoc**nyLoc**nzLoc)*sizeof(float));
    }
    return;
}
                                  
//============================================================================//
/*!
 * @brief Processes each read a chunk of a traveltime table for the given
 *        station and model 
 *
 * @param[in] comm       MPI communicator
 * @param[in] tttfileID  HDF5 traveltime file handle.  The file must be
 *                       opened in read/write mode. 
 * @param[in] station    staion number
 * @param[in] model      earth model number
 * @param[in] iphase     if 1 then this is a P phase.
 *                       if 2 then this is an S phase.
 * @param[in] ix0        processes first global x index in grid
 * @param[in] iy0        processes first global y index in grid
 * @param[in] iz0        processes first global z index in grid
 * @param[in] nxLoc      number of local x grid points.
 * @param[in] nyLoc      number of local y grid points.
 * @param[in] nzLoc      number of local z grid points.
 *
 * @param[out] ttimes    traveltimes read from disk.  the (ix,iy,iz)'th grid
 *                       point is given by
 *                         (iz-1)*nxLoc*nyLoc + (iy-1)*nxLoc + ix
 *                       [nzLoc*nyLoc*nxLoc] 
 *
 * @result 0 indicates success
 * 
 * @author Ben Baker
 *
 * @copyright Apache 2
 *
 */
int eikonal_h5io_readTravelTimes(const MPI_Comm comm,
                                 const hid_t tttFileID,
                                 const int station, const int model,
                                 const int iphase,
                                 const int ix0, const int iy0, const int iz0,
                                 const int nxLoc, const int nyLoc,
                                 const int nzLoc,
                                 float *__restrict__ ttimes)
{
    const char *fcnm = "eikonal_h5io_readTravelTimes\0";
    char dataSetName[512];
    hid_t dataSetID, dataSpace, memSpace, plistID;
    hsize_t dimsLocal[3], block[3];
    hsize_t count[3] = {1, 1, 1};
    hsize_t stride[3] = {1, 1, 1};
    hsize_t offset[3] = {ix0, iy0, iz0}; //{iz0, iy0, ix0};
    herr_t status;
    int i, indx, j, jndx, k, nxMax, nyMax, nzMax;
    bool isP;
    const int rank = 3;
    float *data;
    //------------------------------------------------------------------------//
    //  
    // Make the dataset name and check that it exists
    status = 0;
    isP = true;
    if (iphase == 2){isP = false;}
    eikonal_h5io_setTravelTimeName(model, station, isP, dataSetName);
    if (H5Lexists(tttFileID, dataSetName, H5P_DEFAULT) != 1)
    {   
        printf("%s: Error dataset %s doesn't exist\n", fcnm, dataSetName);
        return -1; 
    }   
    // Set the chunk sizes 
    MPI_Allreduce(&nzLoc, &nzMax, 1, MPI_INTEGER, MPI_MAX, comm);
    MPI_Allreduce(&nyLoc, &nyMax, 1, MPI_INTEGER, MPI_MAX, comm);
    MPI_Allreduce(&nxLoc, &nxMax, 1, MPI_INTEGER, MPI_MAX, comm);
    dimsLocal[0] = nxMax; //nzMax;
    dimsLocal[1] = nyMax; //nyMax;
    dimsLocal[2] = nzMax; //nxMax;
    block[0] = nxMax; //nzMax;
    block[1] = nyMax; //nyMax;
    block[2] = nzMax; //nxMax;
    // Set space
    data = (float *)calloc((size_t) (nxMax*nyMax*nzMax), sizeof(float));
    // Open the memory spaces
    memSpace  = H5Screate_simple(rank, dimsLocal,  NULL);
    // Open the dataset
    dataSetID = H5Dopen2(tttFileID, dataSetName, H5P_DEFAULT);
    // Open the chunked dataset
    dataSpace = H5Dget_space(dataSetID);
    // Define the block to read
    status = H5Sselect_hyperslab(dataSpace, H5S_SELECT_SET, offset,
                                 stride, count, block);
    if (status < 0)
    {
        printf("%s: Error setting slab dimensions\n", fcnm);
        return -1;
    }        
    // Create a property list for collective dataset write
    plistID = H5Pcreate(H5P_DATASET_XFER);
    status = H5Pset_dxpl_mpio(plistID, H5FD_MPIO_COLLECTIVE);
    if (status < 0)
    {
        printf("%s: Error modifying plist\n", fcnm);
        return -1;
    }
    status = H5Dread(dataSetID, H5T_NATIVE_FLOAT, memSpace,
                     dataSpace, plistID, data); 
    if (status < 0)
    {
        printf("%s: Error reading dataset %s\n", fcnm, dataSetName);
        return -1;
    }
    // Copy result
    for (k=0; k<nzLoc; k++)
    {
        for (j=0; j<nyLoc; j++)
        {
            for (i=0; i<nxLoc; i++)
            {
                indx = k*nxMax*nyMax + j*nxMax + i;
                jndx = k*nxLoc*nyLoc + j*nxLoc + i;
                ttimes[jndx] = data[indx];
            }
        }
    }
    // Free space
    free(data);
    status = H5Pclose(plistID);
    if (status < 0)
    {
        printf("%s: Error closing plistID\n", fcnm);
        return -1;
    }
    status = H5Sclose(dataSpace);
    if (status < 0)
    {
        printf("%s: Error closing dataSpace\n", fcnm);
        return -1;
    }
    status = H5Dclose(dataSetID);
    if (status < 0)
    {
        printf("%s: Error closing dataSetID\n", fcnm);
        return -1;
    }
    status = H5Sclose(memSpace);
    if (status < 0)
    {
        printf("%s: Error closing memSpace\n", fcnm);
        return -1;
    }
    MPI_Barrier(comm);
    return 0;
}

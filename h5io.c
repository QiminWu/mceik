#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include <mpi.h>
#include "h5io.h"
#include "os.h"

int eikonal_h5io_setFileName(const char *dirnm, const char *projnm,
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
    strcat(fileName, ".h5\0");
    return 0;
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
 * @brief Initializes the HDF5 file file directory structure such that:
 *
 *        /ModelGeometry
 *          /Domain
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
                            hid_t *tttFileID)
{
    const char *fcnm = "eikonal_h5io_initialize\0";
    MPI_Info info = MPI_INFO_NULL;
    char h5name[PATH_MAX], dataSetName[512], groupName[512];
    hid_t dataSetID, dataSpace, groupID, memSpace, plistID;
    herr_t status;
    size_t blockSize;
    float *data;
    int i, ierr, indx, iphase, ivar, j, k, myid, nxMax, nyMax, nzMax;
    const int rank = 3;
    hsize_t dimsLocal[3], block[3];
    hsize_t count[3] = {1, 1, 1};
    hsize_t dimsGlobal[3] = {nz, ny, nx};
    hsize_t stride[3] = {1, 1, 1};
    hsize_t offset[3] = {iz0, iy0, ix0};
    const bool lsaveRAM = false; // Not supported
    //------------------------------------------------------------------------//
    //
    // Set the scratch file name
    MPI_Comm_rank(comm, &myid);
    status = 0;
    ierr = eikonal_h5io_setFileName(dirnm, projnm, h5name); 
    if (ierr != 0)
    {
        printf("%s: Error setting filename\n", fcnm);
        return -1;
    }
    // Get the chunk sizes 
    MPI_Allreduce(&nzLoc, &nzMax, 1, MPI_INTEGER, MPI_MAX, comm);
    MPI_Allreduce(&nyLoc, &nyMax, 1, MPI_INTEGER, MPI_MAX, comm);
    MPI_Allreduce(&nxLoc, &nxMax, 1, MPI_INTEGER, MPI_MAX, comm);
    dimsLocal[0] = nzMax;
    dimsLocal[1] = nyMax;
    dimsLocal[2] = nxMax;
    block[0] = nzMax;
    block[1] = nyMax;
    block[2] = nxMax;
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
    groupID = H5Gcreate2(*tttFileID, "/Model\0", H5P_DEFAULT,
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
        dataSetID = H5Dcreate(*tttFileID, dataSetName,
                              H5T_NATIVE_FLOAT,
                              dataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Sclose(dataSpace);
        // Open the chunked dataset
        dataSpace = H5Dget_space(dataSetID);
        //globalDataSpace = H5Dget_space(dataSetID);
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
    } // Loop on variables
    MPI_Barrier(comm);
    // Make the traveltime table groups
    groupID = H5Gcreate2(*tttFileID, "/TravelTimeTables\0", H5P_DEFAULT,
                         H5P_DEFAULT, H5P_DEFAULT);
    status = H5Gclose(groupID);
    MPI_Barrier(comm);
    // Make the traveltime table for the model and then the stations.  
    // Since this is a pure initialization phase we will set the space
    // and update throughout the program execution
    data = (float *)calloc((size_t) (nzMax*nyMax*nxMax), sizeof(float)); 
    for (i=0; i<nmodels; i++)
    {
        memset(groupName, 0, sizeof(groupName));
        sprintf(groupName, "/TravelTimeTables/Model_%d", i+1);
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
            memset(groupName, 0, sizeof(groupName));
            sprintf(groupName, "/TravelTimeTables/Model_%d/Station_%d",
                    i+1, k+1);
            groupID = H5Gcreate2(*tttFileID, groupName, H5P_DEFAULT,
                                 H5P_DEFAULT, H5P_DEFAULT);
            status = H5Gclose(groupID);
            if (status < 0)
            {
                printf("%s: Failed to create group %s\n", fcnm, groupName);
                return -1;
            }
            // Create dataspace for dataset
            for (iphase=0; iphase<2; iphase++)
            {
                dataSpace = H5Screate_simple(rank, dimsGlobal, NULL);
                memSpace  = H5Screate_simple(rank, dimsLocal,  NULL);
                memset(dataSetName, 0, sizeof(dataSetName));
                strcpy(dataSetName, groupName);
                if (iphase == 1)
                {
                    strcat(dataSetName, "/PTravelTimeTable\0");
                }
                else
                {
                    strcat(dataSetName, "/STravelTimeTable\0");
                }
                dataSetID = H5Dcreate(*tttFileID, dataSetName,
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
                status = H5Dwrite(dataSetID, H5T_NATIVE_FLOAT, memSpace,
                                  dataSpace, plistID, data);
                if (status < 0)
                {
                    printf("%s: Error writing dataset: %s!\n",
                           fcnm, dataSetName);
                    return -1;
                }
                // Free space
                H5Pclose(plistID);
                H5Sclose(dataSpace);
                H5Dclose(dataSetID);
                H5Sclose(memSpace);
            }
        } // loop on stations
    } // loop on models
    free(data);
    MPI_Barrier(comm);
    return 0;
}
//============================================================================//
/*!
 * @param[in] ix0       processes first global x index in grid
 * @param[in] iy0       processes first global y index in grid
 * @param[in] iz0       processes first global z index in grid
 * @param[in] nx        number of global x grid points
 * @param[in] ny        number of global y grid points
 * @param[in] nz        number of global z grid points
 *
 */
int eikonal_h5io_writeTravelTimes(const MPI_Comm comm,
                             const int station, const int model,
                             const bool isP,
                             const hid_t tttFileID,
                             const int ix0, const int iy0, const int iz0,
                             //const int nx, const int ny, const int nz,
                             const int nxLoc, const int nyLoc, const int nzLoc,
                             const float *__restrict__ ttimes)
{
    const char *fcnm = "eikonal_h5io_writeTravelTimes\0";
    char dataSetName[512];
    hsize_t block[3], count[3], dimsLocal[3], offset[3], stride[3];
    hid_t dataSetID, globalDataSpace, localDataSpace, plistID;
    herr_t status;
    const int rank = 3; // Data is 3D
    // Make the dataset name and check that it exists
    status = 0;
    eikonal_h5io_setTravelTimeName(model, station, isP, dataSetName);
    if (H5Lexists(tttFileID, dataSetName, H5P_DEFAULT) != 1)
    {
        printf("%s: Error dataset %s doesn't exist\n", fcnm, dataSetName);
        return -1;
    }
    // Set the global size and local dataspaces 
    //dimsGlobal[0] = nz;
    //dimsGlobal[1] = ny;
    //dimsGlobal[2] = nx;
    dimsLocal[0] = nzLoc;
    dimsLocal[1] = nyLoc;
    dimsLocal[2] = nxLoc;
    //globalDataSpace = H5Screate_simple(rank, dimsGlobal, NULL);
    localDataSpace  = H5Screate_simple(rank, dimsLocal,  NULL); 
    // Open the chunked dataset 
/*
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    status = H5Pset_chunk(plistID, rank, dimsLocal);
    globalDataSpace = H5Dget_space(dataSetID);
    dataSetID = H5Dcreate(fileID, cDataName, H5T_NATIVE_FLOAT,
                          globalDataSpace, H5P_DEFAULT, plistID,
                          H5P_DEFAULT); 
    H5Pclose(plistID);
    H5Sclose(globalDataSpace); 
*/
    // Open the chunked dataset
    dataSetID = H5Dopen(tttFileID, dataSetName, H5P_DEFAULT);
    globalDataSpace = H5Dget_space(dataSetID);
    // Set the local positions - each process defines dataset in memory
    // and writes it to the hyperslab in the file.  count and block are
    // interchanged so that each process writes 1 continuous block in x,
    // y, and z
    count[0] = 1;
    count[1] = 1;
    count[2] = 1;
    stride[0] = 1; // Unit stride in x
    stride[1] = 1; // Unit stride in y
    stride[2] = 1; // Unit stride in z
    block[0] = dimsLocal[0]; // local x grid points 
    block[1] = dimsLocal[1]; // local y grid points
    block[2] = dimsLocal[2]; // local z grid points
    offset[0] = iz0; // global x start location
    offset[1] = iy0; // global y start location
    offset[2] = ix0; // global z start location
    // Select hyperslab in the file
    status = H5Sselect_hyperslab(globalDataSpace, H5S_SELECT_SET,
                                 offset, stride, count, block);
    // Create property list for collective dataset write
    plistID = H5Pcreate(H5P_DATASET_XFER);
    status = H5Dwrite(dataSetID, H5T_NATIVE_FLOAT, localDataSpace,
                      globalDataSpace, plistID, ttimes);  
    // Close it up
    status += H5Pclose(plistID);
    status += H5Sclose(globalDataSpace);
    status += H5Sclose(localDataSpace);
    status += H5Dclose(dataSetID);
    if (status < 0)
    {
        printf("%s: Failed to write dataset %s\n", fcnm, dataSetName);
        return -1;
    } 
    return 0;
}
//============================================================================//
int eikonal_h5io_readTravelTimes(const hid_t tttFileID,
                                 const int nx, const int ny, const int nz,
                                 const int ix0, const int iy0, const int iz0,
                                 const int nxLoc, const int nyLoc, const int nzLoc,
                                 float *__restrict__ ttimes)
{
    return 0;
}
//============================================================================//

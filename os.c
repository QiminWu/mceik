#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <limits.h>
#include <sys/stat.h>
#include <errno.h>
#include "os.h"

/*! 
 * @brief Tests if dirnm is a directory 
 * 
 * @param[in] dirnm    name of directory to test 
 *
 * @result true -> dirnm is an existing directory 
 *         false -> dirnm is not a directory
 * 
 * @author Ben Baker, ISTI
 *
 */
bool os_path_isdir(const char *dirnm)
{
    struct stat s;
    int err;
    if (dirnm == NULL){return false;}
    if (strlen(dirnm) == 0){return false;}
    err = stat(dirnm, &s);
    // Doesn't exist
    if (err == -1) 
    {   
        if (ENOENT == errno)
        {
            return false;
        }
        // Exists
        else
        {
            return true;
        }
    }   
    // Exists
    else
    {   
        // Test it is a directory
        if (S_ISDIR(s.st_mode))
        {
            return true;
        }
        else
        {
            return false;
        }
    }   
}
//============================================================================//
/*!
 * @brief Recursive directory creation function
 *
 * @param[in] path    directory tree to make
 *
 * @result 0 indicates success
 *
 * @author Ben Baker, ISTI
 *
 */
int os_makedirs(const char *path)
{
    const char *fcnm = "os_makedirs\0";
    char *where, *dname, *work, directory[PATH_MAX];
    int ierr, indx, lenos;
    const char find[] = "/\0";
    //------------------------------------------------------------------------//
    //
    // Errors 
    if (path == NULL){return -1;}
    lenos = strlen(path);
    if (lenos == 0){return -1;}
    if (lenos > PATH_MAX - 2)
    {
        printf("%s: Error directory %s is too long\n", fcnm, path);
        return -1;
    }
    // Already exists
    if (os_path_isdir(path)){return 0;}
    // Initialize
    work = (char *)calloc(lenos+2, sizeof(char));
    strcpy(work, path);
    memset(directory, 0, sizeof(directory));
    dname = work;
    dname[lenos] = '/'; // Try to catch the final case
    where = strpbrk(dname, find);
    while ((where != NULL))
    {
        indx = where - dname;
        lenos = strlen(directory);
        strncat(directory, dname, indx);
        // If directory doesn't exist then make it
        if (!os_path_isdir(directory))
        {
            ierr = os_mkdir(directory);
            if (ierr != 0)
            {
                printf("%s: Error making subdirectory: %s\n",
                       fcnm, directory);
                printf("%s: Error making directory: %s\n",
                        fcnm, path);
                free(dname);
                return -1;
            }
        } // End check on if directory exists
        // Add directory delimiter
        strcat(directory, "/\0");
        dname = dname + indx + 1;
        where = strpbrk(dname, find);
    } // End while
    // Add name of final subdirectory and make it
    strcat(directory, dname);
    if (!os_path_isdir(directory))
    {   
         ierr = os_mkdir(directory);
         if (ierr != 0)
         {
             printf("%s: Error making directory: %s\n", fcnm, directory);
             return -1; 
         }
    }   
    // Free space
    dname = NULL;
    free(work);
    return ierr;
}
//============================================================================//
/*! 
 * Makes a directory named dirnm with full permissions 
 *
 * @param[in] dirnm    name of directory to make
 *
 * @result 0 if success
 *
 * @author Ben Baker, ISTI
 *
 */
int os_mkdir(const char *dirnm)
{
    const char *fcnm = "os_mkdir\0";
    int ierr;
    if (dirnm == NULL){return -1;}
    if (strlen(dirnm) == 0){return -1;}
    ierr = mkdir(dirnm, 0777);
    if (ierr != 0)
    {   
        printf("%s: Error making directory: %s\n", fcnm, dirnm);
        return -1; 
    }   
    return 0;
}
//============================================================================//
/*! 
 * @brief Tests if filenm is a file
 * 
 * @param[in] filenm    name of file to test 
 * 
 * @result  true  -> filenm is an existing file
 *          false -> filenm is not a file
 *
 * @author Ben Baker, ISTI
 *
 */
bool os_path_isfile(const char *filenm)
{
    struct stat info;
    if (filenm == NULL){return false;}
    if (strlen(filenm) == 0){return false;}
    // Doesn't exist
    if (stat(filenm, &info) ==-1)
    {   
        return false;
    }   
    // Exists -> check it is a file
    else
    {   
        if (S_ISREG(info.st_mode))
        {
            return true;
        }
        else
        {
            return false;
        }
    }   
}

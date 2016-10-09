#include <stdbool.h>
#ifndef _os_os_h__
#define _os_os_h__ 1
#ifdef __cplusplus
extern "C" 
{
#endif
/* Determines if pathnm exists */
bool os_path_exists(const char *pathnm);
/* Determines if dirnm is an existing directory */
bool os_path_isdir(const char *dirnm);
/* Determines if filenm is an existing file */
bool os_path_isfile(const char *filenm);
/* Recursively make directories */
int os_makedirs(const char *path);
/* Makes a directory */
int os_mkdir(const char *dirnm);

#ifdef __cplusplus
}
#endif
#endif /* __OS_OS_H__ */

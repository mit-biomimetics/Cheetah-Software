#ifndef LIB_HANDLER_H
#define LIB_HANDLER_H

#include "glob_opts.h"

#ifdef IS_WINDOWS
#include <windows.h>
typedef HINSTANCE soHandle_t;
#else
#include <unistd.h>
#include <dlfcn.h>
typedef void *soHandle_t;
#endif

#ifdef IS_WINDOWS
#define SHAREDLIBEXT "dll"
#elif defined IS_MAC
#define SHAREDLIBEXT "dylib"
#else
#define SHAREDLIBEXT "so"
#endif


// NB: A shared library should be put into the shared library search path:
//  Linux:    LD_LIBRARY_PATH
//  Mac OSX:  DYLD_LIBRARY_PATH
//  Windows:  PATH


/** Loads a dynamically linked library.
 * @param libname The name of the library to load.
 * @return Shared library handle, or OSQP_NULL if failure.
 */
soHandle_t lh_load_lib(const char* libname);


/** Unloads a shared library.
 * @param libhandle Handle of shared library to unload.
 * @return Zero on success, nonzero on failure.
 */
c_int lh_unload_lib(soHandle_t libhandle);


#endif /*LIB_HANDLER_H*/

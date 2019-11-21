#ifndef PARDISOLOADER_H
#define PARDISOLOADER_H

#ifdef __cplusplus
extern "C" {
#endif



/**
 * Tries to load a shared library with Pardiso.
 * Return a failure if the library cannot be loaded or not all Pardiso symbols are found.
 * @param libname The name under which the Pardiso lib can be found, or OSQP_NULL to use a default name (mkl_rt.SHAREDLIBEXT).
 * @return Zero on success, nonzero on failure.
 */
c_int lh_load_pardiso(const char* libname);

/**
 * Unloads the loaded Pardiso shared library.
 * @return Zero on success, nonzero on failure.
 */
c_int lh_unload_pardiso();


#ifdef __cplusplus
}
#endif

#endif /*PARADISOLOADER_H*/

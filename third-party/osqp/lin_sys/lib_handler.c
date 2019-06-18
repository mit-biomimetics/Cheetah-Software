#include "lib_handler.h"
#include <ctype.h> // Needed for tolower functions

#include "constants.h"

soHandle_t lh_load_lib(const char *libName) {
    soHandle_t h = OSQP_NULL;

    if (!libName) {
        #ifdef PRINTING
        c_eprint("no library name given");
        #endif
        return OSQP_NULL;
    }

#ifdef IS_WINDOWS
    h = LoadLibrary (libName);
    if (!h) {
        #ifdef PRINTING
        c_eprint("Windows error while loading dynamic library %s, error = %d",
                libName, (int)GetLastError());
        #endif
    }
#else
    h = dlopen (libName, RTLD_LAZY);
    if (!h) {
        #ifdef PRINTING
        c_eprint("Error while loading dynamic library %s: %s", libName, dlerror());
        #endif
    }
#endif

    return h;
} /* lh_load_lib */


c_int lh_unload_lib (soHandle_t h) {
    c_int rc = 1;

#ifdef IS_WINDOWS
    rc = FreeLibrary (h);
    rc = ! rc;
#else
    rc = dlclose (h);
#endif

    return rc;
} /* LSL_unLoadLib */


#ifdef IS_WINDOWS
typedef FARPROC symtype;
#else
typedef void* symtype;
#endif
/** Loads a symbol from a dynamically linked library.
 * This function is not defined in the header to allow a workaround for the problem that dlsym returns an object instead of a function pointer.
 * However, Windows also needs special care.
 *
 * The method does six attempts to load the symbol. Next to its given name, it also tries variations of lower case and upper case form and with an extra underscore.
 * @param h Handle of dynamically linked library.
 * @param symName Name of the symbol to load.
 * @return A pointer to the symbol, or OSQP_NULL if not found.
 */
symtype lh_load_sym (soHandle_t h, const char *symName) {
    symtype s;
    const char *from;
    char *to;
    const char *tripSym;
    char* err;
    char lcbuf[257];
    char ucbuf[257];
    char ocbuf[257];
    size_t symLen;
    int trip;

    s = OSQP_NULL;
    err = OSQP_NULL;

    /* search in this order:
     *  1. original
     *  2. lower_
     *  3. upper_
     *  4. original_
     *  5. lower
     *  6. upper
     */

    symLen = 0;
    for (trip = 1;  trip <= 6;  trip++) {
        switch (trip) {
        case 1:                             /* original */
            tripSym = symName;
            break;
        case 2:                             /* lower_ */
            for (from = symName, to = lcbuf;  *from;  from++, to++) {
                *to = tolower(*from);
            }
            symLen = from - symName;
            *to++ = '_';
            *to = '\0';
            tripSym = lcbuf;
            break;
        case 3:                             /* upper_ */
            for (from = symName, to = ucbuf;  *from;  from++, to++) {
                *to = toupper(*from);
            }
            *to++ = '_';
            *to = '\0';
            tripSym = ucbuf;
            break;
        case 4:                             /* original_ */
            memcpy (ocbuf, symName, symLen);
            ocbuf[symLen] = '_';
            ocbuf[symLen+1] = '\0';
            tripSym = ocbuf;
            break;
        case 5:                             /* lower */
            lcbuf[symLen] = '\0';
            tripSym = lcbuf;
            break;
        case 6:                             /* upper */
            ucbuf[symLen] = '\0';
            tripSym = ucbuf;
            break;
        default:
            tripSym = symName;
        } /* end switch */
#ifdef IS_WINDOWS
        s = GetProcAddress (h, tripSym);
        if (s) {
            return s;
        } else {
            #ifdef PRINTING
            c_eprint("Cannot find symbol %s in dynamic library, error = %d",
                    symName, (int)GetLastError());
            #endif
        }
#else
        s = dlsym (h, tripSym);
        err = dlerror();  /* we have only one chance; a successive call to dlerror() returns OSQP_NULL */
        if (err) {
            #ifdef PRINTING
            c_eprint("Cannot find symbol %s in dynamic library, error = %s",
                    symName, err);
            #endif
        } else {
            return s;
        }
#endif
    } /* end loop over symbol name variations */

    return OSQP_NULL;
} /* lh_load_sym */

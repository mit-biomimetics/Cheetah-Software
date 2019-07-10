#ifndef _common_types
#define _common_types
#include <stdint.h>
#include <eigen3/Eigen/Dense>

//adding this line adds print statements and sanity checks 
//that are too slow for realtime use.
//#define K_DEBUG

typedef double dbl;
typedef float flt;

//floating point type used whenever possible
typedef float fpt;

//floating point type used when interfacing with MATLAB
typedef double mfp;
typedef int mint;

typedef uint64_t u64;
typedef uint32_t u32;
typedef uint16_t u16;
typedef uint8_t  u8;
typedef int8_t   s8;
typedef int16_t  s16;
typedef int32_t  s32;
typedef int64_t  s64;

#endif

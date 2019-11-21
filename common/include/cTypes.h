/*! @file cTypes.h
 *  @brief Common types that can also be included in C code
 *
 *  This file contains types which are shared between C and C++ code.  The
 * low-level drivers (rt folder) are all in C and everything else is C++.
 * Because this file is included in both C and C++, it can't contain C++ type
 * alias ("using"), namespaces, or templates.
 */

#ifndef PROJECT_CTYPES_H
#define PROJECT_CTYPES_H

#include <stddef.h>  // for size_t
#include <stdint.h>

// short version of the stdint default integer types
typedef uint64_t u64;
typedef uint32_t u32;
typedef uint16_t u16;
typedef uint8_t u8;
typedef int8_t s8;
typedef int16_t s16;
typedef int32_t s32;
typedef int64_t s64;

#endif  // PROJECT_CTYPES_H

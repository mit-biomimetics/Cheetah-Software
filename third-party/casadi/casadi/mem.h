/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


#ifndef CASADI_MEM_H
#define CASADI_MEM_H

#ifdef __cplusplus
extern "C" {
#endif

#ifndef CASADI_STATIC
#include <stdlib.h>
#endif /* CASADI_STATIC */
#include <assert.h>

/* Floating point type */
#ifndef casadi_real
#define casadi_real double
#endif

/* Integer type */
#ifndef casadi_int
#define casadi_int long long int
#endif

/* Function types corresponding to entry points in CasADi's C API */
typedef void (*casadi_signal_t)(void);
typedef casadi_int (*casadi_getint_t)(void);
typedef const casadi_int* (*casadi_sparsity_t)(casadi_int i);
typedef const char* (*casadi_name_t)(casadi_int i);
typedef int (*casadi_work_t)(casadi_int* sz_arg, casadi_int* sz_res,
                             casadi_int* sz_iw, casadi_int* sz_w);
typedef int (*casadi_eval_t)(const casadi_real** arg, casadi_real** res,
                             casadi_int* iw, casadi_real* w, void* mem);

/* Structure to hold meta information about an input or output */
typedef struct {
  const char* name;
  casadi_int nrow;
  casadi_int ncol;
  casadi_int nnz;
  casadi_int numel;
  const casadi_int* colind;
  const casadi_int* row;
} casadi_io;

/* Decompress a sparsity pattern */
inline void casadi_decompress(const int* sp, int* nrow, int* ncol,
                              int* nnz, int* numel,
                              const int** colind, const int** row) {
  if (sp==0) {
    /* Scalar sparsity pattern if sp is null */
    static const casadi_int scalar_colind[2] = {0, 1};
    *ncol = *nrow = 1;
    *nnz = *numel = 1;
    *colind = scalar_colind;
    *row = 0;
  } else {
    /* Decompress */
    *nrow = *sp++;
    *ncol = *sp++;
    *colind = sp;
    *nnz = sp[*ncol];
    *numel = *nrow * *ncol;
    *row = *nnz==*numel ? 0 : sp + *ncol + 1;
  }
}

/* All entry points */
typedef struct {
  casadi_signal_t incref;
  casadi_signal_t decref;
  casadi_getint_t n_in;
  casadi_getint_t n_out;
  casadi_name_t name_in;
  casadi_name_t name_out;
  casadi_sparsity_t sparsity_in;
  casadi_sparsity_t sparsity_out;
  casadi_work_t work;
  casadi_eval_t eval;
} casadi_functions;

/* Memory needed for evaluation */
typedef struct {
  /* Function pointers */
  casadi_functions* f;

  /* Work arrays */
  casadi_int sz_arg, sz_res, sz_iw, sz_w;
  const casadi_real** arg;
  casadi_real** res;
  casadi_int* iw;
  casadi_real* w;
  void* mem;

  /* Meta information */
  casadi_int n_in, n_out;
  casadi_io* in;
  casadi_io* out;
} casadi_mem;

/* Initialize */
inline void casadi_init(casadi_mem* mem, casadi_functions* f) {
  int flag;

  /* Check arguments */
  assert(mem!=0);
  assert(f!=0);

  /* Store function pointers */
  mem->f = f;

  /* Increase reference counter */
  if (f->incref) f->incref();

  /* Number of inputs and outputs */
  mem->n_in = f->n_in ? f->n_in() : 1;
  mem->n_out = f->n_out ? f->n_out() : 1;

  /* Work vector sizes */
  mem->sz_arg=mem->n_in;
  mem->sz_res=mem->n_out;
  mem->sz_iw= mem->sz_w = 0;
  if (f->work) {
    flag = f->work(&mem->sz_arg, &mem->sz_res, &mem->sz_iw, &mem->sz_w);
    assert(flag==0);
  }

  /* TODO: Check out a memory object */
  mem->mem = 0;

  /* No io structs allocated */
  mem->in = 0;
  mem->out = 0;

  /* No work vectors allocated */
  mem->arg = 0;
  mem->res = 0;
  mem->iw = 0;
  mem->w = 0;
}

/* Free claimed static memory */
inline void casadi_deinit(casadi_mem* mem) {
  assert(mem!=0);

  /* TODO: Release a memory object */

  /* Decrease reference counter */
  if (mem->f->decref) mem->f->decref();
}

/* Initialize */
inline void casadi_init_arrays(casadi_mem* mem) {
  casadi_int i;
  assert(mem!=0);
  casadi_functions* f = mem->f;

  /* Input meta data */
  for (i=0; i<mem->n_in; ++i) {
    mem->in[i].name = f->name_in ? f->name_in(i) : 0;
    casadi_decompress(f->sparsity_in ? f->sparsity_in(i) : 0,
                      &mem->in[i].nrow, &mem->in[i].ncol,
                      &mem->in[i].nnz, &mem->in[i].numel,
                      &mem->in[i].colind, &mem->in[i].row);
  }

  /* Output meta data */
  for (i=0; i<mem->n_out; ++i) {
    mem->out[i].name = f->name_out ? f->name_out(i) : 0;
    casadi_decompress(f->sparsity_out ? f->sparsity_out(i) : 0,
                      &mem->out[i].nrow, &mem->out[i].ncol,
                      &mem->out[i].nnz, &mem->out[i].numel,
                      &mem->out[i].colind, &mem->out[i].row);
  }
}

/* Allocate dynamic memory */
#ifndef CASADI_STATIC
inline int casadi_alloc_arrays(casadi_mem* mem) {
  /* Allocate io memory */
  mem->in = (casadi_io*)malloc(mem->n_in*sizeof(casadi_io));
  if (mem->n_in!=0 && mem->in==0) return 1;
  mem->out = (casadi_io*)malloc(mem->n_out*sizeof(casadi_io));
  if (mem->n_out!=0 && mem->out==0) return 1;

  /* Allocate work vectors */
  mem->arg = (const casadi_real**)malloc(mem->sz_arg*sizeof(const casadi_real*));
  if (mem->sz_arg!=0 && mem->arg==0) return 1;
  mem->res = (casadi_real**)malloc(mem->sz_res*sizeof(casadi_real*));
  if (mem->sz_res!=0 && mem->res==0) return 1;
  mem->iw = (casadi_int*)malloc(mem->sz_iw*sizeof(casadi_int));
  if (mem->sz_iw!=0 && mem->iw==0) return 1;
  mem->w = (casadi_real*)malloc(mem->sz_w*sizeof(casadi_real));
  if (mem->sz_w!=0 && mem->w==0) return 1;

  return 0;
}
#endif /* CASADI_STATIC */

/* Free dynamic memory */
#ifndef CASADI_STATIC
inline void casadi_free_arrays(casadi_mem* mem) {
  assert(mem!=0);

  /* Free io meta data */
  if (mem->in) free(mem->in);
  if (mem->out) free(mem->out);

  /* Free work vectors */
  if (mem->arg) free(mem->arg);
  if (mem->res) free(mem->res);
  if (mem->iw) free(mem->iw);
  if (mem->w) free(mem->w);
}
#endif /* CASADI_STATIC */

/* Evaluate */
inline int casadi_eval(casadi_mem* mem) {
  assert(mem!=0);
  return mem->f->eval(mem->arg, mem->res, mem->iw, mem->w, mem->mem);
}

/* Create a memory struct with dynamic memory allocation */
#ifndef CASADI_STATIC
inline casadi_mem* casadi_alloc(casadi_functions* f) {
  int flag;

  /* Allocate struct */
  casadi_mem* mem = (casadi_mem*)malloc(sizeof(casadi_mem));
  assert(mem!=0);

  /* Initialize struct */
  casadi_init(mem, f);

  /* Dynamically allocate arrays */
  flag = casadi_alloc_arrays(mem);
  assert(flag==0);

  /* Initialize allocated arrays */
  casadi_init_arrays(mem);

  return mem;
}
#endif /* CASADI_STATIC */

/* Free memory struct */
#ifndef CASADI_STATIC
inline void casadi_free(casadi_mem* mem) {
  assert(mem!=0);

  /* Free dynamically allocated arrays */
  casadi_free_arrays(mem);

  /* Free claimed static memory */
  casadi_deinit(mem);

  /* Free memory structure */
  free(mem);
}
#endif /* CASADI_STATIC */

#ifdef __cplusplus
}
#endif

#endif // CASADI_MEM_H

#include "util.h"

/***************
* Versioning  *
***************/
const char* osqp_version(void) {
  return OSQP_VERSION;
}

/************************************
* Printing Constants to set Layout *
************************************/
#ifdef PRINTING
# define HEADER_LINE_LEN 65
#endif /* ifdef PRINTING */

/**********************
* Utility Functions  *
**********************/
void c_strcpy(char dest[], const char source[]) {
  int i = 0;

  while (1) {
    dest[i] = source[i];

    if (dest[i] == '\0') break;
    i++;
  }
}

#ifdef PRINTING

static void print_line(void) {
  char  the_line[HEADER_LINE_LEN + 1];
  c_int i;

  for (i = 0; i < HEADER_LINE_LEN; ++i) the_line[i] = '-';
  the_line[HEADER_LINE_LEN] = '\0';
  c_print("%s\n", the_line);
}

void print_header(void) {
  // Different indentation required for windows
# ifdef IS_WINDOWS
# ifndef PYTHON
  c_print("iter  ");
# endif /* ifdef PYTHON */
# else  /* ifdef IS_WINDOWS */
  c_print("iter   ");
# endif /* ifdef IS_WINDOWS */

  // Main information
  c_print("objective    pri res    dua res    rho");
# ifdef PROFILING
  c_print("        time");
# endif /* ifdef PROFILING */
  c_print("\n");
}

void print_setup_header(const OSQPWorkspace *work) {
  OSQPData *data;
  OSQPSettings *settings;
  c_int nnz; // Number of nonzeros in the problem

  data     = work->data;
  settings = work->settings;

  // Number of nonzeros
  nnz = data->P->p[data->P->n] + data->A->p[data->A->n];

  print_line();
  c_print("           OSQP v%s  -  Operator Splitting QP Solver\n"
          "              (c) Bartolomeo Stellato,  Goran Banjac\n"
          "        University of Oxford  -  Stanford University 2018\n",
          OSQP_VERSION);
  print_line();

  // Print variables and constraints
  c_print("problem:  ");
  c_print("variables n = %i, constraints m = %i\n          ",
                                    (int)data->n,
          (int)data->m);
  c_print("nnz(P) + nnz(A) = %i\n", (int)nnz);

  // Print Settings
  c_print("settings: ");
  c_print("linear system solver = %s",
          LINSYS_SOLVER_NAME[settings->linsys_solver]);

  if (work->linsys_solver->nthreads != 1) {
    c_print(" (%d threads)", (int)work->linsys_solver->nthreads);
  }
  c_print(",\n          ");

  c_print("eps_abs = %.1e, eps_rel = %.1e,\n          ",
          settings->eps_abs, settings->eps_rel);
  c_print("eps_prim_inf = %.1e, eps_dual_inf = %.1e,\n          ",
          settings->eps_prim_inf, settings->eps_dual_inf);
  c_print("rho = %.2e ", settings->rho);

  if (settings->adaptive_rho) c_print("(adaptive)");
  c_print(",\n          ");
  c_print("sigma = %.2e, alpha = %.2f, ",
          settings->sigma, settings->alpha);
  c_print("max_iter = %i\n", (int)settings->max_iter);

  if (settings->check_termination) c_print(
      "          check_termination: on (interval %i),\n",
      (int)settings->check_termination);
  else c_print("          check_termination: off,\n");

# ifdef PROFILING
  if (settings->time_limit) c_print("          time_limit: %.2e sec,\n",
                                    settings->time_limit);
# endif /* ifdef PROFILING */

  if (settings->scaling) c_print("          scaling: on, ");
  else c_print("          scaling: off, ");

  if (settings->scaled_termination) c_print("scaled_termination: on\n");
  else c_print("scaled_termination: off\n");

  if (settings->warm_start) c_print("          warm start: on, ");
  else c_print("          warm start: off, ");

  if (settings->polish) c_print("polish: on\n");
  else c_print("polish: off\n");
  c_print("\n");
}

void print_summary(OSQPWorkspace *work) {
  OSQPInfo *info;

  info = work->info;

  c_print("%4i",     (int)info->iter);
  c_print(" %12.4e", info->obj_val);
  c_print("  %9.2e", info->pri_res);
  c_print("  %9.2e", info->dua_res);
  c_print("  %9.2e", work->settings->rho);
# ifdef PROFILING

  if (work->first_run) {
    // total time: setup + solve
    c_print("  %9.2es", info->setup_time + info->solve_time);
  } else {
    // total time: update + solve
    c_print("  %9.2es", info->update_time + info->solve_time);
  }
# endif /* ifdef PROFILING */
  c_print("\n");

  work->summary_printed = 1; // Summary has been printed
}

void print_polish(OSQPWorkspace *work) {
  OSQPInfo *info;

  info = work->info;

  c_print("%4s",     "plsh");
  c_print(" %12.4e", info->obj_val);
  c_print("  %9.2e", info->pri_res);
  c_print("  %9.2e", info->dua_res);

  // Different characters for windows/unix
# ifdef IS_WINDOWS
# ifndef PYTHON
  c_print("  ---------");
# endif /* ifdef PYTHON */
# else  /* ifdef IS_WINDOWS */
  c_print("   --------");
# endif /* ifdef IS_WINDOWS */
# ifdef PROFILING
  if (work->first_run) {
    // total time: setup + solve
    c_print("  %9.2es", info->setup_time + info->solve_time +
            info->polish_time);
  } else {
    // total time: update + solve
    c_print("  %9.2es", info->update_time + info->solve_time +
            info->polish_time);
  }
# endif /* ifdef PROFILING */
  c_print("\n");
}

void print_footer(OSQPInfo *info, c_int polish) {
  c_print("\n"); // Add space after iterations

  c_print("status:               %s\n", info->status);

  if (polish && (info->status_val == OSQP_SOLVED)) {
    if (info->status_polish == 1) {
      c_print("solution polish:      successful\n");
    } else if (info->status_polish < 0) {
      c_print("solution polish:      unsuccessful\n");
    }
  }

  c_print("number of iterations: %i\n", (int)info->iter);

  if ((info->status_val == OSQP_SOLVED) ||
      (info->status_val == OSQP_SOLVED_INACCURATE)) {
    c_print("optimal objective:    %.4f\n", info->obj_val);
  }

# ifdef PROFILING
  c_print("run time:             %.2es\n", info->run_time);
# endif /* ifdef PROFILING */

# if EMBEDDED != 1
  c_print("optimal rho estimate: %.2e\n", info->rho_estimate);
# endif /* if EMBEDDED != 1 */
  c_print("\n");
}

#endif /* End #ifdef PRINTING */


#ifndef EMBEDDED

OSQPSettings* copy_settings(OSQPSettings *settings) {
  OSQPSettings *new = c_malloc(sizeof(OSQPSettings));

  // Copy settings
  memcpy(new, settings, sizeof(OSQPSettings));
  return new;
}

#endif // #ifndef EMBEDDED


/*******************
* Timer Functions *
*******************/

#ifdef PROFILING

// Windows
# ifdef IS_WINDOWS

void osqp_tic(OSQPTimer *t)
{
  QueryPerformanceFrequency(&t->freq);
  QueryPerformanceCounter(&t->tic);
}

c_float osqp_toc(OSQPTimer *t)
{
  QueryPerformanceCounter(&t->toc);
  return (t->toc.QuadPart - t->tic.QuadPart) / (c_float)t->freq.QuadPart;
}

// Mac
# elif defined IS_MAC

void osqp_tic(OSQPTimer *t)
{
  /* read current clock cycles */
  t->tic = mach_absolute_time();
}

c_float osqp_toc(OSQPTimer *t)
{
  uint64_t duration; /* elapsed time in clock cycles*/

  t->toc   = mach_absolute_time();
  duration = t->toc - t->tic;

  /*conversion from clock cycles to nanoseconds*/
  mach_timebase_info(&(t->tinfo));
  duration *= t->tinfo.numer;
  duration /= t->tinfo.denom;

  return (c_float)duration / 1e9;
}

// Linux
# else  /* ifdef IS_WINDOWS */

/* read current time */
void osqp_tic(OSQPTimer *t)
{
  clock_gettime(CLOCK_MONOTONIC, &t->tic);
}

/* return time passed since last call to tic on this timer */
c_float osqp_toc(OSQPTimer *t)
{
  struct timespec temp;

  clock_gettime(CLOCK_MONOTONIC, &t->toc);

  if ((t->toc.tv_nsec - t->tic.tv_nsec) < 0) {
    temp.tv_sec  = t->toc.tv_sec - t->tic.tv_sec - 1;
    temp.tv_nsec = 1e9 + t->toc.tv_nsec - t->tic.tv_nsec;
  } else {
    temp.tv_sec  = t->toc.tv_sec - t->tic.tv_sec;
    temp.tv_nsec = t->toc.tv_nsec - t->tic.tv_nsec;
  }
  return (c_float)temp.tv_sec + (c_float)temp.tv_nsec / 1e9;
}

# endif /* ifdef IS_WINDOWS */

#endif // If Profiling end


/* ==================== DEBUG FUNCTIONS ======================= */



// If debug mode enabled
#ifdef DDEBUG

#ifdef PRINTING

void print_csc_matrix(csc *M, const char *name)
{
  c_int j, i, row_start, row_stop;
  c_int k = 0;

  // Print name
  c_print("%s :\n", name);

  for (j = 0; j < M->n; j++) {
    row_start = M->p[j];
    row_stop  = M->p[j + 1];

    if (row_start == row_stop) continue;
    else {
      for (i = row_start; i < row_stop; i++) {
        c_print("\t[%3u,%3u] = %.3g\n", (int)M->i[i], (int)j, M->x[k++]);
      }
    }
  }
}

void dump_csc_matrix(csc *M, const char *file_name) {
  c_int j, i, row_strt, row_stop;
  c_int k = 0;
  FILE *f = fopen(file_name, "w");

  if (f != NULL) {
    for (j = 0; j < M->n; j++) {
      row_strt = M->p[j];
      row_stop = M->p[j + 1];

      if (row_strt == row_stop) continue;
      else {
        for (i = row_strt; i < row_stop; i++) {
          fprintf(f, "%d\t%d\t%20.18e\n",
                  (int)M->i[i] + 1, (int)j + 1, M->x[k++]);
        }
      }
    }
    fprintf(f, "%d\t%d\t%20.18e\n", (int)M->m, (int)M->n, 0.0);
    fclose(f);
    c_print("File %s successfully written.\n", file_name);
  } else {
    c_eprint("Error during writing file %s.\n", file_name);
  }
}

void print_trip_matrix(csc *M, const char *name)
{
  c_int k = 0;

  // Print name
  c_print("%s :\n", name);

  for (k = 0; k < M->nz; k++) {
    c_print("\t[%3u, %3u] = %.3g\n", (int)M->i[k], (int)M->p[k], M->x[k]);
  }
}

void print_dns_matrix(c_float *M, c_int m, c_int n, const char *name)
{
  c_int i, j;

  c_print("%s : \n\t", name);

  for (i = 0; i < m; i++) {   // Cycle over rows
    for (j = 0; j < n; j++) { // Cycle over columns
      if (j < n - 1)
        // c_print("% 14.12e,  ", M[j*m+i]);
        c_print("% .3g,  ", M[j * m + i]);

      else
        // c_print("% 14.12e;  ", M[j*m+i]);
        c_print("% .3g;  ", M[j * m + i]);
    }

    if (i < m - 1) {
      c_print("\n\t");
    }
  }
  c_print("\n");
}

void print_vec(c_float *v, c_int n, const char *name) {
  print_dns_matrix(v, 1, n, name);
}

void dump_vec(c_float *v, c_int len, const char *file_name) {
  c_int i;
  FILE *f = fopen(file_name, "w");

  if (f != NULL) {
    for (i = 0; i < len; i++) {
      fprintf(f, "%20.18e\n", v[i]);
    }
    fclose(f);
    c_print("File %s successfully written.\n", file_name);
  } else {
    c_print("Error during writing file %s.\n", file_name);
  }
}

void print_vec_int(c_int *x, c_int n, const char *name) {
  c_int i;

  c_print("%s = [", name);

  for (i = 0; i < n; i++) {
    c_print(" %i ", (int)x[i]);
  }
  c_print("]\n");
}

#endif // PRINTING

#endif // DEBUG MODE

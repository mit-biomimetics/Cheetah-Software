/* QDLDL TESTER MODULE */

/* THE CODE FOR MINIMAL UNIT TESTING HAS BEEN TAKEN FROM
   http://www.jera.com/techinfo/jtns/jtn002.html */

#define mu_assert(message, test) \
  do { if (!(test)) return message; } while (0)
#define mu_run_test(test)                   \
  do { char *message = test(); tests_run++; \
       if (message) return message; } while (0)
extern int tests_run;

#define QDLDL_TESTS_TOL 1e-4 // Define tests tolerance

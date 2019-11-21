static char* test_singleton()
{
  //A matrix data
  QDLDL_int Ap[]  = {0, 1};
  QDLDL_int Ai[]  = {0};
  QDLDL_float Ax[] = {0.2};
  QDLDL_int An = 1;

  // RHS and solution to Ax = b
  QDLDL_float b[]    = {2};
  QDLDL_float xsol[] = {10.0};

  //x replaces b during solve
  int status = ldl_factor_solve(An,Ap,Ai,Ax,b);

  mu_assert("Factorisation failed", status >= 0);
  mu_assert("Solve accuracy failed", vec_diff_norm(b,xsol,An) < QDLDL_TESTS_TOL);

  return 0;
}

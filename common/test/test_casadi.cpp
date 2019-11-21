#include "gmock/gmock.h"
#include "gtest/gtest.h"

#define CASADI_UNIT_TEST 0

#if (CASADI_UNIT_TEST)

#include <casadi/casadi.hpp>
#include <ctime>
#include <fstream>
#include <iostream>

#include "Dynamics/DynamicsSimulator.h"
#include "Dynamics/FloatingBaseModel.h"
#include "Dynamics/MiniCheetah.h"
#include "Dynamics/Quadruped.h"
#include <Utilities/Utilities_print.h>


using namespace casadi;
using namespace std;

TEST(casadi, rocket_ipopt) {
  cout << "program started" << endl;

  // Dimensions
  int nu = 20;   // Number of control segments
  int nj = 100;  // Number of integration steps per control segment

  // optimization variable
  SX u = SX::sym("u", nu);  // control

  SX s_0 = 0;  // initial position
  SX v_0 = 0;  // initial speed
  SX m_0 = 1;  // initial mass

  SX dt = 10.0 / (nj * nu);  // time step
  SX alpha = 0.05;           // friction
  SX beta = 0.1;             // fuel consumption rate

  // Trajectory
  vector<SX> s_k, v_k, m_k;

  // Integrate over the interval with Euler forward
  SX s = s_0, v = v_0, m = m_0;
  for (int k = 0; k < nu; ++k) {
    for (int j = 0; j < nj; ++j) {
      s += dt * v;
      v += dt / m * (u(k) - alpha * v * v);
      m += -dt * beta * u(k) * u(k);
    }
    s_k.push_back(s);
    v_k.push_back(v);
    m_k.push_back(m);
  }
  SX s_all = vertcat(s_k), v_all = vertcat(v_k), m_all = vertcat(m_k);

  // Objective function
  SX f = dot(u, u);

  // Terminal constraints
  SX g = vertcat(s, v, v_all);

  // Create the NLP
  SXDict nlp = {{"x", u}, {"f", f}, {"g", g}};

  // Allocate an NLP solver and buffers
  Function solver = nlpsol("solver", "ipopt", nlp);

  // Bounds on g
  vector<double> gmin = {10, 0};
  vector<double> gmax = {10, 0};
  gmin.resize(2 + nu, -numeric_limits<double>::infinity());
  gmax.resize(2 + nu, 1.1);

  // Solve the problem
  DMDict arg = {
      {"lbx", -10}, {"ubx", 10}, {"x0", 0.4}, {"lbg", gmin}, {"ubg", gmax}};
  DMDict res = solver(arg);

  // Print the optimal cost
  double cost(res.at("f"));
  cout << "optimal cost: " << cost << endl;

  // Print the optimal solution
  vector<double> uopt(res.at("x"));
  cout << "optimal control: " << uopt << endl;

  // Get the state trajectory
  Function xfcn("xfcn", {u}, {s_all, v_all, m_all});
  vector<double> sopt, vopt, mopt;
  xfcn({uopt}, {&sopt, &vopt, &mopt});
  cout << "position: " << sopt << endl;
  cout << "velocity: " << vopt << endl;
  cout << "mass:     " << mopt << endl;

  bool b_save_matlab_file(false);
  if (b_save_matlab_file) {
    // Create Matlab script to plot the solution
    ofstream file;
    string filename = "rocket_ipopt_results.m";
    file.open(filename.c_str());
    file << "% Results file from " __FILE__ << endl;
    file << "% Generated " __DATE__ " at " __TIME__ << endl;
    file << endl;
    file << "cost = " << cost << ";" << endl;
    file << "u = " << uopt << ";" << endl;

    // Save results to file
    file << "t = linspace(0,10.0," << nu << ");" << endl;
    file << "s = " << sopt << ";" << endl;
    file << "v = " << vopt << ";" << endl;
    file << "m = " << mopt << ";" << endl;

    // Finalize the results file
    file << endl;
    file << "% Plot the results" << endl;
    file << "figure(1);" << endl;
    file << "clf;" << endl << endl;

    file << "subplot(2,2,1);" << endl;
    file << "plot(t,s);" << endl;
    file << "grid on;" << endl;
    file << "xlabel('time [s]');" << endl;
    file << "ylabel('position [m]');" << endl << endl;

    file << "subplot(2,2,2);" << endl;
    file << "plot(t,v);" << endl;
    file << "grid on;" << endl;
    file << "xlabel('time [s]');" << endl;
    file << "ylabel('velocity [m/s]');" << endl << endl;

    file << "subplot(2,2,3);" << endl;
    file << "plot(t,m);" << endl;
    file << "grid on;" << endl;
    file << "xlabel('time [s]');" << endl;
    file << "ylabel('mass [kg]');" << endl << endl;

    file << "subplot(2,2,4);" << endl;
    file << "plot(t,u);" << endl;
    file << "grid on;" << endl;
    file << "xlabel('time [s]');" << endl;
    file << "ylabel('Thrust [kg m/s^2]');" << endl << endl;

    file.close();
    cout << "Results saved to \"" << filename << "\"" << endl;
  }
}

TEST(casadi, part1) {
  /** Test problem
   *
   *    min x0^2 + x1^2
   *    s.t.    x0 + x1 - 10 = 0
   */

  // Optimization variables
  MX x = MX::sym("x", 2);

  // Objective
  MX f = x(0) * x(0) + x(1) * x(1);

  // Constraints
  MX g = x(0) + x(1) - 10;

  // Create an NLP solver instance
  Function solver = nlpsol("solver", "ipopt", {{"x", x}, {"f", f}, {"g", g}});

  // Generate C code for the NLP functions
  solver.generate_dependencies("nlp.c");

  // Just-in-time compilation?
  bool jit = false;
  if (jit) {
    // Create a new NLP solver instance using just-in-time compilation
    solver = nlpsol("solver", "ipopt", "nlp.c");
  } else {
    // Compile the c-code
    int flag = system("gcc -fPIC -shared -O3 nlp.c -o nlp.so");
    casadi_assert(flag == 0, "Compilation failed");

    // Create a new NLP solver instance from the compiled code
    solver = nlpsol("solver", "ipopt", "nlp.so");
  }

  // Bounds and initial guess
  std::map<std::string, DM> arg, res;
  arg["lbx"] = -DM::inf();
  arg["ubx"] = DM::inf();
  arg["lbg"] = 0;
  arg["ubg"] = DM::inf();
  arg["x0"] = 0;

  // Solve the NLP
  res = solver(arg);

  // Print solution
  cout << "-----" << endl;
  cout << "objective at solution = " << res.at("f") << endl;
  cout << "primal solution = " << res.at("x") << endl;
  cout << "dual solution (x) = " << res.at("lam_x") << endl;
  cout << "dual solution (g) = " << res.at("lam_g") << endl;
}

// dx/dt = f(x,u)
MX f(const MX& x, const MX& u) { return vertcat(x(1), u - x(1)); }

TEST(casadi, race_car) {
  // Car race along a track
  // ----------------------
  // An optimal control problem (OCP),
  // solved with direct multiple-shooting.
  //
  // For more information see: http://labs.casadi.org/OCP

  int N = 100;  // number of control intervals

  auto opti = casadi::Opti();  // Optimization problem

  Slice all;
  // ---- decision variables ---------
  auto X = opti.variable(2, N + 1);  // state trajectory
  auto pos = X(0, all);
  auto speed = X(1, all);
  auto U = opti.variable(1, N);  // control trajectory (throttle)
  auto T = opti.variable();      // final time

  // ---- objective          ---------
  opti.minimize(T);  // race in minimal time

  // ---- dynamic constraints --------
  auto dt = T / N;
  for (int k = 0; k < N; ++k) {
    auto k1 = f(X(all, k), U(all, k));
    auto k2 = f(X(all, k) + dt / 2 * k1, U(all, k));
    auto k3 = f(X(all, k) + dt / 2 * k2, U(all, k));
    auto k4 = f(X(all, k) + dt * k3, U(all, k));
    auto x_next = X(all, k) + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
    opti.subject_to(X(all, k + 1) == x_next);  // close the gaps
  }

  // ---- path constraints -----------
  opti.subject_to(speed <= 1 - sin(2 * pi * pos) / 2);  // track speed limit
  opti.subject_to(0 <= U <= 1);                         // control is limited

  // ---- boundary conditions --------
  opti.subject_to(pos(1) == 0);    // start at position 0 ...
  opti.subject_to(speed(1) == 0);  // ... from stand-still
  opti.subject_to(pos(N) == 1);    // finish line at position 1

  // ---- misc. constraints  ----------
  opti.subject_to(T >= 0);  // Time must be positive

  // ---- initial values for solver ---
  opti.set_initial(speed, 1);
  opti.set_initial(T, 1);

  // ---- solve NLP              ------
  opti.solver("ipopt");     // set numerical backend
  auto sol = opti.solve();  // actual solve

  bool b_matlab_file(false);
  if (b_matlab_file) {
    // Create Matlab script to plot the solution
    ofstream file;
    string filename = "race_car_results.m";
    file.open(filename.c_str());
    file << "% Results file from " __FILE__ << endl;
    file << "% Generated " __DATE__ " at " __TIME__ << endl;
    file << endl;

    // Save results to file
    file << "t = linspace(0," << sol.value(T) << "," << N << "+1);" << endl;
    file << "speed = " << std::vector<double>(sol.value(speed)) << ";" << endl;
    file << "pos = " << std::vector<double>(sol.value(pos)) << ";" << endl;
    file << "U = " << std::vector<double>(sol.value(U)) << ";" << endl;

    file << "figure;" << endl;
    file << "hold on;" << endl;
    file << "plot(t,speed);" << endl;
    file << "plot(t,pos);" << endl;
    file << "plot(t,1-sin(2*pi*pos)/2,'r--');" << endl;
    file << "stairs(t(1:end-1),U,'k');" << endl;
    file << "xlabel('Time [s]');" << endl;
    file << "legend('speed','pos','speed "
            "limit','throttle','Location','northwest');"
         << endl;

    // Have a look at the constraint Jacobian
    jacobian(opti.g(), opti.x()).sparsity().spy_matlab("race_car_jac_g.m");

    file << "figure" << std::endl;
    file << "race_car_jac_g;" << std::endl;
    file << "xlabel('decision variables');" << std::endl;
    file << "ylabel('constraints');" << std::endl;
    file << "print('jac_sp','-dpng');" << std::endl;
  }
}

// dx/dt = f(x,u)
MX simple_dyn(const MX& x, const MX& u, const MX& P, int i) { 
  double I(0.21);
  double M(8.91);
  double g(9.81);

  if( i<15 ){
    return vertcat (vertcat(x(3), x(4), x(5)), 
        vertcat(u(0)/M, u(1)/M -g), ((P(1) - x(1))*u(0) - (P(0) - x(0)) * u(1))/I );  // Velocity
  }else{
    return vertcat (vertcat(x(3), x(4), x(5)), 
        vertcat(u(0)/M, u(1)/M -g), ((P(3) - x(1))*u(0) - (P(2) - x(0)) * u(1))/I );  // Velocity
  }
}

DVec<double> simple_dyn_vec(const DVec<double> & x, 
    const DVec<double>& u, const DVec<double>& P, int i) { 
  double I(0.21);
  double M(8.91);
  double g(9.81);

  DVec<double> ret(6);
  if( i<15 ){
    ret << x(3), x(4), x(5), 
        u(0)/M, u(1)/M -g, ((P(1) - x(1))*u(0) - (P(0) - x(0)) * u(1))/I ;  // Velocity
  }else{
    ret << x(3), x(4), x(5), 
        u(0)/M, u(1)/M -g, ((P(3) - x(1))*u(0) - (P(2) - x(0)) * u(1))/I ;  // Velocity
  }

  return ret;
}

// [x, y, theta, dot_x, dot_y, dot_theta]
TEST(casadi, jump_opt) {
  // Mass Matrix Test
  //FloatingBaseModel<double> cheetah = buildMiniCheetah<double>().buildModel();
  //FBModelState<double> state;
  //state.q = DVec<double>::Zero(cheetah::num_act_joint);
  //state.qd = DVec<double>::Zero(cheetah::num_act_joint);
 
  //for(size_t i(0); i<4; ++i){
    //state.q[1 + 3*i] = -M_PI/2.;
    //state.q[2 + 3*i] = M_PI;
  //}
    //state.q[7] = M_PI/2.;
    //state.q[10] = M_PI/2.;
  //state.bodyPosition.setZero();
  //state.bodyVelocity.setZero();
  //state.bodyOrientation[0] = 1.;
  //state.bodyOrientation[1] = 0.;
  //state.bodyOrientation[2] = 0.;
  //state.bodyOrientation[3] = 0.;

  //cheetah.setState(state);
  //cheetah.forwardKinematics();

  //DMat<double> A = cheetah.massMatrix();
  //std::cout<<A<<std::endl;
  
  // Jump optimization
  // ----------------------
  // An optimal control problem (OCP),
  // solved with direct multiple-shooting.
  //
  // For more information see: http://labs.casadi.org/OCP

  int N = 31;  // number of control intervals
  int N_fc = 15;
  //int N_air = 1;
  //int N_hc = 15;
  auto opti = casadi::Opti();  // Optimization problem

  Slice all;
  // ---- decision variables ---------
  auto X = opti.variable(6, N+1);  // state trajectory
  auto P = opti.variable(4);  // Landing Location
  auto F = opti.variable(2, N);  // Reaction force (15: front, 1: air, 15: hind)
  auto cost = opti.variable();

  DVec<double> X0(6); X0.setZero();
  X0<< 0.0, 0.2, 0.2, 0.5, -0.1, 0.1;

  DVec<double> Xf(6); Xf.setZero();
  //Xf<< 0.1, 0.2, -M_PI/4., 0.1, 0.3, 0.0;
  Xf<< 0.5, 0.23, -M_PI/4., 0.8, 0.8, 0.1;

  cost = (X(0,N) - Xf(0)) * (X(0,N) - Xf(0))
    + (X(1,N) - Xf(1)) * (X(1,N) - Xf(1))
    + (X(2,N) - Xf(2)) * (X(2,N) - Xf(2))
    + (X(3,N) - Xf(3)) * (X(3,N) - Xf(3))
    + (X(4,N) - Xf(4)) * (X(4,N) - Xf(4))
    + (X(5,N) - Xf(5)) * (X(5,N) - Xf(5));
  //
  // ---- objective          ---------
  opti.minimize(cost);  // race in minimal time

  // ---- dynamic constraints --------
  auto dt = 0.01;
  double mu(0.5);

  for (int k = 0; k < N; ++k) {
    auto dX = simple_dyn(X(all, k), F(all, k), P, k);
    auto X_next = X(all, k) + dt * dX;
    opti.subject_to(X(all, k + 1) == X_next);  // close the gaps

    // ---- contact constraints -----------
    opti.subject_to(-mu * F(1,k) <= F(0, k) <= mu* F(1,k)); 
    opti.subject_to(0<=F(1,k) <= 300.); 
    // Kinematics constraints
    opti.subject_to (-0.15 <= P(0) <= 0.15);
    opti.subject_to (-0.15 <= P(2) <= 0.15);
    opti.subject_to (P(1) == 0.);
    opti.subject_to (P(3) == 0.);

  }
  opti.subject_to(F(1, N_fc) == 0.);

  for(int i(0); i<6; ++i){
    opti.subject_to(X(i, 1) == X0(i)); 
  }
  // ---- solve NLP              ------
  opti.solver("ipopt");     // set numerical backend
  auto sol = opti.solve();  // actual solve
  //std::vector<double> result = std::vector<double>(sol.debug.value(X(0, all) ) );
  
  //for(int i(0); i<N; ++i){
    //printf("xb: %f\n", result[i]);
  //}



  for(int k(0); k<N; ++k){
    std::vector<double> x_std_vec = std::vector<double>(sol.value(X(all, k)));
    std::vector<double> x_next_std_vec = std::vector<double>(sol.value(X(all, k+1)));
    std::vector<double> f = std::vector<double>(sol.value(F(all, k)));
    std::vector<double> p = std::vector<double>(sol.value(P(all)));

    DVec<double> x_current(6);
    DVec<double> x_next(6);
    DVec<double> f_eigen(2); 
    f_eigen[0] = f[0];
    f_eigen[1] = f[1];

    DVec<double> p_eigen(4); 
    p_eigen[0] = p[0]; 
    p_eigen[1] = p[1]; 
    p_eigen[2] = p[2]; 
    p_eigen[3] = p[3]; 

    for(int i(0); i<6; ++i){
      x_current[i] = x_std_vec[i];
      x_next[i] = x_next_std_vec[i];
    }
    DVec<double> dx_eigen = simple_dyn_vec(x_current, f_eigen, p_eigen, k);
    DVec<double> x_next_eigen = x_current + dt * dx_eigen;
    EXPECT_TRUE(almostEqual(x_next, x_next_eigen, .0001));
    pretty_print(x_next, std::cout, "x_next_casadi");
    pretty_print(x_next_eigen, std::cout, "x_next_eigen");
    printf("\n");
  }

  bool b_matlab_file(true);
  if (b_matlab_file) {
    ofstream file;
    string filename = "jump_result.m";
    file.open(filename.c_str());
    file << "% Results file from " __FILE__ << endl;
    file << "% Generated " __DATE__ " at " __TIME__ << endl;
    file << endl;

    //file << "t = linspace(0," << 0.31 << "," << N << ");" << endl;
    file << "xb= " << std::vector<double>(sol.value(X(0, all) ) ) << ";" << endl;
    file << "zb= " << std::vector<double>(sol.value(X(1, all))) << ";" << endl;
    file << "theta= " << std::vector<double>(sol.value(X(2, all))) << ";" << endl;
    file << "dot_xb= " << std::vector<double>(sol.value(X(3, all) ) ) << ";" << endl;
    file << "dot_zb= " << std::vector<double>(sol.value(X(4, all))) << ";" << endl;
    file << "dot_theta= " << std::vector<double>(sol.value(X(5, all))) << ";" << endl;

    file << "Fx = " << std::vector<double>(sol.value(F(0, all) ) ) << ";" << endl;
    file << "Fz = " << std::vector<double>(sol.value(F(1, all) ) )<< ";" << endl;

    file << "P = " << std::vector<double>(sol.value(P(all) ) ) << ";" << endl;

   //file << "figure;" << endl;
    //file << "hold on;" << endl;
    //file << "plot(xb, zb);" << endl;

    //jacobian(opti.g(), opti.x()).sparsity().spy_matlab("jump_jac_g.m");

    //file << "figure" << std::endl;
    //file << "jump_jac_g;" << std::endl;
    //file << "xlabel('decision variables');" << std::endl;
    //file << "ylabel('constraints');" << std::endl;
    //file << "print('jac_sp','-dpng');" << std::endl;
    file.close();
  }
}
#endif

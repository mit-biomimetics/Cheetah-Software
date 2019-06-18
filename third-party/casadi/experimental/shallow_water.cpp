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

#include <casadi/casadi.hpp>

#include <iomanip>
#include <ctime>
#include <cstdlib>

using namespace casadi;
using namespace std;

class Tester{
public:
  // Constructor
  Tester(int n, int n_euler, int n_finite_elements, int n_meas) : n_(n), n_euler_(n_euler), n_finite_elements_(n_finite_elements), n_meas_(n_meas){}

  // Perform the modelling
  void model();

  // Simulate to generae measurements
  void simulate(double drag_true, double depth_true);

  // Transscribe as an NLP
  void transcribe(bool single_shooting, bool gauss_newton, bool codegen, bool ipopt_as_qpsol, bool regularize, double reg_threshold);

  // Solve the NLP
  void optimize(double drag_guess, double depth_guess, int& iter_count, double& sol_time, double& drag_est, double& depth_est);

  // Dimensions
  int n_;
  int n_euler_;
  int n_finite_elements_;
  int n_meas_;

  // Initial conditions
  DM u0_;
  DM v0_;
  DM h0_;

  // Discrete time dynamics
  Function f_;

  // Generated measurements
  vector<DM> H_meas_;

  // Height of the splash
  double spheight_;

  // Scaling factors for the parameters
  vector<double> p_scale_;

  /// NLP solver
  Function nlpsol_;
};

void Tester::model(){
  // Physical parameters
  double g = 9.81; // gravity
  double poolwidth = 0.2;
  double sprad = 0.03;
  spheight_ = 0.01;
  double endtime = 1.0;

  // Discretization
  int ntimesteps = n_euler_*n_finite_elements_*n_meas_;
  double dt = endtime/ntimesteps;
  double dx = poolwidth/n_;
  double dy = poolwidth/n_;
  vector<double> x(n_), y(n_);
  for(int i=0; i<n_; ++i){
    x[i] = (i+0.5)*dx;
    y[i] = (i+0.5)*dy;
  }

  // Initial conditions
  u0_ = DM::zeros(n_+1,n_  );
  v0_ = DM::zeros(n_  ,n_+1);
  h0_ = DM::zeros(n_  ,n_  );
  bool any_point_in_domain = false;
  for(int i=0; i<n_; ++i){
    for(int j=0; j<n_; ++j){
      double spdist = sqrt(pow((x[i]-0.04),2.) + pow((y[j]-0.04),2.));
      if(spdist<sprad/3.0){
	h0_(i,j) = spheight_ * cos(3.0*M_PI*spdist/(2.0*sprad));
	any_point_in_domain = true;
      }
    }
  }

  // Make sure that there is at least one point with nonzero initial values
  if(!any_point_in_domain){
    int i_splash = std::min(int(0.04/dx),n_-1);
    int j_splash = std::min(int(0.04/dy),n_-1);
    h0_(i_splash,j_splash) = spheight_;
  }

  // Free parameters (nominal values)
  MX drag_nom = MX::sym("drag_nom");
  MX depth_nom = MX::sym("depth_nom");
  MX p = MX::vertcat({drag_nom, depth_nom});

  // Scaling factors for the parameters
  double drag_scale = 1;
  double depth_scale = 0.01;
  p_scale_.resize(2);
  p_scale_[0] = drag_scale;
  p_scale_[1] = depth_scale;

  // Real parameter values
  MX drag = drag_nom*drag_scale;
  MX depth = depth_nom*depth_scale;

  // The state at a measurement
  MX uk = MX::sym("uk",n_+1, n_);
  MX vk = MX::sym("vk",n_  , n_+1);
  MX hk = MX::sym("hk",n_  , n_);

  // Take one step of the integrator
  MX u = uk;
  MX v = vk;
  MX h = hk;

  // Update u
  MX d1 = -dt*g/dx;
  MX d2 = dt*drag;
  u(Slice(1,n_),Slice()) += d1*(h(Slice(1,n_),Slice())-h(Slice(0,n_-1),Slice())) - d2*u(Slice(1,n_),Slice());

  // Update v
  d1 = -dt*g/dy;
  v(Slice(),Slice(1,n_)) += d1*(h(Slice(),Slice(1,n_))-h(Slice(),Slice(0,n_-1))) - d2*v(Slice(),Slice(1,n_));

  // Update h
  d1 = (-depth*dt)*(1.0/dx);
  d2 = (-depth*dt)*(1.0/dy);
  h += d1*(u(Slice(1,n_+1),Slice())-u(Slice(0,n_),Slice())) + d2*(v(Slice(),Slice(1,n_+1))-v(Slice(),Slice(0,n_)));

  // Create an integrator function
  vector<MX> f_step_in(4);
  f_step_in[0] = p;
  f_step_in[1] = uk;
  f_step_in[2] = vk;
  f_step_in[3] = hk;
  vector<MX> f_step_out(3);
  f_step_out[0] = u;
  f_step_out[1] = v;
  f_step_out[2] = h;
  Function f_step("f_step_mx", f_step_in, f_step_out, {"p", "u0", "v0", "h0"}, {"uf", "vf", "hf"});
  cout << "generated single step dynamics (" << f_step.n_nodes() << " nodes)" << endl;

  // Expand the discrete dynamics?
  if(false){
    f_step = f_step.expand();
    cout << "generated single step dynamics, SX (" << f_step.n_nodes() << " nodes)" << endl;
  }

  // Integrate over one subinterval
  vector<MX> f_in(4);
  MX P = MX::sym("P",2);
  MX Uk = MX::sym("Uk",n_+1, n_);
  MX Vk = MX::sym("Vk",n_  , n_+1);
  MX Hk = MX::sym("Hk",n_  , n_);
  f_in[0] = P;
  f_in[1] = Uk;
  f_in[2] = Vk;
  f_in[3] = Hk;
  vector<MX> f_inter = f_in;
  vector<MX> f_out;
  for(int j=0; j<n_euler_; ++j){
    // Create a call node
    f_out = f_step(f_inter);

    // Save intermediate state
    f_inter[1] = f_out[0];
    f_inter[2] = f_out[1];
    f_inter[3] = f_out[2];
  }

  // Create an integrator function
  f_ = Function("f_mx", f_in, f_out, {"P", "U0", "V0", "H0"}, {"UF", "VF", "HF"});
  cout << "generated discrete dynamics for one finite element (" << f_.n_nodes() << " MX nodes)" << endl;

  // Integrate over the complete interval
  if(n_finite_elements_>1){
    f_in[0] = P;
    f_in[1] = Uk;
    f_in[2] = Vk;
    f_in[3] = Hk;
    f_inter = f_in;
    for(int j=0; j<n_finite_elements_; ++j){
      // Create a call node
      f_out = f_(f_inter);

      // Save intermediate state
      f_inter[1] = f_out[0];
      f_inter[2] = f_out[1];
      f_inter[3] = f_out[2];
    }

    // Create an integrator function
    f_ = Function("f_mx", f_in, f_out, {"P", "U0", "V0", "H0"}, {"UF", "VF", "HF"});
    cout << "generated discrete dynamics for complete interval (" << f_.n_nodes() << " MX nodes)" << endl;
  }

  // Expand the discrete dynamics
  if(false){
    f_ = f_.expand("f_sx");
    cout << "generated discrete dynamics, SX (" << f_.n_nodes() << " nodes)" << endl;
  }
}

void Tester::simulate(double drag_true, double depth_true){

  // Measurements
  H_meas_.reserve(n_meas_);

  // Unscaled parameter values
  vector<double> p_true(2); p_true[0]=drag_true; p_true[1]=depth_true;
  for(int i=0; i<2; ++i){
    p_true[i] /= p_scale_[i];
  }

  // Simulate once to generate "measurements"
  vector<DM> arg = {p_true, u0_, v0_, h0_};
  clock_t time1 = clock();
  for(int k=0; k<n_meas_; ++k){
    vector<DM> res = f_(arg);
    const DM& u = res.at(0);
    const DM& v = res.at(1);
    const DM& h = res.at(2);
    arg.at(1) = u;
    arg.at(2) = v;
    arg.at(3) = h;

    // Save a copy of h
    H_meas_.push_back(h);
  }
  clock_t time2 = clock();
  double t_elapsed = double(time2-time1)/CLOCKS_PER_SEC;
  cout << "measurements generated in " << t_elapsed << " seconds." << endl;
}


void Tester::transcribe(bool single_shooting, bool gauss_newton, bool codegen, bool ipopt_as_qpsol, bool regularize, double reg_threshold){

  // NLP variables
  MX P = MX::sym("P",2);

  // Variables in the lifted NLP
  stringstream ss;

  // Objective function terms
  vector<MX> nlp_fv;
  if(!gauss_newton) nlp_fv.push_back(0);

  // Constraint function terms
  vector<MX> nlp_gv;

  // Generate full-space NLP
  MX U = u0_;
  MX V = v0_;
  MX H = h0_;
  for(int k=0; k<n_meas_; ++k){
    // Take a step
    vector<MX> f_res = f_(vector<MX>{P,U,V,H});
    U = f_res[0];
    V = f_res[1];
    H = f_res[2];

    if(!single_shooting){
      // Lift the heights, initialized with measurements
      H = lift(H, H_meas_[k]);

      // Initialize with initial conditions
      // U = lift(U, u0_);
      // V = lift(V, v0_);
      // H = lift(H, DM::zeros(n_  ,n_));

      // Initialize through simulation
      // U = lift(U, U);
      // V = lift(V, V);
      // H = lift(H, H);
    }

    // Objective function term
    MX H_dev = H-H_meas_[k];
    if(gauss_newton){
      nlp_fv.push_back(vec(H_dev));
    } else {
      nlp_fv.front() += dot(H_dev, H_dev)/2;
    }

    // Add to the constraints
    nlp_gv.push_back(vec(H));
  }

  MXDict nlp = {{"x", P}, {"f", vertcat(nlp_fv)}, {"g", vertcat(nlp_gv)}};
  cout << "Generated single-shooting NLP" << endl;

  // NLP Solver
  Dict opts;
  opts["verbose"] = true;
  opts["regularize"] = regularize;
  opts["codegen"] = codegen;
  opts["reg_threshold"] = reg_threshold;
  opts["max_iter_ls"] = 3;
  opts["beta"] = 0.5;
  opts["merit_start"] = 1e-3;
  //opts["merit_memory"] = 1;
  opts["max_iter"] = 100;
  opts["compiler"] = "shell";
  opts["jit_options"] = Dict{{"compiler", "clang"}, {"flags", vector<string>{"-O2"}}};
  if(gauss_newton){
    opts["hessian_approximation"] = "gauss-newton";
  }

  // Print both of the variables
  opts["name_x"] = vector<string>{"drag", "depth"};
  opts["print_x"] = range(2);

  if(ipopt_as_qpsol){
    opts["qpsol"] = "nlpsol";
    Dict nlp_opts = {{"ipopt.tol", 1e-12}, {"ipopt.print_level", 0}, {"print_time", false}};
    opts["qpsol_options"] = Dict{{"nlpsol", "ipopt"}, {"nlpsol_options", nlp_opts}};
  } else {
    opts["qpsol"] = "qpoases";
    opts["qpsol_options"] = Dict{{"printLevel", "none"}};
  }

  // Create NLP solver instance
  nlpsol_ = nlpsol("nlpsol", "scpgen", nlp, opts);
}

void Tester::optimize(double drag_guess, double depth_guess, int& iter_count, double& sol_time, double& drag_est, double& depth_est){
  cout << "Starting parameter estimation" << endl;

  // Initial guess
  vector<double> p_init(2);
  p_init[0] = drag_guess/p_scale_[0];
  p_init[1] = depth_guess/p_scale_[1];

  // Bounds on the variables
  vector<double> lbu(2), ubu(2);
  lbu.at(0) = 1.0e-1 / p_scale_[0]; // drag positive
  lbu.at(1) = 5.0e-4 / p_scale_[1]; // depth positive
  ubu.at(0) = 100.0 / p_scale_[0]; // max drag
  ubu.at(1) =  0.10 / p_scale_[1]; // max depth

  clock_t time1 = clock();
  map<string, DM> w = {{"x0", p_init},
                       {"lbx", lbu},
                       {"ubx", ubu},
                       {"lbg", -spheight_},
                       {"ubg",  spheight_}};
  w = nlpsol_(w);
  clock_t time2 = clock();

  // Solution statistics
  sol_time = double(time2-time1)/CLOCKS_PER_SEC;
  const vector<double>& x_opt = w.at("x").nonzeros();
  drag_est = x_opt.at(0)*p_scale_[0];
  depth_est = x_opt.at(1)*p_scale_[1];
  iter_count = nlpsol_.stats().at("iter_count");
}

int main(){

  // True parameter values
  double drag_true = 2.0, depth_true = 0.01;

  // Use IPOPT as QP solver (can handle non-convex QPs)
  bool ipopt_as_qpsol = true;

  // Use Gauss-Newton method
  bool gauss_newton = true;

  // Codegen the Lifted Newton functions
  bool codegen = true;

  // Regularize the QP
  bool regularize = true;

  // Smallest allowed eigenvalue for the regularization
  double reg_threshold = 1e-8;

  // Problem size
  //  int  n = 100, n_euler = 100, n_finite_elements = 1, n_meas = 100;
  //int  n = 30, n_euler = 100, n_finite_elements = 1, n_meas = 100; // Paper
  int n = 15, n_euler = 20, n_finite_elements = 1, n_meas = 20;

  // Initial guesses
  vector<double> drag_guess, depth_guess;
  drag_guess.push_back( 2.0); depth_guess.push_back(0.01); // Optimal solution
  drag_guess.push_back( 0.5); depth_guess.push_back(0.01);
  drag_guess.push_back( 5.0); depth_guess.push_back(0.01);
  drag_guess.push_back(15.0); depth_guess.push_back(0.01);
  drag_guess.push_back(30.0); depth_guess.push_back(0.01);
  drag_guess.push_back( 2.0); depth_guess.push_back(0.005);
  drag_guess.push_back( 2.0); depth_guess.push_back(0.02);
  drag_guess.push_back( 2.0); depth_guess.push_back(0.1);
  drag_guess.push_back( 0.2); depth_guess.push_back(0.001);
  drag_guess.push_back( 1.0); depth_guess.push_back(0.005);
  drag_guess.push_back( 4.0); depth_guess.push_back(0.02);
  drag_guess.push_back( 1.0); depth_guess.push_back(0.02);
  drag_guess.push_back(20.0); depth_guess.push_back(0.001);

  // Number of tests
  const int n_tests = drag_guess.size();

  // Number of iterations
  vector<int> iter_count_gn(n_tests,-1);
  vector<int> iter_count_eh(n_tests,-1);

  // Solution time
  vector<double> sol_time_gn(n_tests,-1);
  vector<double> sol_time_eh(n_tests,-1);

  // Estimated drag and depth
  vector<double> drag_est_gn(n_tests,-1);
  vector<double> depth_est_gn(n_tests,-1);
  vector<double> drag_est_eh(n_tests,-1);
  vector<double> depth_est_eh(n_tests,-1);

  // Create a tester object
  Tester t(n,n_euler,n_finite_elements,n_meas);

  // Perform the modelling
  t.model();

  // Optimization parameters
  t.simulate(drag_true, depth_true);

  // For both single and multiple shooting
  for(int sol=0; sol<2; ++sol){

    // Transcribe as an NLP
    bool single_shooting = sol==0;
    t.transcribe(single_shooting, gauss_newton, codegen, ipopt_as_qpsol, regularize, reg_threshold);

    // Run tests
    for(int test=0; test<n_tests; ++test){
      // Print progress
      cout << "test " << test << endl;
      try{
	t.optimize(drag_guess[test],depth_guess[test],
		   sol==0 ? iter_count_gn[test] : iter_count_eh[test],
		   sol==0 ? sol_time_gn[test] : sol_time_eh[test],
		   sol==0 ? drag_est_gn[test] : drag_est_eh[test],
		   sol==0 ? depth_est_gn[test] : depth_est_eh[test]);

	// Estimated drag
      } catch(exception& ex){
	cout << "Test " << test << " failed: " << ex.what() << endl;
      }
    }
  }

  // Tolerance
  double tol=1e-3;

  cout <<
  setw(10) << "drag" <<  "  &" <<
  setw(10) << "depth" << "  &" <<
  setw(10) << "iter_ss" << "  &" <<
  setw(10) << "time_ss" << "  &" <<
  setw(10) << "iter_ms" << "  &" <<
  setw(10) << "time_ms" << "  \\\\ \%" <<
  setw(10) << "edrag_ss" <<
  setw(10) << "edepth_ss" <<
  setw(10) << "edrag_ms" <<
  setw(10) << "edepth_ms" << endl;
  for(int test=0; test<n_tests; ++test){
    cout << setw(10) << drag_guess[test] << "  &";
    cout << setw(10) << depth_guess[test] << "  &";
    if(fabs(drag_est_gn[test]-drag_true) + fabs(depth_est_gn[test]-depth_true)<tol){
      cout << setw(10) << iter_count_gn[test] << "  &";
      cout << setw(10) << sol_time_gn[test] << "  &";
    } else {
      cout << setw(10) << "$\\infty$" << "  &";
      cout << setw(10) << "$\\infty$" << "  &";
    }
    if(fabs(drag_est_eh[test]-drag_true) + fabs(depth_est_eh[test]-depth_true)<tol){
      cout << setw(10) << iter_count_eh[test] << "  &";
      cout << setw(10) << sol_time_eh[test] << "  \\\\ \%";
    } else {
      cout << setw(10) << "$\\infty$" << "  &";
      cout << setw(10) << "$\\infty$" << "  \\\\ \%";
    }
    cout << setw(10) << drag_est_gn[test];
    cout << setw(10) << depth_est_gn[test];
    cout << setw(10) << drag_est_eh[test];
    cout << setw(10) << depth_est_eh[test] << endl;
  }

  return 0;
}

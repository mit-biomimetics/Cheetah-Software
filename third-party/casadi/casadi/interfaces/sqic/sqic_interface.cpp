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


#include "sqic_interface.hpp"
#include "casadi/core/conic.hpp"
#include "casadi/core/mx_function.hpp"
#include "casadi/core/casadi_misc.hpp"

#include "wsqic.hpp"
#include "casadi/interfaces/sqic/resource_sqic.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_CONIC_SQIC_EXPORT
  casadi_register_conic_sqic(Conic::Plugin* plugin) {
    plugin->creator = SqicInterface::creator;
    plugin->name = "sqic";
    plugin->doc = SqicInterface::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &SqicInterface::options_;
    return 0;
  }

  extern "C"
  void CASADI_CONIC_SQIC_EXPORT casadi_load_conic_sqic() {
    Conic::registerPlugin(casadi_register_conic_sqic);
  }

  SqicInterface::SqicInterface(const std::map<std::string, Sparsity>& st) : Conic(st) {
    is_init_ = false;
  }

  SqicInterface::~SqicInterface() {
    sqicDestroy();
  }

  void SqicInterface::evaluate() {
    if (inputs_check_) {
      check_inputs(input(CONIC_LBX).ptr(), input(CONIC_UBX).ptr(),
                  input(CONIC_LBA).ptr(), input(CONIC_UBA).ptr());
    }

    std::copy(input(CONIC_X0)->begin(), input(CONIC_X0)->end(), x_.begin());
    std::fill(x_.begin()+n_, x_.end(), 0);

    std::transform(input(CONIC_LAM_X0)->begin(), input(CONIC_LAM_X0)->end(),
                   rc_.begin(), negate<double>());
    std::fill(rc_.begin()+n_, rc_.end(), 0);

    std::copy(input(CONIC_LBX)->begin(), input(CONIC_LBX)->end(), bl_.begin());
    std::copy(input(CONIC_UBX)->begin(), input(CONIC_UBX)->end(), bu_.begin());

    std::copy(input(CONIC_LBA)->begin(), input(CONIC_LBA)->end(), bl_.begin()+n_);
    std::copy(input(CONIC_UBA)->begin(), input(CONIC_UBA)->end(), bu_.begin()+n_);

    for (casadi_int i=0;i<n_+nc_+1;++i) {
      if (bl_[i]==-std::numeric_limits<double>::infinity()) bl_[i]=-inf_;
      if (bu_[i]==std::numeric_limits<double>::infinity()) bu_[i]=inf_;
    }

    formatA_.setInput(input(CONIC_A), 0);
    formatA_.setInput(input(CONIC_G), 1);
    formatA_.evaluate();

    sqicSolve(&output(CONIC_COST).nonzeros()[0]);

    std::copy(x_.begin(), x_.begin()+n_, output(CONIC_X)->begin());
    std::transform(rc_.begin(), rc_.begin()+n_, output(CONIC_LAM_X)->begin(),
                   negate<double>());
    std::transform(rc_.begin()+n_, rc_.begin()+n_+nc_, output(CONIC_LAM_A)->begin(),
                   negate<double>());

    output(CONIC_COST)[0]+= x_[n_+nc_];
  }

  void SqicInterface::init() {
    // Call the init method of the base class
    Conic::init();

    if (is_init_) sqicDestroy();

    inf_ = 1.0e+20;

    // Allocate data structures for SQIC
    bl_.resize(n_+nc_+1, 0);
    bu_.resize(n_+nc_+1, 0);
    x_.resize(n_+nc_+1, 0);
    hs_.resize(n_+nc_+1, 0);
    hEtype_.resize(n_+nc_+1, 0);
    pi_.resize(nc_+1, 0);
    rc_.resize(n_+nc_+1, 0);

    locH_ = st_[QP_STRUCT_H].colind();
    indH_ = st_[QP_STRUCT_H].row();

    // Fortran indices are one-based
    for (casadi_int i=0;i<indH_.size();++i) indH_[i]+=1;
    for (casadi_int i=0;i<locH_.size();++i) locH_[i]+=1;

    // Sparsity of augmented linear constraint matrix
    Sparsity A_ = vertcat(st_[QP_STRUCT_A], Sparsity::dense(1, n_));
    locA_ = A_.colind();
    indA_ = A_.row();

    // Fortran indices are one-based
    for (casadi_int i=0;i<indA_.size();++i) indA_[i]+=1;
    for (casadi_int i=0;i<locA_.size();++i) locA_[i]+=1;

    // helper functions for augmented linear constraint matrix
    MX a = MX::sym("A", st_[QP_STRUCT_A]);
    MX g = MX::sym("g", n_);
    std::vector<MX> ins;
    ins.push_back(a);
    ins.push_back(g);
    formatA_ = Function("formatA", ins, vertcat(a, g.T()));

    // Set objective row of augmented linear constraints
    bu_[n_+nc_] = inf_;
    bl_[n_+nc_] = -inf_;

    is_init_ = true;

    casadi_int n = n_;
    casadi_int m = nc_+1;

    casadi_int nnzA=formatA_.size_out(0);
    casadi_int nnzH=input(CONIC_H).size();

    std::fill(hEtype_.begin()+n_, hEtype_.end(), 3);

    sqic(&m , &n, &nnzA, &indA_[0], &locA_[0], &formatA_.output().nonzeros()[0], &bl_[0], &bu_[0],
         &hEtype_[0], &hs_[0], &x_[0], &pi_[0], &rc_[0], &nnzH, &indH_[0], &locH_[0],
         &input(CONIC_H).nonzeros()[0]);

  }

  std::map<casadi_int, string> SqicInterface::calc_flagmap() {
    std::map<casadi_int, string> f;

    return f;
  }

  std::map<casadi_int, string> SqicInterface::flagmap = SqicInterface::calc_flagmap();

  void SqicInterface::sqic_error(const string& module, casadi_int flag) {
    // Find the error
    std::map<casadi_int, string>::const_iterator it = flagmap.find(flag);

    stringstream ss;
    if (it == flagmap.end()) {
      ss << "Unknown error (" << flag << ") from module \"" << module << "\".";
    } else {
      ss << "Module \"" << module << "\" returned flag \"" << it->second << "\".";
    }
    ss << " Consult SQIC documentation.";
    casadi_error(ss.str());
  }

  void SqicInterface::generateNativeCode(std::ostream& file) const {
    // Dump the contents of resource_sqic, but filter out the C bind stuff
    std::string resource_sqic_input(resource_sqic);
    std::istringstream stream(resource_sqic_input);
    std::string line;
    while (std::getline(stream, line)) {
      size_t b_i = line.find("bind ( C, ");
      if (b_i!=std::string::npos) {
        file << line.substr(0, b_i) << std::endl;
      } else {
        file << line << std::endl;
      }
    }

    file.precision(std::numeric_limits<double>::digits10+2);
    file << std::scientific; // This is really only to force a decimal dot,
    // would be better if it can be avoided

    file << "program exported" << std::endl;
    file << "  use SQICModule" << std::endl;
    file << "  implicit none" << std::endl;
    file << "  integer(ip)               :: m, n, n_inf, nnH, nnzH, nnzA, nS" << std::endl;


    file << "  real(rp)                  :: Obj" << std::endl;

    file << "  real(rp), allocatable:: bl(:), bu(:), x(:), valA(:), valH(:) , pi(:), rc(:)"
         << std::endl;
    file << "  integer(ip), allocatable:: indA(:), locA(:), indH(:), locH(:), hEtype(:), hs(:)"
         << std::endl;

    casadi_int n = n_;
    casadi_int m = nc_+1;
    casadi_int nnzA=formatA_.size_out(0);
    casadi_int nnzH=input(CONIC_H).size();

    file << "  n = " << n << std::endl;
    file << "  m = " << m << std::endl;
    file << "  nnzA = " << nnzA << std::endl;
    file << "  nnzH = " << nnzH << std::endl;

    file << "  allocate ( bl(n+m), bu(n+m) )" << std::endl;
    file << "  allocate ( hEtype(n+m) )" << std::endl;
    file << "  allocate ( locA(n+1), valA(nnzA), indA(nnzA) )" << std::endl;
    file << "  allocate ( pi(m), rc(n+m), x(n+m) )" << std::endl;
    file << "  allocate ( hs(n+m) )" << std::endl;
    file << "  allocate ( valH(nnzH), locH(n+1), indH(nnzH) )" << std::endl;

    for (casadi_int i=0;i<indA_.size();++i) {
      file << "  indA(" << i +1 << ") = " << indA_[i] << std::endl;
    }
    for (casadi_int i=0;i<locA_.size();++i) {
      file << "  locA(" << i +1 << ") = " << locA_[i] << std::endl;
    }
    for (casadi_int i=0;i<formatA_.size_out(0);++i) {
      file << "  valA(" << i +1 << ") = " << formatA_.output().at(i) << std::endl;
    }
    for (casadi_int i=0;i<bl_.size();++i) {
      file << "  bl(" << i +1 << ") = " << bl_[i] << std::endl;
      file << "  bu(" << i +1 << ") = " << bu_[i] << std::endl;
    }
    for (casadi_int i=0;i<hEtype_.size();++i) {
      file << "  hEtype(" << i +1 << ") = " << hEtype_[i] << std::endl;
    }
    for (casadi_int i=0;i<hs_.size();++i) {
      file << "  hs(" << i +1 << ") = " << hs_[i] << std::endl;
    }
    for (casadi_int i=0;i<indH_.size();++i) {
      file << "  indH(" << i +1 << ") = " << indH_[i] << std::endl;
    }
    for (casadi_int i=0;i<locH_.size();++i) {
      file << "  locH(" << i +1 << ") = " << locH_[i] << std::endl;
    }
    for (casadi_int i=0;i<input(CONIC_H).size();++i) {
      file << "  valH(" << i +1 << ") = " << input(CONIC_H).at(i) << std::endl;
    }
    for (casadi_int i=0;i<input(CONIC_X0).size();++i) {
      file << "  x(" << i +1 << ") = " << input(CONIC_X0).at(i) << std::endl;
    }
    for (casadi_int i=0;i<pi_.size();++i) {
      file << "  pi(" << i +1 << ") = " <<  0 << std::endl; //pi_[i] << std::endl;
    }
    uout() << "lam_x0:::" << input(CONIC_LAM_X0) << std::endl;
    for (casadi_int i=0;i<rc_.size();++i) {
      file << "  rc(" << i +1 << ") = "
           << ((i<input(CONIC_LAM_X0).size()) ? -input(CONIC_LAM_X0).at(i) : 0.0)
           << std::endl;
    }

    file << "  call wsqic (m, n, nnzA, indA, locA, valA, bl, bu, hEtype, "
         << "hs, x, pi, rc, nnzH, indH, locH, valH)" << std::endl;
    /**for (casadi_int i=0;i<input(CONIC_X0).size();++i) {
       file << "  x(" << i +1 << ") = " << input(CONIC_X0).at(i) << std::endl;
       }
       for (casadi_int i=0;i<pi_.size();++i) {
       file << "  pi(" << i +1 << ") = " << pi_[i] << std::endl;
       }
       uout() << "lam_x0:::" << input(CONIC_LAM_X0) << std::endl;
       for (casadi_int i=0;i<rc_.size();++i) {
       file << "  rc(" << i +1 << ") = "
       << ((i<input(CONIC_LAM_X0).size()) ? -input(CONIC_LAM_X0).at(i) : 0.0)
       << std::endl;
       }*/
    /**
       file << "  call sqicSolve(Obj)" << std::endl;
       for (casadi_int i=0;i<input(CONIC_X0).size();++i) {
       file << "  x(" << i +1 << ") = " << input(CONIC_X0).at(i) << std::endl;
       }
       for (casadi_int i=0;i<pi_.size();++i) {
       file << "  pi(" << i +1 << ") = " << pi_[i] << std::endl;
       }
       uout() << "lam_x0:::" << input(CONIC_LAM_X0) << std::endl;
       for (casadi_int i=0;i<rc_.size();++i) {
       file << "  rc(" << i +1 << ") = "
       << ((i<input(CONIC_LAM_X0).size()) ? -input(CONIC_LAM_X0).at(i) : 0.0)
       << std::endl;
       }
    */
    file << "  call sqicSolve(Obj)" << std::endl;
    file << "  deallocate ( bl, bu )" << std::endl;
    file << "  deallocate ( hEtype )" << std::endl;
    file << "  deallocate ( locA, valA, indA )" << std::endl;
    file << "  deallocate ( pi, rc, x )" << std::endl;
    file << "  deallocate ( valH, locH, indH )" << std::endl;
    file << "  call sqicDestroy()" << std::endl;
    file << "end program exported" << std::endl;


  }

} // namespace casadi

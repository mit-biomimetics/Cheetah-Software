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

#include "hpmpc_interface.hpp"
#include <numeric>

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_CONIC_HPMPC_EXPORT
  casadi_register_conic_hpmpc(Conic::Plugin* plugin) {
    plugin->creator = HpmpcInterface::creator;
    plugin->name = "hpmpc";
    plugin->doc = HpmpcInterface::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &HpmpcInterface::options_;
    return 0;
  }

  extern "C"
  void CASADI_CONIC_HPMPC_EXPORT casadi_load_conic_hpmpc() {
    Conic::registerPlugin(casadi_register_conic_hpmpc);
  }

  HpmpcInterface::HpmpcInterface(const std::string& name,
                                     const std::map<std::string, Sparsity>& st)
    : Conic(name, st) {
  }

  HpmpcInterface::~HpmpcInterface() {
    clear_mem();
  }

  Options HpmpcInterface::options_
  = {{&Conic::options_},
     {{"N",
       {OT_INT,
        "OCP horizon"}},
      {"nx",
       {OT_INTVECTOR,
        "Number of states, length N+1"}},
      {"nu",
       {OT_INTVECTOR,
        "Number of controls, length N"}},
      {"ng",
       {OT_INTVECTOR,
        "Number of non-dynamic constraints, length N+1"}},
      {"mu0",
       {OT_DOUBLE,
        "Max element in cost function as estimate of max multiplier"}},
      {"max_iter",
       {OT_INT,
        "Max number of iterations"}},
      {"tol",
       {OT_DOUBLE,
        "Tolerance in the duality measure"}},
      {"warm_start",
       {OT_BOOL,
        "Use warm-starting"}},
      {"inf",
       {OT_DOUBLE,
        "HPMPC cannot handle infinities. Infinities will be replaced by this option's value."}},
      {"target",
       {OT_STRING,
        "hpmpc target"}},
      {"blasfeo_target",
       {OT_STRING,
        "hpmpc target"}}}
  };

  void HpmpcInterface::init(const Dict& opts) {
    Conic::init(opts);

    // Default options
    mu0_ = 1;
    max_iter_ = 1000;
    print_level_ = 1;
    tol_ = 1e-8;
    warm_start_ = false;
    inf_ = 1e6;
    target_ = "C99_4X4";
    blasfeo_target_ = "GENERIC";
    casadi_int struct_cnt=0;
    // Read options
    for (auto&& op : opts) {
      if (op.first=="N") {
        N_ = op.second;
        struct_cnt++;
      } else if (op.first=="nx") {
        nxs_ = op.second;
        struct_cnt++;
      } else if (op.first=="nu") {
        nus_ = op.second;
        struct_cnt++;
      } else if (op.first=="ng") {
        ngs_ = op.second;
        struct_cnt++;
      } else if (op.first=="mu0") {
        mu0_ = op.second;
      } else if (op.first=="max_iter") {
        max_iter_ = op.second;
      } else if (op.first=="tol") {
        tol_ = op.second;
      } else if (op.first=="warm_start") {
        warm_start_ = op.second;
      } else if (op.first=="target") {
        target_ = static_cast<std::string>(op.second);
      } else if (op.first=="blasfeo_target") {
        blasfeo_target_ = static_cast<std::string>(op.second);
      } else if (op.first=="print_level") {
        print_level_ = op.second;
      } else if (op.first=="inf") {
        inf_ = op.second;
      }
    }

    bool detect_structure = struct_cnt==0;
    casadi_assert(struct_cnt==0 || struct_cnt==4,
      "You must either set all of N, nx, nu, ng; "
      "or set none at all (automatic detection).");

    const std::vector<casadi_int>& nx = nxs_;
    const std::vector<casadi_int>& ng = ngs_;
    const std::vector<casadi_int>& nu = nus_;

    if (detect_structure) {
      /* General strategy: look for the xk+1 diagonal part in A
      */

      // Find the right-most column for each row in A -> A_skyline
      // Find the second-to-right-most column -> A_skyline2
      // Find the left-most column -> A_bottomline
      Sparsity AT = A_.T();
      std::vector<casadi_int> A_skyline;
      std::vector<casadi_int> A_skyline2;
      std::vector<casadi_int> A_bottomline;
      for (casadi_int i=0;i<AT.size2();++i) {
        casadi_int pivot = AT.colind()[i+1];
        A_bottomline.push_back(AT.row()[AT.colind()[i]]);
        if (pivot>AT.colind()[i]) {
          A_skyline.push_back(AT.row()[pivot-1]);
          if (pivot>AT.colind()[i]+1) {
            A_skyline2.push_back(AT.row()[pivot-2]);
          } else {
            A_skyline2.push_back(-1);
          }
        } else {
          A_skyline.push_back(-1);
          A_skyline2.push_back(-1);
        }
      }

      /*
      Loop over the right-most columns of A:
      they form the diagonal part due to xk+1 in gap constraints.
      detect when the diagonal pattern is broken -> new stage
      */
      casadi_int pivot = 0; // Current right-most element
      casadi_int start_pivot = pivot; // First right-most element that started the stage
      casadi_int cg = 0; // Counter for non-gap-closing constraints
      for (casadi_int i=0;i<na_;++i) { // Loop over all rows
        bool commit = false; // Set true to jump to the stage
        if (A_skyline[i]>pivot+1) { // Jump to a diagonal in the future
          nus_.push_back(A_skyline[i]-pivot-1); // Size of jump equals number of states
          commit = true;
        } else if (A_skyline[i]==pivot+1) { // Walking the diagonal
          if (A_skyline2[i]<start_pivot) { // Free of below-diagonal entries?
            pivot++;
          } else {
            nus_.push_back(0); // We cannot but conclude that we arrived at a new stage
            commit = true;
          }
        } else { // non-gap-closing constraint detected
          cg++;
        }

        if (commit) {
          nxs_.push_back(pivot-start_pivot+1);
          ngs_.push_back(cg); cg=0;
          start_pivot = A_skyline[i];
          pivot = A_skyline[i];
        }
      }
      nxs_.push_back(pivot-start_pivot+1);

      // Correction for k==0
      nxs_[0] = A_skyline[0];
      nus_[0] = 0;
      ngs_.erase(ngs_.begin());
      casadi_int cN=0;
      for (casadi_int i=na_-1;i>=0;--i) {
        if (A_bottomline[i]<start_pivot) break;
        cN++;
      }
      ngs_.push_back(cg-cN);
      ngs_.push_back(cN);

      N_ = nus_.size();
      if (verbose_) {
        casadi_message("Detected structure: N " + str(N_) + ", nx " + str(nx) + ", "
          "nu " + str(nu) + ", ng " + str(ng) + ".");
      }
    }

    casadi_assert_dev(nx.size()==N_+1);
    casadi_assert_dev(nu.size()==N_);
    casadi_assert_dev(ng.size()==N_+1);

    casadi_assert(nx_ == std::accumulate(nx.begin(), nx.end(), 0) +
      std::accumulate(nu.begin(), nu.end(), 0),
      "sum(nx)+sum(nu) = must equal total size of variables (" + str(nx_) + "). "
      "Structure is: N " + str(N_) + ", nx " + str(nx) + ", "
      "nu " + str(nu) + ", ng " + str(ng) + ".");
    casadi_assert(na_ == std::accumulate(nx.begin()+1, nx.end(), 0) +
      std::accumulate(ng.begin(), ng.end(), 0),
      "sum(nx+1)+sum(ng) = must equal total size of constraints (" + str(na_) + "). "
      "Structure is: N " + str(N_) + ", nx " + str(nx) + ", "
      "nu " + str(nu) + ", ng " + str(ng) + ".");
    // Load library HPMPC when applicable
    std::string searchpath;

#ifdef HPMPC_DLOPEN

    std::string libname = "casadi_hpmpc_" + target_ + "_" + blasfeo_target_;
    DL_HANDLE_TYPE handle = load_library(libname, searchpath, true);

    std::string work_size_name = "hpmpc_d_ip_ocp_hard_tv_work_space_size_bytes";

#ifdef _WIN32
    hpmpc_d_ip_ocp_hard_tv_work_space_size_bytes =
      (Work_size)GetProcAddress(handle, TEXT(work_size_name.c_str()));
#else // _WIN32
    // Reset error
    dlerror();

    // Load creator
    hpmpc_d_ip_ocp_hard_tv_work_space_size_bytes = (Work_size)dlsym(handle, work_size_name.c_str());
#endif // _WIN32

    casadi_assert(hpmpc_d_ip_ocp_hard_tv_work_space_size_bytes!=nullptr,
      "HPMPC interface: symbol \"" + work_size_name + "\" found in " + searchpath + ".");

    std::string ocp_solve_name = "fortran_order_d_ip_ocp_hard_tv";

#ifdef _WIN32
    fortran_order_d_ip_ocp_hard_tv =
      (Ocp_solve)GetProcAddress(handle, TEXT(ocp_solve_name.c_str()));
#else // _WIN32
    // Reset error
    dlerror();

    // Load creator
    fortran_order_d_ip_ocp_hard_tv = (Ocp_solve)dlsym(handle, ocp_solve_name.c_str());
#endif // _WIN32

    casadi_assert(fortran_order_d_ip_ocp_hard_tv!=nullptr,
      "HPMPC interface: symbol \"" + ocp_solve_name + "\" found in " + searchpath + ".");
#endif

    /* Disassemble A input into:
       A B I
       C D
           A B I
           C D
    */
    casadi_int offset_r = 0, offset_c = 0;
    for (casadi_int k=0;k<N_;++k) { // Loop over blocks
      A_blocks.push_back({offset_r,        offset_c,            nx[k+1], nx[k]});
      B_blocks.push_back({offset_r,        offset_c+nx[k],      nx[k+1], nu[k]});
      C_blocks.push_back({offset_r+nx[k+1], offset_c,           ng[k], nx[k]});
      D_blocks.push_back({offset_r+nx[k+1], offset_c+nx[k],     ng[k], nu[k]});

      offset_c+= nx[k]+nu[k];
      if (k+1<N_)
        I_blocks.push_back({offset_r, offset_c, nx[k+1], nx[k+1]});
      else
        I_blocks.push_back({offset_r, offset_c, nx[k+1], nx[k+1]});
      offset_r+= nx[k+1]+ng[k];
    }

    C_blocks.push_back({offset_r, offset_c,      ng[N_], nx[N_]});

    Asp_ = blocksparsity(na_, nx_, A_blocks);
    Bsp_ = blocksparsity(na_, nx_, B_blocks);
    Csp_ = blocksparsity(na_, nx_, C_blocks);
    Dsp_ = blocksparsity(na_, nx_, D_blocks);
    Isp_ = blocksparsity(na_, nx_, I_blocks, true);

    Sparsity total = Asp_ + Bsp_ + Csp_ + Dsp_ + Isp_;
    casadi_assert((A_ + total).nnz() == total.nnz(),
      "HPMPC: specified structure of A does not correspond to what the interface can handle. "
      "Structure is: N " + str(N_) + ", nx " + str(nx) + ", nu " + str(nu) + ", "
      "ng " + str(ng) + ".");
    casadi_assert_dev(total.nnz() == Asp_.nnz() + Bsp_.nnz() + Csp_.nnz() + Dsp_.nnz()
                      + Isp_.nnz());

    /* Disassemble H input into:
       Q S'
       S R
           Q S'
           S R

       Multiply by 2
    */
    casadi_int offset = 0;
    for (casadi_int k=0;k<N_;++k) { // Loop over blocks
      R_blocks.push_back({offset+nx[k], offset+nx[k],       nu[k], nu[k]});
      S_blocks.push_back({offset+nx[k], offset,             nu[k], nx[k]});
      Q_blocks.push_back({offset,       offset,             nx[k], nx[k]});
      offset+= nx[k]+nu[k];
    }
    Q_blocks.push_back({offset,         offset,       nx[N_], nx[N_]});

    Rsp_ = blocksparsity(nx_, nx_, R_blocks);
    Ssp_ = blocksparsity(nx_, nx_, S_blocks);
    Qsp_ = blocksparsity(nx_, nx_, Q_blocks);

    total = Rsp_ + Ssp_ + Qsp_ + Ssp_.T();
    casadi_assert((H_ + total).nnz() == total.nnz(),
      "HPMPC: specified structure of H does not correspond to what the interface can handle. "
      "Structure is: N " + str(N_) + ", nx " + str(nx) + ", nu " + str(nu) + ", "
      "ng " + str(ng) + ".");
    casadi_assert_dev(total.nnz() == Rsp_.nnz() + 2*Ssp_.nnz() + Qsp_.nnz());

    /* Disassemble LBA/UBA input into:
       b
       lg/ug

       b
       lg/ug
    */
    offset = 0;

    for (casadi_int k=0;k<N_;++k) {
      b_blocks.push_back({offset,   0, nx[k+1], 1}); offset+= nx[k+1];
      lug_blocks.push_back({offset, 0, ng[k], 1}); offset+= ng[k];
    }
    lug_blocks.push_back({offset, 0, ng[N_], 1});

    bsp_ = blocksparsity(na_, 1, b_blocks);
    lugsp_ = blocksparsity(na_, 1, lug_blocks);
    total = bsp_ + lugsp_;
    casadi_assert_dev(total.nnz() == bsp_.nnz() + lugsp_.nnz());
    casadi_assert_dev(total.nnz() == na_);

    /* Disassemble G/X0 input into:
       r/u
       q/x

       r/u
       q/x
    */
    offset = 0;

    for (casadi_int k=0;k<N_;++k) {
      x_blocks.push_back({offset, 0, nx[k], 1}); offset+= nx[k];
      u_blocks.push_back({offset, 0, nu[k], 1}); offset+= nu[k];
    }
    x_blocks.push_back({offset, 0, nx[N_], 1});

    usp_ = blocksparsity(nx_, 1, u_blocks);
    xsp_ = blocksparsity(nx_, 1, x_blocks);
    total = usp_ + xsp_;
    casadi_assert_dev(total.nnz() == usp_.nnz() + xsp_.nnz());
    casadi_assert_dev(total.nnz() == nx_);

    std::vector< Block > theirs_u_blocks, theirs_x_blocks;
    offset = 0;

    for (casadi_int k=0;k<N_;++k) {
      theirs_u_blocks.push_back({offset, 0, nu[k], 1}); offset+= nu[k];
      theirs_x_blocks.push_back({offset, 0, nx[k], 1}); offset+= nx[k];
    }
    theirs_x_blocks.push_back({offset, 0, nx[N_], 1});

    theirs_usp_ = blocksparsity(nx_, 1, theirs_u_blocks);
    theirs_xsp_ = blocksparsity(nx_, 1, theirs_x_blocks);
    total = theirs_usp_ + theirs_xsp_;
    casadi_assert_dev(total.nnz() == theirs_usp_.nnz() + theirs_xsp_.nnz());
    casadi_assert_dev(total.nnz() == nx_);

    offset = 0;
    std::vector< Block > lamg_gap_blocks;
    for (casadi_int k=0;k<N_;++k) {
      lamg_gap_blocks.push_back({offset,       0, nx[k+1], 1});offset+= nx[k+1] + ng[k];
    }
    lamg_gapsp_ = blocksparsity(na_, 1, lamg_gap_blocks);
    lamg_csp_ = lamg_gapsp_.pattern_inverse();

    offset = 0;

    for (casadi_int k=0;k<N_;++k) {
      lam_ul_blocks.push_back({offset, 0, nu[k], 1}); offset+= nu[k];
      lam_xl_blocks.push_back({offset, 0, nx[k], 1}); offset+= nx[k];
      lam_uu_blocks.push_back({offset, 0, nu[k], 1}); offset+= nu[k];
      lam_xu_blocks.push_back({offset, 0, nx[k], 1}); offset+= nx[k];
      lam_cl_blocks.push_back({offset, 0, ng[k], 1}); offset+= ng[k];
      lam_cu_blocks.push_back({offset, 0, ng[k], 1}); offset+= ng[k];
    }
    lam_xl_blocks.push_back({offset, 0, nx[N_], 1}); offset+= nx[N_];
    lam_xu_blocks.push_back({offset, 0, nx[N_], 1}); offset+= nx[N_];
    lam_cl_blocks.push_back({offset, 0, ng[N_], 1}); offset+= ng[N_];
    lam_cu_blocks.push_back({offset, 0, ng[N_], 1}); offset+= ng[N_];

    lam_ulsp_ = blocksparsity(offset, 1, lam_ul_blocks);
    lam_uusp_ = blocksparsity(offset, 1, lam_uu_blocks);
    lam_xlsp_ = blocksparsity(offset, 1, lam_xl_blocks);
    lam_xusp_ = blocksparsity(offset, 1, lam_xu_blocks);
    lam_clsp_ = blocksparsity(offset, 1, lam_cl_blocks);
    lam_cusp_ = blocksparsity(offset, 1, lam_cu_blocks);

    pisp_ = Sparsity::dense(std::accumulate(nx.begin()+1, nx.end(), 0), 1);

    total = lam_ulsp_ + lam_uusp_ + lam_xlsp_ + lam_xusp_ + lam_clsp_ + lam_cusp_;
    casadi_assert_dev(total.nnz() == lam_ulsp_.nnz() + lam_uusp_.nnz() + lam_xlsp_.nnz() +
      lam_xusp_.nnz() + lam_clsp_.nnz() + lam_cusp_.nnz());
    casadi_assert_dev(total.nnz() == offset);

    theirs_Xsp_ = Sparsity::dense(std::accumulate(nx.begin(), nx.end(), 0), 1);
    theirs_Usp_ = Sparsity::dense(std::accumulate(nu.begin(), nu.end(), 0), 1);

  }


  int HpmpcInterface::init_mem(void* mem) const {
    if (Conic::init_mem(mem)) return 1;
    auto m = static_cast<HpmpcMemory*>(mem);

    init_vector(m->nx, nxs_);
    init_vector(m->nu, nus_); m->nu.push_back(0);
    init_vector(m->ng, ngs_);

    const std::vector<int>& nx = m->nx;
    const std::vector<int>& nu = m->nu;
    const std::vector<int>& ng = m->ng;
    const std::vector<int>& nb = m->nb;

    m->nb.resize(N_+1);
    casadi_int offset = 0;
    for (casadi_int k=0;k<N_;++k) m->nb[k] = nx[k]+nu[k];
    m->nb[N_] = nx[N_];

    m->A.resize(Asp_.nnz());
    m->B.resize(Bsp_.nnz());
    m->C.resize(Csp_.nnz());
    m->D.resize(Dsp_.nnz());
    m->R.resize(Rsp_.nnz());
    m->I.resize(Isp_.nnz());
    m->S.resize(Ssp_.nnz());
    m->Q.resize(Qsp_.nnz());
    m->b.resize(bsp_.nnz());
    m->b2.resize(bsp_.nnz());
    m->x.resize(xsp_.nnz());m->q.resize(xsp_.nnz());
    m->u.resize(usp_.nnz());m->r.resize(usp_.nnz());
    m->lg.resize(lugsp_.nnz());
    m->ug.resize(lugsp_.nnz());
    m->lb.resize(nx_);m->ub.resize(nx_);
    m->pi.resize(pisp_.nnz());

    offset = 0;
    for (casadi_int k=0;k<N_+1;++k) offset+=m->nx[k]+nu[k];
    m->hidxb.resize(offset);

    offset = 0;
    for (casadi_int k=0;k<N_;++k) offset+=ng[k]+nx[k]+nu[k];
    offset+=ng[N_]+nx[N_];
    m->lam.resize(2*offset);

    // Allocate double* work vectors
    blockptr(m->As, m->A, A_blocks);
    blockptr(m->Bs, m->B, B_blocks);
    blockptr(m->Cs, m->C, C_blocks);
    blockptr(m->Ds, m->D, D_blocks);
    blockptr(m->Is, m->I, I_blocks, true);
    blockptr(m->Rs, m->R, R_blocks);
    blockptr(m->Qs, m->Q, Q_blocks);
    blockptr(m->Ss, m->S, S_blocks);
    blockptr(m->us, m->u, u_blocks);
    blockptr(m->xs, m->x, x_blocks);
    blockptr(m->rs, m->r, u_blocks);
    blockptr(m->qs, m->q, x_blocks);
    blockptr(m->lgs, m->lg, lug_blocks);
    blockptr(m->ugs, m->ug, lug_blocks);
    blockptr(m->bs, m->b, b_blocks);

    m->pis.resize(N_);
    offset = 0;
    for (casadi_int k=0;k<N_;++k) {
      m->pis[k] = get_ptr(m->pi)+offset;
      offset+=nx[k+1];
    }

    m->lbs.resize(N_+1);
    m->ubs.resize(N_+1);
    offset = 0;
    for (casadi_int k=0;k<N_+1;++k) {
      m->lbs[k] = get_ptr(m->lb)+offset;
      m->ubs[k] = get_ptr(m->ub)+offset;
      offset+=nu[k]+nx[k];
    }

    m->lams.resize(N_+1);
    offset = 0;
    for (casadi_int k=0;k<N_+1;++k) {
      m->lams[k] = get_ptr(m->lam)+offset;
      offset+=2*(ng[k]+nb[k]);
    }

    m->hidxbs.resize(N_+1);
    offset = 0;
    for (casadi_int k=0;k<N_+1;++k) {
      m->hidxbs[k] = get_ptr(m->hidxb)+offset;
      for (casadi_int i=0;i<m->nb[k];++i) m->hidxbs[k][i] = i;
      offset+=nb[k];
    }

    // Allocate extra workspace as per HPMPC request
    int workspace_size = hpmpc_d_ip_ocp_hard_tv_work_space_size_bytes(
      N_, get_ptr(m->nx), get_ptr(m->nu), get_ptr(m->nb), get_ptr(m->hidxbs), get_ptr(m->ng), N_);
    m->workspace.resize(workspace_size);
    m->stats.resize(max_iter_*5);

    m->pv.resize(2*(nx_+na_));

    m->res.resize(4);

    m->fstats["preprocessing"]  = FStats();
    m->fstats["solver"]         = FStats();
    m->fstats["postprocessing"] = FStats();
    return 0;
  }

  void HpmpcInterface::mproject(double factor, const double* x, const casadi_int* sp_x,
                                double* y, const casadi_int* sp_y, double* w) {
    casadi_int ncol_y = sp_y[1];
    const casadi_int *colind_y = sp_y+2;
    casadi_project(x, sp_x, y, sp_y, w);
    casadi_scal(colind_y[ncol_y], factor, y);
  }

  void HpmpcInterface::dense_transfer(double factor, const double* x,
                                      const casadi_int* sp_x, double* y,
                                      const casadi_int* sp_y, double* w) {
    CASADI_PREFIX(sparsify)(x, w, sp_x, false);
    casadi_int nrow_y = sp_y[0];
    casadi_int ncol_y = sp_y[1];
    const casadi_int *colind_y = sp_y+2, *row_y = sp_y + 2 + ncol_y+1;
    /* Loop over columns of y */
    casadi_int i, el;
    for (i=0; i<ncol_y; ++i) {
      for (el=colind_y[i]; el<colind_y[i+1]; ++el) y[nrow_y*i + row_y[el]] += factor*(*w++);
    }
  }

  int HpmpcInterface::
  eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    if (inputs_check_) {
      check_inputs(arg[CONIC_LBX], arg[CONIC_UBX], arg[CONIC_LBA], arg[CONIC_UBA]);
    }

    auto m = static_cast<HpmpcMemory*>(mem);
    // Statistics
    for (auto&& s : m->fstats) s.second.reset();
    m->fstats.at("preprocessing").tic();

    double* pv =  get_ptr(m->pv);

    // Dissect A matrix
    casadi_project(arg[CONIC_A], A_, get_ptr(m->A), Asp_, pv);
    casadi_project(arg[CONIC_A], A_, get_ptr(m->B), Bsp_, pv);
    casadi_project(arg[CONIC_A], A_, get_ptr(m->C), Csp_, pv);
    casadi_project(arg[CONIC_A], A_, get_ptr(m->D), Dsp_, pv);
    casadi_project(arg[CONIC_A], A_, get_ptr(m->I), Isp_, pv);

    // Dissect H matrix; definition of HPMPC lacks a factor 2
    mproject(0.5, arg[CONIC_H], H_, get_ptr(m->R), Rsp_, pv);
    mproject(0.5, arg[CONIC_H], H_, get_ptr(m->S), Ssp_, pv);
    mproject(0.5, arg[CONIC_H], H_, get_ptr(m->Q), Qsp_, pv);

    // Dissect LBA/UBA
    mproject(-1.0, arg[CONIC_LBA], sparsity_in_.at(CONIC_LBA), get_ptr(m->b), bsp_, pv);
    mproject(-1.0, arg[CONIC_UBA], sparsity_in_.at(CONIC_UBA), get_ptr(m->b2), bsp_, pv);
    casadi_assert_dev(std::equal(m->b.begin(), m->b.end(), m->b2.begin()));
    casadi_project(arg[CONIC_LBA], sparsity_in_.at(CONIC_LBA), get_ptr(m->lg), lugsp_, pv);
    casadi_project(arg[CONIC_UBA], sparsity_in_.at(CONIC_UBA), get_ptr(m->ug), lugsp_, pv);

    // Dissect LBX/UBX input
    std::fill(m->lb.begin(), m->lb.end(), 0);
    std::fill(m->ub.begin(), m->ub.end(), 0);

    dense_transfer(1.0, arg[CONIC_LBX], xsp_, get_ptr(m->lb), theirs_xsp_, pv);
    dense_transfer(1.0, arg[CONIC_UBX], xsp_, get_ptr(m->ub), theirs_xsp_, pv);
    dense_transfer(1.0, arg[CONIC_LBX], usp_, get_ptr(m->lb), theirs_usp_, pv);
    dense_transfer(1.0, arg[CONIC_UBX], usp_, get_ptr(m->ub), theirs_usp_, pv);

    // Dissect G
    mproject(0.5, arg[CONIC_G], sparsity_in_.at(CONIC_G), get_ptr(m->r), usp_, pv);
    mproject(0.5, arg[CONIC_G], sparsity_in_.at(CONIC_G), get_ptr(m->q), xsp_, pv);

    // Dissect X0
    casadi_project(arg[CONIC_X0], sparsity_in_.at(CONIC_X0), get_ptr(m->u), usp_, pv);
    casadi_project(arg[CONIC_X0], sparsity_in_.at(CONIC_X0), get_ptr(m->x), xsp_, pv);

    m->iter_count = -1;

    // Deal with non-unity I block
    for (casadi_int k=0;k<N_;++k) {
      casadi_int n_row = m->nx[k+1];
      for (casadi_int i=0;i<n_row;++i) {
        double f = -1/m->Is[k][i];
        m->bs[k][i]*=f;
        for (casadi_int j=0;j<m->nx[k];++j) m->As[k][i+j*n_row]*=f;
        for (casadi_int j=0;j<m->nu[k];++j) m->Bs[k][i+j*n_row]*=f;
      }
    }

    // replace infinities
    for (casadi_int i=0;i<m->lb.size();++i) {
      if (m->lb[i]==-std::numeric_limits<double>::infinity()) m->lb[i] = -inf_;
    }
    for (casadi_int i=0;i<m->ub.size();++i) {
      if (m->ub[i]==std::numeric_limits<double>::infinity()) m->ub[i] = inf_;
    }
    for (casadi_int i=0;i<m->lg.size();++i) {
      if (m->lg[i]==-std::numeric_limits<double>::infinity()) m->lg[i] = -inf_;
    }
    for (casadi_int i=0;i<m->ug.size();++i) {
      if (m->ug[i]==std::numeric_limits<double>::infinity()) m->ug[i] = inf_;
    }

    m->fstats.at("preprocessing").toc();
    m->fstats.at("solver").tic();


    std::fill(m->pi.begin(), m->pi.end(), 0);
    std::fill(m->lam.begin(), m->lam.end(), 0);

    if (arg[CONIC_LAM_A0]) {
      dense_transfer(0.5, arg[CONIC_LAM_A0], lamg_gapsp_, get_ptr(m->pi), pisp_, pv);
      // Deal with non-unity I block
      for (casadi_int k=0;k<N_;++k) {
        casadi_int n_row = m->nx[k+1];
        for (casadi_int i=0;i<n_row;++i) {
          double f = -m->Is[k][i];
          m->pis[k][i]*=f;
        }
      }

      dense_transfer(0.5, arg[CONIC_LAM_A0], lamg_csp_, get_ptr(m->lam), lam_cusp_, pv);
    }

    if (arg[CONIC_LAM_X0]) {
      dense_transfer(0.5, arg[CONIC_LAM_X0], usp_, get_ptr(m->lam), lam_uusp_, pv);
      dense_transfer(0.5, arg[CONIC_LAM_X0], xsp_, get_ptr(m->lam), lam_xusp_, pv);
    }

    m->return_status =
      fortran_order_d_ip_ocp_hard_tv(&m->iter_count, max_iter_, mu0_, tol_, N_, get_ptr(m->nx),
      get_ptr(m->nu), get_ptr(m->nb), get_ptr(m->hidxbs), get_ptr(m->ng), N_, warm_start_,
      get_ptr(m->As), get_ptr(m->Bs), get_ptr(m->bs), get_ptr(m->Qs), get_ptr(m->Ss),
      get_ptr(m->Rs), get_ptr(m->qs), get_ptr(m->rs), get_ptr(m->lbs), get_ptr(m->ubs),
      get_ptr(m->Cs), get_ptr(m->Ds), get_ptr(m->lgs), get_ptr(m->ugs), get_ptr(m->xs),
      get_ptr(m->us), get_ptr(m->pis), get_ptr(m->lams), get_ptr(m->res),
      get_ptr(m->workspace), get_ptr(m->stats));

    m->fstats.at("solver").toc();
    m->fstats.at("postprocessing").tic();
    if (print_level_>0) {
      uout() << "HPMPC finished after " << m->iter_count << " iterations." << std::endl;
      uout() << "return status: " << m->return_status << std::endl;
      uout() << "residuals: " << m->res << std::endl;
    }

    std::fill(res[CONIC_X], res[CONIC_X]+nx_, 0);
    dense_transfer(1.0, get_ptr(m->x), theirs_Xsp_, res[CONIC_X], xsp_, pv);
    dense_transfer(1.0, get_ptr(m->u), theirs_Usp_, res[CONIC_X], usp_, pv);

    std::fill(res[CONIC_LAM_X], res[CONIC_LAM_X]+nx_, 0);
    std::fill(res[CONIC_LAM_A], res[CONIC_LAM_A]+na_, 0);

    // Deal with non-unity I block
    for (casadi_int k=0;k<N_;++k) {
      casadi_int n_row = m->nx[k+1];
      for (casadi_int i=0;i<n_row;++i) {
        double f = -1/m->Is[k][i];
        m->pis[k][i]*=f;
      }
    }

    dense_transfer(2.0, get_ptr(m->pi), pisp_, res[CONIC_LAM_A], lamg_gapsp_, pv);
    dense_transfer(2.0, get_ptr(m->lam), lam_cusp_, res[CONIC_LAM_A], lamg_csp_, pv);
    dense_transfer(-2.0, get_ptr(m->lam), lam_clsp_, res[CONIC_LAM_A], lamg_csp_, pv);

    dense_transfer(-2.0, get_ptr(m->lam), lam_ulsp_, res[CONIC_LAM_X], usp_, pv);
    dense_transfer(2.0, get_ptr(m->lam), lam_uusp_, res[CONIC_LAM_X], usp_, pv);
    dense_transfer(-2.0, get_ptr(m->lam), lam_xlsp_, res[CONIC_LAM_X], xsp_,  pv);
    dense_transfer(2.0, get_ptr(m->lam), lam_xusp_, res[CONIC_LAM_X], xsp_,  pv);

    // Construct f
    double f = casadi_dot(nx_, arg[CONIC_G], res[CONIC_X]);
    f += 0.5*casadi_bilin(arg[CONIC_H], H_, res[CONIC_X], res[CONIC_X]);

    if (res[CONIC_COST]) res[CONIC_COST][0] = f;

    m->fstats.at("postprocessing").toc();

    // Show statistics
    if (print_time_)  print_fstats(static_cast<ConicMemory*>(mem));
    return 0;
  }

  Dict HpmpcInterface::get_stats(void* mem) const {
    Dict stats = Conic::get_stats(mem);
    auto m = static_cast<HpmpcMemory*>(mem);
    stats["return_status"] = m->return_status;
    stats["iter_count"] = m->iter_count;
    stats["res"] = m->res;
    return stats;
  }

  HpmpcMemory::HpmpcMemory() {
  }

  HpmpcMemory::~HpmpcMemory() {

  }

  Sparsity HpmpcInterface::blocksparsity(casadi_int rows, casadi_int cols,
      const std::vector<Block>& blocks, bool eye) {
    DM r(rows, cols);
    for (auto && b : blocks) {
      if (eye) {
        r(range(b.offset_r, b.offset_r+b.rows),
          range(b.offset_c, b.offset_c+b.cols)) = DM::eye(b.rows);
        casadi_assert_dev(b.rows==b.cols);
      } else {
        r(range(b.offset_r, b.offset_r+b.rows),
        range(b.offset_c, b.offset_c+b.cols)) = DM::zeros(b.rows, b.cols);
      }
    }
    return r.sparsity();
  }
  void HpmpcInterface::blockptr(std::vector<double *>& vs, std::vector<double>& v,
      const std::vector<Block>& blocks, bool eye) {
    casadi_int N = blocks.size();
    vs.resize(N);
    casadi_int offset=0;
    for (casadi_int k=0;k<N;++k) {
      vs[k] = get_ptr(v)+offset;
      if (eye) {
        casadi_assert_dev(blocks[k].rows==blocks[k].cols);
        offset+=blocks[k].rows;
      } else {
        offset+=blocks[k].rows*blocks[k].cols;
      }
    }
  }
} // namespace casadi

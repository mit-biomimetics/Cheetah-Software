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

#ifndef CASADI_SX_HPP
#define CASADI_SX_HPP

/// \cond INTERNAL

namespace casadi {
  template<> CASADI_EXPORT bool Matrix<SXElem>::__nonzero__() const;
  template<> inline std::string matrixName<SXElem>() { return "SX"; }
  template<> SX SX::_sym(const std::string& name, const Sparsity& sp);
  template<> bool SX::is_regular() const;
  template<> bool SX::is_smooth() const;
  template<> bool SX::is_leaf() const;
  template<> bool SX::is_commutative() const;
  template<> bool SX::is_symbolic() const;
  template<> bool SX::is_valid_input() const;
  template<> bool SX::has_duplicates() const;
  template<> bool SX::is_op(casadi_int op) const;
  template<> casadi_int SX::op() const;
  template<> void SX::reset_input() const;
  template<> SX SX::dep(casadi_int ch) const;
  template<> casadi_int SX::n_dep() const;
  template<> std::string SX::name() const;
  template<> void SX::set_max_depth(casadi_int eq_depth);
  template<> casadi_int SX::get_max_depth();
  template<> casadi_int SX::element_hash() const;
  template<> void SX::expand(const SX& f, SX& weights, SX& terms);
  template<> SX SX::pw_const(const SX& t, const SX& tval, const SX& val);
  template<> SX SX::pw_lin(const SX& t, const SX &tval, const SX &val);
  template<> SX SX::gauss_quadrature(const SX& f, const SX &x, const SX &a,
                                     const SX &b, casadi_int order,
                                     const SX& w);
  template<> SX SX::simplify(const SX& x);
  template<> SX SX::substitute(const SX& ex, const SX& v, const SX& vdef);
  template<> std::vector<SX > SX::substitute(const std::vector<SX >& ex,
                                                const std::vector<SX >& v,
                                                const std::vector<SX >& vdef);
  template<> void SX::substitute_inplace(const std::vector<SX >& v,
                                        std::vector<SX >& vdef,
                                        std::vector<SX >& ex,
                                        bool reverse);
  template<> bool SX::depends_on(const SX &x, const SX &arg);
  template<> std::vector<SX > SX::symvar(const SX &x);
  template<> SX SX::jacobian(const SX &f, const SX &x, const Dict& opts);
  template<> SX SX::hessian(const SX &f, const SX &x);
  template<> SX SX::hessian(const SX &f, const SX &x, SX &g);
  template<> std::vector<std::vector<SX> >
  SX::forward(const std::vector<SX> &ex, const std::vector<SX> &arg,
          const std::vector<std::vector<SX> > &v, const Dict& opts);
  template<> std::vector<std::vector<SX> >
  SX::reverse(const std::vector<SX> &ex, const std::vector<SX> &arg,
          const std::vector<std::vector<SX> > &v, const Dict& opts);
  template<> std::vector<bool> SX::which_depends(const SX &expr, const SX &var,
                                                  casadi_int order, bool tr);
  template<> SX SX::taylor(const SX& f, const SX& x, const SX& a, casadi_int order);
  template<> SX SX::mtaylor(const SX& f, const SX& x, const SX& a, casadi_int order);
  template<> SX SX::mtaylor(const SX& f, const SX& x, const SX& a, casadi_int order,
                            const std::vector<casadi_int>& order_contributions);
  template<> casadi_int SX::n_nodes(const SX& x);
  template<> std::string
  SX::print_operator(const SX& x, const std::vector<std::string>& args);
  template<> void SX::shared(std::vector<SX >& ex,
                                    std::vector<SX >& v,
                                    std::vector<SX >& vdef,
                                    const std::string& v_prefix,
                                    const std::string& v_suffix);
  template<> SX SX::poly_coeff(const SX& f, const SX& x);
  template<> SX SX::poly_roots(const SX& p);
  template<> SX SX::eig_symbolic(const SX& m);
  template<> void SX::print_split(std::vector<std::string>& nz,
                                 std::vector<std::string>& inter) const;

  template<> std::vector<SX> SX::get_input(const Function& f);
  template<> std::vector<SX> SX::get_free(const Function& f);

  // Templates instantiated in matrix.cpp
#ifndef CASADI_MATRIX_CPP
  extern template class Matrix<double>;
  extern template class Matrix<casadi_int>;
  extern template class Matrix<SXElem>;
#endif // CASADI_MATRIX_CPP
} // namespace casadi

/// \endcond

#endif // CASADI_SX_HPP

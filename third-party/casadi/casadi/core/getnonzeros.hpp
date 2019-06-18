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


#ifndef CASADI_GETNONZEROS_HPP
#define CASADI_GETNONZEROS_HPP

#include "mx_node.hpp"
#include <map>
#include <stack>

/// \cond INTERNAL

namespace casadi {
  /** \brief Get nonzeros of a matrix
      \author Joel Andersson
      \date 2013
  */
  class CASADI_EXPORT GetNonzeros : public MXNode {
  public:

    ///@{
    /// Create functions
    static MX create(const Sparsity& sp, const MX& x, const std::vector<casadi_int>& nz);
    static MX create(const Sparsity& sp, const MX& x, const Slice& s);
    static MX create(const Sparsity& sp, const MX& x, const Slice& inner, const Slice& outer);
    ///@}

    /// Constructor
    GetNonzeros(const Sparsity& sp, const MX& x);

    /// Destructor
    ~GetNonzeros() override {}

    /** \brief  Evaluate symbolically (MX) */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief Calculate forward mode directional derivatives */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /// Get an IM representation of a GetNonzeros or SetNonzeros node
    Matrix<casadi_int> mapping() const override;

    /// Get all the nonzeros
    virtual std::vector<casadi_int> all() const = 0;

    /** \brief Get the operation */
    casadi_int op() const override { return OP_GETNONZEROS;}

    /// Get the nonzeros of matrix
    MX get_nzref(const Sparsity& sp, const std::vector<casadi_int>& nz) const override;
  };

  class CASADI_EXPORT GetNonzerosVector : public GetNonzeros {
  public:
    /// Constructor
    GetNonzerosVector(const Sparsity& sp, const MX& x,
                      const std::vector<casadi_int>& nz) : GetNonzeros(sp, x), nz_(nz) {}

    /// Destructor
    ~GetNonzerosVector() override {}

    /// Get all the nonzeros
    std::vector<casadi_int> all() const override { return nz_;}

    /** \brief  Propagate sparsity forward */
    int sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Propagate sparsity backwards */
    int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /// Evaluate the function (template)
    template<typename T>
    int eval_gen(const T* const* arg, T* const* res, casadi_int* iw, T* w) const;

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /// Evaluate the function symbolically (SX)
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const override;

    /** \brief  Print expression */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief Generate code for the operation */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res) const override;

    /** \brief Check if two nodes are equivalent up to a given depth */
    bool is_equal(const MXNode* node, casadi_int depth) const override;

    /** Obtain information about node */
    Dict info() const override { return {{"nz", nz_}}; }

    /// Operation sequence
    std::vector<casadi_int> nz_;
  };

  // Specialization of the above when nz_ is a Slice
  class CASADI_EXPORT GetNonzerosSlice : public GetNonzeros {
  public:

    /// Constructor
    GetNonzerosSlice(const Sparsity& sp, const MX& x, const Slice& s) : GetNonzeros(sp, x), s_(s) {}

    /// Destructor
    ~GetNonzerosSlice() override {}

    /// Get all the nonzeros
    std::vector<casadi_int> all() const override { return s_.all(s_.stop);}

    /** \brief  Propagate sparsity forward */
    int sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Propagate sparsity backwards */
    int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /// Evaluate the function (template)
    template<typename T>
    int eval_gen(const T* const* arg, T* const* res, casadi_int* iw, T* w) const;

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /// Evaluate the function symbolically (SX)
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const override;

    /** \brief  Print expression */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief Generate code for the operation */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res) const override;

    /** \brief Check if two nodes are equivalent up to a given depth */
    bool is_equal(const MXNode* node, casadi_int depth) const override;

    /** Obtain information about node */
    Dict info() const override { return {{"slice", s_.info()}}; }

    // Data member
    Slice s_;
  };

  // Specialization of the above when nz_ is a nested Slice
  class CASADI_EXPORT GetNonzerosSlice2 : public GetNonzeros {
  public:

    /// Constructor
    GetNonzerosSlice2(const Sparsity& sp, const MX& x, const Slice& inner,
                      const Slice& outer) : GetNonzeros(sp, x), inner_(inner), outer_(outer) {}

    /// Destructor
    ~GetNonzerosSlice2() override {}

    /// Get all the nonzeros
    std::vector<casadi_int> all() const override { return inner_.all(outer_, outer_.stop);}

    /** \brief  Propagate sparsity forward */
    int sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Propagate sparsity backwards */
    int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /// Evaluate the function (template)
    template<typename T>
    int eval_gen(const T* const* arg, T* const* res, casadi_int* iw, T* w) const;

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /// Evaluate the function symbolically (SX)
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const override;

    /** \brief  Print expression */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief Generate code for the operation */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res) const override;

    /** \brief Check if two nodes are equivalent up to a given depth */
    bool is_equal(const MXNode* node, casadi_int depth) const override;

    /** Obtain information about node */
    Dict info() const override { return {{"inner", inner_.info()}, {"outer", outer_.info()}}; }

    // Data members
    Slice inner_, outer_;
  };


} // namespace casadi
/// \endcond

#endif // CASADI_GETNONZEROS_HPP

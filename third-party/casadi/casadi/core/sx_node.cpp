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


#include "sx_node.hpp"
#include <limits>
#include <stack>

using namespace std;
namespace casadi {

  SXNode::SXNode() {
    count = 0;
    temp = 0;
  }

  SXNode::~SXNode() {
    #ifdef WITH_REFCOUNT_WARNINGS
    // Make sure that this is there are no scalar expressions pointing to it when it is destroyed
    if (count!=0) {
      // Note that casadi_assert_warning cannot be used in destructors
      std::cerr << "Reference counting failure." <<
                   "Possible cause: Circular dependency in user code." << std::endl;
    }
    #endif // WITH_REFCOUNT_WARNINGS
  }

  double SXNode::to_double() const {
    return numeric_limits<double>::quiet_NaN();
  }

  casadi_int SXNode::to_int() const {
    casadi_error("to_int not defined for " + class_name());
  }

  bool SXNode::is_equal(const SXNode* node, casadi_int depth) const {
    return false;
  }

  const std::string& SXNode::name() const {
    casadi_error("'name' not defined for " + class_name());
  }

  const SXElem& SXNode::dep(casadi_int i) const {
    casadi_error("'dep' not defined for " + class_name());
  }

  SXElem& SXNode::dep(casadi_int i) {
    casadi_error("'dep' not defined for " + class_name());
  }

  void SXNode::disp(std::ostream& stream, bool more) const {
    // Find out which noded can be inlined
    std::map<const SXNode*, casadi_int> nodeind;
    can_inline(nodeind);

    // Print expression
    vector<string> intermed;
    string s = print_compact(nodeind, intermed);

    // Print intermediate expressions
    for (casadi_int i=0; i<intermed.size(); ++i)
      stream << "@" << (i+1) << "=" << intermed[i] << ", ";

    // Print this
    stream << s;
  }

  bool SXNode::marked() const {
    return temp<0;
  }

  void SXNode::mark() const {
    temp = -temp-1;
  }

  void SXNode::can_inline(std::map<const SXNode*, casadi_int>& nodeind) const {
    // Add or mark node in map
    std::map<const SXNode*, casadi_int>::iterator it=nodeind.find(this);
    if (it==nodeind.end()) {
      // First time encountered, mark inlined
      nodeind.insert(it, make_pair(this, 0));

      // Handle dependencies with recursion
      for (casadi_int i=0; i<n_dep(); ++i) {
        dep(i)->can_inline(nodeind);
      }
    } else if (it->second==0 && op()!=OP_PARAMETER) {
      // Node encountered before, do not inline (except if symbolic primitive)
      it->second = -1;
    }
  }

  std::string SXNode::print_compact(std::map<const SXNode*, casadi_int>& nodeind,
                                   std::vector<std::string>& intermed) const {
    // Get reference to node index
    casadi_int& ind = nodeind[this];

    // If positive, already in intermediate expressions
    if (ind>0) {
      stringstream ss;
      ss << "@" << ind;
      return ss.str();
    }

    // Get expressions for dependencies
    std::string arg[2];
    for (casadi_int i=0; i<n_dep(); ++i) {
      arg[i] = dep(i)->print_compact(nodeind, intermed);
    }

    // Get expression for this
    string s = print(arg[0], arg[1]);

    // Decide what to do with the expression
    if (ind==0) {
      // Inline expression
      return s;
    } else {
      // Add to list of intermediate expressions and return reference
      intermed.push_back(s);
      ind = intermed.size(); // For subsequent references
      stringstream ss;
      ss << "@" << ind;
      return ss.str();
    }
  }

  void SXNode::safe_delete(SXNode* n) {
    // Quick return if more owners
    if (n->count>0) return;
    // Delete straight away if it doesn't have any dependencies
    if (!n->n_dep()) {
      delete n;
      return;
    }
    // Stack of expressions to be deleted
    std::stack<SXNode*> deletion_stack;
    // Add the node to the deletion stack
    deletion_stack.push(n);
    // Process stack
    while (!deletion_stack.empty()) {
      // Top element
      SXNode *t = deletion_stack.top();
      // Check if the top element has dependencies with dependencies
      bool added_to_stack = false;
      for (casadi_int c2=0; c2<t->n_dep(); ++c2) { // for all dependencies of the dependency
        // Get the node of the dependency of the top element
        // and remove it from the smart pointer
        SXNode *n2 = t->dep(c2).assignNoDelete(casadi_limits<SXElem>::nan);
        // Check if this is the only reference to the element
        if (n2->count == 0) {
          // Check if unary or binary
          if (!n2->n_dep()) {
            // Delete straight away if not binary
            delete n2;
          } else {
            // Add to deletion stack
            deletion_stack.push(n2);
            added_to_stack = true;
          }
        }
      }
      // Delete and pop from stack if nothing added to the stack
      if (!added_to_stack) {
        delete deletion_stack.top();
        deletion_stack.pop();
      }
    }
  }

  casadi_int SXNode::eq_depth_ = 1;

} // namespace casadi

/*------------------------------------------------------------------------*/
/*! \file int_polynomial.hpp
    \brief Simplified polynomial with integer coefficients and Boolean algebra
           Implements x*x = x property for nonlinear polynomials.
           All substitution and simplification done here, then converted to
           multimod Polynomial only for specification reduction.

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/

#ifndef TALISMAN_SRC_INT_POLYNOMIAL_HPP
#define TALISMAN_SRC_INT_POLYNOMIAL_HPP

#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "gate.hpp"
#include "polynomial.hpp"

namespace talisman {

// Simplified monomial for integer polynomial
// (Same structure as regular Monomial)
using IntMonomial = Monomial;

// Integer polynomial: monomials with int32_t coefficients
// Property: x*x = x (idempotent Boolean algebra)
class IntPolynomial {
  public:
  using IntTermMap = std::map<IntMonomial, int32_t>;

  IntTermMap terms;

  IntPolynomial() = default;
  explicit IntPolynomial(Gate* gate);
  explicit IntPolynomial(int32_t coeff);

  // Add a term: if monomial exists, add coefficient; otherwise insert
  void add_term(const IntMonomial& mono, int32_t coeff);

  // Multiply by another polynomial
  IntPolynomial operator*(const IntPolynomial& other) const;

  // Multiply in place
  void mul_inplace(const IntPolynomial& other);

  // Check if polynomial is zero
  bool is_zero() const { return terms.empty(); }

  // Get size (number of terms)
  size_t get_size() const { return terms.size(); }

  // Print polynomial
  void print() const;
  std::string to_string() const;

  // Convert to multimodular Polynomial for specification reduction
  Polynomial to_multimod() const;

  // Extract substitution expression: if this = x - expr, returns expr
  IntPolynomial extract_substitution_expr(Gate* x) const;

  // Create gate polynomial: g - child0_poly * child1_poly
  static IntPolynomial get_gate_poly(Gate* g);
};

}  // namespace talisman

#endif  // TALISMAN_SRC_INT_POLYNOMIAL_HPP

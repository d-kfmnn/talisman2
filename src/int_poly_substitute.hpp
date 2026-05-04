/*------------------------------------------------------------------------*/
/*! \file int_poly_substitute.hpp
    \brief Gate polynomial substitution using IntPolynomial

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/

#ifndef TALISMAN_SRC_INT_POLY_SUBSTITUTE_HPP
#define TALISMAN_SRC_INT_POLY_SUBSTITUTE_HPP

#include "gate.hpp"
#include "int_polynomial.hpp"

namespace talisman {

// Helper: Substitute all occurrences of gate n in a polynomial with its expression
// Performs the full substitution and returns the result
// The substitution_poly should be extracted from n's polynomial using
// extract_substitution_expr(n, n_poly) when n_poly = n - expr
IntPolynomial substitute_gate(Gate* n, const IntPolynomial& parent_poly,
                              const IntPolynomial& substitution_expr);

}  // namespace talisman

#endif  // TALISMAN_SRC_INT_POLY_SUBSTITUTE_HPP

/*------------------------------------------------------------------------*/
/*! \file int_polynomial.cpp
    \brief Implementation of IntPolynomial with integer coefficients
           and Boolean algebra properties (x*x = x)

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/

#include "int_polynomial.hpp"

#include <algorithm>
#include <iostream>
#include <sstream>

using namespace talisman;

namespace talisman {

IntPolynomial::IntPolynomial(Gate* gate) {
  IntMonomial mono;
  mono.vars.push_back(gate);
  terms[mono] = 1;
}

IntPolynomial::IntPolynomial(int32_t coeff) {
  if(coeff != 0) {
    IntMonomial mono;  // Empty monomial = constant
    terms[mono] = coeff;
  }
}

void IntPolynomial::add_term(const IntMonomial& mono, int32_t coeff) {
  if(coeff == 0)
    return;

  auto it = terms.find(mono);
  if(it == terms.end()) {
    terms[mono] = coeff;
  } else {
    it->second += coeff;
    // Remove term if coefficient becomes zero
    if(it->second == 0) {
      terms.erase(it);
    }
  }
}

IntPolynomial IntPolynomial::operator*(const IntPolynomial& other) const {
  IntPolynomial result;

  for(const auto& [mono1, coeff1] : terms) {
    for(const auto& [mono2, coeff2] : other.terms) {
      // Merge monomials: concatenate variables
      IntMonomial new_mono;
      new_mono.vars.insert(new_mono.vars.end(), mono1.vars.begin(),
                           mono1.vars.end());
      new_mono.vars.insert(new_mono.vars.end(), mono2.vars.begin(),
                           mono2.vars.end());

      // Sort and apply Boolean algebra: x*x = x (remove duplicates)
      new_mono.sort_vars();
      auto last = std::unique(new_mono.vars.begin(), new_mono.vars.end());
      new_mono.vars.erase(last, new_mono.vars.end());

      // Multiply coefficients
      int32_t new_coeff = coeff1 * coeff2;

      result.add_term(new_mono, new_coeff);
    }
  }

  return result;
}

void IntPolynomial::mul_inplace(const IntPolynomial& other) {
  IntPolynomial result = *this * other;
  terms = std::move(result.terms);
}

void IntPolynomial::print() const {
  std::cout << to_string() << std::endl;
}

std::string IntPolynomial::to_string() const {
  std::stringstream ss;
  bool first = true;

  for(const auto& [mono, coeff] : terms) {
    if(!first) {
      if(coeff > 0)
        ss << " + ";
      else
        ss << " - ";
    }

    ss << "{" << std::abs(coeff) << "}";

    for(const auto& gate : mono.vars) {
      ss << "*" << gate->get_var_name();
    }

    first = false;
  }

  if(first)
    ss << "0";

  return ss.str();
}

Polynomial IntPolynomial::to_multimod() const {
  Polynomial result;

  for(const auto& [mono, coeff] : terms) {
    // Create a CoeffVec from the integer coefficient
    CoeffVec cv;
    cv.set_constant(coeff);
    result.add_term(mono, cv);
  }

  return result;
}

IntPolynomial IntPolynomial::extract_substitution_expr(Gate* x) const {
  IntPolynomial substitution;

  for(const auto& [mono, coeff] : terms) {
    // Skip the monomial that is just gate x itself
    if(mono.vars.size() == 1 && mono.vars[0] == x) {
      continue;
    }

    // Add all other terms with negated coefficient
    substitution.add_term(mono, -coeff);
  }

  return substitution;
}

// Helper: Convert a literal (aiger format) to an IntPolynomial
// Positive lit = gate(lit), Negative lit = 1 - gate(lit)
IntPolynomial literal_to_int_poly(unsigned lit);

IntPolynomial IntPolynomial::get_gate_poly(Gate* g) {
  if(!g)
    throw std::runtime_error("IntPolynomial::get_gate_poly: null gate");

  // Create gate polynomial: g
  IntPolynomial gout(g);

  // Child polynomials respecting signedness/complementation
  IntPolynomial c0 = literal_to_int_poly(g->get_ch0());
  IntPolynomial c1 = literal_to_int_poly(g->get_ch1());

  // Multiply children: c0 * c1
  IntPolynomial prod = c0 * c1;

  // Gate polynomial: g - c0*c1
  // We need to subtract prod from gout
  IntPolynomial result = gout;
  for(const auto& [mono, coeff] : prod.terms) {
    result.add_term(mono, -coeff);
  }

  return result;
}

// Helper: Convert a literal (aiger format) to an IntPolynomial
IntPolynomial literal_to_int_poly(unsigned lit) {
  if(lit == 0)
    return IntPolynomial(0);  // Constant 1
  if(lit == 1)
    return IntPolynomial(1);  // Constant -1 (logical NOT of constant 1)

  unsigned node_id = lit / 2;
  Gate* g = gate(node_id * 2);  // Get gate from node id
  if(!g)
    throw std::runtime_error("literal_to_int_poly: invalid gate");

  IntPolynomial result(g);
  
  if(lit & 1) {
    // Complemented literal: NOT g = 1 - g
    // Create the negation: 1 - g
    IntPolynomial negated;
    
    // First add constant 1
    IntMonomial const_mono;
    negated.add_term(const_mono, 1);
    
    // Then subtract all terms from g
    for(const auto& [mono, coeff] : result.terms) {
      negated.add_term(mono, -coeff);
    }
    
    result = negated;
  }

  return result;
}

}  // namespace talisman

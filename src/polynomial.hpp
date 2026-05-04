// Polynomial.hpp

#ifndef POLYNOMIAL_HPP
#define POLYNOMIAL_HPP

#include <algorithm>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "coeff_vec.hpp"
#include "gate.hpp"
#include "linpolynomial.hpp"

namespace talisman {

// In Monomial.hpp
struct Monomial {
  std::vector<Gate*> vars;

  Monomial() = default;
  Monomial(std::vector<Gate*> v);

  // Sort variables (decreasing level, tie by name)
  void sort_vars();

  // Comparison operator for std::map ordering
  bool operator<(const Monomial& other) const;

  void print() const;

  // Extract remainder: all variables except the leading one
  Monomial get_remainder() const;
};

struct MonomialHash {
  std::size_t operator()(const Monomial& m) const noexcept {
    std::size_t h = 0;
    for(const auto& g : m.vars) {
      if(g) {
        h ^= std::hash<const void*>{}(g) + 0x9e3779b97f4a7c15ULL + (h << 6)
             + (h >> 2);
      }
    }
    return h;
  }
};

class Polynomial {
  public:
  using MultimodTermMap = std::map<Monomial, CoeffVec>;

  MultimodTermMap multimod_terms;

  Polynomial() = default;
  explicit Polynomial(Gate* gate);
  explicit Polynomial(int coeff);

  // Convert from LinPolynomial
  static Polynomial from_linpoly(const LinPolynomial& linpoly);

  // Multimodular operations
  void add_term(const Monomial& mono, int32_t coeff, int idx);
  void add_term(const Monomial& mono, CoeffVec& coeff);
  void addmul_term(const Monomial& mono,
                   const CoeffVec& multiplier,
                   const CoeffVec& coeff);
  // Multiply monomial by both a monomial remainder and a coefficient multiplier
  void addmul_term(const Monomial& mono,
                   const Monomial& remainder,
                   const CoeffVec& multiplier,
                   const CoeffVec& coeff);
  void reduce_by_multimod(const Polynomial& p2);

  void mul_inplace(const Polynomial& other);

  bool is_zero() const;

  size_t get_size() const {
    return multimod_terms.size();
  }

  std::string to_string() const;
  void print() const;
  void print_multimod() const;

  bool check_leading_coeff();
  // Check if any coefficient equals 2^NN modulo all primes
  void check_coeffs(bool remove = false);

  const Monomial& get_lm() const {
    if(multimod_terms.empty())
      throw std::runtime_error("Polynomial has no terms");
    return multimod_terms.begin()->first;
  }

  // Printing
  friend std::ostream& operator<<(std::ostream& os, const Polynomial& p) {
    os << p.to_string();
    return os;
  }

  // Operators needed for gate polynomial generation (g - c0*c1)
  friend Polynomial operator*(Polynomial a, const Polynomial& b);
  friend Polynomial operator-(Polynomial a, const Polynomial& b);

  // Polynomial reduction: reduce this polynomial by dividing by p2
  // p2's leading term must contain only a single variable
  void reduce_by(const Polynomial& p2);

  // Build nonlinear gate polynomial for AND-like gate semantics:
  // get_gate_poly(g) = g - child0_poly * child1_poly
  // where child_poly(lit) = x if lit is positive edge; (1 - x) if lit is
  // complemented; and constants: lit==1 -> 1, lit==0 -> 0. If multimod is true,
  // creates a multimodular polynomial with CoeffVec coefficients.
  static Polynomial get_gate_poly(Gate* g);

  // Helper to convert an AIG literal to a polynomial
  static Polynomial literal_to_poly(unsigned lit);

  // Reduce gate polynomial by children that have only one parent
  // Iteratively reduces until no more reductions are possible
  void reduce_single_children();
};

void
test_polynomial_operations();

void translate_all_stored_gate_polys();

}// namespace talisman

#endif// POLYNOMIAL_HPP

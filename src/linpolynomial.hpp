/*------------------------------------------------------------------------*/
/*! \file linpolynomial.hpp
    \brief contains the class LinPolynomial and further functions to
    organize the linear polynomial, such as initializing the polynomial
  constraints

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#ifndef TALISMAN_SRC_LINPOLYNOMIAL_H_
#define TALISMAN_SRC_LINPOLYNOMIAL_H_
/*------------------------------------------------------------------------*/
#include <map>
#include <vector>

#include "aig.hpp"
#include "coeff_vec.hpp"
#include "fast_mod.hpp"
#include "gate.hpp"
/*------------------------------------------------------------------------*/

namespace talisman {

// Forward declaration: reduce a value modulo (p,c) into the same
// representative used for stored coefficients.
int32_t
mod(int64_t value, int32_t p, int32_t c);

struct GateLevelComparator {
  bool operator()(const Gate* a, const Gate* b) const {
    // Constants (nullptr) come last (lowest level)
    if(a == nullptr)
      return false;
    if(b == nullptr)
      return true;

    // Sort descending by level
    return a->get_level() > b->get_level();
  }
};

// images of +-k*2^NN (k = 1..1024) modulo all PRIMES
extern std::vector<CoeffVec> two_pow_NN_imgs;

void
init_two_pow_NN_imgs();

/*------------------------------------------------------------------------*/

class LinPolynomial {
  private:
  // coeffs are 32-bit integers, mapped to Gate* as keys
  std::map<Gate*, int32_t, GateLevelComparator> terms;
  static int32_t p;// Prime modulus for the field, if c= 0 we calculate mod p,
                   // else mod 2^p-c
  static int32_t c;// Prime modulus for the field
  static int32_t two_pow_nn_mod;// cached value of 2^NN modulo current modulus

  int32_t mod(int64_t value) const;

  public:
  bool multimod = false;
  std::map<Gate*, CoeffVec, GateLevelComparator> multimod_terms;

  LinPolynomial();
  LinPolynomial(std::vector<int32_t> coeffs,
                std::vector<Gate*> gates,
                bool raw = false);
  ~LinPolynomial();

  static void set_modulus(int32_t p_, int32_t c_ = 0) {
    if(p_ < 0)
      die(100, "Modulus must be positive");
    if(c_ != 0 && c_ != 1)
      die(100, "c must be 0 or 1 for Mersenne prime");
    p = p_;
    c = c_;
    // Precompute 2^NN modulo the chosen modulus and cache it
    if(p == 0) {
      two_pow_nn_mod = 0;
    } else {
      // Prefer FLINT's optimized 2^exp mod p implementation when a modulus
      // is configured. Fall back to local logic only when FLINT isn't
      // available (p==0 handled above).
      nmod_t nmod;
      nmod_init(&nmod, static_cast<ulong>(p));
      uint32_t raw_u
        = flint_nmod_2_pow_ui(static_cast<unsigned long>(NN), nmod, false);
      int32_t raw = static_cast<int32_t>(raw_u);
      // Normalize to the same representative used by the free `mod` helper
      two_pow_nn_mod = ::talisman::mod(static_cast<int64_t>(raw), p, c);
    }
  };// Set the prime modulus

  static int32_t get_two_pow_nn_mod() { return two_pow_nn_mod; }

  static int32_t get_modulus() { return c == 0 ? p : ((int32_t)2 << p) - c; }
  static int32_t get_mersenne_constant() { return c; }

  unsigned get_size() const { return terms.size(); }
  std::vector<std::pair<Gate*, int32_t>> get_terms() const;
  std::vector<Gate*> get_gates() const;
  std::vector<int32_t> get_literals() const;
  int32_t get_constant()
    const;// Get the constant term (coefficient of null gate)
  std::vector<int32_t> get_coeffs() const;
  std::vector<int64_t> get_coeffs_64() const;

  int minimize_coeffs(int p);

  int64_t get_max_eval() const;// Get the maximum value of coefficients

  int32_t get_coeff(Gate* gate) const;
  std::pair<Gate*, int32_t> get_leading_term() const;

  void print() const;
  void print_multimod() const;

  bool check_leading_coeff();
  void check_coeffs(bool remove = false);

  void add_term(Gate* gate, int64_t coeff);
  void add_term(Gate* gate, int32_t coeff, int idx);
  void add_term(Gate* gate, CoeffVec& coeff);
  void addmul_term(Gate* gate,
                   const CoeffVec& multiplier,
                   const CoeffVec& coeff);

  void add_polynomial(const LinPolynomial& other);
  void subtract_polynomial(const LinPolynomial& other);
  void multiply_by_constant(int32_t scalar);
  // Reduce all coefficients modulo the configured modulus and remove zero
  // terms
  void reduce_mod();
  void my_mod(int p);

  // void normalize();
  void reduce_by(const LinPolynomial& other);

  LinPolynomial operator+(const LinPolynomial& other) const;
  LinPolynomial operator-(const LinPolynomial& other) const;
  LinPolynomial operator*(int32_t scalar) const;
  friend LinPolynomial operator*(int32_t scalar, const LinPolynomial& poly);
  // Return a new polynomial with all coefficients reduced modulo `modulus`.
  // LinPolynomial operator%(int32_t modulus) const;
};

int32_t
mod(int64_t value, int32_t p, int32_t c);

void
test_linpolynomial();

// When true (default), coefficients equal to 2^NN (and their power-of-two
// multiples) modulo the configured modulus are removed during reduction.
extern bool enable_two_pow_nn_pruning;

}// namespace talisman

#endif// TALISMAN_SRC_LINPOLYNOMIAL_H_

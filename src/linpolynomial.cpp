/*------------------------------------------------------------------------*/
/*! \file linpolynomial.cpp
    \brief contains functions necessary to initialize the linear polynomial

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#include "linpolynomial.hpp"

#include <algorithm>// for std::find
#include <iostream>
#include <nmod.h>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "aig.hpp"
#include "fast_mod.hpp"
#include "signal_statistics.hpp"

using namespace talisman;

namespace talisman {

/*------------------------------------------------------------------------*/
int32_t LinPolynomial::p = 0;             // default: no modulus
int32_t LinPolynomial::c = 0;             // default: no modulus
int32_t LinPolynomial::two_pow_nn_mod = 0;// cached 2^NN modulo current modulus
/*------------------------------------------------------------------------*/
LinPolynomial::LinPolynomial() {
  // empty constructor
}
/*------------------------------------------------------------------------*/
LinPolynomial::LinPolynomial(std::vector<int32_t> coeffs,
                             std::vector<Gate*> gates,
                             bool raw) {
  if(coeffs.size() != gates.size()) {
    die(100, "Coeffs and gates must have the same size");
    return;
  }

  // check that all gates are distinct
  for(size_t i = 0; i < gates.size(); ++i) {
    for(size_t j = i + 1; j < gates.size(); ++j) {
      if(gates[i] == gates[j]) {
        die(100, "Gates must be distinct");
      }
    }
  }

  multimod = false;

  for(size_t i = 0; i < coeffs.size(); ++i) {
    if(coeffs[i] != 0) {
      if(raw)
        terms[gates[i]] = coeffs[i];
      else
        terms[gates[i]] = mod(coeffs[i]);
    }
  }
}
/*------------------------------------------------------------------------*/
LinPolynomial::~LinPolynomial() {
  terms.clear();
}
/*------------------------------------------------------------------------*/

std::vector<int32_t>
LinPolynomial::get_coeffs() const {
  std::vector<int32_t> coeffs;
  for(const auto& term : terms) {
    coeffs.push_back(term.second);
  }
  return coeffs;
}
/*------------------------------------------------------------------------*/
std::vector<int64_t>
LinPolynomial::get_coeffs_64() const {
  std::vector<int64_t> coeffs;

  for(const auto& term : terms) {
    if(term.first != nullptr) {
      coeffs.push_back(static_cast<int64_t>(term.second));
    }
  }
  return coeffs;
}

/*------------------------------------------------------------------------*/
thread_local std::vector<slong> coeffs;

int
LinPolynomial::minimize_coeffs(int p) {
  int n = terms.size();
  coeffs.clear();
  coeffs.reserve(n);

  int64_t best_sum = 0;
  // Copy coefficients to vector in map order
  for(const auto& [key, value] : terms) {
    slong c = value >= 0 ? value : value + p;
    coeffs.push_back(c);
    best_sum += std::abs(value);
  }

  if(best_sum < 16)
    return best_sum;

  slong best_k = 1;

  // try scalar multiples
  for(slong k = 2; k <= p / 2; k++) {
    int64_t new_sum = 0;
    for(int i = 0; i < n; i++) {
      slong c = (coeffs[i] * k) % p;
      new_sum += c < p / 2 ? c : p - c;
    }
    if(new_sum < best_sum) {
      best_sum = new_sum;
      best_k = k;
    }
  }

  if(best_k > 1) {
    // update the map
    int i = 0;
    for(auto& [key, value] : terms) {
      slong c = (coeffs[i++] * best_k) % p;
      value = c < p / 2 ? c : c - p;
    }
  }
  return best_sum;
}
/*------------------------------------------------------------------------*/

std::vector<std::pair<Gate*, int32_t>>
LinPolynomial::get_terms() const {
  std::vector<std::pair<Gate*, int32_t>> term_list(terms.begin(), terms.end());
  return term_list;
}
/*------------------------------------------------------------------------*/
int32_t
LinPolynomial::get_coeff(Gate* gate) const {
  if(terms.find(gate) != terms.end()) {
    return terms.at(gate);
  }
  return 0;
}
/*------------------------------------------------------------------------*/
std::pair<Gate*, int32_t>
LinPolynomial::get_leading_term() const {
  return *terms.begin();
}
/*------------------------------------------------------------------------*/
std::vector<Gate*>
LinPolynomial::get_gates() const {
  std::vector<Gate*> gates;
  for(const auto& term : terms) {
    gates.push_back(term.first);
  }
  return gates;
}
/*------------------------------------------------------------------------*/
std::vector<int32_t>
LinPolynomial::get_literals() const {
  std::vector<int32_t> literals;
  for(const auto& term : terms) {
    if(term.first != nullptr) {
      int var = term.first->get_var();
      literals.push_back(var >> 1);
    }
  }
  return literals;
}
/*------------------------------------------------------------------------*/
int32_t
LinPolynomial::get_constant() const {
  // The constant term is the coefficient of the null gate (if it exists)
  auto it = terms.find(nullptr);
  if(it != terms.end()) {
    return it->second;
  }
  return 0;// No constant term
}

/*------------------------------------------------------------------------*/
int64_t
LinPolynomial::get_max_eval() const {
  // Find the maximum absolute value of coefficients
  int64_t max_eval = 0;
  for(const auto& term : terms) {
    int32_t coeff = term.second;

    if(__builtin_add_overflow(max_eval, coeff, &max_eval))
      die(7, "addition overflow");
    // max_eval += coeff;
  }
  return max_eval;
}
/*------------------------------------------------------------------------*/
std::vector<CoeffVec> two_pow_NN_imgs;

void
init_two_pow_NN_imgs() {
  // 2^NN
  two_pow_NN_imgs.emplace_back(NN);
  // -2^NN
  two_pow_NN_imgs.emplace_back(NN, -1);

  for(ulong k = 2; k <= 1024; ++k) {
    // k*2^NN
    two_pow_NN_imgs.push_back(two_pow_NN_imgs[0].scalar_mul(k));
    // -k*2^NN
    two_pow_NN_imgs.push_back(two_pow_NN_imgs[1].scalar_mul(k));
  }
}
/*------------------------------------------------------------------------*/
bool
LinPolynomial::check_leading_coeff() {

  auto lc = multimod_terms.begin()->second;

  if(std::any_of(two_pow_NN_imgs.begin(),
                 two_pow_NN_imgs.end(),
                 [&lc](const CoeffVec& v) { return lc == v; })) {

    multimod_terms.erase(multimod_terms.begin());
    terms.erase(terms.begin());
    return true;
  }
  return false;
}
/*------------------------------------------------------------------------*/

void
LinPolynomial::check_coeffs(bool remove) {

  if(!remove) {
    for(const auto& term : multimod_terms)
      if(std::any_of(two_pow_NN_imgs.begin(),
                     two_pow_NN_imgs.end(),
                     [&term](const CoeffVec& v) { return term.second == v; }))
        die(101, "spec coeff equals 2^NN modulo all primes for some gate");
  } else {
    std::vector<Gate*> to_erase;
    for(const auto& term : multimod_terms) {
      if(std::any_of(two_pow_NN_imgs.begin(),
                     two_pow_NN_imgs.end(),
                     [&term](const CoeffVec& v) { return term.second == v; }))
        to_erase.push_back(term.first);
    }
    for(Gate* g : to_erase) {
      if(g != nullptr)
        warning("Removing gate %s because coeff is multiple of 2^NN",
                g->get_var_name().c_str());
      else {
        warning("Removing constant term because coeff is multiple of 2^NN");
      }
      multimod_terms.erase(g);
      terms.erase(g);
    }
  }
}

/*------------------------------------------------------------------------*/
void
LinPolynomial::print() const {
  if(verbose < 1)
    return;// Only print when verbose >= 1

  std::vector<std::pair<Gate*, int32_t>> term_list(terms.begin(), terms.end());

  bool first = true;
  for(const auto& term : term_list) {
    Gate* gate = term.first;
    int32_t coeff = term.second;

    if(coeff == 0)
      continue;

    // Print sign

    if(coeff < 0) {
      std::cout << " - ";
      coeff = -coeff;// Make it positive for printing
    } else if(!first) {
      std::cout << " + ";
    }

    if((coeff != 1 && coeff != -1) || gate == nullptr) {
      std::cout << coeff;
      if(gate != nullptr)
        std::cout << "*";
    } else if(coeff == -1 && gate != nullptr) {
      std::cout << " - ";
    }

    // Print gate name or constant
    if(gate != nullptr) {
      std::cout << gate->get_var_name();
    }

    first = false;
  }

  if(first) {
    std::cout << "0";// If nothing printed, it's zero
  }

  std::cout << std::endl;
}
/*------------------------------------------------------------------------*/

void
LinPolynomial::print_multimod() const {
  if(verbose < 1)
    return;// Only print when verbose >= 1

  std::vector<std::pair<Gate*, CoeffVec>> term_list(multimod_terms.begin(),
                                                    multimod_terms.end());

  bool first = true;
  for(const auto& term : term_list) {
    Gate* gate = term.first;
    CoeffVec coeff = term.second;

    if(coeff.is_zero())
      continue;

    if(!first) {
      std::cout << " + ";
    }

    coeff.print();

    if(gate != nullptr) {
      std::cout << "*";
      std::cout << gate->get_var_name();
    }

    first = false;
  }

  if(first) {
    std::cout << "0";// If nothing printed, it's zero
  }

  std::cout << std::endl;
}
/*------------------------------------------------------------------------*/
// Polynomial arithmetic
/*------------------------------------------------------------------------*/
int32_t
mod(int64_t value,
    int32_t p,
    int32_t c) {// TODO change this to algorithms of fast_mod.hpp?
  if(p == 0 && c == 0)
    return value;// no mod
  if(value == 0)
    return 0;

  if(c == 0) {
    // Case 1: Standard mod p
    int32_t r = (value % p);
    if(r > p / 2) {
      if(__builtin_sub_overflow(r, p, &r))
        die(7, "subtraction overflow");
      // r -= p;
    }
    if(r < -p / 2) {
      if(__builtin_add_overflow(r, p, &r))
        die(7, "addition overflow");
      // r -= p;
    }

    return r;
  }

  if(c == 1) {
    // Case 2: Mersenne prime: mod (2^p - 1)
    uint64_t m = (1ULL << p) - 1;
    uint64_t r = static_cast<uint64_t>(value);
    while(r > m) {
      r = (r & m) + (r >> p);
    }
    if(r == m)
      r = 0;
    return static_cast<int32_t>(r);
  }

  assert(false && "Unsupported case: y must be 0 or 1");
  return 0;
}
/*------------------------------------------------------------------------*/
int32_t
LinPolynomial::mod(int64_t value) const {
  return ::mod(value, p, c);
}

/*------------------------------------------------------------------------*/
void
LinPolynomial::add_term(Gate* gate, int64_t coeff) {

  // coeff = mod(coeff);// We assume coeff is already in the correct field
  if(coeff == 0) {
    return;
  }

  auto it = terms.find(gate);
  if(it != terms.end()) {
    int64_t new_coeff = it->second + coeff;
    if(new_coeff == 0) {
      terms.erase(it);
    } else {
      it->second = new_coeff;
    }
  } else {
    terms[gate] = coeff;
  }
}
/*------------------------------------------------------------------------*/
// add term to multimod poly at index idx
void
LinPolynomial::add_term(Gate* gate, int32_t coeff, int idx) {

  multimod = true;

  if(coeff == 0)
    return;

  multimod_terms[gate].set_at_idx(coeff, idx);
  terms[gate] = 67;
}
/*------------------------------------------------------------------------*/

void
LinPolynomial::add_term(Gate* gate, CoeffVec& coeff) {

  if(coeff.is_zero())
    return;

  multimod_terms.emplace(gate, std::move(coeff));
  terms[gate] = 42;
}
/*------------------------------------------------------------------------*/
void
LinPolynomial::addmul_term(Gate* gate,
                           const CoeffVec& multiplier,
                           const CoeffVec& coeff) {

  if(coeff.is_zero())
    return;

  auto it = multimod_terms.find(gate);
  if(it != multimod_terms.end()) {
    CoeffVec& my_coeff = it->second;
    my_coeff.addmul(multiplier, coeff);

    if(my_coeff.is_zero()) {
      multimod_terms.erase(it);
      terms.erase(gate);
    }
  } else {
    multimod_terms.emplace(gate, multiplier * coeff);
    terms[gate] = 21;
  }
}
/*------------------------------------------------------------------------*/
void
LinPolynomial::add_polynomial(const LinPolynomial& other) {
  for(const auto& term : other.terms) {
    add_term(term.first, term.second);
  }
}
/*------------------------------------------------------------------------*/
void
LinPolynomial::subtract_polynomial(const LinPolynomial& other) {
  for(const auto& term : other.terms) {
    add_term(term.first, -term.second);
  }
}
/*------------------------------------------------------------------------*/
void
LinPolynomial::multiply_by_constant(int32_t scalar) {
  if(scalar == 0) {
    terms.clear();
    return;
  }

  for(auto& term : terms)
    term.second *= scalar;
}
/*------------------------------------------------------------------------*/

void
LinPolynomial::my_mod(int p) {
  std::vector<Gate*> to_erase;
  for(auto& term : terms) {
    term.second %= p;
    if(term.second == 0)
      to_erase.push_back(term.first);
  }

  for(Gate* g : to_erase) {
    terms.erase(g);
  }
}
/*------------------------------------------------------------------------*/
// Reduce all coefficients modulo the configured modulus and remove zero
// coefficients from the polynomial.
void
LinPolynomial::reduce_mod() {
  std::vector<Gate*> to_erase;
  for(auto& term : terms) {
    int32_t reduced = mod(term.second);
    if(reduced == 0) {
      to_erase.push_back(term.first);
    } else {
      term.second = reduced;
    }
  }

  for(Gate* g : to_erase) {
    terms.erase(g);
  }
}
/*------------------------------------------------------------------------*/
// Operation overloading
/*------------------------------------------------------------------------*/
LinPolynomial
LinPolynomial::operator+(const LinPolynomial& other) const {
  LinPolynomial result(*this);
  result.add_polynomial(other);
  return result;
}
/*------------------------------------------------------------------------*/
LinPolynomial
LinPolynomial::operator-(const LinPolynomial& other) const {
  LinPolynomial result(*this);
  result.subtract_polynomial(other);
  return result;
}
/*------------------------------------------------------------------------*/
LinPolynomial
LinPolynomial::operator*(int32_t scalar) const {
  LinPolynomial result(*this);
  result.multiply_by_constant(scalar);
  return result;
}
/*------------------------------------------------------------------------*/
LinPolynomial
operator*(int32_t scalar, const LinPolynomial& poly) {
  return poly * scalar;
}
/*------------------------------------------------------------------------*/
std::vector<Gate*> to_erase;
void
LinPolynomial::reduce_by(const LinPolynomial& other) {

  auto my_lead = *multimod_terms.begin();
  auto other_lead = other.get_leading_term();

  assert(my_lead.first != nullptr && "Leading term cannot be constant");
  assert(other_lead.first != nullptr && "Leading term cannot be constant");
  if(my_lead.first != other_lead.first)
    return;

  lin_steps++;

  CoeffVec my_coeff = my_lead.second;
  CoeffVec other_coeff;

  if(other.multimod)
    other_coeff = other.multimod_terms.begin()->second;
  else
    other_coeff.set_constant(other_lead.second);

  CoeffVec multiplier = my_coeff.negdiv(other_coeff);

  if(other.multimod) {
    for(const auto& term : other.multimod_terms) {
      other_coeff = term.second;
      addmul_term(term.first, multiplier, other_coeff);
    }
  } else {
    for(const auto& term : other.terms) {
      other_coeff.set_constant(term.second);
      addmul_term(term.first, multiplier, other_coeff);
    }
  }
  // check_coeffs(true);

  if(other.multimod)
    multimod = true;

  if(verbose > 0) {
    msg("CURRENT POLYNOMIAL:");
    print_multimod();
    std::cout << std::endl;
  }
}
/*------------------------------------------------------------------------*/

void
test_linpolynomial() {
  // Example usage of LinPolynomial
  LinPolynomial::set_modulus(5);// Set modulus to Mersenne prime 2^3 - 1

  LinPolynomial p1;
  p1.add_term(gates[num_gates - 2], 1);
  p1.add_term(gates[3], 3);
  p1.add_term(gates[10], -3);
  p1.add_term(0, -3);
  std::cout << "Polynomial p1:      ";
  p1.print();

  // now from vector
  std::vector<int32_t> coeffs = { 128, 1, 3, -1 };
  std::vector<Gate*> gates2 = { gates[num_gates - 3], gates[10], gates[3], 0 };
  LinPolynomial p2(coeffs, gates2);
  std::cout << "Polynomial p2:      ";
  p2.print();

  LinPolynomial p3 = p1 + p2;
  std::cout << "Polynomial p1 + p2: ";
  p3.print();
  LinPolynomial p4 = p1 - p2;
  std::cout << "Polynomial p1 - p2: ";
  p4.print();
  LinPolynomial p5 = p1 * 2;
  std::cout << "Polynomial p1 * 2:  ";
  p5.print();
  LinPolynomial p6 = 3 * p1;
  std::cout << "Polynomial 3 * p1:  ";
  p6.print();
  std::pair<Gate*, int32_t> leading_term = p6.get_leading_term();
  std::cout << "    Leading term: " << leading_term.first->get_var_name()
            << " with coeff " << leading_term.second << std::endl;
  int32_t coeff = p6.get_coeff(gates[10]);
  std::cout << "    Coeff of gate " << gates[10]->get_var_name() << " is "
            << coeff << std::endl;
  coeff = p6.get_coeff(gates[1]);
  std::cout << "    Coeff of gate " << gates[1]->get_var_name() << " is "
            << coeff << std::endl;
}

}// namespace talisman

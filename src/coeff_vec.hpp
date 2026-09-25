/*------------------------------------------------------------------------*/
/*! \file coeff_vec.hpp

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#ifndef TALISMAN_SRC_COEFF_VEC_H_
#define TALISMAN_SRC_COEFF_VEC_H_
/*------------------------------------------------------------------------*/
#include <algorithm>
#include <iostream>
#include <vector>

#include <nmod.h>

namespace talisman {

extern std::vector<int32_t> PRIMES;
inline std::vector<nmod_t> moduli;

class CoeffVec {

  public:
  std::vector<int32_t> entries;

  static void set_moduli() {
    moduli.clear();
    for(size_t i = 0; i < PRIMES.size(); i++) {
      nmod_t mod;
      nmod_init(&mod, PRIMES[i]);
      moduli.push_back(mod);
    }
  }

  CoeffVec()
    : entries(PRIMES.size(), 0) {}

  explicit CoeffVec(std::vector<int32_t>& vec)
    : entries(std::move(vec)) {}

  // init with 2^e but take entry i mod moduli[i]
  explicit CoeffVec(int32_t e, int sign = 1)
    : entries(PRIMES.size()) {
    for(std::size_t i = 0; i < entries.size(); ++i) {
      entries[i] = static_cast<int32_t>(
        nmod_pow_ui(2UL, static_cast<ulong>(e), moduli[i]));
      entries[i] = sign > 0 ? entries[i] : moduli[i].n - entries[i];
    }
  }

  void set_constant(int32_t c) {
    std::fill(entries.begin(), entries.end(), c);
    if(c < 0) {
      for(size_t i = 0; i < moduli.size(); i++)
        entries[i] += moduli[i].n;
    }
  }

  void set_at_idx(int32_t c, int i) {
    entries[i] = c > 0 ? c : c + moduli[i].n;
  }

  void print() const {
    std::cout << "{";
    for(std::size_t i = 0; i < entries.size(); ++i) {
      std::cout << entries[i];
      if(i + 1 < entries.size())
        std::cout << ", ";
    }
    std::cout << "}";
  }

  bool is_zero() const {
    return std::all_of(
      entries.begin(), entries.end(), [](int32_t x) { return x == 0; });
  }
  bool is_nonzero() const {
    return std::all_of(
      entries.begin(), entries.end(), [](int32_t x) { return x != 0; });
  }

  bool operator==(const CoeffVec& other) const {
    return entries == other.entries;
  }
  bool operator!=(const CoeffVec& other) const { return !(*this == other); }

  CoeffVec operator*(const CoeffVec& other) const {
    CoeffVec res;
    for(size_t i = 0; i < moduli.size(); ++i)
      res.entries[i] = nmod_mul(entries[i], other.entries[i], moduli[i]);
    return res;
  }

  void add_inplace(const CoeffVec& other) {
    for(size_t i = 0; i < moduli.size(); ++i)
      entries[i] = nmod_add(entries[i], other.entries[i], moduli[i]);
  }

  void addmul(const CoeffVec& multiplier, const CoeffVec& other) {
    for(size_t i = 0; i < moduli.size(); i++)
      NMOD_ADDMUL(
        entries[i], multiplier.entries[i], other.entries[i], moduli[i]);
  }

  CoeffVec negate() const {
    CoeffVec res;
    for(size_t i = 0; i < moduli.size(); i++)
      res.entries[i] = nmod_neg(entries[i], moduli[i]);
    return res;
  }

  CoeffVec negdiv(const CoeffVec& other) const {
    CoeffVec res;
    for(size_t i = 0; i < moduli.size(); i++)
      res.entries[i] = nmod_neg(
        nmod_div(entries[i], other.entries[i], moduli[i]), moduli[i]);
    return res;
  }

  CoeffVec operator/(const CoeffVec& other) const {
    CoeffVec res;
    for(size_t i = 0; i < moduli.size(); i++)
      res.entries[i] = nmod_div(entries[i], other.entries[i], moduli[i]);
    return res;
  }

  CoeffVec scalar_mul(ulong c) const {
    CoeffVec res;
    for(size_t i = 0; i < moduli.size(); i++)
      res.entries[i] = nmod_mul(entries[i], c, moduli[i]);
    return res;
  }
};

}// namespace talisman

#endif// TALISMAN_SRC_COEFF_VEC_H_

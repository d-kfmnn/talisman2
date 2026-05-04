#pragma once

#include <cstddef>
#include <cstdint>
// FLINT nmod helpers (provides nmod_t and nmod_pow_ui)
#include <flint/nmod.h>

// The operations in this file are inspired by:
//
//   - [1]: https://arxiv.org/pdf/2008.08654

namespace talisman {

[[nodiscard]] constexpr inline int64_t
mersenne_mod(int64_t v) noexcept {
  // Algorithm 4 of [1]
  //
  // Mersenne Prime: 2^31-1
  // b = 31
  const int64_t v_prime = v + 1;
  const int64_t z = ((v_prime >> 31) + v_prime) >> 31;
  return (v + z) & 2'147'483'647;// 2^31-1
}

[[nodiscard]] constexpr inline int64_t
pseudo_mersenne_mod(int64_t v) noexcept {
  // Algorithm 5 & 6 of [1]
  // using the pseudo Mersenne prime p = 2^24-3 = 2^24 - (2^1 + 1)
  // using m = 3 (then it works for all v <= p^2, which should be
  // enough

  const int64_t v_prime = v + 3;
  int64_t z = v_prime >> 24;
  z = ((z << 1) + z + v_prime) >> 24;// z = (3*z + v') >> 24
  z = ((z << 1) + z + v_prime) >> 24;
  z = (v + (z << 1) + z) & 16'777'215;// z = (v + 3*z) & (2^24 - 1)

  return z;
}

template<uint32_t p>
[[nodiscard]] constexpr inline uint32_t
v_mod_tmpl(uint32_t v) noexcept {
  // Algorithm 5 of [1]
  const int64_t c = 2'147'483'648 - p;// 2^31 - p = c; b = 31
  const int64_t v_prime = v + c;
  int64_t z = v_prime >> 31;
  z = (z * c + v_prime) >> 31;
  z = (z * c + v_prime) >> 31;
  return (v - z * p);
}

[[nodiscard]] constexpr inline uint32_t
mult_mod_2_31_1(uint32_t a, uint32_t b, uint32_t x) noexcept {
  // Algorithm 4 of [1], merged with multiply and addition
  //
  // Mersenne Prime: 2^31-1
  // b = 31
  const uint64_t ax = (uint64_t)a * (uint64_t)x;
  const uint64_t ax_b_1 = ax + b + 1ul;
  const uint64_t z = ((ax_b_1 >> 31) + ax_b_1) >> 31;
  return (ax_b_1 + z - 1) & 2'147'483'647ul;// 2^31-1
}

[[nodiscard]] constexpr inline uint32_t
mult_mod_p(uint32_t a, uint32_t b, uint32_t x, uint32_t p) noexcept {
  // Algorithm 5 of [1] together with multiplication.
  return ((uint64_t)a * (uint64_t)x + (uint64_t)b) % p;
}

/*
static int32_t
modpow(int32_t a, int32_t e, int32_t p) {
  int64_t base = a % p;
  int64_t result = 1;

  while(e > 0) {
    if(e & 1)
      result = (result * base) % p;
    base = (base * base) % p;
    e >>= 1;
  }

  return static_cast<int32_t>(result);
}*/



/*-------------------------------------------------------------------------*/

// Wrapper around FLINT's nmod_2_pow_ui to compute 2^exp mod p using
// FLINT's optimized implementation. Returns a value in [0, p-1]. If
// `is_negative` is true the additive inverse modulo p is returned
// (i.e., p - value), with 0 mapping to 0.
static inline uint32_t
flint_nmod_2_pow_ui(unsigned long exp,
                    const nmod_t& mod,
                    bool is_negative = false) {
  // Use FLINT's general modular exponentiation with base = 2
  uint32_t result = static_cast<uint32_t>(nmod_pow_ui(2UL, exp, mod));
  if(is_negative)
    result = (result == 0) ? 0 : static_cast<uint32_t>(mod.n - result);
  return result;
}

/*------------------------------------------------------------------------*/

}

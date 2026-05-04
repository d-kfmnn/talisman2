/*------------------------------------------------------------------------*/
/*! \file specpolynomial.hpp
    \brief contains the class SpecPolynomial and further functions to
    organize the specification polynomial, such as initializing the polynomial constraints

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#ifndef TALISMAN_SRC_SPECPOLYNOMIAL_H_
#define TALISMAN_SRC_SPECPOLYNOMIAL_H_
/*------------------------------------------------------------------------*/
#include <map>
#include <vector>

#include "gate.hpp"
#include "fast_mod.hpp"
#include "linpolynomial.hpp"
/*------------------------------------------------------------------------*/

namespace talisman {
  int32_t mod_pow2_exponent(int32_t k, int32_t n, int32_t c, bool is_negative = false);
  LinPolynomial gen_mult_spec_poly();
}  // namespace talisman

#endif  // TALISMAN_SRC_SPECPOLYNOMIAL_H_  

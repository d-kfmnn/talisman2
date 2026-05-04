/*------------------------------------------------------------------------*/
/*! \file specpolynomial.cpp
    \brief contains functions necessary to initialize the specification
  polynomial

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#include "specpoly.hpp"

#include <algorithm>
#include <iostream>
#include <random>

#include "coeff_vec.hpp"
#include "fast_mod.hpp"

using namespace talisman;

namespace talisman {

/*-------------------------------------------------------------------------*/
// FLINT-backed variant: generate the same specification polynomial but use
// FLINT's `nmod_2_pow_ui` via our wrapper `flint_nmod_2_pow_ui` so the two
// implementations can be compared.
LinPolynomial
gen_mult_spec_poly() {
  LinPolynomial spec_poly;

  // Add terms for the main gates (outputs)
  for(int i = MM - 1; i >= 0; i--) {
    Gate* g = gates[i + M - 1];

    assert(g->get_output());

    CoeffVec coeff(i);
    spec_poly.add_term(g, coeff);
  }

  // Add terms for extension gates (partial products)
  unsigned k = 1;
  for(int i = NN / 2 - 1; i >= 0; i--) {
    for(int j = NN / 2 - 1; j >= 0; j--) {
      Gate* ext = gates[num_gates - k];
      k++;
      assert(ext->is_extension());

      Gate* use = ext;
      for(int gind = 0; gind < num_gates; gind++) {
        Gate* cand = gates[gind];
        if(cand->get_dist() > 2)
          continue;
        if(cand == ext)
          continue;
        if(cand->get_ch0() != ext->get_ch0()
           && cand->get_ch0() != ext->get_ch1())
          continue;
        if(cand->get_ch1() != ext->get_ch0()
           && cand->get_ch1() != ext->get_ch1())
          continue;
        use = cand;
        break;
      }

      CoeffVec coeff(i + j, -1);
      spec_poly.add_term(use, coeff);
    }
  }

  // Check if any coefficient equals 2^NN modulo the current modulus and die
  // spec_poly.check_coeffs();

  return spec_poly;
}

}// namespace talisman

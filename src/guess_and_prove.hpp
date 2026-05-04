/*------------------------------------------------------------------------*/
/*! \file guess_and_prove.hpp
    \brief contains the functions to prove the guess and prove the guess

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#ifndef TALISMAN_SRC_GUESS_AND_PROVE_H_
#define TALISMAN_SRC_GUESS_AND_PROVE_H_
/*------------------------------------------------------------------------*/
#include <list>

#include "guess.hpp"
#include "prove.hpp"


namespace talisman {

void aggregate_thread_times(int32_t prime);
void print_individual_thread_times();

bool guess_and_prove(LinPolynomial& poly,int sc_depth_init);
bool
guess_and_prove_test(LinPolynomial& poly,
                     int32_t prime,
                     std::vector<int>& reverse_fsa_nodes,
                     std::vector<int>& bondary_nodes);
}  // namespace talisman

#endif  // TALISMAN_SRC_PROVE_H_

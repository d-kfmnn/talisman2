/*------------------------------------------------------------------------*/
/*! \file guess.hpp
    \brief contains the functions to guess

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#ifndef TALISMAN_SRC_GUESS_H_
#define TALISMAN_SRC_GUESS_H_
/*------------------------------------------------------------------------*/
#include <chrono>
#include <list>

#include "specpoly.hpp"
#include "subcircuit.hpp"

namespace talisman {

extern std::chrono::duration<double> total_la_time;
extern std::chrono::duration<double> total_repair_time;

void
cleanup();

std::pair<bool, std::vector<LinPolynomial>>
guess_naive(Subcircuit& sc, nmod_mat_t mat_p);

bool
linear_polies_from_truth_table(Subcircuit& sc,
                               std::vector<LinPolynomial>& polies);

std::pair<bool, std::vector<LinPolynomial>>
repair_guess(Subcircuit& sc,
             std::vector<std::vector<int>>& samples,
             nmod_mat_t mat);

}// namespace talisman

#endif// TALISMAN_SRC_GUESS_H_

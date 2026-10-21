/*------------------------------------------------------------------------*/
/*! \file candidates.hpp
    \brief structural detection of half-/full-adders and the reverse
           "free simplification area" boundary, computed directly from
           the parsed AIG (no external EDA library required)

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#ifndef TALISMAN_SRC_CANDIDATES_H_
#define TALISMAN_SRC_CANDIDATES_H_
/*------------------------------------------------------------------------*/
#include "parser.hpp"

#include <vector>
/*------------------------------------------------------------------------*/

namespace talisman {

/**
    Detects half-adders, full-adders (and, for API compatibility, an always
    empty list of leftover XOR structures -- affine cuts already cover that
    case) directly on the AIG produced by parse_aig().

    For each AND node with a k<=3-leaf cut (k=2 for a half-adder, k=3 for a
    full-adder) whose function matches XOR(k) paired with AND2/MAJ3 sharing
    the same leaves, a candidate is emitted with the leaves, the internal
    (body) nodes connecting leaves to the two roots, and the input/output
    polarities needed to state the arithmetic identity
        2*carry + sum = a + b (+ c)
    as a linear relation.

    Caches the accepted candidate pairs so a following call to
    get_kcover_nodes() can reuse them instead of recomputing.

    @return false if the AIG has not been parsed yet, true otherwise
*/
bool detect_ha_fa_xor(std::vector<HalfAdder>& half_adders,
                       std::vector<FullAdder>& full_adders,
                       std::vector<XorRemain>& xor_list);

/**
    Computes the reverse "free simplification area" (FSA): starting from
    each primary output, finds the first half-/full-adder pair reachable
    and grows a local cone bounded by those adders' inputs and by AIG
    level. Reuses the adder pairs cached by detect_ha_fa_xor() if present.

    @param reverse_fsa_nodes all node ids in the computed area (border+middle)
    @param boundary_nodes    the border (leaf) nodes of that area

    @return false if the AIG has not been parsed yet, true otherwise
*/
bool get_kcover_nodes(std::vector<int>& reverse_fsa_nodes,
                       std::vector<int>& boundary_nodes);

/**
    Drops the adder-pair cache populated by detect_ha_fa_xor(). Call once
    the results of get_kcover_nodes() are no longer needed.
*/
void clear_detection_cache();

/**
    Prints basic AIG statistics (PI/PO/AND counts), gated on 'verbose'.
*/
void print_aig_stats();

}  // namespace talisman
#endif  // TALISMAN_SRC_CANDIDATES_H_

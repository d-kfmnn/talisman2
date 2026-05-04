/*------------------------------------------------------------------------*/
/*! \file parser.h
    \brief contains functions necessary to parse the AIG

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#ifndef TALISMAN_SRC_PARSER_H_
#define TALISMAN_SRC_PARSER_H_
/*------------------------------------------------------------------------*/
#include "aig.hpp"
#include <cstdint>
/*------------------------------------------------------------------------*/

namespace talisman {
extern bool mult_spec;   // /< true if multiplier spec is used
extern bool miter_spec;  // /< true if miter spec is used

/**
    Reads the input aiger given in the file called input_name to the aiger 'model'
    using the parserer function of aiger.h

    @param input_name char * ame of input file
*/
struct HalfAdder {
  int id;
  std::vector<int> nodes;            // parsed from "nodes"
  std::vector<uint8_t> tt;           // parsed from "tt" (hex string array -> bytes)
  std::vector<int> input_phase;      // parsed from "input_phase"
  std::vector<int> output_phase;     // parsed from "output_phase"
};

struct FullAdder {
  int id;
  std::vector<int> nodes;            // parsed from "nodes"    
  std::vector<uint8_t> tt;           // parsed from "tt" (hex string array -> bytes)
  std::vector<int> input_phase;      // parsed from "input_phase"
  std::vector<int> output_phase;     // parsed from "output_phase"
};

// XOR remaining structure, parsed from "xor_remaining_phase"
struct XorRemain {
  int id;
  std::vector<int> nodes;            // parsed from "nodes" or "input_nodes" (see implementation)
  std::vector<uint8_t> tt;           // parsed from "tt" (hex string -> single byte or array)
  std::vector<int> input_phase;      // parsed from "input_phase"
  std::vector<int> output_phase;     // parsed from "output_phase" (single-element array)
};

struct LinearCut {
  int root;                    // root node id
  std::vector<int> leaves;     // leaf ids of the cut
  std::vector<int> nodes;      // internal nodes between leaves and root
  uint64_t truth;              // 6-input truth table (affine) stored as 64-bit word

  // ANF polynomial representation:
  // For affine function: root = constant_coeff + sum(linear_coeffs[i] * leaves[i])
  // This stores: [constant_coeff, linear_coeff_0, linear_coeff_1, ...]
  // where linear_coeff_i corresponds to leaves[i]
  // Coefficients are in {0, 1} for XOR-based (mod 2) representation
  int constant_coeff;                  // constant term (0 or 1)
  std::vector<int> linear_coeffs;      // coefficients for each leaf variable (0 or 1)
};

void parse_aig(const char *input_name);

bool read_adder_json(
  const char* filename,
  std::vector<int>& out_reverse_fsa_nodes,
  std::vector<int>& out_bondary_nodes);

// Extended to also parse xor_remaining_phase into xor_list.
bool read_faha_json(
    const char* filename,
    std::vector<HalfAdder>& half_adders,
    std::vector<FullAdder>& full_adders,
    std::vector<XorRemain>& xor_list);

}  // namespace talisman

#endif  // TALISMAN_SRC_PARSER_H_

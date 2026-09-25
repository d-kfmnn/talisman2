/*------------------------------------------------------------------------*/
/*! \file gate.h
    \brief contains the class Gate and further functions to
    organize the gate structure, such as initializing the gate constraints

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#ifndef TALISMAN_SRC_GATE_H_
#define TALISMAN_SRC_GATE_H_
/*------------------------------------------------------------------------*/
#include <list>
#include <map>
#include <vector>
#include <algorithm>
#include <queue>
#include <set>
#include <string>
#include <memory>

#include "parser.hpp"
/*------------------------------------------------------------------------*/
extern int add_var;
extern int max_dist;
extern int fsa_size;

// / Counts the number of gates
extern unsigned size_gates;

namespace talisman {

class LinPolynomial;
class Polynomial;
class IntPolynomial;

/*------------------------------------------------------------------------*/

// New function to evaluate bitwise logic on 64-bit integers
typedef unsigned long long sample_t;
struct LingateInfo {
  int id = -1;
  int type = 0;                     // 1=ha, 2=fa （当前只用到这两类作为root标注）
  std::vector<int> input_phase;     // JSON: input_phase
  std::vector<int> output_phase;    // JSON: output_phase
};
class Gate {
  private:
  unsigned var = 0;   // variable id
  std::string varname;// variable name

  unsigned ch0 = 1;          // child 0 (default 1 = true)
  unsigned ch1 = 1;          // child 1 (default 1 = true)
  std::set<unsigned> parents;// set of parent gates

  int carry = 0; // is carry in gate

  bool elim = false;   // elimination flag

  bool input = false;     // is input gate
  bool output = false;    // is output gate
  bool aig_output = false;// is AIG output
  bool fsa = false;
  bool fsa_inp = false;
  // A node may belong to multiple FA/HA instances; store list of (id,type)
  // where type = 1 (HA) or 2 (FA). Input-mark lists are stored similarly.
  std::vector<std::pair<int,int>> faha_ids;          // (id, type)
  std::vector<std::pair<int,int>> faha_inp_ids;      // (id, type)
  std::vector<std::pair<int,int>> faha_xor_ids;      // (id, type)
  std::vector<std::pair<int,int>> faha_xor_inp_ids;  // (id, type)
  std::vector<std::pair<int,int>> faha_maj_ids;      // (id, type)
  std::vector<std::pair<int,int>> faha_maj_inp_ids;  // (id, type)
  int xor_inp_id = -1;
  int xor_id = -1;        // id of XOR structure this gate belongs to (XorRemain)
  int xor_root_id = -1;
  int xor_root_inp_id = -1;
  int faha_type = 0;     // 1=ha, 2=fa
  bool extension = false;// extension flag
  unsigned distance = 0; // distance to inputs
  unsigned level = 0;    // unique leevel in the AIG
  unsigned orig_level = 0; // original level before any reassignment
  int xor_gate = 0;      //
  sample_t value = 0;    // value of the gate, unknown if not set

  int cone_id = -1; // cone id for output cones
  // optional linear polynomial associated with this gate
  std::unique_ptr<LinPolynomial> gatepoly;
  std::unique_ptr<Polynomial> gatepoly_nonlinear;
  // optional integer polynomial associated with this gate (for decomposition)
  std::unique_ptr<IntPolynomial> gatepoly_int;
  public:
  // Constructor
  Gate(int n_, std::string name_, bool input_ = false, bool output_ = false);

  // Destructor
  ~Gate();

  // Accessors for variable id
  unsigned get_var() const { return var; }
  const std::string& get_var_name() const { return varname; }

  // Accessors for input/output flags
  bool get_input() const { return input; }
  bool get_output() const { return output; }

  // Accessors for AIG output
  bool get_aig_output() const { return aig_output; }
  void mark_aig_output() { aig_output = true; }

  // Set/get extension flag
  void set_ext() { extension = true; }
  bool is_extension() const { return extension; }

  // Distance accessors
  unsigned get_dist() const { return distance; }
  void set_dist(unsigned l) { distance = l; }

  // Level accessors
  unsigned get_level() const { return level; }
  void set_level(unsigned l) { level = l; }

  unsigned get_orig_level() const { return orig_level; }
  void set_orig_level(unsigned l) { orig_level = l; }

  // Getter for ch0
  unsigned get_ch0() const { return ch0; }

  // Setter for ch0
  void set_ch0(unsigned val) { ch0 = val; }

  // Getter for ch1
  unsigned get_ch1() const { return ch1; }

  // Setter for ch1
  void set_ch1(unsigned val) { ch1 = val; }

  bool is_fsa() const { return fsa; }
  void set_fsa() { fsa = true; }
  void reset_fsa() { fsa = false; }

  bool is_fsa_inp() const { return fsa_inp; }
  void set_fsa_inp() { fsa_inp = true; }
  int get_xor_id() const { return xor_id; }
  void set_xor_id(int id) { xor_id = id; }
  int get_xor_inp_id() const { return xor_inp_id; }
  void set_xor_inp_id(int id) { xor_inp_id = id; }
    
  void reset_xor_id() { xor_id = -1; }

  int get_xor_root_id() const { return xor_root_id; }
  void set_xor_root_id(int id) { xor_root_id = id; }
  int get_xor_root_inp_id() const { return xor_root_inp_id; }
  void set_xor_root_inp_id(int id) { xor_root_inp_id = id; }

  void reset_xor_root_id() { xor_root_id = -1; }

  
  // Returns whether this gate belongs to any FA/HA instance
  bool has_any_faha_id() const { return !faha_ids.empty(); }
  // Contains check for a specific FA/HA id (ignores type)
  bool has_faha_id(int id) const {
    if (id == -1) return false;
    for (const auto &p : faha_ids) if (p.first == id) return true;
    return false;
  }
  // Contains check for specific id+type
  bool has_faha_id_with_type(int id, int type) const {
    if (id == -1) return false;
    for (const auto &p : faha_ids) if (p.first == id && p.second == type) return true;
    return false;
  }
  // Append an FA/HA id (avoid duplicates). Type may be set later.
  void set_faha_id(int id) {
    if(id == -1) return;
    if(!has_faha_id(id)) faha_ids.emplace_back(id, 0);
  }
  // Append or set id with explicit type
  void set_faha_id_with_type(int id, int type) {
    if(id == -1) return;
    for (auto &p : faha_ids) if (p.first == id) { p.second = type; return; }
    faha_ids.emplace_back(id, type);
  }
  void reset_faha_id() { faha_ids.clear(); }

  const std::vector<std::pair<int,int>>& get_faha_ids() const { return faha_ids; }

  // Return the first recorded type for backward compatibility, or -1
  int get_faha_type() const { return faha_ids.empty() ? -1 : faha_ids.front().second; }
  void set_faha_type(int type) {
    if (faha_ids.empty()) return;
    faha_ids.front().second = type;
  }
  // Return the type recorded for a specific id, or -1 if not present
  int get_faha_type_for_id(int id) const {
    for (const auto &p : faha_ids) if (p.first == id) return p.second;
    return -1;
  }

  // Returns whether this gate is marked as an input to any FA/HA instance
  bool has_any_faha_inp_id() const { return !faha_inp_ids.empty(); }
  // Contains check for a specific FA/HA input id (ignores type)
  bool has_faha_inp_id(int id) const {
    if (id == -1) return false;
    for (const auto &p : faha_inp_ids) if (p.first == id) return true;
    return false;
  }
  bool has_faha_inp_id_with_type(int id, int type) const {
    if (id == -1) return false;
    for (const auto &p : faha_inp_ids) if (p.first == id && p.second == type) return true;
    return false;
  }
  void set_faha_inp_id(int id) {
    if(id == -1) return;
    if(!has_faha_inp_id(id)) faha_inp_ids.emplace_back(id, 0);
  }
  void set_faha_inp_id_with_type(int id, int type) {
    if(id == -1) return;
    for (auto &p : faha_inp_ids) if (p.first == id) { p.second = type; return; }
    faha_inp_ids.emplace_back(id, type);
  }
  const std::vector<std::pair<int,int>>& get_faha_inp_ids() const { return faha_inp_ids; }


  // XOR-related membership checks
  bool has_any_faha_xor_id() const { return !faha_xor_ids.empty(); }
  bool has_faha_xor_id(int id) const {
    if(id == -1) return false;
    for (const auto &p : faha_xor_ids) if (p.first == id) return true;
    return false;
  }
  bool has_faha_xor_id_with_type(int id, int type) const {
    if(id == -1) return false;
    for (const auto &p : faha_xor_ids) if (p.first == id && p.second == type) return true;
    return false;
  }
  void set_faha_xor_id(int id) {
    if(id == -1) return;
    if(!has_faha_xor_id(id)) faha_xor_ids.emplace_back(id, 0);
  }
  void set_faha_xor_id_with_type(int id, int type) {
    if(id == -1) return;
    for (auto &p : faha_xor_ids) if (p.first == id) { p.second = type; return; }
    faha_xor_ids.emplace_back(id, type);
  }
  bool has_any_faha_xor_inp_id() const { return !faha_xor_inp_ids.empty(); }
  bool has_faha_xor_inp_id(int id) const {
    if(id == -1) return false;
    for (const auto &p : faha_xor_inp_ids) if (p.first == id) return true;
    return false;
  }
  bool has_faha_xor_inp_id_with_type(int id, int type) const {
    if(id == -1) return false;
    for (const auto &p : faha_xor_inp_ids) if (p.first == id && p.second == type) return true;
    return false;
  }
  void set_faha_xor_inp_id(int id) {
    if(id == -1) return;
    if(!has_faha_xor_inp_id(id)) faha_xor_inp_ids.emplace_back(id, 0);
  }
  void set_faha_xor_inp_id_with_type(int id, int type) {
    if(id == -1) return;
    for (auto &p : faha_xor_inp_ids) if (p.first == id) { p.second = type; return; }
    faha_xor_inp_ids.emplace_back(id, type);
  }
  const std::vector<std::pair<int,int>>& get_faha_xor_ids() const { return faha_xor_ids; }
  const std::vector<std::pair<int,int>>& get_faha_xor_inp_ids() const { return faha_xor_inp_ids; }

  // Majority-related membership checks
  bool has_any_faha_maj_id() const { return !faha_maj_ids.empty(); }
  bool has_faha_maj_id(int id) const {
    if(id == -1) return false;
    for (const auto &p : faha_maj_ids) if (p.first == id) return true;
    return false;
  }
  bool has_faha_maj_id_with_type(int id, int type) const {
    if(id == -1) return false;
    for (const auto &p : faha_maj_ids) if (p.first == id && p.second == type) return true;
    return false;
  }
  void set_faha_maj_id(int id) {
    if(id == -1) return;
    if(!has_faha_maj_id(id)) faha_maj_ids.emplace_back(id, 0);
  }
  void set_faha_maj_id_with_type(int id, int type) {
    if(id == -1) return;
    for (auto &p : faha_maj_ids) if (p.first == id) { p.second = type; return; }
    faha_maj_ids.emplace_back(id, type);
  }
  bool has_any_faha_maj_inp_id() const { return !faha_maj_inp_ids.empty(); }
  bool has_faha_maj_inp_id(int id) const {
    if(id == -1) return false;
    for (const auto &p : faha_maj_inp_ids) if (p.first == id) return true;
    return false;
  }
  bool has_faha_maj_inp_id_with_type(int id, int type) const {
    if(id == -1) return false;
    for (const auto &p : faha_maj_inp_ids) if (p.first == id && p.second == type) return true;
    return false;
  }
  void set_faha_maj_inp_id(int id) {
    if(id == -1) return;
    if(!has_faha_maj_inp_id(id)) faha_maj_inp_ids.emplace_back(id, 0);
  }
  void set_faha_maj_inp_id_with_type(int id, int type) {
    if(id == -1) return;
    for (auto &p : faha_maj_inp_ids) if (p.first == id) { p.second = type; return; }
    faha_maj_inp_ids.emplace_back(id, type);
  }
  const std::vector<std::pair<int,int>>& get_faha_maj_ids() const { return faha_maj_ids; }
  const std::vector<std::pair<int,int>>& get_faha_maj_inp_ids() const { return faha_maj_inp_ids; }

  // set _input/output phases and truth table for ha & ha



  int get_cone_id() const { return cone_id; }
  void set_cone_id(int id) { cone_id = id; }
  // Get AIG parents and modify
  const std::set<unsigned>& get_parents() const { return parents; }
  void add_parent(unsigned parent) { parents.insert(parent); }
  void remove_parent(unsigned parent) { parents.erase(parent); }
  size_t parent_size() const { return parents.size(); }

  // Value accessors
  sample_t get_value() const { return value; }
  void set_value(sample_t val) { value = val; }
  void reset_value() { value = 0; }
  void eval_gate();

  // LinPolynomial accessors
  bool has_gatepoly() const;
  const LinPolynomial* get_gatepoly() const;
  LinPolynomial* get_gatepoly_mut();
  void set_gatepoly(const LinPolynomial& lp);
  void clear_gatepoly();

  // Polynomial accessors
  bool has_gatepoly_nonlinear() const;
  const Polynomial* get_gatepoly_nonlinear() const;
  Polynomial* get_gatepoly_nonlinear_mut();
  void set_gatepoly_nonlinear(const Polynomial& poly);
  void clear_gatepoly_nonlinear();

  // IntPolynomial accessors (for decomposition with integer coefficients)
  bool has_gatepoly_int() const { return gatepoly_int != nullptr; }
  const IntPolynomial* get_gatepoly_int() const { return gatepoly_int.get(); }
  IntPolynomial* get_gatepoly_int_mut() { return gatepoly_int.get(); }
  void set_gatepoly_int(const IntPolynomial& poly);
  void clear_gatepoly_int();

// carry gate accessors
  int get_carry_gate() const { return carry; }
  void set_carry_gate(int val) { carry = val; }
  void inc_carry_gate() { carry++; }
  void dec_carry_gate() { carry--; }

  // elimination flag accessors
  bool is_elim() const { return elim; }
  void set_elim() { elim = true; }
  void reset_elim() { elim = false; }


  // Get CNF encoding of the gate
  std::vector<int> get_cnf_encoding();

  // Get PB encoding of the gate
  std::vector<std::tuple<int, int, std::string>> get_pb_encoding();
};

/*------------------------------------------------------------------------*/

// Global pointer array to all gates
extern Gate** gates;
extern std::vector<std::vector<Gate*>> gates_by_dist;

// Number of gates
extern unsigned num_gates;

// Initialize gates with extension flag
void
init_gates(int ext);

bool mark_fahaxor(const std::vector<FullAdder>& fa_list, int fa_flag);
bool mark_fahaxor(const std::vector<HalfAdder>& ha_list, int fa_flag);
bool mark_fahaxor(const std::vector<XorRemain>& xor_list);

// Build and cache LinPolynomials from affine (linear) cuts
// The polynomial represents: root = constant_coeff XOR sum(linear_coeffs[i] * leaves[i])
// which translates to: root - constant_coeff - sum(linear_coeffs[i] * leaves[i]) = 0
void build_linpoly_from_linear_cuts(const std::vector<LinearCut>& cuts, int32_t prime);

bool
mark_fsa(std::vector<int> nodes, std::vector<int> inp_nodes);

// Return pointer to gate with given aiger literal
Gate*
gate(int lit);


void reset_levels_to_orig();

void adjust_cones_by_merging();

static inline Gate* getGateByNodeId(int node_id) {
  if (node_id < 1) return nullptr;
  return gate(node_id * 2);
}
// Delete all gates and free memory
void
delete_gates();

// Delete all gates and reset the global gate state for a new run
void
reset_gates();

void decompose_all_gate_polynomials(Gate * border);

void reassign_levels_by_cone_order();

std::vector<std::vector<Gate*>> get_gates_by_slice();

std::vector<Gate*> get_gates_by_decreasing_level();

void print_cones();

// Forward declaration for IntPolynomial
class IntPolynomial;

// Substitute gate n in parent polynomial with its expression
// Uses integer coefficients with Boolean algebra (x*x = x)
IntPolynomial substitute_gate(Gate* n, const IntPolynomial& parent_poly,
                              const IntPolynomial& substitution_expr);

}// namespace talisman

#endif// TALISMAN_SRC_GATE_H_

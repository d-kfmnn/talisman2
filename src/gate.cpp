/*------------------------------------------------------------------------*/
/*! \file gate.cpp
    \brief contains the class Gate and further functions to
    organize the gate structure, such as initializing the gate constraints

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#include "gate.hpp"

#include "guess.hpp"
#include "int_polynomial.hpp"
#include "linpoly_cache.hpp"
#include "linpolynomial.hpp"
#include "polynomial.hpp"
#include "subcircuit.hpp"
#include "thread_pool.hpp"
#include <algorithm>
#include <atomic>
#include <bitset>
#include <iomanip>
#include <list>
#include <string>
#include <utility>
using namespace talisman;

namespace talisman {
/*------------------------------------------------------------------------*/
// Global variables
int add_var = 0;
unsigned max_dist = 0;
int fsa_size = 0;
/*------------------------------------------------------------------------*/
// ERROR CODES:
static int err_allocate = 91;// failed to allocate gates

/*------------------------------------------------------------------------*/
Gate::Gate(int n_, std::string name_, bool input_, bool output_)
  : input(input_)
  , output(output_) {
  var = n_;
  varname = name_;
}

/*------------------------------------------------------------------------*/
Gate::~Gate() {}

bool
Gate::has_gatepoly() const {
  return gatepoly != nullptr;
}

LinPolynomial*
Gate::get_gatepoly_mut() {
  return gatepoly.get();
}

const LinPolynomial*
Gate::get_gatepoly() const {
  return gatepoly.get();
}

void
Gate::set_gatepoly(const LinPolynomial& lp) {
  gatepoly = std::make_unique<LinPolynomial>(lp);
}

bool
Gate::has_gatepoly_nonlinear() const {
  return gatepoly_nonlinear != nullptr;
}

Polynomial*
Gate::get_gatepoly_nonlinear_mut() {
  return gatepoly_nonlinear.get();
}

const Polynomial*
Gate::get_gatepoly_nonlinear() const {
  return gatepoly_nonlinear.get();
}

void
Gate::set_gatepoly_nonlinear(const Polynomial& poly) {
  gatepoly_nonlinear = std::make_unique<Polynomial>(poly);
}

void
Gate::clear_gatepoly_nonlinear() {
  gatepoly_nonlinear.reset();
}

void
Gate::set_gatepoly_int(const IntPolynomial& poly) {
  gatepoly_int = std::make_unique<IntPolynomial>(poly);
}

void
Gate::clear_gatepoly_int() {
  gatepoly_int.reset();
}

void
Gate::clear_gatepoly() {
  gatepoly.reset();
}

/*------------------------------------------------------------------------*/

sample_t
eval_from_ints(sample_t a, sample_t b, unsigned l, unsigned r) {
  sample_t left = (l & 1) ? ~a : a;
  sample_t right = (r & 1) ? ~b : b;
  return left & right;
}

/*------------------------------------------------------------------------*/
void
Gate::eval_gate() {
  // only one child (second is true)
  if(ch1 == 1) {
    Gate* l = gate(ch0);

    value = (ch0 & 1) ? ~l->get_value() : l->get_value();

    LOGDBG("evaluating gate %s with child %u (%u=%s) to value %s",
           varname.c_str(),
           ch0,
           l->get_var(),
           std::bitset<64>(l->get_value()).to_string().c_str(),
           std::bitset<64>(value).to_string().c_str());
    return;
  }

  Gate* l = gate(ch0);
  sample_t l_value = l->get_value();
  Gate* r = gate(ch1);
  sample_t r_value = r->get_value();

  value = eval_from_ints(l_value, r_value, ch0, ch1);

  LOGDBG("evaluating gate %s with children %u (l) and %u (r):\n      l = %s\n  "
         "    r = %s\n  value = %s",
         varname.c_str(),
         ch0,
         ch1,
         std::bitset<64>(l->get_value()).to_string().c_str(),
         std::bitset<64>(r->get_value()).to_string().c_str(),
         std::bitset<64>(value).to_string().c_str());
}

/*------------------------------------------------------------------------*/
std::vector<int>
Gate::get_cnf_encoding() {
  if(input) {
    warning("input gate %s has no CNF encoding", varname.c_str());
    return {};
  }

  int lhs = var >> 1;
  std::vector<int> cnf;

  // Helper to convert a child into its CNF literal
  auto to_literal = [](int var) {
    int lit = var >> 1;
    return (var & 1) ? -lit : lit;
  };

  int lit0 = to_literal(ch0);
  int lit1 = to_literal(ch1);

  if(ch1 == 1) {
    // Single-input gate (literal)
    cnf.insert(cnf.end(), { -lhs, lit0, 0, lhs, -lit0, 0 });

    LOGDBG("CNF encoding of gate %s with child %u: ", varname.c_str(), ch0);
    LOGDBG("  (%d, %d)", -lhs, lit0);
    LOGDBG("  (%d, %d)", lhs, -lit0);
  } else {
    // Two-input AND gate
    cnf.insert(cnf.end(),
               { lhs, -lit0, -lit1, 0, -lhs, lit0, 0, -lhs, lit1, 0 });

    LOGDBG("CNF encoding of gate %s with children %u and %u: ",
           varname.c_str(),
           ch0,
           ch1);
    LOGDBG("  (%d, %d, %d)", lhs, -lit0, -lit1);
    LOGDBG("  (%d, %d)", -lhs, lit0);
    LOGDBG("  (%d, %d)", -lhs, lit1);
  }

  return cnf;
}
/*------------------------------------------------------------------------*/
std::vector<std::tuple<int, int, std::string>>
Gate::get_pb_encoding() {
  /*
   * Translate out <-> (in0 & in1) into
   *  out - in0       <= 0
   *  out - in1       <= 0
   *  in0 + in1 - out <= 1
   *
   *  or if only one input
   *  out - in <= 0
   *  in - out <= 0
   */
  if(input) {
    warning("input gate %s has no pb encoding", varname.c_str());
    return {};
  }

  std::string out = varname;
  std::vector<std::tuple<int, int, std::string>> pb;

  // Helper to get sign of child
  auto sign = [](int var) { return (var & 1) ? -1 : 1; };

  if(ch1 == 1) {
    // Single-input gate (literal)
    std::string in = gate(ch0)->get_var_name();

    pb.insert(pb.end(),
              { std::make_tuple(1, 1, out),
                std::make_tuple(-1, sign(ch0), in),
                std::make_tuple(0, 0, ""),
                std::make_tuple(1, sign(ch0), in),
                std::make_tuple(-1, 1, out),
                std::make_tuple(0, 0, "") });

    std::string log = sign(ch0) > 0 ? in : "(1 - " + in + ")";
    LOGDBG("PB encoding of gate %s with child %u: ", varname.c_str(), ch0);
    LOGDBG("  %s - %s <= 0", out.c_str(), log.c_str());
    LOGDBG("  %s - %s <= 0", log.c_str(), out.c_str());
  } else {
    // Two-input AND gate
    std::string in0 = gate(ch0)->get_var_name();
    std::string in1 = gate(ch1)->get_var_name();

    pb.insert(pb.end(),
              {
                std::make_tuple(1, 1, out),
                std::make_tuple(-1, sign(ch0), in0),
                std::make_tuple(0, 0, ""),
                std::make_tuple(1, 1, out),
                std::make_tuple(-1, sign(ch1), in1),
                std::make_tuple(0, 0, ""),
                std::make_tuple(1, sign(ch0), in0),
                std::make_tuple(1, sign(ch1), in1),
                std::make_tuple(-1, 1, out),
                std::make_tuple(1, 0, ""),

              });

    std::string log0 = sign(ch0) > 0 ? in0 : "(1 - " + in0 + ")";
    std::string log1 = sign(ch1) > 0 ? in1 : "(1 - " + in1 + ")";

    LOGDBG("PB encoding of gate %s with children %u and %u: ",
           varname.c_str(),
           ch0,
           ch1);
    LOGDBG("  %s - %s <= 0", out.c_str(), log0.c_str());
    LOGDBG("  %s - %s <= 0", out.c_str(), log1.c_str());
    LOGDBG("  %s + %s - %s <= 1", log0.c_str(), log1.c_str(), out.c_str());
  }

  return pb;
}

/*------------------------------------------------------------------------*/
Gate** gates = nullptr;
unsigned num_gates = 0;
unsigned size_gates = 0;

/*------------------------------------------------------------------------*/

static void
mark_aig_outputs() {
  for(unsigned i = 0; i < MM; i++) {
    unsigned lit = slit(i);
    if(lit < 2)
      continue;
    Gate* n = gate(lit);
    n->mark_aig_output();

    LOGDBG("marking %s as AIG output", n->get_var_name().c_str());
  }
}
/*------------------------------------------------------------------------*/

static void
allocate_gates(int expected_extension) {
  unsigned aiger;
  // num of gates + num of outputs + expected extension - index correction
  num_gates = M + MM + expected_extension - 1;

  msg("allocating %i gates", num_gates);
  gates = new Gate*[num_gates];
  size_gates = num_gates;

  if(!gates)
    die(err_allocate, "failed to allocate gates");

  if(mult_spec) {
    // inputs a
    for(unsigned i = a0; i <= al; i += ainc) {
      aiger = 2 * (i + 1);

      assert(is_model_input(aiger));

      std::string name = std::string("a") + std::to_string((i - a0) / ainc);
      gates[i] = new Gate(aiger, name, 1);

      LOGDBG("input %s at %u", name.c_str(), i);
    }

    // inputs b
    for(unsigned i = b0; i <= bl; i += binc) {
      aiger = 2 * (i + 1);
      assert(is_model_input(aiger));

      std::string name = std::string("b") + std::to_string((i - b0) / binc);
      gates[i] = new Gate(aiger, name, 1);

      LOGDBG("input %s at %u", name.c_str(), i);
    }
  } else {
    for(unsigned i = 0; i < NN; i++) {
      aiger = 2 * (i + 1);
      assert(is_model_input(aiger));

      std::string name = std::string("i") + std::to_string(i);
      gates[i] = new Gate(aiger, name, 1);

      LOGDBG("input %s at %u", name.c_str(), i);
    }
  }

  // internal gates
  for(unsigned i = NN; i < M - 1; i++) {
    aiger = 2 * (i + 1);
    assert(is_model_and(aiger));

    std::string name = std::string("l") + std::to_string(aiger);
    gates[i] = new Gate(aiger, name, 0);

    LOGDBG("internal gate %s at %u", name.c_str(), i);
  }

  // output s
  int index = 0;
  for(unsigned i = M - 1; i < M - 1 + MM; i++) {
    aiger = 2 * (i + 1);
    std::string name = std::string("s") + std::to_string(index++);
    gates[i] = new Gate(aiger, name, 0, 1);

    LOGDBG("output %s at %u", name.c_str(), i);
  }

  // extensions
  index = 0;
  for(unsigned i = M - 1 + MM; i < M + MM + expected_extension - 1; i++) {
    aiger = 2 * (i + 1);
    std::string name = std::string("t") + std::to_string(index++);
    gates[i] = new Gate(aiger, name, 0, 0);

    gates[i]->set_ext();
    // extension nodes are not part of the original AIG, distance = 0
    gates[i]->set_dist(0);
    LOGDBG("extension %s at %u", name.c_str(), i);
  }

  msg("finished allocation");
}

/*------------------------------------------------------------------------*/
static void
set_parents_and_children() {
  unsigned pp = 0;
  for(unsigned i = M - 1; i >= NN; i--) {
    Gate* g = gates[i];

    assert(!g->get_input());
    aiger_and* and1 = is_model_and(g->get_var());
    if(!and1)
      continue;

    unsigned l = and1->rhs0, r = and1->rhs1;
    g->set_ch0(l);
    g->set_ch1(r);

    if(gate(l)->get_input() && gate(r)->get_input())
      pp++;

    if(l > 1)
      gate(l)->add_parent(g->get_var());
    if(r > 1)
      gate(r)->add_parent(g->get_var());

    LOGDBG("internal gate %s at %u with children %s (%u) and %s (%u)",
           g->get_var_name().c_str(),
           i,
           gate(l)->get_var_name().c_str(),
           l,
           gate(r)->get_var_name().c_str(),
           r);
  }
  // set children for outputs
  for(unsigned i = 0; i < MM; i++) {
    Gate* g = gates[i + M - 1];
    assert(g->get_output());
    unsigned lit = slit(i);
    if(lit < 2)
      continue;// special case if output is a constant

    g->set_ch0(lit);

    gate(lit)->add_parent(g->get_var());

    LOGDBG("output %s at %u with child %s",
           g->get_var_name().c_str(),
           i + M - 1,
           gate(lit)->get_var_name().c_str());
  }

  if(mult_spec) {
    int indx = M - 1 + MM;

    for(unsigned i = a0; i <= al; i += ainc) {
      for(unsigned j = b0; j <= bl; j += binc) {
        Gate* g = gates[indx];
        unsigned l = gates[i]->get_var();
        unsigned r = gates[j]->get_var();
        g->set_ch0(l);
        g->set_ch1(r);

        if(l > 1)
          gate(l)->add_parent(g->get_var());
        if(r > 1)
          gate(r)->add_parent(g->get_var());

        LOGDBG("extension gate %s at %u with children %s and %s",
               g->get_var_name().c_str(),
               indx,
               gate(l)->get_var_name().c_str(),
               gate(r)->get_var_name().c_str());

        indx++;
      }
    }
  }

  if(pp != NN / 2 * NN / 2) {
    booth = 1;
    msg("detected booth multiplier");
  }
}
/*------------------------------------------------------------------------*/
static void
print_node_statistics() {
  for(unsigned i = 0; i < num_gates; i++) {
    Gate* g = gates[i];
    if(!g)
      continue;

    std::string parents_str;
    for(const auto& p : g->get_parents()) {
      parents_str += gate(p)->get_var_name() + " ";
    }

    msg("gate %s at %u: parents: (%s), children: (%s, %s), distance: %d, "
        "level: %d, extension: %d",
        g->get_var_name().c_str(),
        i,
        parents_str.c_str(),
        gate(g->get_ch0()) ? gate(g->get_ch0())->get_var_name().c_str() : "1",
        gate(g->get_ch1()) ? gate(g->get_ch1())->get_var_name().c_str() : "1",
        g->get_dist(),
        g->get_level(),
        g->is_extension());
  }
}

/*------------------------------------------------------------------------*/
static void
set_distance_to_input() {
  max_dist = 0;
  for(unsigned i = 0; i < num_gates; i++) {
    Gate* g = gates[i];
    if(!g)
      continue;

    if(g->get_input()) {
      g->set_dist(0);

      LOGDBG("input %s at %u has distance 0", g->get_var_name().c_str(), i);

      continue;
    }

    unsigned l = g->get_ch0();
    unsigned r = g->get_ch1();

    if(g->get_output()) {
      // output gates should only have one child
      assert(r == 1);
      if(l < 2) {
        // if the only child is a constant, distance is 0
        g->set_dist(0);
        LOGDBG("output gate %s at %u has distance 0 (child is constant)",
               g->get_var_name().c_str(),
               i);

        continue;
      }

      Gate* l_gate = gate(l);
      unsigned dist_l = l_gate->get_dist();
      g->set_dist(dist_l + 1);
      if(dist_l + 1 > max_dist) {
        max_dist = dist_l + 1;
      }
      LOGDBG("output gate %s at %u has distance %d",
             g->get_var_name().c_str(),
             i,
             g->get_dist());

      continue;
    }

    if((l < 2 || r < 2)) {
      // if one child is constant, distance is 0
      g->set_dist(0);
      LOGDBG("gate %s at %u has distance 0 (one child is constant)",
             g->get_var_name().c_str(),
             i);

      continue;
    }
    // both children are not constants, calculate distance
    assert(l > 1 && r > 1);// both children should be gates
    Gate *l_gate = gate(l), *r_gate = gate(r);
    unsigned dist_l = l_gate->get_dist();
    unsigned dist_r = r_gate->get_dist();
    unsigned dist_n = dist_l > dist_r ? dist_l + 1 : dist_r + 1;

    g->set_dist(dist_n);

    if(dist_n > max_dist)
      max_dist++;

    LOGDBG("gate %s at %u has distance %d",
           g->get_var_name().c_str(),
           i,
           g->get_dist());
  }

  msg("maximum distance to input: %d", max_dist);
}
/*------------------------------------------------------------------------*/
std::vector<std::vector<Gate*>> gates_by_dist;

static void
set_levels() {
  // Create a vector of buckets for gates by distance (1-based indexing)
  gates_by_dist.resize(max_dist + 1);

  // Bucket gates by their distance
  for(unsigned i = 0; i < num_gates; i++) {
    Gate* g = gates[i];
    unsigned dist = g->get_dist();
    assert(dist <= max_dist);
    gates_by_dist[dist].push_back(g);
  }

  // Assign levels sequentially
  int level = 1;
  for(unsigned dist = 0; dist <= max_dist; ++dist) {
    for(auto* gate : gates_by_dist[dist]) {
      if(gate->get_output())
        continue;// skip outputs

      gate->set_level(level++);
      LOGDBG("gate %s at dist %u has level %d",
             gate->get_var_name().c_str(),
             dist,
             gate->get_level());
    }
  }

  for(unsigned i = M - 1; i < M - 1 + MM; i++) {
    Gate* g = gates[i];
    g->set_level(level++);
    LOGDBG(
      "output gate %s has level %d", g->get_var_name().c_str(), g->get_level());
  }
}
/*----------------------------------------------------------------------*/
/**
 * Helper for marking XOR-only instances.
 * It marks all nodes of the XOR structure with the XOR id.
 */
static void
mark_one_xor_instance(const XorRemain& xr, int num_inputs = 2) {
  const auto& nodes = xr.nodes;
  const int size = static_cast<int>(nodes.size());
  if(size == 0) {
    warning("xor_remain id=%d has empty nodes, skip", xr.id);
    return;
  }
  if(num_inputs <= 0 || size <= num_inputs) {
    warning("xor_remain id=%d nodes size=%d <= num_inputs=%d, skip",
            xr.id,
            size,
            num_inputs);
    return;
  }

  // First node is the XOR root.
  int root_node_id = nodes[0];
  Gate* root_gate = getGateByNodeId(root_node_id);
  if(!root_gate) {
    warning(
      "xor_remain id=%d root node %d has no gate, skip", xr.id, root_node_id);
    return;
  }
  //
  const int body_limit = size - num_inputs;// boundary between body and inputs

  for(int i = 0; i < size; ++i) {
    int node_id = nodes[i];
    Gate* g = getGateByNodeId(node_id);
    if(!g)
      continue;

    if(i >= body_limit) {
      // Last num_inputs nodes are treated as XOR inputs.
      g->set_xor_inp_id(xr.id);
    } else {
      // All other nodes belong to this XOR structure as internal/body nodes.
      g->set_xor_id(xr.id);
    }
  }

  // Additionally mark the XOR root.
  root_gate->set_xor_root_id(xr.id);
}
/*------------------------------------------------------------------------*/
/**
 * Helper for marking one HA/FA instance.
 * It marks input nodes (first num_inputs nodes) and internal nodes,
 * and additionally marks XOR and MAJ roots (last two nodes in the list).
 */
template<typename AdderT>
static void
mark_one_faha_instance(const AdderT& adder,
                       int num_inputs,
                       int type_label,// 1 = ha, 2 = fa
                       const char* kind_label) {
  const auto& nodes = adder.nodes;
  const int size = (int)nodes.size();
  if(size < 2) {
    warning("%s id=%d nodes size < 2, skip", kind_label, adder.id);
    return;
  }

  // Find the minimum cone_id among all nodes
  unsigned min_cone_id = UINT_MAX;
  for(int i = num_inputs; i < size; ++i) {
    int node_id = nodes[i];
    Gate* g = getGateByNodeId(node_id);
    if(!g)
      continue;
    if(g->get_cone_id() < min_cone_id) {
      min_cone_id = g->get_cone_id();
    }
  }

  // Mark all nodes and set them to the same minimum cone_id
  for(int i = 0; i < size; ++i) {
    int node_id = nodes[i];
    Gate* g = getGateByNodeId(node_id);
    if(!g)
      continue;
    if(i < num_inputs) {
      g->set_faha_inp_id_with_type(adder.id, type_label);
    } else {
      g->set_faha_id_with_type(adder.id, type_label);
      // Set all nodes (except inputs) to the same minimum cone_id
      if(min_cone_id != UINT_MAX && type_label == 2) {
        g->set_cone_id(min_cone_id);
      }
    }
  }

  int xorRoot = nodes[size - 2];
  int majRoot = nodes[size - 1];
  Gate* gx = getGateByNodeId(xorRoot);
  Gate* gm = getGateByNodeId(majRoot);

  auto attach_info_to_root = [&](Gate* gr, const char* root_name) {
    if(!gr)
      return;
    gr->set_faha_id_with_type(adder.id, type_label);
    if(root_name && std::string(root_name) == "xor")
      gr->set_faha_xor_id_with_type(adder.id, type_label);
    if(root_name && std::string(root_name) == "maj")
      gr->set_faha_maj_id_with_type(adder.id, type_label);
    // Also set cone_id for roots
    if(min_cone_id != UINT_MAX && type_label == 2) {
      gr->set_cone_id(min_cone_id);
    }
  };

  attach_info_to_root(gx, "xor");
  attach_info_to_root(gm, "maj");
}

/*------------------------------------------------------------------------*/
/**
 * Original build_linpoly_from_lingates, renamed to
 * build_linpoly_from_lingates_hafa. It builds linear polynomials for both XOR
 * and MAJ parts of HA/FA.
 */
template<typename AdderT>
void
build_linpoly_from_lingates_hafa(int info_type, const AdderT& adder) {
  const auto& nodes = adder.nodes;
  const int size = static_cast<int>(nodes.size());
  if(!(info_type == 1 || info_type == 2))
    return;// 1 = HA, 2 = FA
  if(size < 2)
    return;

  Gate* xor_root = getGateByNodeId(nodes[size - 2]);
  Gate* maj_root = getGateByNodeId(nodes[size - 1]);

  const auto& in_phase = adder.input_phase;
  const auto& out_phase = adder.output_phase;

  int K = (info_type == 1) ? 2 : 3;
  if(in_phase.size() == 2)
    K = 2;
  if(in_phase.size() == 3)
    K = 3;
  if(size < K + 2)
    return;

  std::vector<Gate*> inputs;
  inputs.reserve(K);
  for(int i = 0; i < K; ++i)
    inputs.push_back(getGateByNodeId(nodes[i]));

  // Helper: add variable with phase to polynomial
  auto add_phase_var
    = [](LinPolynomial& lp, Gate* g, int phase, int32_t coeff) {
        if(!g)
          return;
        if(!phase) {
          lp.add_term(g, coeff);
        } else {
          // phase 1 means we encode (1 - g)
          lp.add_term(nullptr, coeff);
          lp.add_term(g, -coeff);
        }
      };

  // Output phase convention
  int phase_xor = out_phase.size() > 0 ? static_cast<int>(out_phase[0]) : 0;
  int phase_maj = out_phase.size() > 1 ? static_cast<int>(out_phase[1]) : 0;

  // Build input forms according to input phases
  LinPolynomial p3, p4, p5;
  for(int i = 0; i < K; ++i) {
    int ph = in_phase.size() > static_cast<size_t>(i)
               ? static_cast<int>(in_phase[i])
               : 0;
    if(i == 0)
      add_phase_var(p3, inputs[0], ph, +1);
    if(i == 1)
      add_phase_var(p4, inputs[1], ph, +1);
    if(i == 2)
      add_phase_var(p5, inputs[2], ph, +1);
  }

  // MAJ polynomial
  if(maj_root) {
    LinPolynomial p1, p2;
    add_phase_var(p1, maj_root, phase_maj, +1);
    add_phase_var(p2, xor_root, phase_xor, +1);

    LinPolynomial p_maj = (p1 * 2) + p2 - p3 - p4;
    if(info_type == 2)
      p_maj = p_maj - p5;// FA has 3 inputs

    Gate* g = p_maj.get_leading_term().first;
    if(g)
      linpoly_cache_add(g, p_maj);
  }
}

/*------------------------------------------------------------------------*/
/**
 * New builder for XOR-only structures coming from XorRemain.
 * The logic of the XOR polynomial is consistent with the XOR part
 * of build_linpoly_from_lingates_hafa.
 *
 * Assumptions on XorRemain:
 *  - nodes: list of input nodes (usually input_nodes from JSON).
 *  - tt: encoded truth table (1 byte), currently not used here.
 *  - input_phase: phase for each input node.
 *  - output_phase: single-element array, phase for XOR output.
 *
 * We currently assume the last node in 'nodes' is the XOR output root.
 * If your JSON provides a dedicated 'xor_root', you can change it here.
 */
void
build_linpoly_from_lingates_xor(const XorRemain& xr) {
  const auto& nodes = xr.nodes;
  const int K = static_cast<int>(nodes.size());
  // For a 2-input XOR we at least need: 1 root + 2 inputs
  if(K < 3)
    return;

  // According to the convention:
  //   nodes[0]     : XOR root
  //   nodes[K-2]   : input 0
  //   nodes[K-1]   : input 1
  Gate* xor_root = getGateByNodeId(nodes[0]);
  if(!xor_root)
    return;
  msg("building XOR linpoly for xor_remain id=%d, root=%s, cone_id=%d, NN-1=%d",
      xr.id,
      xor_root->get_var_name().c_str(),
      xor_root->get_cone_id(),
      NN - 1);
  if(xor_root->get_cone_id() != (unsigned int)(NN - 1))
    return;

  const auto& in_phase = xr.input_phase;  // size >= 2 is ideal
  const auto& out_phase = xr.output_phase;// size >= 1 is ideal

  int phase_xor = !out_phase.empty() ? static_cast<int>(out_phase[0]) : 0;

  // Collect input gates: last two nodes
  std::vector<Gate*> inputs;
  inputs.reserve(2);
  Gate* in0 = getGateByNodeId(nodes[K - 2]);
  Gate* in1 = getGateByNodeId(nodes[K - 1]);
  if(!in0 || !in1)
    return;
  inputs.push_back(in0);
  inputs.push_back(in1);

  auto add_phase_var
    = [](LinPolynomial& lp, Gate* g, int phase, int32_t coeff) {
        if(!g)
          return;
        if(!phase) {
          // un-negated:  coeff * g
          lp.add_term(g, coeff);
        } else {
          // negated literal: coeff * (1 - g) = coeff*1 - coeff*g
          lp.add_term(nullptr, coeff);// constant term
          lp.add_term(g, -coeff);     // -coeff * g
        }
      };

  // p_xor = (phase_xor ? (1 - xor_root) : xor_root)
  LinPolynomial p_xor;
  add_phase_var(p_xor, xor_root, phase_xor, +1);

  // rhs_xor = sum of phased inputs
  LinPolynomial rhs_xor;
  for(int i = 0; i < 2; ++i) {
    int ph = (in_phase.size() > static_cast<size_t>(i))
               ? static_cast<int>(in_phase[i])
               : 0;
    LinPolynomial pi;
    add_phase_var(pi, inputs[i], ph, +1);
    rhs_xor = rhs_xor + pi;
  }

  // p = p_xor - rhs_xor
  LinPolynomial p = p_xor - rhs_xor;

  Gate* glead = p.get_leading_term().first;
  p.print();
  // glead->set_gatepoly(p);
  if(glead)
    linpoly_cache_add(glead, p);
}

/*------------------------------------------------------------------------*/
/**
 * mark_fahaxor for FullAdder list.
 * This marks HA/FA structure and builds corresponding XOR/MAJ polynomials.
 */
bool
mark_fahaxor(const std::vector<FullAdder>& fa_list, int /*fa_flag*/) {
  constexpr int num_inputs = 3;

  // Phase 1: Serial marking (modifies Gate objects)
  for(const auto& f : fa_list) {
    mark_one_faha_instance(f, num_inputs, 2, "full adder");
    build_linpoly_from_lingates_hafa(2, f);
  }

  return true;
}

/*------------------------------------------------------------------------*/
/**
 * mark_fahaxor for HalfAdder list.
 * This marks HA/FA structure and builds corresponding XOR/MAJ polynomials.
 */
bool
mark_fahaxor(const std::vector<HalfAdder>& ha_list, int /*fa_flag*/) {
  constexpr int num_inputs = 2;

  // Phase 1: Serial marking (modifies Gate objects)
  for(const auto& h : ha_list) {
    mark_one_faha_instance(h, num_inputs, 1, "half adder");
    build_linpoly_from_lingates_hafa(1, h);
  }

  return true;
}

/*----------------------------------------------------------------------*/
/**
 * mark_fahaxor for XOR-only structures.
 * It marks XOR instances and builds XOR polynomials.
 */
bool
mark_fahaxor(const std::vector<XorRemain>& xor_list) {
  for(const auto& xr : xor_list) {
    mark_one_xor_instance(xr);
    build_linpoly_from_lingates_xor(xr);
  }

  return true;
}

/*----------------------------------------------------------------------*/
/**
 * Build LinPolynomials from affine (linear) cuts and add to cache.
 *
 * For an affine function with k inputs, the ANF is:
 *   root = constant_coeff XOR (coeff_0 * leaf_0) XOR ... XOR (coeff_{k-1} *
 * leaf_{k-1})
 *
 * This is a XOR relation in the Boolean domain. To convert to integer domain
 * for check_linpoly_cache verification, we need to handle this carefully.
 *
 * For a 2-input XOR: s = a XOR b
 *   Integer relation: s = a + b - 2*ab (but this is non-linear)
 *   Alternative: We use the subcircuit's truth table to derive the correct
 *   linear polynomial via linear_polies_from_truth_table.
 *
 * This function builds a Subcircuit from the LinearCut and uses
 * linear_polies_from_truth_table to derive the correct linear polynomials
 * that will pass check_linpoly_cache verification.
 */
void
build_linpoly_from_linear_cuts(const std::vector<LinearCut>& cuts,
                               int32_t prime) {
  int count = 0;
  for(const auto& cut : cuts) {
    Gate* root_gate = getGateByNodeId(cut.root);
    if(!root_gate)
      continue;

    // Skip if no internal nodes or leaves
    if(cut.nodes.empty() || cut.leaves.empty())
      continue;

    // Skip cuts with too many inputs (truth table would be too large)
    if(cut.leaves.size() > 6)
      continue;

    // Build Subcircuit from cut.nodes and cut.leaves
    Subcircuit sc;

    // Build nodes list from cut.nodes (internal nodes between leaves and
    // root)
    for(int node_id : cut.nodes) {
      Gate* g = getGateByNodeId(node_id);
      if(g) {
        sc.nodes.push_back(g);
      }
    }
    // Add root to nodes if not already included
    if(std::find(sc.nodes.begin(), sc.nodes.end(), root_gate)
       == sc.nodes.end()) {
      sc.nodes.push_back(root_gate);
    }

    // Build inputs list from cut.leaves
    for(int leaf_id : cut.leaves) {
      Gate* g = getGateByNodeId(leaf_id);
      if(g) {
        sc.inputs.push_back(g);
      }
    }

    // Sort nodes by descending level (required by subcircuit evaluation)
    std::sort(
      sc.nodes.begin(), sc.nodes.end(), [](const Gate* a, const Gate* b) {
        return a->get_level() > b->get_level();
      });

    // Skip if subcircuit is invalid
    if(sc.inputs.empty() || sc.nodes.empty())
      continue;

    // Print cut info
    if(verbose >= 1) {
      std::cout << "  [LinearCut] root=" << cut.root << " leaves=[";
      for(size_t j = 0; j < cut.leaves.size(); j++) {
        if(j)
          std::cout << ",";
        std::cout << cut.leaves[j];
      }
      std::cout << "] nodes=[";
      for(size_t j = 0; j < cut.nodes.size(); j++) {
        if(j)
          std::cout << ",";
        std::cout << cut.nodes[j];
      }
      std::cout << "] tt=0x" << std::hex << cut.truth << std::dec << std::endl;

      std::cout << "    Subcircuit: " << sc.nodes.size() << " nodes, "
                << sc.inputs.size() << " inputs" << std::endl;
    }

    // Use linear_polies_from_truth_table to derive correct linear polynomials
    std::vector<LinPolynomial> polies;
    bool has_candidate
      = false;// linear_polies_from_truth_table(sc, prime, polies);

    if(verbose >= 1) {
      std::cout << "    linear_polies_from_truth_table: has_candidate="
                << (has_candidate ? "true" : "false")
                << ", polies.size()=" << polies.size() << std::endl;
    }

    // Add each polynomial to the cache
    for(size_t p_idx = 0; p_idx < polies.size(); p_idx++) {
      LinPolynomial& poly = polies[p_idx];
      Gate* glead = poly.get_leading_term().first;

      if(glead) {
        if(verbose >= 1) {
          std::cout << "    polies[" << p_idx
                    << "] (lead=" << glead->get_var_name() << "): ";
          poly.print();
          std::cout << std::endl;
        }

        // Add to linpoly cache
        linpoly_cache_add(glead, poly);
        count++;
      }
    }

    // Clean up after linear_polies_from_truth_table
    // cleanup();
  }

  msg("Built %d LinPolynomials from %zu affine cuts", count, cuts.size());
}

/*----------------------------------------------------------------------*/
void
reset_levels_to_orig() {
  for(unsigned i = 0; i < num_gates; i++) {
    Gate* g = gates[i];
    if(!g)
      continue;
    if(g->get_orig_level() != 0) {
      g->set_level(g->get_orig_level());
      LOGDBG("reset gate %s to original level %u",
             g->get_var_name().c_str(),
             g->get_level());
    }
  }
}
/*----------------------------------------------------------------------*/

void
adjust_cones_by_merging() {
  int adjustments = 0;

  for(unsigned i = 0; i < num_gates; i++) {
    Gate* g = gates[i];
    if(!g)
      continue;

    // Skip extension gates
    if(g->is_extension())
      continue;

    // Skip input gates
    if(g->get_input())
      continue;

    int current_cone_id = g->get_cone_id();
    if(current_cone_id == -1)
      continue;

    unsigned ch0_lit = g->get_ch0();
    unsigned ch1_lit = g->get_ch1();

    // Collect cone IDs of non-input children
    int child0_cone = -1;
    int child1_cone = -1;
    bool has_child0 = false;
    bool has_child1 = false;

    // Check child 0
    if(ch0_lit > 1) {
      Gate* c0 = gate(ch0_lit);
      if(c0 && !c0->get_input()) {
        child0_cone = c0->get_cone_id();
        has_child0 = true;
      }
    }

    // Check child 1
    if(ch1_lit > 1) {
      Gate* c1 = gate(ch1_lit);
      if(c1 && !c1->get_input()) {
        child1_cone = c1->get_cone_id();
        has_child1 = true;
      }
    }

    // If both children are non-input gates and both have smaller cone_ids
    if(has_child0 && has_child1) {
      if(child0_cone != -1 && child1_cone != -1 && child0_cone < current_cone_id
         && child1_cone < current_cone_id) {
        // Update to the larger of the two children cone_ids
        int new_cone_id
          = (child0_cone > child1_cone) ? child0_cone : child1_cone;
        g->set_cone_id(new_cone_id);
        adjustments++;

        LOGDBG("adjusted cone_id for gate %s from %d to %d (children: %d, %d)",
               g->get_var_name().c_str(),
               current_cone_id,
               new_cone_id,
               child0_cone,
               child1_cone);
      }
    }
  }

  msg("adjusted cone_ids for %d gates by merging with children", adjustments);
}
/*----------------------------------------------------------------------*/

static void
reassign_levels_by_cone_id() {
  // Collect only FSA-marked gates
  std::vector<Gate*> fsa_gates;
  for(unsigned i = 0; i < num_gates; i++) {
    Gate* g = gates[i];
    if(!g || !g->is_fsa())
      continue;
    if(g->get_dist() < 7)
      continue;
    fsa_gates.push_back(g);
  }

  if(fsa_gates.empty())
    return;

  // Sort by decreasing cone_id, then by decreasing level
  std::sort(
    fsa_gates.begin(), fsa_gates.end(), [](const Gate* a, const Gate* b) {
      int cone_a = a->get_cone_id();
      int cone_b = b->get_cone_id();
      if(cone_a != cone_b)
        return cone_a > cone_b;// decreasing cone_id
      return a->get_level()
             > b->get_level();// decreasing level within same cone
    });

  // Assign levels in decreasing order: num_gates + num_fsa_gates - i
  unsigned num_fsa = fsa_gates.size();
  for(unsigned i = 0; i < num_fsa; i++) {
    unsigned new_level = num_gates + num_fsa - i;
    fsa_gates[i]->set_orig_level(fsa_gates[i]->get_level());
    fsa_gates[i]->set_level(new_level);
    LOGDBG("reassigned FSA gate %s to level %u (cone_id=%d)",
           fsa_gates[i]->get_var_name().c_str(),
           new_level,
           fsa_gates[i]->get_cone_id());
  }

  msg("reassigned levels for %u FSA gates by cone_id (levels %u to %u)",
      num_fsa,
      num_gates + 1,
      num_gates + num_fsa);
}

bool
mark_fsa(std::vector<int> nodes, std::vector<int> inp_nodes) {
  for(auto& node_id : inp_nodes) {
    if(node_id < 1)
      continue;
    Gate* g = gate(node_id * 2);
    g->set_fsa_inp();
    g->reset_fsa();
    // mark FSA input nodes
  }

  for(auto& node_id : nodes) {
    if(node_id < 1)
      continue;
    Gate* g = gate(node_id * 2);
    if(g->is_fsa_inp())
      continue;
    g->set_fsa();
    fsa_size++;
    // mark FSA internal nodes
  }

  if(fsa_size > 10000) {
    reassign_levels_by_cone_id();
  }

  return true;
}

/*------------------------------------------------------------------------*/
static void
mark_coi_internal(int i, Gate* g) {
  if(g->get_cone_id() != -1)
    return;
  if(g->get_input())
    return;
  g->set_cone_id(i);

  unsigned l = g->get_ch0();
  unsigned r = g->get_ch1();
  if(l > 1)
    mark_coi_internal(i, gate(l));
  if(r > 1)
    mark_coi_internal(i, gate(r));
}

/*------------------------------------------------------------------------*/
static void
mark_coi() {
  for(unsigned i = M - 1; i < M - 1 + MM; i++) {
    Gate* g = gates[i];
    mark_coi_internal(i + 1 - M, g);
  }
}

//*------------------------------------------------------------------------*/
void
init_gates(int ext) {
  allocate_gates(ext);
  mark_aig_outputs();
  set_parents_and_children();
  set_distance_to_input();
  mark_coi();
  set_levels();

  if(verbose >= 4) {
    print_node_statistics();
  }
}
/*------------------------------------------------------------------------*/

Gate*
gate(int lit) {
  if(lit <= 0)
    return gates[M - lit - 1];
  if(lit < 2)
    return 0;
  return gates[lit / 2 - 1];
}

/*------------------------------------------------------------------------*/

void
delete_gates() {
  for(unsigned i = 0; i < num_gates; i++) {
    if(gates[i])
      delete(gates[i]);
  }
  delete[] gates;
}

// set once decompose_all_gate_polynomials() ran
static bool decomposed = false;

/*------------------------------------------------------------------------*/
void
reset_gates() {
  if(gates)
    delete_gates();
  gates = nullptr;
  num_gates = 0;
  size_gates = 0;
  gates_by_dist.clear();
  add_var = 0;
  max_dist = 0;
  fsa_size = 0;
  decomposed = false;
}

void
reassign_levels_by_cone_order() {
  // find max cone id
  int max_cone = -1;
  for(unsigned i = 0; i < num_gates; ++i) {
    Gate* g = gates[i];
    if(!g)
      continue;
    int cid = g->get_cone_id();
    if(cid > max_cone)
      max_cone = cid;
  }

  unsigned level = 1;
  // iterate cones from -1 up to max_cone
  for(int cid = -1; cid <= max_cone; ++cid) {
    // collect gates in this cone with their indices
    std::vector<std::pair<int, Gate*>> vec;
    vec.reserve(64);
    for(int i = 0; i < (int)num_gates; ++i) {
      Gate* g = gates[i];
      if(!g)
        continue;
      if(g->get_cone_id() == cid)
        vec.emplace_back(i, g);
    }
    if(vec.empty())
      continue;

    // sort by index ascending (smallest gate first)
    std::sort(vec.begin(), vec.end(), [](const auto& a, const auto& b) {
      return a.first < b.first;
    });

    // assign increasing levels
    for(const auto& p : vec) {
      Gate* g = p.second;
      g->set_orig_level(g->get_level());
      g->set_level(level++);
      LOGDBG("reassigned gate %s (cone %d, idx=%d) to level %u",
             g->get_var_name().c_str(),
             cid,
             p.first,
             g->get_level());
    }
  }

  msg("reassigned levels by cone order (-1..%d), assigned %u levels",
      max_cone,
      level - 1);
}

std::vector<std::vector<Gate*>>
get_gates_by_slice() {
  // Find max slice/cone id
  int max_slice = -1;
  for(unsigned i = 0; i < num_gates; ++i) {
    Gate* g = gates[i];
    if(!g)
      continue;
    int cid = g->get_cone_id();
    if(cid > max_slice)
      max_slice = cid;
  }

  // Create a vector of vectors, one per slice
  std::vector<std::vector<Gate*>> slices(max_slice + 1);

  // Populate each slice with its gates
  for(int i = num_gates - 1; i >= 0; --i) {
    Gate* g = gates[i];
    if(!g)
      continue;
    int cid = g->get_cone_id();
    if(cid >= 0 && cid <= max_slice) {
      slices[cid].push_back(g);
    }
  }

  return slices;
}

std::vector<Gate*>
get_gates_by_decreasing_level() {
  std::vector<Gate*> result;

  int max_cone = -1;
  for(unsigned i = 0; i < num_gates; ++i) {
    Gate* g = gates[i];
    if(!g)
      continue;
    int cid = g->get_cone_id();
    if(cid > max_cone)
      max_cone = cid;
  }

  if(max_cone < 0)
    return result;

  // iterate cones in descending order (include cone -1 as special case)
  for(int cid = max_cone; cid >= 0; --cid) {
    // iterate gates in reverse index order so gates inside cone are reversed
    for(int i = (int)num_gates - 1; i >= 0; --i) {
      Gate* g = gates[i];
      if(!g)
        continue;
      if(g->get_cone_id() == cid) {
        result.push_back(g);
      }
    }
  }

  return result;
}

void
print_cones() {
  std::vector<Gate*> gates_in_order = get_gates_by_decreasing_level();

  int max_cone = -1;
  for(unsigned i = 0; i < num_gates; ++i) {
    Gate* g = gates[i];
    if(!g)
      continue;
    int cid = g->get_cone_id();
    if(cid > max_cone)
      max_cone = cid;
  }

  if(max_cone < 0)
    return;

  // iterate cones in descending order (include cone -1 as special case)
  for(int cid = max_cone; cid >= 0; --cid) {
    msg("slice %i", cid);
    unsigned prev_level = UINT_MAX;
    // iterate through gates in order and filter by cone_id
    for(Gate* g : gates_in_order) {
      if(g->get_cone_id() == cid) {
        unsigned lvl = g->get_level();
        if(prev_level != UINT_MAX && lvl >= prev_level) {
          msg("cone %i: levels not decreasing at gate %s (level=%u >= prev=%u)",
              cid,
              g->get_var_name().c_str(),
              lvl,
              prev_level);
        }
        msg("%s", g->get_var_name().c_str(), lvl);
        prev_level = lvl;
      }
    }
    msg("");
  }
}

// Helper: Extract the substitution expression for gate n from its polynomial
// If n_poly = n - expr, returns expr (with positive coefficients)
Polynomial
extract_substitution_expr(Gate* n, const Polynomial& n_poly) {
  Polynomial substitution;

  LOGDBG("Extracting substitution for gate %s from polynomial",
         n->get_var_name().c_str());

  for(const auto& [mono, coeff] : n_poly.multimod_terms) {
    // Skip the monomial that is just gate n itself
    if(mono.vars.size() == 1 && mono.vars[0] == n) {
      LOGDBG("  Skipping self term: %s", n->get_var_name().c_str());
      continue;
    }

    // Add all other terms with negated coefficient
    CoeffVec neg_coeff = coeff.negate();
    substitution.add_term(mono, neg_coeff);

    LOGDBG("  Adding term with vars.size=%zu", mono.vars.size());
  }

  return substitution;
}

// Helper: Substitute gate n with its expression in a single monomial
// Returns the resulting polynomial from expanding mono * substitution
Polynomial
substitute_gate_in_monomial(Gate* n,
                            const Monomial& mono,
                            const CoeffVec& coeff,
                            const Polynomial& substitution) {
  Polynomial result;

  // Build remainder monomial: all vars in mono except first occurrence of n
  Monomial remainder;
  bool found_n = false;
  for(Gate* g : mono.vars) {
    if(g == n && !found_n) {
      found_n = true;
      LOGDBG("    Replacing occurrence of %s in monomial",
             n->get_var_name().c_str());
      continue;
    }
    remainder.vars.push_back(g);
  }

  if(!found_n) {
    // This monomial doesn't contain n, shouldn't happen

    CoeffVec coeff_copy = coeff;
    result.add_term(mono, coeff_copy);
    return result;
  }

  // Multiply each term in substitution by remainder and original coefficient
  for(const auto& [sub_mono, sub_coeff] : substitution.multimod_terms) {
    Monomial new_mono = sub_mono;

    // Append remainder variables
    new_mono.vars.insert(
      new_mono.vars.end(), remainder.vars.begin(), remainder.vars.end());

    // Sort and remove duplicates (Boolean algebra: x*x = x)
    new_mono.sort_vars();
    auto last = std::unique(new_mono.vars.begin(), new_mono.vars.end());
    new_mono.vars.erase(last, new_mono.vars.end());

    // Compute final coefficient: sub_coeff * original_coeff
    CoeffVec final_coeff = sub_coeff * coeff;
    result.add_term(new_mono, final_coeff);
  }

  return result;
}

// Helper: Substitute all occurrences of gate n in a polynomial with its
// expression
Polynomial
substitute_gate_in_polynomial(Gate* n,
                              const Polynomial& parent_poly,
                              const Polynomial& substitution) {
  Polynomial result;
  int subst_count = 0;

  LOGDBG("Substituting gate %s in parent polynomial",
         n->get_var_name().c_str());

  for(const auto& [mono, coeff] : parent_poly.multimod_terms) {
    // Check if this monomial contains gate n
    bool contains_n = false;
    for(Gate* g : mono.vars) {
      if(g == n) {
        contains_n = true;
        break;
      }
    }

    if(!contains_n) {
      // Monomial doesn't contain n, keep as-is
      CoeffVec coeff_copy = coeff;
      result.add_term(mono, coeff_copy);

    } else {
      // Monomial contains n: substitute it
      subst_count++;
      Polynomial expanded
        = substitute_gate_in_monomial(n, mono, coeff, substitution);

      // Add all expanded terms to result
      for(const auto& [exp_mono, exp_coeff] : expanded.multimod_terms) {
        CoeffVec exp_coeff_copy = exp_coeff;
        result.add_term(exp_mono, exp_coeff_copy);
      }
    }
  }

  LOGDBG("Substitution complete: %d monomials", subst_count);

  return result;
}

/*------------------------------------------------------------------------*/
// Integer polynomial substitution functions
// These use simplified integer coefficients with Boolean algebra (x*x = x)
/*------------------------------------------------------------------------*/

// Helper: Substitute gate n with its expression in a single monomial
// Returns the resulting polynomial from expanding mono * substitution
IntPolynomial
int_substitute_gate_in_monomial(Gate* n,
                                const IntMonomial& mono,
                                int32_t coeff,
                                const IntPolynomial& substitution) {
  IntPolynomial result;

  // Build remainder monomial: all vars in mono except first occurrence of n
  IntMonomial remainder;
  bool found_n = false;
  for(Gate* g : mono.vars) {
    if(g == n && !found_n) {
      found_n = true;
      LOGDBG("    Replacing occurrence of %s in monomial",
             n->get_var_name().c_str());
      continue;
    }
    remainder.vars.push_back(g);
  }

  if(!found_n) {
    // This monomial doesn't contain n, shouldn't happen
    result.add_term(mono, coeff);
    return result;
  }

  // Multiply each term in substitution by remainder and original coefficient
  for(const auto& [sub_mono, sub_coeff] : substitution.terms) {
    IntMonomial new_mono = sub_mono;

    // Append remainder variables
    new_mono.vars.insert(
      new_mono.vars.end(), remainder.vars.begin(), remainder.vars.end());

    // Sort and remove duplicates (Boolean algebra: x*x = x)
    new_mono.sort_vars();
    auto last = std::unique(new_mono.vars.begin(), new_mono.vars.end());
    new_mono.vars.erase(last, new_mono.vars.end());

    // Compute final coefficient: sub_coeff * original_coeff
    int32_t final_coeff = sub_coeff * coeff;
    result.add_term(new_mono, final_coeff);
  }

  return result;
}

// Helper: Substitute all occurrences of gate n in a polynomial with its
// expression
IntPolynomial
int_substitute_gate_in_polynomial(Gate* n,
                                  const IntPolynomial& parent_poly,
                                  const IntPolynomial& substitution) {
  IntPolynomial result;
  int subst_count = 0;

  LOGDBG("Substituting gate %s in parent polynomial",
         n->get_var_name().c_str());

  for(const auto& [mono, coeff] : parent_poly.terms) {
    // Check if this monomial contains gate n
    bool contains_n = false;
    for(Gate* g : mono.vars) {
      if(g == n) {
        contains_n = true;
        break;
      }
    }

    if(!contains_n) {
      // Monomial doesn't contain n, keep as-is
      result.add_term(mono, coeff);
    } else {
      // Monomial contains n: substitute it
      subst_count++;
      IntPolynomial expanded
        = int_substitute_gate_in_monomial(n, mono, coeff, substitution);

      // Add all expanded terms to result
      for(const auto& [exp_mono, exp_coeff] : expanded.terms) {
        result.add_term(exp_mono, exp_coeff);
      }
    }
  }

  LOGDBG("Substitution complete: %d monomials", subst_count);

  return result;
}

// Public API: Substitute gate n in parent polynomial with its expression
IntPolynomial
substitute_gate(Gate* n,
                const IntPolynomial& parent_poly,
                const IntPolynomial& substitution_expr) {
  return int_substitute_gate_in_polynomial(n, parent_poly, substitution_expr);
}

void
decompose_all_gate_polynomials(Gate* border) {
  if(decomposed)
    return;
  decomposed = true;

  if(verbose > 0)
    msg("eliminate single occs");

  auto slices = get_gates_by_slice();
  int begin = NN - 1;

  int counter = 0;
  int total_rounds = 0;

  // Iterate through slices from begin down to 0
  for(int slice_idx = begin; slice_idx >= 0 && slice_idx < (int)slices.size();
      slice_idx--) {
    if(verbose > 0)
      msg("Processing slice %d with %zu gates",
          slice_idx,
          slices[slice_idx].size());

    bool change = true;
    int round = 0;

    while(change) {
      change = false;
      round++;
      total_rounds++;

      LOGDBG("  Round %d: %zu gates remaining in slice",
             round,
             slices[slice_idx].size());

      for(auto it = slices[slice_idx].begin(); it != slices[slice_idx].end();) {
        Gate* n = *it;

        if(n->get_level() > border->get_level()) {
          ++it;
          continue;
        }

        // Skip invalid or non-decomposable gates
        if(!n) {
          ++it;
          continue;
        }
        if(n->get_input()) {
          ++it;
          continue;
        }
        if(n->is_extension()) {
          ++it;
          continue;
        }
        if(n->is_elim()) {
          ++it;
          continue;
        }
        if(n->parent_size() != 1) {
          LOGDBG("    Gate %s has %d parents, skipping",
                 n->get_var_name().c_str(),
                 n->parent_size());
          ++it;
          continue;
        }
        if(n->get_carry_gate() > 0) {
          ++it;
          continue;
        }

        // Get the single parent
        const auto& parents = n->get_parents();
        if(parents.empty()) {
          ++it;
          continue;
        }
        Gate* parent = gate(*parents.begin());
        if(!parent) {
          ++it;
          continue;
        }

        if(verbose > 0)
          msg("  Decomposing %s (parent: %s)",
              n->get_var_name().c_str(),
              parent->get_var_name().c_str());

        try {
          // Get parent's polynomial (using cached IntPolynomial if available)
          IntPolynomial parent_poly;
          if(parent->has_gatepoly_int()) {
            parent_poly = *(parent->get_gatepoly_int());
          } else {
            // Only compute IntPolynomial if children are valid
            if(parent->get_ch0() > 1 || parent->get_ch1() > 1) {
              parent_poly = IntPolynomial::get_gate_poly(parent);
              parent->set_gatepoly_int(parent_poly);
            } else {
              // Fall back to empty polynomial if children not set
              parent_poly = IntPolynomial(parent);
              parent->set_gatepoly_int(parent_poly);
            }
          }

          // Get n's polynomial (using cached IntPolynomial if available)
          IntPolynomial n_poly;
          if(n->has_gatepoly_int()) {
            n_poly = *(n->get_gatepoly_int());
          } else {
            // Only compute IntPolynomial if children are valid
            if(n->get_ch0() > 1 || n->get_ch1() > 1) {
              n_poly = IntPolynomial::get_gate_poly(n);
              n->set_gatepoly_int(n_poly);
            } else {
              // Fall back to polynomial with just the gate itself
              n_poly = IntPolynomial(n);
              n->set_gatepoly_int(n_poly);
            }
          }

          // Extract substitution: if n_poly = n - expr, get expr
          IntPolynomial n_substitution = n_poly.extract_substitution_expr(n);

          // Substitute n in parent_poly with n_substitution (using integer
          // arithmetic)
          IntPolynomial new_parent_poly
            = substitute_gate(n, parent_poly, n_substitution);

          // Store result back on parent (cache as IntPolynomial)
          parent->set_gatepoly_int(new_parent_poly);

          // Also convert to multimod for compatibility with other code paths
          Polynomial multimod_result = new_parent_poly.to_multimod();

          parent->set_gatepoly_nonlinear(multimod_result);

          // Update parent/child relationships
          unsigned parent_id = parent->get_var();
          unsigned n_id = n->get_var();

          // Disconnect parent from n
          n->remove_parent(parent_id);
          n->set_elim();

          // Connect parent to n's children
          for(unsigned lit : { n->get_ch0(), n->get_ch1() }) {
            if(lit <= 1)
              continue;
            Gate* n_child = gate(lit);
            if(!n_child)
              continue;

            n_child->remove_parent(n_id);
            if(n_child->get_parents().find(parent_id)
               == n_child->get_parents().end()) {
              n_child->add_parent(parent_id);
            }
          }

          counter++;
          change = true;

          // Remove from slice and advance
          it = slices[slice_idx].erase(it);
        } catch(const std::exception& e) {
          LOGDBG(
            "Error decomposing %s: %s", n->get_var_name().c_str(), e.what());
          msg("Warning: Failed to decompose %s: %s",
              n->get_var_name().c_str(),
              e.what());
          ++it;
        }
      }
    }

    if(verbose > 0)
      msg("  Slice %d: %d rounds, %zu gates remaining",
          slice_idx,
          round,
          slices[slice_idx].size());
  }

  if(verbose > 0)
    msg("decomposed %i variables in %d total rounds", counter, total_rounds);
}

/*------------------------------------------------------------------------*/

}// namespace talisman

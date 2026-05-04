/*------------------------------------------------------------------------*/
/*! \file parser.cpp
    \brief contains functions necessary to parse the AIG

  Part of TalisMan
*/
/*------------------------------------------------------------------------*/
#include "parser.hpp"

#include <algorithm>
#include <cstdint>
#include <cctype>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>

/*------------------------------------------------------------------------*/
// ERROR CODES:
static int err_parsing      = 20;   // general parsing error
static int err_latches      = 21;   // cannot handle latches
static int err_no_inputs    = 22;   // no inputs
static int err_odd_inputs   = 23;   // odd inputs
static int err_no_outputs   = 24;   // no outputs
static int err_odd_outputs  = 25;   // odd outputs
static int err_wrong_outputs= 26;   // wrong number of outputs

/*------------------------------------------------------------------------*/
using namespace talisman;

namespace talisman {

bool mult_spec  = false;  ///< true if multiplier spec is used
bool miter_spec = false;  ///< true if miter spec is used

/*------------------------------------------------------------------------*/
// Helper: match a literal with an AND gate in the AIG model
static bool match_and(unsigned lhs, unsigned rhs0, unsigned rhs1) {
  if (lhs == aiger_false)
    return 0;
  if (aiger_sign(lhs))
    return 0;
  assert(lhs != aiger_true);

  aiger_and* and1 = is_model_and(lhs);
  if (!and1)
    return 0;
  if (and1->rhs0 == rhs0 && and1->rhs1 == rhs1)
    return 1;
  if (and1->rhs0 == rhs1 && and1->rhs1 == rhs0)
    return 1;
  return 0;
}

/*------------------------------------------------------------------------*/
// Determine input ordering for multiplier specifications.
static void determine_multiplier_input_order() {
  unsigned s0 = 0, sl = NN - 1;
  if (match_and(slit(0), get_model_inputs_lit(0), get_model_inputs_lit(1))) {
    a0   = 0;      al   = NN - 2; ainc = 2;
    b0   = 1;      bl   = NN - 1; binc = 2;
    msg("assuming ordering as in BTOR generated benchmarks");
  } else {
    a0   = 0;            al   = NN / 2 - 1; ainc = 1;
    b0   = NN / 2;       bl   = NN - 1;     binc = 1;
    msg("assuming ordering as in the ABC generated or AOKI benchmarks");
  }
  if (verbose >= 2) {
    if (NN == 2) {
      msg("a[0] = input[%d]", a0);
      msg("b[0] = input[%d]", b0);
      msg("s[0] = output[%d]", s0);
    } else if (NN == 4) {
      msg("(a[0], a[1]) =(input[%d], input[%d])", a0, al);
      msg("(b[0], b[1]) =(input[%d], input[%d])", b0, bl);
      msg("(s[0], ..., s[3]) =(output[%d], ..., output[%d])", s0, sl);
    } else if (NN == 6) {
      msg("(a[0], a[1], a[2]) =(input[%d], input[%d], input[%d])",
          a0, a0 + ainc, al);
      msg("(b[0], b[1], b[2]) =(input[%d], input[%d], input[%d])",
          b0, b0 + binc, bl);
      msg("(s[0], ..., s[5]) =(output[%d], ..., output[%d])", s0, sl);
    } else {
      msg("(a[0], a[1], ..., a[%d]) =(input[%d], input[%d], ..., input[%d])",
          NN / 2 - 1, a0, a0 + ainc, al);
      msg("(b[0], b[1], ..., b[%d]) =(input[%d], input[%d], ..., input[%d])",
          NN / 2 - 1, b0, b0 + binc, bl);
      msg("(s[0], ..., s[%d]) =(output[%d], ..., output[%d])", NN - 1, s0, sl);
    }
  }
}

/*------------------------------------------------------------------------*/
// Initialize AIG structure and perform consistency checks.
static void init_aiger_with_checks() {
  if (get_model_num_latches())
    die(err_latches, "can not handle latches");
  if (!get_model_num_inputs())
    die(err_no_inputs, "no inputs");
  if (!get_model_num_outputs())
    die(err_no_outputs, "no outputs");

  if (mult_spec) {
    if (get_model_num_inputs() & 1)
      die(err_odd_inputs,
          "odd number of inputs not allowed in multiplier input");
    if (get_model_num_outputs() & 1)
      die(err_odd_outputs,
          "odd number of outputs not allowed in multiplier input");
    if (get_model_num_outputs() != get_model_num_inputs())
      die(err_wrong_outputs,
          "number of inputs and outputs does not match in multiplier input");
  } else if (miter_spec) {
    if (get_model_num_outputs() != 1)
      die(err_wrong_outputs, "too many outputs in miter input");
  }

  M  = get_model_maxvar() + 1;
  NN = get_model_num_inputs();
  MM = get_model_num_outputs();

  msg("MILOA %u %u %u %u %u",
      get_model_maxvar(),
      get_model_num_inputs(),
      get_model_num_latches(),
      get_model_num_outputs(),
      get_model_num_ands());
}

/*------------------------------------------------------------------------*/
// Top-level AIG parsing entry.
void parse_aig(const char* input_name) {
  init_aig_parsing();

  msg("reading '%s'", input_name);
  const char* err = aiger_open_and_read_to_model(input_name);
  if (err)
    die(err_parsing, "error parsing '%s': %s", input_name, err);

  init_aiger_with_checks();
  if (mult_spec)
    determine_multiplier_input_order();
  msg("finished parsing '%s'", input_name);
}

/*------------------------------------------------------------------------*/
/*                     JSON helpers and adder.json parsing                */
/*------------------------------------------------------------------------*/

// Parse a simple integer array from a bracketed string: "[1, 2, 3]"
static std::vector<int> parse_int_array_from_brackets(const std::string& s) {
  std::vector<int> out;
  size_t i = 0, n = s.size();

  // Skip leading whitespace
  while (i < n && std::isspace((unsigned char)s[i])) ++i;

  // Optional '[' at the beginning
  if (i < n && s[i] == '[') ++i;

  while (i < n) {
    // Skip whitespace
    while (i < n && std::isspace((unsigned char)s[i])) ++i;
    if (i >= n) break;
    if (s[i] == ']') { ++i; break; }

    // Optional '+' sign (we do not support negative numbers here)
    if (s[i] == '+') ++i;
    if (i >= n || !std::isdigit((unsigned char)s[i])) {
      // Skip unexpected characters until next ',' or ']'
      while (i < n && s[i] != ',' && s[i] != ']') ++i;
      if (i < n && s[i] == ',') { ++i; continue; }
      if (i < n && s[i] == ']') { ++i; break; }
      continue;
    }

    int val = 0;
    while (i < n && std::isdigit((unsigned char)s[i])) {
      val = val * 10 + (s[i] - '0');
      ++i;
    }
    out.push_back(val);

    // Skip whitespace and comma
    while (i < n && std::isspace((unsigned char)s[i])) ++i;
    if (i < n && s[i] == ',') { ++i; continue; }
    if (i < n && s[i] == ']') { ++i; break; }
  }
  return out;
}

// Read a very simple JSON file that contains two integer arrays
//   "reverse fsa nodes": [...]
//   "bondary nodes":     [...]
bool read_adder_json(const char* filename,
                     std::vector<int>& out_reverse_fsa_nodes,
                     std::vector<int>& out_bondary_nodes) {
  out_reverse_fsa_nodes.clear();
  out_bondary_nodes.clear();

  std::ifstream fin(filename, std::ios::in | std::ios::binary);
  if (!fin) {
    warning("Cannot open JSON file '%s'.", filename);
    return false;
  }

  std::ostringstream oss;
  oss << fin.rdbuf();
  std::string content = oss.str();
  fin.close();

  // Helper to find and parse a key's array value
  auto find_key_array = [&](const char* key, std::vector<int>& out) -> bool {
    std::string pat = std::string("\"") + key + "\"";
    size_t pos = content.find(pat);
    if (pos == std::string::npos) return false;

    pos = content.find(':', pos + pat.size());
    if (pos == std::string::npos) return false;

    size_t lb = content.find('[', pos + 1);
    if (lb == std::string::npos) return false;

    size_t rb = content.find(']', lb + 1);
    if (rb == std::string::npos) return false;

    std::string arr = content.substr(lb, rb - lb + 1);
    out = parse_int_array_from_brackets(arr);
    return true;
  };

  bool ok1 = find_key_array("reverse fsa nodes", out_reverse_fsa_nodes);
  bool ok2 = find_key_array("bondary nodes", out_bondary_nodes);

  if (!ok1 && !ok2) {
    warning("JSON file '%s' does not contain expected keys.", filename);
    return false;
  }

  // Sort and deduplicate results
  auto sort_uniq = [](std::vector<int>& v) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
  };
  sort_uniq(out_reverse_fsa_nodes);
  sort_uniq(out_bondary_nodes);

  msg("Parsed JSON '%s': reverse fsa nodes = %zu, bondary nodes = %zu",
      filename,
      out_reverse_fsa_nodes.size(),
      out_bondary_nodes.size());
  return true;
}

/*------------------------------------------------------------------------*/
/*                         Generic JSON utilities                         */
/*------------------------------------------------------------------------*/

// Simple whitespace checker
static inline bool is_ws(char c) {
  return c == ' ' || c == '\t' || c == '\n' || c == '\r';
}

// In content, locate the array value of a given key and
// return its span [lb, rb] including the brackets.
static bool find_key_array_span(const std::string& content,
                                const std::string& key,
                                size_t& lb_out,
                                size_t& rb_out) {
  std::string pat = std::string("\"") + key + "\"";
  size_t pos = content.find(pat);
  if (pos == std::string::npos) return false;

  pos = content.find(':', pos + pat.size());
  if (pos == std::string::npos) return false;

  size_t lb = content.find('[', pos + 1);
  if (lb == std::string::npos) return false;

  int depth = 0;
  size_t rb = std::string::npos;
  for (size_t i = lb; i < content.size(); ++i) {
    if (content[i] == '[')      ++depth;
    else if (content[i] == ']') --depth;
    if (depth == 0) { rb = i; break; }
  }
  if (rb == std::string::npos) return false;

  lb_out = lb;
  rb_out = rb;
  return true;
}

// Given an array text "[ { .. }, { .. }, ... ]", slice out each "{...}" object.
static std::vector<std::string> slice_objects_from_array(const std::string& arr) {
  std::vector<std::string> objs;

  size_t p = 0;
  while (p < arr.size() && is_ws(arr[p])) ++p;
  if (p >= arr.size() || arr[p] != '[') return objs;

  ++p; // skip '['
  while (p < arr.size()) {
    // Skip whitespace and commas
    while (p < arr.size() && (is_ws(arr[p]) || arr[p] == ',')) ++p;
    if (p >= arr.size()) break;
    if (arr[p] == ']') break;

    if (arr[p] != '{') {
      // Not an object start, skip one char
      ++p;
      continue;
    }

    // Match object braces
    int depth = 0;
    size_t start = p;
    size_t end   = std::string::npos;
    for (size_t i = p; i < arr.size(); ++i) {
      if (arr[i] == '{')      ++depth;
      else if (arr[i] == '}') --depth;
      if (depth == 0) { end = i; break; }
    }
    if (end == std::string::npos) break;

    objs.emplace_back(arr.substr(start, end - start + 1));
    p = end + 1;
  }

  return objs;
}

// Parse "id": <int> from an object text
static bool parse_id_from_object(const std::string& obj, int& id_out) {
  std::string key = "\"id\"";
  size_t pos = obj.find(key);
  if (pos == std::string::npos) return false;

  pos = obj.find(':', pos + key.size());
  if (pos == std::string::npos) return false;

  size_t i = pos + 1;
  while (i < obj.size() && is_ws(obj[i])) ++i;

  bool neg = false;
  if (i < obj.size() && obj[i] == '-') {
    neg = true;
    ++i;
  }

  long val = 0;
  bool any = false;
  while (i < obj.size() && std::isdigit((unsigned char)obj[i])) {
    val = val * 10 + (obj[i] - '0');
    ++i;
    any = true;
  }
  if (!any) return false;

  id_out = neg ? (int)-val : (int)val;
  return true;
}

// From an object text, parse an integer array field, e.g. "nodes": [ ... ]
static bool parse_int_array_from_object(const std::string& obj,
                                        const std::string& key,
                                        std::vector<int>& out) {
  out.clear();
  std::string pat = std::string("\"") + key + "\"";
  size_t pos = obj.find(pat);
  if (pos == std::string::npos) return false;

  pos = obj.find(':', pos + pat.size());
  if (pos == std::string::npos) return false;

  size_t lb = obj.find('[', pos + 1);
  if (lb == std::string::npos) return false;

  int depth = 0;
  size_t rb = std::string::npos;
  for (size_t i = lb; i < obj.size(); ++i) {
    if (obj[i] == '[')      ++depth;
    else if (obj[i] == ']') --depth;
    if (depth == 0) { rb = i; break; }
  }
  if (rb == std::string::npos) return false;

  std::string inner = obj.substr(lb + 1, rb - (lb + 1));

  // Parse comma-separated integers
  size_t p = 0;
  while (p < inner.size()) {
    while (p < inner.size() && (is_ws(inner[p]) || inner[p] == ',')) ++p;
    if (p >= inner.size()) break;

    bool neg = false;
    if (inner[p] == '-') { neg = true; ++p; }

    long val = 0;
    bool any = false;
    while (p < inner.size() &&
           std::isdigit((unsigned char)inner[p])) {
      val = val * 10 + (inner[p] - '0');
      ++p;
      any = true;
    }
    if (any) out.push_back(neg ? (int)-val : (int)val);

    while (p < inner.size() && inner[p] != ',') ++p;
  }
  return true;
}

// Parse a hex string array field, e.g. "tt": ["0x69","0x8E"]
// into a vector<uint8_t>.
static bool parse_hex_tt_from_object(const std::string& obj,
                                     const std::string& key,
                                     std::vector<uint8_t>& out) {
  out.clear();
  std::string pat = std::string("\"") + key + "\"";
  size_t pos = obj.find(pat);
  if (pos == std::string::npos) return false;

  pos = obj.find(':', pos + pat.size());
  if (pos == std::string::npos) return false;

  size_t lb = obj.find('[', pos + 1);
  if (lb == std::string::npos) return false;

  int depth = 0;
  size_t rb = std::string::npos;
  for (size_t i = lb; i < obj.size(); ++i) {
    if (obj[i] == '[')      ++depth;
    else if (obj[i] == ']') --depth;
    if (depth == 0) { rb = i; break; }
  }
  if (rb == std::string::npos) return false;

  std::string inner = obj.substr(lb + 1, rb - (lb + 1));
  size_t p = 0;

  // Extract tokens inside quotes
  while (p < inner.size()) {
    while (p < inner.size() && (is_ws(inner[p]) || inner[p] == ',')) ++p;
    if (p >= inner.size()) break;
    if (inner[p] != '\"') {
      ++p;
      continue;
    }
    size_t q = inner.find('\"', p + 1);
    if (q == std::string::npos) break;

    std::string tok = inner.substr(p + 1, q - (p + 1));  // e.g. "0x69"
    const char* s = tok.c_str();
    if (tok.size() >= 2 && s[0] == '0' &&
        (s[1] == 'x' || s[1] == 'X')) {
      s += 2;
    }

    unsigned val = 0;
    for (; *s; ++s) {
      char c = *s;
      unsigned d = 0;
      if (c >= '0' && c <= '9')      d = c - '0';
      else if (c >= 'a' && c <= 'f') d = c - 'a' + 10;
      else if (c >= 'A' && c <= 'F') d = c - 'A' + 10;
      else                           d = 0;
      val = (val << 4) | d;
    }
    out.push_back((uint8_t)(val & 0xFF));

    p = q + 1;
  }

  return !out.empty();
}

// Parse "tt": "0x09" style field (single hex string), used in xor_remaining_phase.
static bool parse_tt_single_hex_from_object(const std::string& obj,
                                            const std::string& key,
                                            std::vector<uint8_t>& out) {
  out.clear();
  std::string pat = std::string("\"") + key + "\"";
  size_t pos = obj.find(pat);
  if (pos == std::string::npos) return false;

  pos = obj.find(':', pos + pat.size());
  if (pos == std::string::npos) return false;

  size_t i = pos + 1;
  while (i < obj.size() && is_ws(obj[i])) ++i;
  if (i >= obj.size()) return false;

  // Try direct string style: "0x.."
  if (obj[i] == '\"') {
    size_t q = obj.find('\"', i + 1);
    if (q == std::string::npos) return false;

    std::string tok = obj.substr(i + 1, q - (i + 1));
    const char* s = tok.c_str();
    if (tok.size() >= 2 && s[0] == '0' &&
        (s[1] == 'x' || s[1] == 'X')) {
      s += 2;
    }

    unsigned val = 0;
    for (; *s; ++s) {
      char c = *s;
      unsigned d = 0;
      if (c >= '0' && c <= '9')      d = c - '0';
      else if (c >= 'a' && c <= 'f') d = c - 'a' + 10;
      else if (c >= 'A' && c <= 'F') d = c - 'A' + 10;
      else                           d = 0;
      val = (val << 4) | d;
    }
    out.push_back((uint8_t)(val & 0xFF));
    return true;
  }

  // Fallback: array-style "tt": ["0x.."] (kept for backward compatibility).
  return parse_hex_tt_from_object(obj, key, out);
}

/*------------------------------------------------------------------------*/
/*               Parsing of FA/HA/XOR structures from JSON                */
/*------------------------------------------------------------------------*/

// Parse "half adders" array into vector<HalfAdder>.
static bool parse_half_adders_block(const std::string& content,
                                    std::vector<HalfAdder>& half_adders) {
  size_t lb, rb;
  if (!find_key_array_span(content, "half adders", lb, rb)) return false;

  std::string arr = content.substr(lb, rb - lb + 1);
  auto objs       = slice_objects_from_array(arr);

  for (const auto& o : objs) {
    HalfAdder h;
    if (!parse_id_from_object(o, h.id)) continue;

    std::vector<int> nodes;
    if (parse_int_array_from_object(o, "nodes", nodes)) {
      h.nodes = std::move(nodes);
    }

    half_adders.push_back(std::move(h));
  }
  return true;
}

// Parse "full adders" array into vector<FullAdder>.
static bool parse_full_adders_block(const std::string& content,
                                    std::vector<FullAdder>& full_adders) {
  size_t lb, rb;
  if (!find_key_array_span(content, "full adders", lb, rb)) return false;

  std::string arr = content.substr(lb, rb - lb + 1);
  auto objs       = slice_objects_from_array(arr);

  for (const auto& o : objs) {
    FullAdder f;
    if (!parse_id_from_object(o, f.id)) continue;

    std::vector<int> nodes;
    if (parse_int_array_from_object(o, "nodes", nodes)) {
      f.nodes = std::move(nodes);
    }

    full_adders.push_back(std::move(f));
  }
  return true;
}

/*---------------- XOR: xor_remaining_xors + xor_remaining_phase ---------*/

// Step 1: parse "xor_remaining_xors" -> id + nodes only.
// This mirrors parse_half_adders_block style.
static bool parse_xor_remaining_xors_block(const std::string& content,
                                           std::vector<XorRemain>& xor_list) {
  xor_list.clear();

  size_t lb, rb;
  if (!find_key_array_span(content, "xor_remaining_xors", lb, rb))
    return false;

  std::string arr = content.substr(lb, rb - lb + 1);
  auto objs       = slice_objects_from_array(arr);

  for (const auto& o : objs) {
    XorRemain xr;
    xr.id = -1;

    // id
    if (!parse_id_from_object(o, xr.id)) continue;

    // nodes: full XOR cone nodes
    std::vector<int> nodes;
    if (parse_int_array_from_object(o, "nodes", nodes)) {
      xr.nodes = std::move(nodes);
    }

    // tt / input_phase / output_phase will be filled later
    xor_list.push_back(std::move(xr));
  }

  return true;
}

// Step 2: parse "xor_remaining_phase" as a temporary list.
// Only parse fields that exist in struct XorRemain.
static bool parse_xor_remaining_phase_block(const std::string& content,
                                            std::vector<XorRemain>& xor_phase) {
  xor_phase.clear();

  size_t lb, rb;
  if (!find_key_array_span(content, "xor_remaining_phase", lb, rb))
    return false;

  std::string arr = content.substr(lb, rb - lb + 1);
  auto objs       = slice_objects_from_array(arr);

  for (const auto& o : objs) {
    XorRemain xr;
    xr.id = -1;

    if (!parse_id_from_object(o, xr.id)) continue;

    // tt: single hex like "0x09"
    parse_tt_single_hex_from_object(o, "tt", xr.tt);

    // input_phase: [ ... ]
    parse_int_array_from_object(o, "input_phase", xr.input_phase);

    // output_phase: [ ... ]
    parse_int_array_from_object(o, "output_phase", xr.output_phase);

    // NOTE: fields like "input_nodes", "group_last_two_nodes",
    // "check", "xor_root", "n_inputs" are intentionally ignored
    // because they are not part of struct XorRemain you provided.

    xor_phase.push_back(std::move(xr));
  }

  return true;
}

/*------------------------------------------------------------------------*/
/*                       Merge tt/phase information                       */
/*------------------------------------------------------------------------*/

// Merge "ha_tt_phase" into half_adders (by id).
static bool merge_ha_tt_phase(const std::string& content,
                              std::vector<HalfAdder>& half_adders) {
  size_t lb, rb;
  if (!find_key_array_span(content, "ha_tt_phase", lb, rb)) return false;

  std::string arr = content.substr(lb, rb - lb + 1);
  auto objs       = slice_objects_from_array(arr);

  // Build index by id for fast lookup.
  std::unordered_map<int, size_t> id2pos;
  id2pos.reserve(half_adders.size() * 2);
  for (size_t i = 0; i < half_adders.size(); ++i)
    id2pos[half_adders[i].id] = i;

  for (const auto& o : objs) {
    int id = -1;
    if (!parse_id_from_object(o, id)) continue;

    auto it = id2pos.find(id);
    if (it == id2pos.end()) {
      std::fprintf(stderr,
                   "warning: ha_tt_phase id=%d not found in 'half adders'\n",
                   id);
      continue;
    }

    auto& h = half_adders[it->second];

    std::vector<uint8_t> tt;
    if (parse_hex_tt_from_object(o, "tt", tt))
      h.tt = std::move(tt);

    std::vector<int> inph, outph;
    if (parse_int_array_from_object(o, "input_phase", inph))
      h.input_phase = std::move(inph);
    if (parse_int_array_from_object(o, "output_phase", outph))
      h.output_phase = std::move(outph);
  }

  return true;
}

// Merge "fa_tt_phase" into full_adders (by id).
static bool merge_fa_tt_phase(const std::string& content,
                              std::vector<FullAdder>& full_adders) {
  size_t lb, rb;
  if (!find_key_array_span(content, "fa_tt_phase", lb, rb)) return false;

  std::string arr = content.substr(lb, rb - lb + 1);
  auto objs       = slice_objects_from_array(arr);

  std::unordered_map<int, size_t> id2pos;
  id2pos.reserve(full_adders.size() * 2);
  for (size_t i = 0; i < full_adders.size(); ++i)
    id2pos[full_adders[i].id] = i;

  for (const auto& o : objs) {
    int id = -1;
    if (!parse_id_from_object(o, id)) continue;

    auto it = id2pos.find(id);
    if (it == id2pos.end()) {
      std::fprintf(stderr,
                   "warning: fa_tt_phase id=%d not found in 'full adders'\n",
                   id);
      continue;
    }

    auto& f = full_adders[it->second];

    std::vector<uint8_t> tt;
    if (parse_hex_tt_from_object(o, "tt", tt))
      f.tt = std::move(tt);

    std::vector<int> inph, outph;
    if (parse_int_array_from_object(o, "input_phase", inph))
      f.input_phase = std::move(inph);
    if (parse_int_array_from_object(o, "output_phase", outph))
      f.output_phase = std::move(outph);
  }

  return true;
}

// Merge xor_remaining_phase info into xor_list (from xor_remaining_xors) by id.
static bool merge_xor_tt_phase(const std::string& content,
                               std::vector<XorRemain>& xor_list) {
  std::vector<XorRemain> xor_phase;
  if (!parse_xor_remaining_phase_block(content, xor_phase))
    return false;

  // Build index by id for xor_list
  std::unordered_map<int, size_t> id2pos;
  id2pos.reserve(xor_list.size() * 2);
  for (size_t i = 0; i < xor_list.size(); ++i)
    id2pos[xor_list[i].id] = i;

  for (const auto& xp : xor_phase) {
    auto it = id2pos.find(xp.id);
    if (it == id2pos.end()) {
      std::fprintf(stderr,
                   "warning: xor_remaining_phase id=%d not found in xor_remaining_xors\n",
                   xp.id);
      continue;
    }

    XorRemain& xr = xor_list[it->second];

    // Merge only the fields that exist in struct XorRemain
    if (!xp.tt.empty())           xr.tt          = xp.tt;
    if (!xp.input_phase.empty())  xr.input_phase = xp.input_phase;
    if (!xp.output_phase.empty()) xr.output_phase= xp.output_phase;
  }

  return true;
}

/*------------------------------------------------------------------------*/
/*                            read_faha_json                              */
/*------------------------------------------------------------------------*/

// Read FA/HA/XOR metadata from JSON file (e.g., faha.json).
bool read_faha_json(const char* filename,
                    std::vector<HalfAdder>&  half_adders,
                    std::vector<FullAdder>&  full_adders,
                    std::vector<XorRemain>&  xor_list) {
  half_adders.clear();
  full_adders.clear();
  xor_list.clear();

  std::ifstream fin(filename, std::ios::in | std::ios::binary);
  if (!fin) {
    warning("Cannot open JSON file '%s'.", filename);
    return false;
  }

  std::ostringstream oss;
  oss << fin.rdbuf();
  std::string content = oss.str();
  fin.close();

  bool ok_full = parse_full_adders_block(content, full_adders);
  bool ok_half = parse_half_adders_block(content, half_adders);
  bool ok_xor  = parse_xor_remaining_xors_block(content, xor_list);

  if (!ok_half && !ok_full && !ok_xor) {
    warning("JSON file '%s' does not contain 'half adders', "
            "'full adders', or 'xor_remaining_xors'.",
            filename);
  }

  // Merge phase/tt information for HA, FA and XOR structures.
  merge_ha_tt_phase(content, half_adders);
  merge_fa_tt_phase(content, full_adders);
  merge_xor_tt_phase(content, xor_list);

  msg("Parsed JSON '%s': half adders = %zu, full adders = %zu, xor_remaining = %zu",
      filename,
      half_adders.size(),
      full_adders.size(),
      xor_list.size());
  std::fflush(stdout);
  // ----------------------------------------------------------

  return true;
}

/*------------------------------------------------------------------------*/

}  // namespace talisman
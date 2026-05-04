/*------------------------------------------------------------------------*/
/*! \file talisman.cpp
    \brief main file of our tool TalisMan

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#define VERSION "2.0"
/*------------------------------------------------------------------------*/
// / Manual of TalisMan, will be printed with command line '-h'
static const char* USAGE
  = "\n"
    "Usage: talisman <input-file> <spec-mode> [options]\n"
    "\n"
    "Positional arguments:\n"
    "  <input-file>            File containing the AIG graph\n"
    "  <spec_file>             Custom spec file (use when not using built-in "
    "specs)\n"
    "\n"
    "Specification modes (pick one):\n"
    "  -miter-spec             Use the miter specification\n"
    "  -mult-spec              Use the multiplier specification\n"
    "  -assert-spec            Use the assertion specification\n"
    "\n"  
    "General options:\n"
    "  -h, --help              Show this help message and exit\n"
    "  -v <int>                Verbosity level (0..4)\n"
    "  --custom-linalg <0|1>   Use custom routines for linear algebra (default "
    "1)\n"
    "  -primes <0|1>           Primes to use (0: structured, 1: random)\n"
    "  --reduction <0|1|2>     Reduction strategy (0:hybrid, 1:linear, "
    "2:nonlinear) (default 0)\n"
    "  --sc-depth <int>         Subcircuit depth (levels, default 2)\n"
    "  --no-trust-abc           Do not trust linear polynomial cache from ABC (perform "
    "SAT "
    "checks)\n"
    "  --uniform-sampling <0|1> Disable uniform sampling (default 1)\n"
    "\n"
    "Examples:\n"
    "  talisman input.aig -mult-spec\n"
    "  talisman input.aig spec.txt\n"
    "  talisman input.aig -mult-spec --no-trust-abc\n"
    "";
/*------------------------------------------------------------------------*/
#include <algorithm>// for strcmp
#include <cctype>   // for std::isdigit
#include <cstdlib>  // for std::atoi
#include <cstring>  // for strcmp 
#include <ctime>
#include <iomanip>

#include "abc_interface.hpp"
#include "coeff_vec.hpp"
#include "gate.hpp"
#include "guess_and_prove.hpp"
#include "linear_algebra.hpp"
#include "linpoly_cache.hpp"
#include "linpolynomial.hpp"
#include "parser.hpp"
#include "signal_statistics.hpp"
#include "specpoly.hpp"
#include "subcircuit.hpp"
#include "thread_pool.hpp"

using namespace talisman;

/*------------------------------------------------------------------------*/
// / Name of the input file
const char* input_file = nullptr;
const char* spec_file = nullptr;
// If true, trust ABC-produced linear polynomial cache and skip SAT checks
bool trust_abc = true;
/*------------------------------------------------------------------------*/
// ERROR CODES:

static int err_no_file = 100;     // no input file given
static int err_spec_sel = 101;    // mode has already been selected/not selected
static int err_ukwn_option = 102; // unknown option given
static int err_no_spec_file = 103;// no spec file provided
static int err_too_many_arguments = 104;// too many arguments given

/*------------------------------------------------------------------------*/

enum class SpecType { None, Miter, Mult, Assert, File };
SpecType spec_type = SpecType::None;

int subcircuit_depth = 2;// Default depth for subcircuit traversal
int prime = 101;

// Define linpolynomial pruning flag (default enabled)
namespace talisman {
bool custom_linalg = true;
bool uniform_sampling = true;
int reduction_strategy = 0;
int primes = 0;

using thread_pool = std::unique_ptr<BS::thread_pool<BS::none>>;

std::vector<int32_t> PRIMES;

thread_pool pool = std::make_unique<BS::thread_pool<BS::none>>(PRIMES.size());
}

/*------------------------------------------------------------------------*/
static bool
is_number(const std::string& s) {
  return !s.empty() && std::all_of(s.begin(), s.end(), ::isdigit);
}
/*------------------------------------------------------------------------*/

// Parses a non-negative integer from argv
static int
parse_non_negative_int(const char* flag, const char* value_str) {
  if(!is_number(value_str))
    die(123, "%s needs to be followed by a non-negative integer", flag);
  int value = std::stoi(value_str);
  if(value < 0)
    die(123, "%s needs to be followed by a non-negative integer", flag);
  return value;
}

/*------------------------------------------------------------------------*/
static void
parse_arguments(int argc, char** argv) {
  for(int i = 1; i < argc; i++) {
    std::string arg = argv[i];

    if(arg == "-h" || arg == "--help") {
      fputs(USAGE, stdout);
      fflush(stdout);
      exit(0);
    } else if(arg == "-v" && i + 1 < argc) {
      verbose = parse_non_negative_int("-v", argv[++i]);
    } else if(arg == "--sc-depth" && i + 1 < argc) {
      subcircuit_depth = parse_non_negative_int("--sc-depth", argv[++i]);
    } else if(arg == "--primes" && i + 1 < argc) {
      primes = parse_non_negative_int("--primes", argv[++i]);
    } else if(arg == "-miter-spec") {
      if(spec_type != SpecType::None)
        die(err_spec_sel, "only one specification is allowed (try '-h')");
      spec_type = SpecType::Miter;
    } else if(arg == "-mult-spec") {
      if(spec_type != SpecType::None)
        die(err_spec_sel, "only one specification is allowed (try '-h')");
      spec_type = SpecType::Mult;
    } else if(arg == "-assert-spec") {
      if(spec_type != SpecType::None)
        die(err_spec_sel, "only one specification is allowed (try '-h')");
      spec_type = SpecType::Assert;
    } else if(arg == "--custom-linalg") {
      if(i + 1 < argc && is_number(argv[i + 1])) {
        int v = parse_non_negative_int("--custom-linalg", argv[++i]);
        talisman::custom_linalg = (v != 0);
      }
    } else if(arg == "--reduction") {
      if(i + 1 < argc && is_number(argv[i + 1])) {
        int v = parse_non_negative_int("--reduction", argv[++i]);
        talisman::reduction_strategy = v;
      }
    } else if(arg == "--no-trust-abc") {
      trust_abc = false;
    } else if(arg == "--uniform-sampling") {
      if(i + 1 < argc && is_number(argv[i + 1])) {
        int v = parse_non_negative_int("--uniform-sampling", argv[++i]);
        talisman::uniform_sampling = (v != 0);
      }
    } else if(arg.size() > 0 && arg[0] == '-') {
      die(err_ukwn_option, "unknown option '%s' (try '-h')", argv[i]);
    } else if(!input_file) {
      input_file = argv[i];
    } else if(spec_type == SpecType::None && !spec_file) {
      spec_file = argv[i];
      spec_type = SpecType::File;
    } else {
      die(
        err_too_many_arguments, "too many arguments '%s' (try '-h')", argv[i]);
    }
  }

  if(!input_file)
    die(err_no_file, "no input file given (try '-h')");
  if(spec_type == SpecType::File && !spec_file)
    die(err_no_spec_file, "no spec file provided (try '-h')");

  mult_spec = (spec_type == SpecType::Mult);
  miter_spec = (spec_type == SpecType::Miter);
}
/*------------------------------------------------------------------------*/
static void
print_settings() {
  print_hline();
  msg("SETTINGS");
  msg("");

  switch(spec_type) {
    case SpecType::Miter:
      msg("specification: miter");
      break;
    case SpecType::Mult:
      msg("specification: unsigned multiplier");
      break;
    case SpecType::Assert:
      msg("specification: assertion");
      break;
    case SpecType::File:
      msg("spec from file %s will be used", spec_file);
      break;
    default:
      msg("specification: none");
      break;
  }

  if(spec_type == SpecType::None) {
    warning("No specification selected. ");
    warning("Circuit will be processed without verification.");
  }

  // print selected subcircuit depth
  msg("subcircuit depth: %d", subcircuit_depth);
  print_hline();
}

/*------------------------------------------------------------------------*/
int
main(int argc, char** argv) {
  msg("TalisMan Version " VERSION);
  msg("");
  msg("Copyright(C) 2025 Daniela Kaufmann, TU Wien");
  msg("                  Clemens Hofstadler, JKU Linz");
  msg("");

  if(argc < 2) {
    fputs(USAGE, stdout);
    fflush(stdout);
    return 0;
  }

  parse_arguments(argc, argv);

  SignalGuard sg;
  parse_aig(input_file);

  init_gates(mult_spec ? NN * NN / 4 : 0);
  print_hline();

  linpoly_cache_init();

  if(primes == 0 and NN <= 128)
    talisman::PRIMES
      = { 65519, 65521, 65537, 65539, 65543, 65551, 65557, 65563 };
  else if(primes == 1 and NN <= 128)
    talisman::PRIMES
      = { 91373, 79601, 91583, 127837, 46639, 57719, 111497, 43103 };
  if(primes == 0 and NN > 128)
    talisman::PRIMES
      = { 65497, 65519, 65521, 65537, 65539, 65543, 65551, 65557,
          65563, 65579, 65581, 65587, 65599, 65609, 65617, 65629 };
  else if(primes == 1 and NN > 128)
    talisman::PRIMES
      = { 76907,  120539, 91513, 45841,  46399,  96601, 72253,  42743,
          122147, 79231,  73961, 119027, 128659, 49727, 119297, 78989 };

  num_primes = talisman::PRIMES.size();

  CoeffVec::set_moduli();
  init_two_pow_NN_imgs();

  msg("Primes: %s", talisman::primes == 0 ? "structured" : "random");
  msg("Reduction strategy: %s",
      talisman::reduction_strategy == 0   ? "hybrid"
      : talisman::reduction_strategy == 1 ? "linear"
                                          : "nonlinear");
  msg("Using ABC internal detection for adders");
  if(talisman::reduction_strategy < 2) {
    msg("Custom linear algebra: %s",
        talisman::custom_linalg ? "enabled" : "disabled");
    msg("Using %s sampling", talisman::uniform_sampling ? "uniform" : "random");
  }

  // Detection of FA/HA/XOR and FSA nodes using ABC static library
  auto abc_start = std::chrono::high_resolution_clock::now();
  try {
    // Initialize ABC
    abc_interface::abc_init();

    // Load AIG file
    Gia_Man_t* pGia = abc_interface::load_aig(input_file);

    if(verbose >= 2) {
      abc_interface::print_gia_stats(pGia);
    }

    // Detect HA/FA/XOR first to populate cache.vAdds for get_kcover_nodes
    std::vector<HalfAdder> ha_list;
    std::vector<FullAdder> fa_list;
    // NOTE: xor_list is kept for API compatibility but not used since affine
    // cuts cover it
    std::vector<XorRemain> xor_list;
    if(abc_interface::detect_ha_fa_xor(pGia, ha_list, fa_list, xor_list)) {
      if(verbose >= 2) {
        // XOR-remaining count is no longer relevant since we use affine cuts
        // instead
        msg("ABC detected FA count: %zu, HA count: %zu",
            fa_list.size(),
            ha_list.size());
      }

      if(!fa_list.empty()) {
        mark_fahaxor(fa_list, 1);
      }
      if(!ha_list.empty()) {
        mark_fahaxor(ha_list, 0);
      }
      // NOTE: xor_remain processing is commented out because affine cuts
      // already cover the linear relations from xor_remaining gates.
      // if(!xor_list.empty()) {
      //   mark_fahaxor(xor_list);
      // }

      if(trust_abc) {
        trust_linpoly_cache();
      } else {
        check_linpoly_cache(subcircuit_depth, prime);
      }
    } else {
      warning("ABC HA/FA/XOR detection failed. Proceeding without it.");
    }

    // Get K-Cover nodes (reverse_fsa_nodes and boundary_nodes)
    // Now uses cached vAdds from detect_ha_fa_xor via get_detection_cache()
    std::vector<int> reverse_fsa_nodes;
    std::vector<int> boundary_nodes;
    int k = 0;
    if(abc_interface::get_kcover_nodes(
         pGia, k, reverse_fsa_nodes, boundary_nodes)) {
      if(verbose >= 2) {
        msg("ABC detected %zu reverse fsa nodes and %zu boundary nodes.",
            reverse_fsa_nodes.size(),
            boundary_nodes.size());
      }
      mark_fsa(reverse_fsa_nodes, boundary_nodes);
    } else {
      warning("ABC K-Cover detection failed. Proceeding without it.");
    }

    // Cleanup ABC
    abc_interface::clear_detection_cache();// Clear cached vAdds/vXors
    abc_interface::free_gia(pGia);
    abc_interface::abc_cleanup();

  } catch(const std::exception& e) {
    warning("ABC Error: %s. Proceeding without adder detection.", e.what());
  }
  // testing cache linear polynomial
  if(verbose >= 2) {
    print_linpoly_cache();
  }
  auto abc_end = std::chrono::high_resolution_clock::now();
  abc_time = std::chrono::duration<double>(abc_end - abc_start).count();

  LinPolynomial spec_poly = gen_mult_spec_poly();

  bool result = guess_and_prove(spec_poly, subcircuit_depth);

  if(!result) {
    msg("Guess and Prove failed.");
  } else {
    msg("Guess and Prove succeeded.");
  }

  // linpoly_cache_clear();
  // reset_aig_parsing();
  // reset_all_signal_handlers();

  print_statistics();
  return result ? 0 : 1;
}

/*------------------------------------------------------------------------*/
/*! \file main.cpp
    \brief command line front end of our tool TalisMan

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
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
    "  --trust-abc             Trust linear polynomial cache from ABC (skip "
    "SAT "
    "checks)\n"
    "  --uniform-sampling <0|1> Disable uniform sampling (default 1)\n"
    "\n"
    "Examples:\n"
    "  talisman input.aig -mult-spec\n"
    "  talisman input.aig spec.txt\n"
    "  talisman input.aig -mult-spec --trust-abc\n"
    "";
/*------------------------------------------------------------------------*/
#include <algorithm>// for strcmp
#include <cctype>   // for std::isdigit
#include <cstdio>
#include <string>

#include "talisman/talisman.hpp"

#include "signal_statistics.hpp"

using namespace talisman;

/*------------------------------------------------------------------------*/
// ERROR CODES:

static int err_no_file = 100;     // no input file given
static int err_spec_sel = 101;    // mode has already been selected/not selected
static int err_ukwn_option = 102; // unknown option given
static int err_no_spec_file = 103;// no spec file provided
static int err_too_many_arguments = 104;// too many arguments given

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
// Returns false if the program should stop right away (e.g. after '-h')
static bool
parse_arguments(int argc, char** argv, Options& opts) {
  for(int i = 1; i < argc; i++) {
    std::string arg = argv[i];

    if(arg == "-h" || arg == "--help") {
      fputs(USAGE, stdout);
      fflush(stdout);
      return false;
    } else if(arg == "-v" && i + 1 < argc) {
      opts.verbose = parse_non_negative_int("-v", argv[++i]);
    } else if(arg == "--sc-depth" && i + 1 < argc) {
      opts.sc_depth = parse_non_negative_int("--sc-depth", argv[++i]);
    } else if(arg == "--primes" && i + 1 < argc) {
      opts.primes = parse_non_negative_int("--primes", argv[++i]);
    } else if(arg == "-miter-spec") {
      if(opts.spec != SpecMode::None)
        die(err_spec_sel, "only one specification is allowed (try '-h')");
      opts.spec = SpecMode::Miter;
    } else if(arg == "-mult-spec") {
      if(opts.spec != SpecMode::None)
        die(err_spec_sel, "only one specification is allowed (try '-h')");
      opts.spec = SpecMode::Mult;
    } else if(arg == "-assert-spec") {
      if(opts.spec != SpecMode::None)
        die(err_spec_sel, "only one specification is allowed (try '-h')");
      opts.spec = SpecMode::Assert;
    } else if(arg == "--custom-linalg") {
      if(i + 1 < argc && is_number(argv[i + 1])) {
        int v = parse_non_negative_int("--custom-linalg", argv[++i]);
        opts.custom_linalg = (v != 0);
      }
    } else if(arg == "--reduction") {
      if(i + 1 < argc && is_number(argv[i + 1])) {
        int v = parse_non_negative_int("--reduction", argv[++i]);
        opts.reduction = v;
      }
    } else if(arg == "--trust-abc") {
      opts.trust_abc = true;
    } else if(arg == "--uniform-sampling") {
      if(i + 1 < argc && is_number(argv[i + 1])) {
        int v = parse_non_negative_int("--uniform-sampling", argv[++i]);
        opts.uniform_sampling = (v != 0);
      }
    } else if(arg.size() > 0 && arg[0] == '-') {
      die(err_ukwn_option, "unknown option '%s' (try '-h')", argv[i]);
    } else if(opts.input_file.empty()) {
      opts.input_file = argv[i];
    } else if(opts.spec == SpecMode::None && opts.spec_file.empty()) {
      opts.spec_file = argv[i];
      opts.spec = SpecMode::File;
    } else {
      die(
        err_too_many_arguments, "too many arguments '%s' (try '-h')", argv[i]);
    }
  }

  if(opts.input_file.empty())
    die(err_no_file, "no input file given (try '-h')");
  if(opts.spec == SpecMode::File && opts.spec_file.empty())
    die(err_no_spec_file, "no spec file provided (try '-h')");

  return true;
}
/*------------------------------------------------------------------------*/
[[maybe_unused]] static void
print_settings(const Options& opts) {
  print_hline();
  msg("SETTINGS");
  msg("");

  switch(opts.spec) {
    case SpecMode::Miter:
      msg("specification: miter");
      break;
    case SpecMode::Mult:
      msg("specification: unsigned multiplier");
      break;
    case SpecMode::Assert:
      msg("specification: assertion");
      break;
    case SpecMode::File:
      msg("spec from file %s will be used", opts.spec_file.c_str());
      break;
    default:
      msg("specification: none");
      break;
  }

  if(opts.spec == SpecMode::None) {
    warning("No specification selected. ");
    warning("Circuit will be processed without verification.");
  }

  // print selected subcircuit depth
  msg("subcircuit depth: %d", opts.sc_depth);
  print_hline();
}

/*------------------------------------------------------------------------*/
int
main(int argc, char** argv) {
  msg("TalisMan Version %s", version());
  msg("");
  msg("Copyright(C) 2025 Daniela Kaufmann, TU Wien");
  msg("                  Clemens Hofstadler, JKU Linz");
  msg("");

  if(argc < 2) {
    fputs(USAGE, stdout);
    fflush(stdout);
    return 0;
  }

  Options opts;
  try {
    if(!parse_arguments(argc, argv, opts))
      return 0;
  } catch(const Error& e) {
    return e.code;
  }

  Hooks hooks;
  hooks.install_signal_handlers = true;

  Result res = run(opts, hooks);
  switch(res.status) {
    case Status::Verified:
      return 0;
    case Status::Error:
      return res.error_code;
    default:
      return 1;
  }
}

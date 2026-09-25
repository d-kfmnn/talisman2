/*------------------------------------------------------------------------*/
/*! \file signal_statistics.cpp
    \brief used to handle signals, messages and statistics

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/

#include "signal_statistics.hpp"

#include <sys/resource.h>
#include <unistd.h>

#include <atomic>
#include <chrono>
#include <csignal>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string_view>

#include "guess_and_prove.hpp"

// ANSI color codes
static const char* RED = "\033[31m";
[[maybe_unused]] static const char* GREEN = "\033[32m";
static const char* YELLOW = "\033[33m";
static const char* RESET = "\033[0m";

/*------------------------------------------------------------------------*/
namespace talisman {

// Global variable
int verbose = 1;
bool booth = 0;

double abc_time = 0;
double extract_sc_time = 0;
double sample_time = 0;
std::atomic<double> guess_time{0.0};
double guess_noparallel_time = 0;
std::atomic<double> prove_time{0.0};
std::atomic<double> repair_time{0.0};
double nonlin_time = 0;
double simplify_time = 0;
double lin_time = 0;
size_t abc_polies = 0;
size_t lin_steps = 0;
size_t nonlin_steps = 0;
size_t num_repair_loop = 0;
size_t max_num_repair_loop = 0;
size_t num_guesses = 0;
size_t num_proves = 0;
size_t num_correct_guesses = 0;
int num_primes = 0;

constexpr std::string_view prefix = "[talisman] ";

enum class Signal { Unknown, INT, SEGV, ABRT, TERM };

const char*
signal_name(int sig) {
  switch(sig) {
    case SIGINT:
      return "SIGINT";
    case SIGSEGV:
      return "SIGSEGV";
    case SIGABRT:
      return "SIGABRT";
    case SIGTERM:
      return "SIGTERM";
    default:
      return "SIGUNKNOWN";
  }
}

using __sighandler_t = void (*)(int);

// Static member definitions:
__sighandler_t SignalGuard::original_SIGINT = nullptr;
__sighandler_t SignalGuard::original_SIGSEGV = nullptr;
__sighandler_t SignalGuard::original_SIGABRT = nullptr;
__sighandler_t SignalGuard::original_SIGTERM = nullptr;

// Constructor:
SignalGuard::SignalGuard() {
  install();
}

// Destructor:
SignalGuard::~SignalGuard() {
  reset_all_signal_handlers();
}

// Install handlers and save originals:
void
SignalGuard::install() {
  original_SIGINT = std::signal(SIGINT, catch_signal);
  original_SIGSEGV = std::signal(SIGSEGV, catch_signal);
  original_SIGABRT = std::signal(SIGABRT, catch_signal);
  original_SIGTERM = std::signal(SIGTERM, catch_signal);
}

// Signal handler function:
void
SignalGuard::catch_signal(int sig) {
  std::cerr << "c\nc caught signal '" << signal_name(sig) << "' (" << sig
            << ")\nc\n";
  std::cerr << "c\nc raising signal '" << signal_name(sig) << "' (" << sig
            << ") again\n";
  reset_all_signal_handlers();
  std::fflush(stdout);
  std::fflush(stderr);
  std::raise(sig);// re-raise to default handler
}

// Reset all signals to original handlers:
void
reset_all_signal_handlers() {
  std::signal(SIGINT, SignalGuard::original_SIGINT);
  std::signal(SIGSEGV, SignalGuard::original_SIGSEGV);
  std::signal(SIGABRT, SignalGuard::original_SIGABRT);
  std::signal(SIGTERM, SignalGuard::original_SIGTERM);
}

// ------------------------------------------------------------------
// Messages and Error Handling
// ------------------------------------------------------------------

void
msg(const char* fmt, ...) {

  std::cout << prefix;
  va_list args;
  va_start(args, fmt);
  vprintf(fmt, args);
  va_end(args);
  std::cout << '\n';
  std::fflush(stdout);
}
// ------------------------------------------------------------------
void
warning(const char* fmt, ...) {
  std::cout << YELLOW;// Set text color to green
  std::cout << prefix;
  std::cout << "WARNING: ";
  va_list args;
  va_start(args, fmt);
  vprintf(fmt, args);
  va_end(args);
  std::cout << '\n';
  std::cout << RESET;// Reset color after printing
  std::fflush(stdout);
}
// ------------------------------------------------------------------

void
msg_nl(const char* fmt, ...) {
  if(verbose < 1)
    return;// Only print when verbose >= 1

  std::cout << prefix;
  va_list args;
  va_start(args, fmt);
  vprintf(fmt, args);
  va_end(args);
  std::fflush(stdout);
}
// ------------------------------------------------------------------
void
print_hline() {
  if(verbose < 1)
    return;// Only print when verbose >= 1

  std::cout << prefix
            << "-------------------------------------------------------\n";
  std::fflush(stdout);
}

// ------------------------------------------------------------------
[[noreturn]]
void
die(int error_code, const char* fmt, ...) {
  va_list args;
  va_start(args, fmt);
  va_list args_len;
  va_copy(args_len, args);
  std::string message(std::vsnprintf(nullptr, 0, fmt, args_len), '\0');
  va_end(args_len);
  std::vsnprintf(message.data(), message.size() + 1, fmt, args);
  va_end(args);

  std::cout << RED;// Set text color to green
  const char* prefix = "[ERROR] ";
  std::fflush(stdout);
  std::cerr << "*** " << prefix << "error code " << error_code << '\n';
  std::cerr << "*** " << prefix << message << '\n';
  std::cout << RESET << '\n';// Reset color after printing
  std::fflush(stderr);
  throw Error{ error_code, message };
}

// ------------------------------------------------------------------
// Timing & Memory Usage
// ------------------------------------------------------------------

auto start_time = std::chrono::steady_clock::now();
// ------------------------------------------------------------------

double
process_time() {
  struct rusage u;
  if(getrusage(RUSAGE_SELF, &u))
    return 0.0;
  double user = u.ru_utime.tv_sec + u.ru_utime.tv_usec / 1e6;
  double sys = u.ru_stime.tv_sec + u.ru_stime.tv_usec / 1e6;
  return user + sys;
}
// ------------------------------------------------------------------
static size_t
maximum_resident_set_size() {
  struct rusage u;
  if(getrusage(RUSAGE_SELF, &u))
    return 0;
  return static_cast<size_t>(u.ru_maxrss) << 10;// in bytes
}
// ------------------------------------------------------------------
static double
percent(double a, double b) {
  return b ? (100.0 * a / b) : 0.0;
}

static double
average(double a, double b) {
  return b ? (a / b) : 0.0;
}
// ------------------------------------------------------------------
void
print_statistics() {

  // nonlin includes simplify time, so we need to subtract it to get the pure nonlin time
  nonlin_time = nonlin_time - simplify_time;
  
  size_t max_rss = maximum_resident_set_size();
  auto now = std::chrono::steady_clock::now();
  std::chrono::duration<double> wall = now - start_time;
  double wall_time = wall.count();
  double cpu_time = process_time();

  std::cout << prefix << "=======================================================\n";
  std::cout << prefix << "STATISTICS\n";
  std::cout << prefix << "=======================================================\n";
  std::cout << prefix << "\n";

  // Resource Usage Section
  std::cout << prefix << "--- Resource Usage ---\n";
  std::cout << prefix
            << "  Maximum Resident Set Size: " << (max_rss / (1024.0 * 1024.0))
            << " MB\n";
  std::cout << prefix << "  Wall (elapsed) Time:       " << wall_time << " seconds\n";
  std::cout << prefix << "  CPU (user+sys) Time:       " << cpu_time << " seconds\n";
  if(wall_time > 0.0) {
    std::cout << prefix << "  CPU Utilization:           " << (100.0 * cpu_time / wall_time)
              << " %\n";
  }
  std::cout << prefix << "\n";

  // Time Breakdown Section
  std::cout << prefix << "--- Time Breakdown ---\n";
  std::cout << prefix << "  ABC Time:                   " << abc_time << " seconds\n";
  std::cout << prefix << "  Extract Subcircuit Time:    " << extract_sc_time << " seconds\n";
  std::cout << prefix << "  Sample Time:                " << sample_time << " seconds\n";
  std::cout << prefix << "  Guess Time (avg/thread):    " << guess_time.load() / num_primes + guess_noparallel_time << " seconds\n";
  std::cout << prefix << "  Prove Time (avg/thread):    " << prove_time.load() / num_primes << " seconds\n";
  std::cout << prefix << "  Repair Time (avg/thread):   " << repair_time.load() / num_primes << " seconds\n";
  std::cout << prefix << "  Linear Rew. Time:           " << lin_time << " seconds\n";
  std::cout << prefix << "  Simplify Time:              " << simplify_time << " seconds\n";
  std::cout << prefix << "  Nonlin. Rew. Time:          " << nonlin_time << " seconds\n";
  std::cout << prefix << "  ----------------------------------------\n";
  std::cout << prefix << "  Sum of Times:               " 
            << (abc_time + extract_sc_time + sample_time + (guess_time.load() + prove_time.load() + repair_time.load())/num_primes + simplify_time + guess_noparallel_time + lin_time + nonlin_time) 
            << " seconds\n";
  std::cout << prefix << "\n";

  // Rewriting Steps Section
  std::cout << prefix << "--- Rewriting Steps ---\n";
  std::cout << prefix << "  Linear Rewriting Steps:    " << lin_steps << "\n";
  std::cout << prefix << "  Nonlinear Rewriting Steps: " << nonlin_steps << "\n";
  std::cout << prefix << "  ABC Polies:                " << abc_polies << "\n";
  std::cout << prefix << "\n";

  // Algorithm Statistics Section
  std::cout << prefix << "--- Algorithm Statistics ---\n";
  std::cout << prefix << "  Nr. Guesses:               " << num_guesses 
            << " (avg per thread: " << (num_primes ? static_cast<double>(num_guesses) / num_primes : 0.0) << ")\n";
  std::cout << prefix << "  Nr. Correct Guesses:       " << num_correct_guesses 
            << " (avg per thread: " << (num_primes ? static_cast<double>(num_correct_guesses) / num_primes : 0.0) << ")\n";
  std::cout << prefix << "  Nr. Proofs:                " << num_proves 
            << " (avg per thread: " << (num_primes ? static_cast<double>(num_proves) / num_primes : 0.0) << ")\n";
  std::cout << prefix << "  Nr. Repair Loops:          " << num_repair_loop 
            << " (avg per thread: " << (num_primes ? static_cast<double>(num_repair_loop) / num_primes : 0.0) << ")\n";
  std::cout << prefix << "  Max Repair Loops/Circuit:  " << max_num_repair_loop << "\n";
  std::cout << prefix << "\n";
  
  // Print individual thread times
  print_individual_thread_times();
  
  std::cout << prefix << "=======================================================\n";

  std::fflush(stdout);
}
// ------------------------------------------------------------------
void
reset_statistics() {
  abc_time = 0;
  extract_sc_time = 0;
  sample_time = 0;
  guess_time = 0.0;
  guess_noparallel_time = 0;
  prove_time = 0.0;
  repair_time = 0.0;
  nonlin_time = 0;
  simplify_time = 0;
  lin_time = 0;
  abc_polies = 0;
  lin_steps = 0;
  nonlin_steps = 0;
  num_repair_loop = 0;
  max_num_repair_loop = 0;
  num_guesses = 0;
  num_proves = 0;
  num_correct_guesses = 0;
  num_primes = 0;
  reset_individual_thread_times();
  start_time = std::chrono::steady_clock::now();
}

}// namespace talisman

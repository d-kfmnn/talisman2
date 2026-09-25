/*------------------------------------------------------------------------*/
/*! \file signal_statistics.h
    \brief used to handle signals, messages and statistics

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#ifndef TALISMAN_SRC_SIGNAL_STATISTICS_H_
#define TALISMAN_SRC_SIGNAL_STATISTICS_H_
/*------------------------------------------------------------------------*/
#include <atomic>
#include <cstdint>
#include <signal.h>
#include <stdarg.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <vector>

#include <iostream>
#include <string>
/*------------------------------------------------------------------------*/
#define LOGDBG(fmt, ...) \
  if(verbose >= 4)       \
  msg(fmt, __VA_ARGS__)

namespace talisman {

// / Level of output verbosity, ranges from 0 to 4
extern int verbose;
extern bool booth;

extern double abc_time;
extern double extract_sc_time;
extern double sample_time;
extern std::atomic<double> guess_time;
extern double guess_noparallel_time;
extern std::atomic<double> prove_time;
extern std::atomic<double> repair_time;
extern double simplify_time;
extern double nonlin_time;
extern double lin_time;
extern size_t abc_polies;
extern size_t lin_steps;
extern size_t nonlin_steps;
extern size_t num_repair_loop;
extern size_t max_num_repair_loop;
extern size_t num_guesses;
extern size_t num_proves;
extern size_t num_correct_guesses;
extern int num_primes;

class SignalGuard {
  public:
  SignalGuard();
  ~SignalGuard();

  void install();

  static void catch_signal(int sig);

  friend void reset_all_signal_handlers();

  private:
  using __sighandler_t = void (*)(int);

  static __sighandler_t original_SIGINT;
  static __sighandler_t original_SIGSEGV;
  static __sighandler_t original_SIGABRT;
  static __sighandler_t original_SIGTERM;
};

/**
    Resets all signal handlers
*/
void
reset_all_signal_handlers();

/*------------------------------------------------------------------------*/

// Thrown by die(). Deliberately not derived from std::exception so that
// generic `catch(const std::exception&)` handlers do not swallow it.
struct Error {
  int code;
  std::string message;
};

/**
    Prints an error message to stderr and throws Error, which aborts the
    current run

    @param char* fmt message
    @param int error_code
*/
[[noreturn]]
void
die(int error_code, const char* fmt, ...);

/**
    Prints a message to stdout

    @param char* fmt message
*/
void
msg(const char* fmt, ...);
void
msg_nl(const char* fmt, ...);

void
warning(const char* fmt, ...);

void
print_hline();
/*------------------------------------------------------------------------*/
// Timing & Memory Usage
// ------------------------------------------------------------------

/**
    Determines the used process time
*/
double
process_time();

/**
    Print statistics of maximum memory and used process time depending on
    selected modus
*/
void
print_statistics();

/**
    Resets all statistics counters and the start time
*/
void
reset_statistics();

}// namespace talisman

#endif// TALISMAN_SRC_SIGNAL_STATISTICS_H_

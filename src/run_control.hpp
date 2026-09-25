/*------------------------------------------------------------------------*/
/*! \file run_control.hpp
    \brief termination of a run and error propagation out of worker tasks

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#ifndef TALISMAN_SRC_RUN_CONTROL_H_
#define TALISMAN_SRC_RUN_CONTROL_H_
/*------------------------------------------------------------------------*/
#include <atomic>
#include <exception>
#include <functional>
#include <mutex>
#include <utility>

#include "signal_statistics.hpp"

namespace talisman {

// Thrown by check_terminate(). Like Error, deliberately not derived from
// std::exception so that generic handlers do not swallow it.
struct Interrupted {};

/**
    Installs the callback polled by terminate_requested(), empty disables it
*/
void
set_terminate_callback(std::function<bool()> cb);

bool
has_terminate_callback();

/**
    True once the terminate callback returned true (sticky until the next
    set_terminate_callback). Thread-safe.
*/
bool
terminate_requested();

/**
    Throws Interrupted if termination was requested
*/
void
check_terminate();

/*------------------------------------------------------------------------*/
// BS::thread_pool silently drops exceptions of detached tasks. Wrap tasks
// with TaskErrors::wrap so that die() (Error) or an interruption on a worker
// thread is rethrown on the thread waiting for the pool. Once a task failed,
// the remaining tasks of the batch are skipped.
class TaskErrors {
  public:
  template<typename F>
  auto wrap(F task) {
    return [this, task = std::move(task)]() mutable {
      if(failed.load(std::memory_order_relaxed))
        return;
      try {
        task();
      } catch(const Error&) {
        record(std::current_exception());
      } catch(const Interrupted&) {
        record(std::current_exception());
      }
    };
  }

  // Rethrows the first recorded exception. Call after the pool was waited for.
  void rethrow_if_any();

  private:
  void record(std::exception_ptr e);

  std::atomic<bool> failed{ false };
  std::mutex mutex;
  std::exception_ptr first;
};

}// namespace talisman

#endif// TALISMAN_SRC_RUN_CONTROL_H_

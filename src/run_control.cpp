/*------------------------------------------------------------------------*/
/*! \file run_control.cpp
    \brief termination of a run and error propagation out of worker tasks

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#include "run_control.hpp"

namespace talisman {

static std::function<bool()> terminate_cb;
static std::atomic<bool> terminated{ false };

/*------------------------------------------------------------------------*/
void
set_terminate_callback(std::function<bool()> cb) {
  terminate_cb = std::move(cb);
  terminated = false;
}
/*------------------------------------------------------------------------*/
bool
has_terminate_callback() {
  return static_cast<bool>(terminate_cb);
}
/*------------------------------------------------------------------------*/
bool
terminate_requested() {
  if(terminated.load(std::memory_order_relaxed))
    return true;
  if(terminate_cb && terminate_cb()) {
    terminated.store(true, std::memory_order_relaxed);
    return true;
  }
  return false;
}
/*------------------------------------------------------------------------*/
void
check_terminate() {
  if(terminate_requested())
    throw Interrupted();
}
/*------------------------------------------------------------------------*/
void
TaskErrors::record(std::exception_ptr e) {
  std::lock_guard<std::mutex> lock(mutex);
  failed.store(true, std::memory_order_relaxed);
  if(!first)
    first = std::move(e);
}
/*------------------------------------------------------------------------*/
void
TaskErrors::rethrow_if_any() {
  std::exception_ptr e;
  {
    std::lock_guard<std::mutex> lock(mutex);
    std::swap(e, first);
    failed.store(false, std::memory_order_relaxed);
  }
  if(e)
    std::rethrow_exception(e);
}

}// namespace talisman

#pragma once
// runner_timeout.hpp — POSIX setitimer wrapper for per-entry timeouts.
//
// Spec: .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Runner_Robustness.md (F7.5.A3).
//
// Why setitimer + SIGALRM rather than threads:
//   - portable POSIX, no thread-cancel undefined behaviour
//   - signal handler can set a sig_atomic_t flag AND call
//     std::atomic<bool>::store (lock-free, async-signal-safe on all
//     platforms we target) on cas::symbolic::CASContext::interrupt()
//   - cleaner termination than pthread_kill
//
// One global timer (this is a single-threaded CLI runner). We expose:
//   set_timeout_target(ctx*)  — register CAS context to interrupt
//   start_entry_timer(sec)    — arm timer
//   stop_entry_timer()        — disarm timer
//   timed_out()               — was SIGALRM raised since last arm?

#include "cas/symbolic.hpp"

#include <atomic>
#include <csignal>
#include <sys/time.h>

namespace cas::golden {

// Lifted global state. Set by main.cpp once at startup.
inline std::atomic<cas::symbolic::CASContext*>& timeout_target() {
    static std::atomic<cas::symbolic::CASContext*> t{nullptr};
    return t;
}

inline volatile std::sig_atomic_t& timed_out_flag() {
    static volatile std::sig_atomic_t f = 0;
    return f;
}

// SIGALRM handler. Async-signal-safe constraints:
//   - sig_atomic_t volatile write: OK
//   - std::atomic<T*>::load(memory_order_relaxed): OK (lock-free)
//   - CASContext::interrupt() sets std::atomic_bool: OK (lock-free)
//   No I/O, no malloc, no locks.
inline void alarm_handler(int) {
    timed_out_flag() = 1;
    if (auto* ctx = timeout_target().load(std::memory_order_relaxed)) {
        ctx->interrupt();
    }
}

inline void install_alarm_handler() {
    struct sigaction sa{};
    sa.sa_handler = &alarm_handler;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = 0;  // restartable not needed: we want syscalls to interrupt
    sigaction(SIGALRM, &sa, nullptr);
}

// Arm ITIMER_REAL to fire once after `seconds`. Zero means disarm.
inline void set_real_timer_seconds(unsigned int seconds) {
    struct itimerval timer{};
    timer.it_value.tv_sec = static_cast<time_t>(seconds);
    timer.it_value.tv_usec = 0;
    // it_interval = 0 → one-shot
    setitimer(ITIMER_REAL, &timer, nullptr);
}

inline void start_entry_timer(unsigned int seconds) {
    timed_out_flag() = 0;
    if (auto* ctx = timeout_target().load(std::memory_order_relaxed)) {
        ctx->clear_interrupt();
    }
    set_real_timer_seconds(seconds);
}

inline void stop_entry_timer() {
    set_real_timer_seconds(0);
}

inline bool entry_timed_out() {
    return timed_out_flag() != 0;
}

}  // namespace cas::golden

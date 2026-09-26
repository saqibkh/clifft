#pragma once

// Only the separately built research executor includes this instrumentation.
#include <cassert>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <span>

namespace wan_profile {
using Clock = std::chrono::steady_clock;
struct Counter {
    uint64_t nanoseconds = 0;
    uint64_t visits = 0;
};
extern std::span<Counter> counters;
inline void record(size_t action, Clock::time_point start) noexcept {
    const auto end = Clock::now();
    assert(action < counters.size());
    counters[action].nanoseconds +=
        std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    ++counters[action].visits;
}
}  // namespace wan_profile

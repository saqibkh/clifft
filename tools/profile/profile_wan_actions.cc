// Source-attributed timing for one precompiled scalar CPU plan.
#include "clifft/circuit/parser.h"
#include "clifft/frontend/frontend.h"
#include "clifft/optimizer/active_width_schedule_pass.h"
#include "clifft/optimizer/pass_factory.h"
#include "clifft/sampling/executable_plan.h"
#include "clifft/sampling/executor.h"
#include "clifft/sampling/planner.h"

#include "wan_action_timer.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace wan_profile {
std::span<Counter> counters;
}

int main(int argc, char** argv) {
    if (argc != 8) {
        std::cerr << "Usage: profiler circuit schedule postselect shots seed output_prefix refs\n";
        return 2;
    }
    const std::string schedule = argv[2];
    const bool postselect = std::stoi(argv[3]) != 0;
    const size_t shots = std::stoull(argv[4]);
    const uint64_t seed = std::stoull(argv[5]);
    const std::string output = argv[6];
    std::ifstream input(argv[1]);
    if (!input)
        throw std::runtime_error("Cannot open circuit");
    const std::string text{std::istreambuf_iterator<char>(input), {}};
    auto hir = clifft::trace(clifft::parse(text));
    auto passes = clifft::default_hir_pass_manager();
    passes.run(hir);
    if (schedule != "default") {
        clifft::ActiveWidthScheduleOptions options;
        if (schedule == "unlimited" || schedule == "nosink") {
            options.search_budget = std::nullopt;
        } else if (schedule != "scheduled") {
            throw std::invalid_argument("Unknown schedule");
        }
        options.sink_neutral_rotations = schedule != "nosink";
        clifft::ActiveWidthSchedulePass scheduler(options);
        scheduler.run(hir);
    }
    std::ifstream references(argv[7]);
    std::string detector_bits, observable_bits;
    if (!std::getline(references, detector_bits) || !std::getline(references, observable_bits)) {
        throw std::runtime_error("Missing reference parities");
    }
    const auto bits = [](const std::string& value) {
        std::vector<uint8_t> result;
        for (char bit : value) {
            if (bit != '0' && bit != '1')
                throw std::invalid_argument("Invalid parity");
            result.push_back(static_cast<uint8_t>(bit - '0'));
        }
        return result;
    };
    const auto detectors = bits(detector_bits);
    const auto observables = bits(observable_bits);
    const std::vector<uint8_t> mask(hir.num_detectors, postselect ? 1 : 0);
    const auto plan = clifft::sampling::plan_sampling(hir, {.postselection_mask = mask,
                                                            .expected_detectors = detectors,
                                                            .expected_observables = observables,
                                                            .retain_source_map = true});
    if (plan.peak_active_width > 24)
        throw std::runtime_error("Dense worker exceeds bound");
    clifft::sampling::ExecutablePlan program(plan);
    std::ofstream(output + ".plan.txt") << program.inspect();
    std::vector<wan_profile::Counter> current(program.num_actions());
    std::vector<wan_profile::Counter> accepted(program.num_actions());
    std::vector<wan_profile::Counter> rejected(program.num_actions());
    wan_profile::counters = current;
    clifft::sampling::Executor executor(program, seed);
    uint64_t record_hash = 14695981039346656037ULL;
    uint64_t accepted_count = 0;
    uint64_t error_count = 0;
    uint64_t total_ns = 0;
    uint64_t accepted_ns = 0;
    uint64_t rejected_ns = 0;
    std::ofstream shot_output(output + ".shots.tsv");
    shot_output << "shot\taccepted\tnanoseconds\n";
    for (size_t shot = 0; shot < shots; ++shot) {
        std::fill(current.begin(), current.end(), wan_profile::Counter{});
        const auto start = wan_profile::Clock::now();
        executor.run_shot();
        const auto elapsed =
            std::chrono::duration_cast<std::chrono::nanoseconds>(wan_profile::Clock::now() - start)
                .count();
        total_ns += elapsed;
        const bool keep = !executor.discarded() &&
                          std::ranges::all_of(executor.detectors(), [](uint8_t d) { return !d; });
        accepted_count += keep;
        if (keep) {
            accepted_ns += elapsed;
            error_count += std::ranges::any_of(executor.observables(), [](uint8_t o) { return o; });
        } else {
            rejected_ns += elapsed;
        }
        auto& target = keep ? accepted : rejected;
        for (size_t action = 0; action < current.size(); ++action) {
            target[action].nanoseconds += current[action].nanoseconds;
            target[action].visits += current[action].visits;
        }
        for (uint8_t record : executor.visible_records()) {
            record_hash = (record_hash ^ record) * 1099511628211ULL;
        }
        shot_output << shot << '\t' << keep << '\t' << elapsed << '\n';
    }
    std::ofstream actions(output + ".actions.tsv");
    actions << "action\taccepted_ns\trejected_ns\taccepted_visits\trejected_visits\tactive_width"
               "\tsource_lines\tinstruction\n";
    for (size_t action = 0; action < current.size(); ++action) {
        actions << action << '\t' << accepted[action].nanoseconds << '\t'
                << rejected[action].nanoseconds << '\t' << accepted[action].visits << '\t'
                << rejected[action].visits << '\t';
        std::set<uint32_t> source_lines;
        const auto range = program.action_plan_range(action);
        if (!range || !plan.source_map)
            throw std::runtime_error("Missing source provenance");
        actions << plan.actions.at(range->begin).active_before << '\t';
        for (size_t i = range->begin; i < range->end; ++i) {
            const auto lines = plan.source_map->lines_for(i);
            source_lines.insert(lines.begin(), lines.end());
        }
        bool first = true;
        for (uint32_t line : source_lines) {
            if (!first)
                actions << ',';
            actions << line;
            first = false;
        }
        actions << '\t' << program.inspect_action(action) << '\n';
    }
    std::ofstream summary(output + ".json");
    summary << "{\"shots\":" << shots << ",\"seed\":" << seed << ",\"accepted\":" << accepted_count
            << ",\"errors\":" << error_count << ",\"nanoseconds\":" << total_ns
            << ",\"accepted_ns\":" << accepted_ns << ",\"rejected_ns\":" << rejected_ns
            << ",\"record_hash\":\"" << record_hash
            << "\",\"peak_active_width\":" << plan.peak_active_width
            << ",\"actions\":" << program.num_actions() << "}\n";
}

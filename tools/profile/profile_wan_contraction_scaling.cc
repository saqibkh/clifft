// Numerical contraction sweeps, excluding physical fault/record binding.
#include "gadget_contraction_kernel.h"

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

int main(int argc, char** argv) {
    if (argc != 3)
        throw std::invalid_argument("usage: kernel directory repeats");
    const std::string directory = argv[1];
    const size_t repeats = std::stoull(argv[2]);
    std::ifstream metadata_file(directory + "/metadata.txt");
    gadget_study::Reader metadata{metadata_file};
    const size_t rank = metadata.size();
    const size_t leaves = metadata.size();
    std::vector<gadget_study::Contraction> plans;
    size_t scratch = 0, product = 0, lookups = 0;
    for (size_t k = 0; k <= rank; ++k) {
        std::ifstream plan_file(directory + "/marginal_" + std::to_string(k) + ".txt");
        gadget_study::Reader reader{plan_file};
        plans.emplace_back(reader, rank + k, 2 * leaves, 4);
        scratch = std::max(scratch, plans.back().scratch_size());
        product = std::max(product, plans.back().product_size());
        lookups += plans.back().lookup_bytes();
    }
    gadget_study::ContractionWorkspace workspace;
    workspace.prepare(scratch, product);
    std::vector<gadget_study::Complex> local(8 * leaves, 1);
    gadget_study::Complex checksum = 0;
    auto sweep = [&](size_t iteration) {
        // Changing a factor exercises fresh table values and blocks memoization.
        local[1] = (iteration & 1) ? gadget_study::Complex(0, 1) : gadget_study::Complex(1, 0);
        local[4 * leaves + 1] = std::conj(local[1]);
        for (const auto& plan : plans)
            checksum += plan.evaluate(local.data(), workspace);
    };
    sweep(0);
    const auto start = std::chrono::steady_clock::now();
    for (size_t iteration = 0; iteration < repeats; ++iteration)
        sweep(iteration);
    const double seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::cout << std::setprecision(17) << "{\"seconds\":" << seconds << ",\"sweeps\":" << repeats
              << ",\"marginal_plans\":" << plans.size() << ",\"lookup_bytes\":" << lookups
              << ",\"workspace_bytes\":" << 16 * (scratch + product)
              << ",\"checksum_real\":" << checksum.real()
              << ",\"checksum_imag\":" << checksum.imag() << "}\n";
}

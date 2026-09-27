// Fixed-plan arithmetic for fault-conditioned folded instrument validation.
#include "gadget_contraction_kernel.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>

namespace {
using namespace gadget_study;

struct Instrument {
    size_t rank, leaves;
    std::unique_ptr<Contraction> amplitude;
    std::vector<Contraction> marginals;
    ContractionWorkspace workspace;
    std::vector<Complex> local;

    explicit Instrument(const std::filesystem::path& path) {
        std::ifstream meta(path / "metadata.txt");
        Reader reader{meta};
        rank = reader.size();
        leaves = reader.size();
        if (leaves < rank)
            throw std::invalid_argument("missing Fourier factors");
        std::ifstream source(path / "amplitude.txt");
        Reader amp{source};
        amplitude = std::make_unique<Contraction>(amp, rank, leaves, 4);
        size_t scratch = amplitude->scratch_size(), product = amplitude->product_size();
        for (size_t k = 0; k <= rank; ++k) {
            std::ifstream input(path / ("marginal_" + std::to_string(k) + ".txt"));
            Reader plan{input};
            marginals.emplace_back(plan, rank + k, 2 * leaves, 4);
            scratch = std::max(scratch, marginals.back().scratch_size());
            product = std::max(product, marginals.back().product_size());
        }
        workspace.prepare(scratch, product);
        local.resize(8 * leaves);
    }

    void bind(const Complex* values, const uint8_t* syndrome, Complex* out,
              bool conjugate) noexcept {
        for (size_t i = 0; i < 4 * leaves; ++i)
            out[i] = conjugate ? std::conj(values[i]) : values[i];
        for (size_t i = 0; i < rank; ++i)
            if (syndrome[i])
                out[4 * (leaves - rank + i) + 1] *= -1;
    }

    double marginal(const Complex* coefficients, const Complex* values, const uint8_t* syndrome,
                    size_t measured) noexcept {
        assert(measured <= rank);
        double result = 0;
        for (size_t logical = 0; logical < 2; ++logical) {
            for (size_t i = 0; i < 4; ++i) {
                const size_t left = 4 * logical + i;
                bind(values + left * 4 * leaves, syndrome, local.data(), false);
                for (size_t j = 0; j <= i; ++j) {
                    const size_t right = 4 * logical + j;
                    bind(values + right * 4 * leaves, syndrome, local.data() + 4 * leaves, true);
                    result += (double(i == j ? 1 : 2) * coefficients[left] *
                               std::conj(coefficients[right]) *
                               marginals[measured].evaluate(local.data(), workspace))
                                  .real();
                }
            }
        }
        assert(result > -1e-10);
        return std::max(0.0, result);
    }

    void amplitudes(const Complex* coefficients, const Complex* values, const uint8_t* syndrome,
                    Complex* result) noexcept {
        for (size_t logical = 0; logical < 2; ++logical) {
            result[logical] = 0;
            for (size_t i = 0; i < 4; ++i) {
                const size_t term = 4 * logical + i;
                bind(values + term * 4 * leaves, syndrome, local.data(), false);
                result[logical] +=
                    coefficients[term] * amplitude->evaluate(local.data(), workspace);
            }
        }
    }
};
}  // namespace

extern "C" {
void* sahay_create(const char* path) {
    try {
        return new Instrument(path);
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return nullptr;
    }
}
void sahay_destroy(void* context) {
    delete static_cast<Instrument*>(context);
}
double sahay_marginal(void* context, const double* coefficients, const double* values,
                      const uint8_t* syndrome, size_t measured) {
    return static_cast<Instrument*>(context)->marginal(
        reinterpret_cast<const Complex*>(coefficients), reinterpret_cast<const Complex*>(values),
        syndrome, measured);
}
void sahay_amplitudes(void* context, const double* coefficients, const double* values,
                      const uint8_t* syndrome, double* result) {
    static_cast<Instrument*>(context)->amplitudes(reinterpret_cast<const Complex*>(coefficients),
                                                  reinterpret_cast<const Complex*>(values),
                                                  syndrome, reinterpret_cast<Complex*>(result));
}
}

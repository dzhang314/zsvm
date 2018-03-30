#ifndef DZNL_PERMUTATION_HPP_INCLUDED
#define DZNL_PERMUTATION_HPP_INCLUDED

// C++ standard library headers
#include <algorithm> // for std::next_permutation
#include <cstddef> // for std::size_t
#include <stdexcept> // for std::invalid_argument
#include <vector> // for std::vector

// Project-specific headers
#include "Particle.hpp" // for zsvm::Particle

namespace dznl {

    bool is_invariant_permutation(const std::vector<zsvm::Particle> &particles,
                                  const std::vector<std::size_t> &permutation) {
        const std::size_t n = particles.size();
        if (permutation.size() != n) {
            throw std::invalid_argument(
                    "is_invariant_permutation received particle and "
                    "permutation vectors of different sizes");
        }
        for (std::size_t i = 0; i < n; ++i) {
            if (particles[i].type_id != particles[permutation[i]].type_id) {
                return false;
            }
        }
        return true;
    }

    std::size_t count_changes(const std::vector<zsvm::Particle> &particles,
                              const std::vector<std::size_t> &permutation) {
        const std::size_t n = particles.size();
        if (permutation.size() != n) {
            throw std::invalid_argument(
                    "is_invariant_permutation received particle and "
                    "permutation vectors of different sizes");
        }
        std::size_t count = 0;
        for (std::size_t i = 0; i < n; ++i) {
            if (particles[i].type_id != particles[permutation[i]].type_id) {
                ++count;
            }
        }
        return count;
    }

    std::size_t count_inversions(const std::vector<std::size_t> &permutation) {
        const std::size_t n = permutation.size();
        std::size_t count = 0;
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = i + 1; j < n; ++j) {
                if (permutation[i] > permutation[j]) { ++count; }
            }
        }
        return count;
    }

    std::vector<std::vector<std::size_t>> invariant_permutations(
            const std::vector<zsvm::Particle> &particles) {
        const std::size_t n = particles.size();
        std::vector<std::size_t> permutation(n);
        for (std::size_t i = 0; i < n; ++i) { permutation[i] = i; }
        std::vector<std::vector<size_t>> result;
        do {
            if (is_invariant_permutation(particles, permutation)) {
                result.push_back(permutation);
            }
        } while (std::next_permutation(permutation.begin(), permutation.end()));
        return result;
    }

} // namespace dznl

#endif // DZNL_PERMUTATION_HPP_INCLUDED

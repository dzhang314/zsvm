#ifndef DZNL_PERMUTATION_HPP_INCLUDED
#define DZNL_PERMUTATION_HPP_INCLUDED

// C++ standard library headers
#include <cstddef> // for std::size_t
#include <vector> // for std::vector

// Project-specific headers
#include "Particle.hpp" // for zsvm::Spin

namespace dznl {

    bool is_invariant_permutation(const std::vector<int> &items,
                                  const std::vector<std::size_t> &permutation);

    std::size_t count_changes(const std::vector<zsvm::Spin> &items,
                              const std::vector<std::size_t> &permutation);

    std::size_t count_inversions(const std::vector<std::size_t> &permutation);

    std::vector<std::vector<std::size_t>> invariant_permutations(
            const std::vector<int> &items);

} // namespace dznl

#endif // DZNL_PERMUTATION_HPP_INCLUDED

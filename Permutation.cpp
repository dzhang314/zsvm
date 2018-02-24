#include "Permutation.hpp"

// C++ standard library headers
#include <algorithm> // for td::next_permutation
#include <cstddef> // for std::size_t
#include <stdexcept> // for std::invalid_argument


bool dznl::is_invariant_permutation(
        const std::vector<int> &items,
        const std::vector<std::size_t> &permutation) {
    const std::size_t n = items.size();
    if (permutation.size() != n) {
        throw std::invalid_argument(
                "is_invariant_permutation received item and "
                        "permutation vectors of different sizes");
    }
    for (std::size_t i = 0; i < n; ++i) {
        if (items[i] != items[permutation[i]]) { return false; }
    }
    return true;
}


std::size_t dznl::count_changes(const std::vector<zsvm::Spin> &items,
                                const std::vector<std::size_t> &permutation) {
    const std::size_t n = items.size();
    if (permutation.size() != n) {
        throw std::invalid_argument(
                "is_invariant_permutation received item and "
                        "permutation vectors of different sizes");
    }
    std::size_t count = 0;
    for (std::size_t i = 0; i < n; ++i) {
        if (!(items[i] == items[permutation[i]])) { ++count; }
    }
    return count;
}


std::size_t dznl::count_inversions(
        const std::vector<std::size_t> &permutation) {
    const std::size_t n = permutation.size();
    std::size_t count = 0;
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            if (permutation[i] > permutation[j]) { ++count; }
        }
    }
    return count;
}


std::vector<std::vector<std::size_t>> dznl::invariant_permutations(
        const std::vector<int> &items) {
    const std::size_t n = items.size();
    std::vector<std::size_t> permutation(n);
    for (std::size_t i = 0; i < n; ++i) { permutation[i] = i; }
    std::vector<std::vector<size_t>> result;
    do {
        if (is_invariant_permutation(items, permutation)) {
            result.push_back(permutation);
        }
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return result;
}

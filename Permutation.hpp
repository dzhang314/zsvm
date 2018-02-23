#ifndef DZNL_PERMUTATION_HPP_INCLUDED
#define DZNL_PERMUTATION_HPP_INCLUDED

#include <algorithm>
#include <cstddef>
#include <vector>
#include <stdexcept>

namespace dznl {

    template <typename T>
    bool is_invariant_permutation(const std::vector<T> &items,
                                  const std::vector<std::size_t> &permutation) {
        const std::size_t n = items.size();
        if (permutation.size() != n) {
            throw std::invalid_argument(
                    "is_invariant_permutation received item and "
                            "permutation vectors of different sizes");
        }
        for (std::size_t i = 0; i < n; ++i) {
            if (!(items[i] == items[permutation[i]])) {
                return false;
            }
        }
        return true;
    }

    template <typename T>
    std::vector<std::vector<std::size_t>> invariant_permutations(
            const std::vector<T> &items) {
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

} // namespace dznl

#endif // DZNL_PERMUTATION_HPP_INCLUDED

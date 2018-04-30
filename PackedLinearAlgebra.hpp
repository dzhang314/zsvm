#ifndef ZSVM_PACKED_LINEAR_ALGEBRA_HPP_INCLUDED
#define ZSVM_PACKED_LINEAR_ALGEBRA_HPP_INCLUDED

// C++ standard library headers
#include <cmath> // for std::sqrt

// Project-specific headers
#include "PackedLinearAlgebraImpl.hpp"

template <typename T>
T half_inverse_pow(const T &x, long long int n) {
    using std::sqrt;
    using std::pow;
    switch (n) {
        case 0:
            return 1;
        case 1:
            return 1 / sqrt(x);
        case 2:
            return 1 / x;
        case 3:
            return 1 / (x * sqrt(x));
        case 4:
            return 1 / (x * x);
        case 5:
            return 1 / (x * x * sqrt(x));
        case 6:
            return 1 / (x * x * x);
        case 7:
            return 1 / (x * x * x * sqrt(x));
        default:
            return 1 / sqrt(pow(x, n));
    }
}

// CLion does not understand how to format function type alias declarations.
// @formatter:off

template <typename T>
using packed_determinant_inverse_function = T (*)(
        const T *__restrict__ x, T *__restrict__ y);

template <typename T>
using packed_kinetic_trace_function = T (*)(
        const T *__restrict__ a, const T *__restrict__ b,
        const T *__restrict__ c, const T *__restrict__ m);

template <typename T>
using packed_quadratic_form_function = T (*)(
        const T *__restrict__ x, const T *__restrict__ v);

template <typename T>
using packed_permutation_conjugate_function = void (*)(
        const T *__restrict__ x, const T *__restrict__ p,
        T *__restrict__ y);

// @formatter:on

template <typename T>
struct PackedDeterminantInverse {
    static const packed_determinant_inverse_function<T> *functions() {
        static const packed_determinant_inverse_function<T> funcs[] = {
                packed_determinant_inverse_1<T>,
                packed_determinant_inverse_2<T>,
                packed_determinant_inverse_3<T>,
                packed_determinant_inverse_4<T>,
                packed_determinant_inverse_5<T>,
        };
        return funcs;
    }
};

template <typename T>
struct PackedKineticTrace {
    static const packed_kinetic_trace_function<T> *functions() {
        static const packed_kinetic_trace_function<T> funcs[] = {
                packed_kinetic_trace_1<T>,
                packed_kinetic_trace_2<T>,
                packed_kinetic_trace_3<T>,
                packed_kinetic_trace_4<T>,
                packed_kinetic_trace_5<T>,
        };
        return funcs;
    }
};

template <typename T>
struct PackedQuadraticForm {
    static const packed_quadratic_form_function<T> *functions() {
        static const packed_quadratic_form_function<T> funcs[] = {
                packed_quadratic_form_1<T>,
                packed_quadratic_form_2<T>,
                packed_quadratic_form_3<T>,
                packed_quadratic_form_4<T>,
                packed_quadratic_form_5<T>,
        };
        return funcs;
    }
};

template <typename T>
struct PackedPermutationConjugate {
    static const packed_permutation_conjugate_function<T> *functions() {
        static const packed_permutation_conjugate_function<T> funcs[] = {
                packed_permutation_conjugate_1<T>,
                packed_permutation_conjugate_2<T>,
                packed_permutation_conjugate_3<T>,
                packed_permutation_conjugate_4<T>,
                packed_permutation_conjugate_5<T>,
        };
        return funcs;
    }
};

#endif // ZSVM_PACKED_LINEAR_ALGEBRA_HPP_INCLUDED

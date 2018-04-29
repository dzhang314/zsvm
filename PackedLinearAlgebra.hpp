#ifndef ZSVM_PACKED_LINEAR_ALGEBRA_HPP_INCLUDED
#define ZSVM_PACKED_LINEAR_ALGEBRA_HPP_INCLUDED

#include "PackedLinearAlgebraImpl.hpp"

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
const packed_determinant_inverse_function<T> PACKED_DETERMINANT_INVERSE[] = {
        packed_determinant_inverse_1<T>,
        packed_determinant_inverse_2<T>,
        packed_determinant_inverse_3<T>,
        packed_determinant_inverse_4<T>,
        packed_determinant_inverse_5<T>,
};

template <typename T>
const packed_kinetic_trace_function<T> PACKED_KINETIC_TRACE[] = {
        packed_kinetic_trace_1<T>,
        packed_kinetic_trace_2<T>,
        packed_kinetic_trace_3<T>,
        packed_kinetic_trace_4<T>,
        packed_kinetic_trace_5<T>,
};

template <typename T>
const packed_quadratic_form_function<T> PACKED_QUADRATIC_FORM[] = {
        packed_quadratic_form_1<T>,
        packed_quadratic_form_2<T>,
        packed_quadratic_form_3<T>,
        packed_quadratic_form_4<T>,
        packed_quadratic_form_5<T>,
};

template <typename T>
const packed_permutation_conjugate_function<T> PACKED_PERMUTATION_CONJUGATE[] = {
        packed_permutation_conjugate_1<T>,
        packed_permutation_conjugate_2<T>,
        packed_permutation_conjugate_3<T>,
        packed_permutation_conjugate_4<T>,
        packed_permutation_conjugate_5<T>,
};

#endif // ZSVM_PACKED_LINEAR_ALGEBRA_HPP_INCLUDED

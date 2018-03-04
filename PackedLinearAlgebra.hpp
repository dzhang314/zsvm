#ifndef ZSVM_PACKED_LINEAR_ALGEBRA_HPP_INCLUDED
#define ZSVM_PACKED_LINEAR_ALGEBRA_HPP_INCLUDED

double packed_determinant_inverse_4(
        const double *__restrict__ x, double *__restrict__ y);

double packed_kinetic_trace_4(
        const double *__restrict__ a, const double *__restrict__ b,
        const double *__restrict__ c, const double *__restrict__ m);

double packed_quadratic_form_4(
        const double *__restrict__ x, const double *__restrict__ v);

void packed_permutation_conjugate_4(
        const double *__restrict__ x, const double *__restrict__ p,
        double *__restrict__ y);

#endif // ZSVM_PACKED_LINEAR_ALGEBRA_HPP_INCLUDED

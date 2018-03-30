#ifndef ZSVM_PACKED_LINEAR_ALGEBRA_HPP_INCLUDED
#define ZSVM_PACKED_LINEAR_ALGEBRA_HPP_INCLUDED

typedef double (*packed_determinant_inverse_function_t)(
        const double *__restrict__ x, double *__restrict__ y);

typedef double (*packed_kinetic_trace_function_t)(
        const double *__restrict__ a, const double *__restrict__ b,
        const double *__restrict__ c, const double *__restrict__ m);

typedef double (*packed_quadratic_form_function_t)(
        const double *__restrict__ x, const double *__restrict__ v);

typedef void (*packed_permutation_conjugate_function_t)(
        const double *__restrict__ x, const double *__restrict__ p,
        double *__restrict__ y);

extern const packed_determinant_inverse_function_t
        PACKED_DETERMINANT_INVERSE[];

extern const packed_kinetic_trace_function_t
        PACKED_KINETIC_TRACE[];

extern const packed_quadratic_form_function_t
        PACKED_QUADRATIC_FORM[];

extern const packed_permutation_conjugate_function_t
        PACKED_PERMUTATION_CONJUGATE[];

double packed_determinant_inverse_1(
        const double *__restrict__ x, double *__restrict__ y) noexcept;

double packed_kinetic_trace_1(
        const double *__restrict__ a, const double *__restrict__ b,
        const double *__restrict__ c, const double *__restrict__ m) noexcept;

double packed_quadratic_form_1(
        const double *__restrict__ x, const double *__restrict__ v) noexcept;

void packed_permutation_conjugate_1(
        const double *__restrict__ x, const double *__restrict__ p,
        double *__restrict__ y) noexcept;

double packed_determinant_inverse_2(
        const double *__restrict__ x, double *__restrict__ y) noexcept;

double packed_kinetic_trace_2(
        const double *__restrict__ a, const double *__restrict__ b,
        const double *__restrict__ c, const double *__restrict__ m) noexcept;

double packed_quadratic_form_2(
        const double *__restrict__ x, const double *__restrict__ v) noexcept;

void packed_permutation_conjugate_2(
        const double *__restrict__ x, const double *__restrict__ p,
        double *__restrict__ y) noexcept;

double packed_determinant_inverse_3(
        const double *__restrict__ x, double *__restrict__ y) noexcept;

double packed_kinetic_trace_3(
        const double *__restrict__ a, const double *__restrict__ b,
        const double *__restrict__ c, const double *__restrict__ m) noexcept;

double packed_quadratic_form_3(
        const double *__restrict__ x, const double *__restrict__ v) noexcept;

void packed_permutation_conjugate_3(
        const double *__restrict__ x, const double *__restrict__ p,
        double *__restrict__ y) noexcept;

double packed_determinant_inverse_4(
        const double *__restrict__ x, double *__restrict__ y) noexcept;

double packed_kinetic_trace_4(
        const double *__restrict__ a, const double *__restrict__ b,
        const double *__restrict__ c, const double *__restrict__ m) noexcept;

double packed_quadratic_form_4(
        const double *__restrict__ x, const double *__restrict__ v) noexcept;

void packed_permutation_conjugate_4(
        const double *__restrict__ x, const double *__restrict__ p,
        double *__restrict__ y) noexcept;

double packed_determinant_inverse_5(
        const double *__restrict__ x, double *__restrict__ y) noexcept;

double packed_kinetic_trace_5(
        const double *__restrict__ a, const double *__restrict__ b,
        const double *__restrict__ c, const double *__restrict__ m) noexcept;

double packed_quadratic_form_5(
        const double *__restrict__ x, const double *__restrict__ v) noexcept;

void packed_permutation_conjugate_5(
        const double *__restrict__ x, const double *__restrict__ p,
        double *__restrict__ y) noexcept;

#endif // ZSVM_PACKED_LINEAR_ALGEBRA_HPP_INCLUDED

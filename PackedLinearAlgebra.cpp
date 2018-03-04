#include "PackedLinearAlgebra.hpp"


double packed_determinant_inverse_4(
        const double *__restrict__ x, double *__restrict__ y) {
    const double det =
            x[4] * x[4] * x[6] * x[6]
            - x[2] * x[5] * x[6] * x[6]
            - 2 * x[3] * x[4] * x[6] * x[7]
            + 2 * x[1] * x[5] * x[6] * x[7]
            + x[3] * x[3] * x[7] * x[7]
            - x[0] * x[5] * x[7] * x[7]
            + 2 * x[2] * x[3] * x[6] * x[8]
            - 2 * x[1] * x[4] * x[6] * x[8]
            - 2 * x[1] * x[3] * x[7] * x[8]
            + 2 * x[0] * x[4] * x[7] * x[8]
            + x[1] * x[1] * x[8] * x[8]
            - x[0] * x[2] * x[8] * x[8]
            - x[2] * x[3] * x[3] * x[9]
            + 2 * x[1] * x[3] * x[4] * x[9]
            - x[0] * x[4] * x[4] * x[9]
            - x[1] * x[1] * x[5] * x[9]
            + x[0] * x[2] * x[5] * x[9];
    y[0] = (-x[5] * x[7] * x[7]
            + 2 * x[4] * x[7] * x[8]
            - x[2] * x[8] * x[8]
            - x[4] * x[4] * x[9]
            + x[2] * x[5] * x[9]) / det;
    y[1] = (x[5] * x[6] * x[7]
            - x[4] * x[6] * x[8]
            - x[3] * x[7] * x[8]
            + x[1] * x[8] * x[8]
            + x[3] * x[4] * x[9]
            - x[1] * x[5] * x[9]) / det;
    y[2] = (-x[5] * x[6] * x[6]
            + 2 * x[3] * x[6] * x[8]
            - x[0] * x[8] * x[8]
            - x[3] * x[3] * x[9]
            + x[0] * x[5] * x[9]) / det;
    y[3] = (-x[4] * x[6] * x[7]
            + x[3] * x[7] * x[7]
            + x[2] * x[6] * x[8]
            - x[1] * x[7] * x[8]
            - x[2] * x[3] * x[9]
            + x[1] * x[4] * x[9]) / det;
    y[4] = (x[4] * x[6] * x[6]
            - x[3] * x[6] * x[7]
            - x[1] * x[6] * x[8]
            + x[0] * x[7] * x[8]
            + x[1] * x[3] * x[9]
            - x[0] * x[4] * x[9]) / det;
    y[5] = (-x[2] * x[6] * x[6]
            + 2 * x[1] * x[6] * x[7]
            - x[0] * x[7] * x[7]
            - x[1] * x[1] * x[9]
            + x[0] * x[2] * x[9]) / det;
    y[6] = (x[4] * x[4] * x[6]
            - x[2] * x[5] * x[6]
            - x[3] * x[4] * x[7]
            + x[1] * x[5] * x[7]
            + x[2] * x[3] * x[8]
            - x[1] * x[4] * x[8]) / det;
    y[7] = (-x[3] * x[4] * x[6]
            + x[1] * x[5] * x[6]
            + x[3] * x[3] * x[7]
            - x[0] * x[5] * x[7]
            - x[1] * x[3] * x[8]
            + x[0] * x[4] * x[8]) / det;
    y[8] = (x[2] * x[3] * x[6]
            - x[1] * x[4] * x[6]
            - x[1] * x[3] * x[7]
            + x[0] * x[4] * x[7]
            + x[1] * x[1] * x[8]
            - x[0] * x[2] * x[8]) / det;
    y[9] = (-x[2] * x[3] * x[3]
            + 2 * x[1] * x[3] * x[4]
            - x[0] * x[4] * x[4]
            - x[1] * x[1] * x[5]
            + x[0] * x[2] * x[5]) / det;
    return det;
}


double packed_kinetic_trace_4(
        const double *__restrict__ a, const double *__restrict__ b,
        const double *__restrict__ c, const double *__restrict__ m) {
    return a[0] * b[0] * c[0] * m[0]
           + a[1] * b[0] * c[1] * m[0]
           + a[0] * b[1] * c[1] * m[0]
           + a[1] * b[1] * c[2] * m[0]
           + a[3] * b[0] * c[3] * m[0]
           + a[0] * b[3] * c[3] * m[0]
           + a[3] * b[1] * c[4] * m[0]
           + a[1] * b[3] * c[4] * m[0]
           + a[3] * b[3] * c[5] * m[0]
           + a[6] * b[0] * c[6] * m[0]
           + a[0] * b[6] * c[6] * m[0]
           + a[6] * b[1] * c[7] * m[0]
           + a[1] * b[6] * c[7] * m[0]
           + a[6] * b[3] * c[8] * m[0]
           + a[3] * b[6] * c[8] * m[0]
           + a[6] * b[6] * c[9] * m[0]
           + a[1] * b[1] * c[0] * m[1]
           + a[2] * b[1] * c[1] * m[1]
           + a[1] * b[2] * c[1] * m[1]
           + a[2] * b[2] * c[2] * m[1]
           + a[4] * b[1] * c[3] * m[1]
           + a[1] * b[4] * c[3] * m[1]
           + a[4] * b[2] * c[4] * m[1]
           + a[2] * b[4] * c[4] * m[1]
           + a[4] * b[4] * c[5] * m[1]
           + a[7] * b[1] * c[6] * m[1]
           + a[1] * b[7] * c[6] * m[1]
           + a[7] * b[2] * c[7] * m[1]
           + a[2] * b[7] * c[7] * m[1]
           + a[7] * b[4] * c[8] * m[1]
           + a[4] * b[7] * c[8] * m[1]
           + a[7] * b[7] * c[9] * m[1]
           + a[3] * b[3] * c[0] * m[2]
           + a[4] * b[3] * c[1] * m[2]
           + a[3] * b[4] * c[1] * m[2]
           + a[4] * b[4] * c[2] * m[2]
           + a[5] * b[3] * c[3] * m[2]
           + a[3] * b[5] * c[3] * m[2]
           + a[5] * b[4] * c[4] * m[2]
           + a[4] * b[5] * c[4] * m[2]
           + a[5] * b[5] * c[5] * m[2]
           + a[8] * b[3] * c[6] * m[2]
           + a[3] * b[8] * c[6] * m[2]
           + a[8] * b[4] * c[7] * m[2]
           + a[4] * b[8] * c[7] * m[2]
           + a[8] * b[5] * c[8] * m[2]
           + a[5] * b[8] * c[8] * m[2]
           + a[8] * b[8] * c[9] * m[2]
           + a[6] * b[6] * c[0] * m[3]
           + a[7] * b[6] * c[1] * m[3]
           + a[6] * b[7] * c[1] * m[3]
           + a[7] * b[7] * c[2] * m[3]
           + a[8] * b[6] * c[3] * m[3]
           + a[6] * b[8] * c[3] * m[3]
           + a[8] * b[7] * c[4] * m[3]
           + a[7] * b[8] * c[4] * m[3]
           + a[8] * b[8] * c[5] * m[3]
           + a[9] * b[6] * c[6] * m[3]
           + a[6] * b[9] * c[6] * m[3]
           + a[9] * b[7] * c[7] * m[3]
           + a[7] * b[9] * c[7] * m[3]
           + a[9] * b[8] * c[8] * m[3]
           + a[8] * b[9] * c[8] * m[3]
           + a[9] * b[9] * c[9] * m[3];
}


double packed_quadratic_form_4(
        const double *__restrict__ x, const double *__restrict__ v) {
    return v[0] * v[0] * x[0]
           + 2 * v[0] * v[1] * x[1]
           + v[1] * v[1] * x[2]
           + 2 * v[0] * v[2] * x[3]
           + 2 * v[1] * v[2] * x[4]
           + v[2] * v[2] * x[5]
           + 2 * v[0] * v[3] * x[6]
           + 2 * v[1] * v[3] * x[7]
           + 2 * v[2] * v[3] * x[8]
           + v[3] * v[3] * x[9];
}


void packed_permutation_conjugate_4(
        const double *__restrict__ x, const double *__restrict__ p,
        double *__restrict__ y) {
    y[0] = p[0] * p[0] * x[0]
           + 2 * p[0] * p[1] * x[1]
           + p[1] * p[1] * x[2]
           + 2 * p[0] * p[2] * x[3]
           + 2 * p[1] * p[2] * x[4]
           + p[2] * p[2] * x[5]
           + 2 * p[0] * p[3] * x[6]
           + 2 * p[1] * p[3] * x[7]
           + 2 * p[2] * p[3] * x[8]
           + p[3] * p[3] * x[9];
    y[1] = p[0] * p[4] * x[0]
           + p[1] * p[4] * x[1]
           + p[0] * p[5] * x[1]
           + p[1] * p[5] * x[2]
           + p[2] * p[4] * x[3]
           + p[0] * p[6] * x[3]
           + p[2] * p[5] * x[4]
           + p[1] * p[6] * x[4]
           + p[2] * p[6] * x[5]
           + p[3] * p[4] * x[6]
           + p[0] * p[7] * x[6]
           + p[3] * p[5] * x[7]
           + p[1] * p[7] * x[7]
           + p[3] * p[6] * x[8]
           + p[2] * p[7] * x[8]
           + p[3] * p[7] * x[9];
    y[2] = p[4] * p[4] * x[0]
           + 2 * p[4] * p[5] * x[1]
           + p[5] * p[5] * x[2]
           + 2 * p[4] * p[6] * x[3]
           + 2 * p[5] * p[6] * x[4]
           + p[6] * p[6] * x[5]
           + 2 * p[4] * p[7] * x[6]
           + 2 * p[5] * p[7] * x[7]
           + 2 * p[6] * p[7] * x[8]
           + p[7] * p[7] * x[9];
    y[3] = p[0] * p[8] * x[0]
           + p[1] * p[8] * x[1]
           + p[0] * p[9] * x[1]
           + p[1] * p[9] * x[2]
           + p[2] * p[8] * x[3]
           + p[0] * p[10] * x[3]
           + p[2] * p[9] * x[4]
           + p[1] * p[10] * x[4]
           + p[2] * p[10] * x[5]
           + p[3] * p[8] * x[6]
           + p[0] * p[11] * x[6]
           + p[3] * p[9] * x[7]
           + p[1] * p[11] * x[7]
           + p[3] * p[10] * x[8]
           + p[2] * p[11] * x[8]
           + p[3] * p[11] * x[9];
    y[4] = p[4] * p[8] * x[0]
           + p[5] * p[8] * x[1]
           + p[4] * p[9] * x[1]
           + p[5] * p[9] * x[2]
           + p[6] * p[8] * x[3]
           + p[4] * p[10] * x[3]
           + p[6] * p[9] * x[4]
           + p[5] * p[10] * x[4]
           + p[6] * p[10] * x[5]
           + p[7] * p[8] * x[6]
           + p[4] * p[11] * x[6]
           + p[7] * p[9] * x[7]
           + p[5] * p[11] * x[7]
           + p[7] * p[10] * x[8]
           + p[6] * p[11] * x[8]
           + p[7] * p[11] * x[9];
    y[5] = p[8] * p[8] * x[0]
           + 2 * p[8] * p[9] * x[1]
           + p[9] * p[9] * x[2]
           + 2 * p[8] * p[10] * x[3]
           + 2 * p[9] * p[10] * x[4]
           + p[10] * p[10] * x[5]
           + 2 * p[8] * p[11] * x[6]
           + 2 * p[9] * p[11] * x[7]
           + 2 * p[10] * p[11] * x[8]
           + p[11] * p[11] * x[9];
    y[6] = p[0] * p[12] * x[0]
           + p[1] * p[12] * x[1]
           + p[0] * p[13] * x[1]
           + p[1] * p[13] * x[2]
           + p[2] * p[12] * x[3]
           + p[0] * p[14] * x[3]
           + p[2] * p[13] * x[4]
           + p[1] * p[14] * x[4]
           + p[2] * p[14] * x[5]
           + p[3] * p[12] * x[6]
           + p[0] * p[15] * x[6]
           + p[3] * p[13] * x[7]
           + p[1] * p[15] * x[7]
           + p[3] * p[14] * x[8]
           + p[2] * p[15] * x[8]
           + p[3] * p[15] * x[9];
    y[7] = p[4] * p[12] * x[0]
           + p[5] * p[12] * x[1]
           + p[4] * p[13] * x[1]
           + p[5] * p[13] * x[2]
           + p[6] * p[12] * x[3]
           + p[4] * p[14] * x[3]
           + p[6] * p[13] * x[4]
           + p[5] * p[14] * x[4]
           + p[6] * p[14] * x[5]
           + p[7] * p[12] * x[6]
           + p[4] * p[15] * x[6]
           + p[7] * p[13] * x[7]
           + p[5] * p[15] * x[7]
           + p[7] * p[14] * x[8]
           + p[6] * p[15] * x[8]
           + p[7] * p[15] * x[9];
    y[8] = p[8] * p[12] * x[0]
           + p[9] * p[12] * x[1]
           + p[8] * p[13] * x[1]
           + p[9] * p[13] * x[2]
           + p[10] * p[12] * x[3]
           + p[8] * p[14] * x[3]
           + p[10] * p[13] * x[4]
           + p[9] * p[14] * x[4]
           + p[10] * p[14] * x[5]
           + p[11] * p[12] * x[6]
           + p[8] * p[15] * x[6]
           + p[11] * p[13] * x[7]
           + p[9] * p[15] * x[7]
           + p[11] * p[14] * x[8]
           + p[10] * p[15] * x[8]
           + p[11] * p[15] * x[9];
    y[9] = p[12] * p[12] * x[0]
           + 2 * p[12] * p[13] * x[1]
           + p[13] * p[13] * x[2]
           + 2 * p[12] * p[14] * x[3]
           + 2 * p[13] * p[14] * x[4]
           + p[14] * p[14] * x[5]
           + 2 * p[12] * p[15] * x[6]
           + 2 * p[13] * p[15] * x[7]
           + 2 * p[14] * p[15] * x[8]
           + p[15] * p[15] * x[9];
}

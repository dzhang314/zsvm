#ifndef ZSVM_REAL_VARIATIONAL_SOLVER_HPP_INCLUDED
#define ZSVM_REAL_VARIATIONAL_SOLVER_HPP_INCLUDED

// C++ standard library headers
#include <cmath> // for std::sqrt, std::pow, std::abs
#include <cstddef> // for std::size_t
#include <limits>

// Eigen linear algebra library headers
#include <Eigen/Core>
#include <Eigen/Eigenvalues>

namespace zsvm {

    template <typename T>
    class RealVariationalSolver {

    private: // =============================================== MEMBER VARIABLES

        typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> MatrixXT;
        typedef Eigen::Matrix<T, Eigen::Dynamic, 1> VectorXT;
        typedef Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> ArrayXT;

        std::size_t size;
        MatrixXT overlap_matrix;
        MatrixXT hamiltonian_matrix;
        VectorXT eigenvalues;
        MatrixXT eigenvectors;
        bool clean;

    public: // ===================================================== CONSTRUCTOR

        RealVariationalSolver() : size(0), clean(false) {}

    public: // ======================================================== MUTATORS

        void set_basis_size_conservative(std::size_t basis_size) {
            size = basis_size;
            overlap_matrix.conservativeResize(basis_size, basis_size);
            hamiltonian_matrix.conservativeResize(basis_size, basis_size);
            clean = false;
        }

        void set_basis_size_destructive(std::size_t basis_size) {
            size = basis_size;
            overlap_matrix.resize(basis_size, basis_size);
            hamiltonian_matrix.resize(basis_size, basis_size);
            clean = false;
        }

        T &overlap_matrix_element(std::size_t i, std::size_t j) {
            clean = false;
            return overlap_matrix(i, j);
        }

        T &hamiltonian_matrix_element(std::size_t i, std::size_t j) {
            clean = false;
            return hamiltonian_matrix(i, j);
        }

        void compute_eigenstates() {
            if (!clean) {
                Eigen::GeneralizedSelfAdjointEigenSolver eigen_solver(
                        hamiltonian_matrix, overlap_matrix,
                        Eigen::ComputeEigenvectors | Eigen::Ax_lBx);
                eigenvalues = eigen_solver.eigenvalues();
                eigenvectors = eigen_solver.eigenvectors();
                clean = true;
            }
        }

    public: // ======================================================== ACCESORS

        bool empty() const { return (size == 0); }

        T get_eigenvalue(std::size_t i) {
            if (!clean) { compute_eigenstates(); }
            return eigenvalues(i);
        }

    private: // ========================== EIGENVALUE COMPUTATION HELPER METHODS

        T secular_objective_function(const T &x, const T &t,
                                     const ArrayXT &beta) const {
            return t - x - (beta / (eigenvalues.array() - x)).sum();
        }

        void find_lower_bracketing_interval(
                T &lower_bound, T &upper_bound,
                const T &strict_upper_bound,
                const T &t, const ArrayXT &beta) const {
            using std::ldexp;
            const T one(1);
            int k = 0;
            lower_bound = strict_upper_bound - ldexp(one, -k);
            upper_bound = strict_upper_bound - ldexp(one, -k - 1);
            T lower_objective = secular_objective_function(
                    lower_bound, t, beta);
            T upper_objective = secular_objective_function(
                    upper_bound, t, beta);
            while (true) {
                if (lower_objective > 0 && upper_objective > 0) {
                    ++k;
                    lower_bound = upper_bound;
                    lower_objective = upper_objective;
                    upper_bound = strict_upper_bound - ldexp(one, -k - 1);
                    upper_objective = secular_objective_function(
                            upper_bound, t, beta);
                } else if (lower_objective < 0 && upper_objective < 0) {
                    --k;
                    upper_bound = lower_bound;
                    upper_objective = lower_objective;
                    lower_bound = strict_upper_bound - ldexp(one, -k);
                    lower_objective = secular_objective_function(
                            lower_bound, t, beta);
                } else {
                    break;
                }
            }
        }

        T solve_secular_equation_bisection(
                T lower_bound, T upper_bound,
                const T &t, const ArrayXT &beta) const {
            using std::abs;
            const T tolerance = 64 * std::numeric_limits<T>::epsilon();
            while (true) {
                const T midpoint = (lower_bound + upper_bound) / 2;
                const T relative_difference = abs(
                        (upper_bound - lower_bound) / midpoint);
                if (relative_difference < tolerance) { return midpoint; }
                const T midpoint_objective = secular_objective_function(
                        midpoint, t, beta);
                if (midpoint_objective == 0) {
                    return midpoint;
                } else if (midpoint_objective > 0) {
                    lower_bound = midpoint;
                } else if (midpoint_objective < 0) {
                    upper_bound = midpoint;
                } else {
                    return std::numeric_limits<T>::quiet_NaN();
                }
            }
        }

    public: // ============================= FAST EIGENVALUE COMPUTATION METHODS

        T minimum_augmented_eigenvalue(
                const T *new_overlap_column,
                const T *new_hamiltonian_column) {
            using std::sqrt;
            if (!clean) { compute_eigenstates(); }
            VectorXT overlap_vector(size);
            for (std::size_t i = 0; i < size; ++i) {
                overlap_vector(i) = new_overlap_column[i];
            }
            VectorXT hamiltonian_vector(size);
            for (std::size_t i = 0; i < size; ++i) {
                hamiltonian_vector(i) = new_hamiltonian_column[i];
            }
            const VectorXT alpha = eigenvectors.transpose() * overlap_vector;
            const T norm_factor = 1 / sqrt(
                    new_overlap_column[size] - alpha.squaredNorm());
            const VectorXT phi = -norm_factor * (eigenvectors * alpha);
            const VectorXT psi = hamiltonian_matrix * phi +
                                 norm_factor * hamiltonian_vector;
            const ArrayXT beta =
                    (eigenvectors.transpose() * psi).array().square();
            const T t = phi.dot(psi) + norm_factor * (
                    hamiltonian_vector.dot(phi) +
                    norm_factor * new_hamiltonian_column[size]);
            T lower_bound, upper_bound;
            find_lower_bracketing_interval(
                    lower_bound, upper_bound, eigenvalues[0], t, beta);
            return solve_secular_equation_bisection(
                    lower_bound, upper_bound, t, beta);
        }

    }; // class RealVariationalSolver

} // namespace zsvm

#endif // ZSVM_REAL_VARIATIONAL_SOLVER_HPP_INCLUDED

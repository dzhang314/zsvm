#ifndef ZSVM_REAL_VARIATIONAL_SOLVER_HPP_INCLUDED
#define ZSVM_REAL_VARIATIONAL_SOLVER_HPP_INCLUDED

// C++ standard library headers
#include <cstddef> // for std::size_t

// Eigen linear algebra library headers
#include <Eigen/Core> // for Eigen::MatrixXd

namespace zsvm {

    class RealVariationalSolver {

    private: // =============================================== MEMBER VARIABLES

        std::size_t size;
        Eigen::MatrixXd overlap_matrix;
        Eigen::MatrixXd hamiltonian_matrix;
        Eigen::VectorXd eigenvalues;
        Eigen::MatrixXd eigenvectors;
        bool clean;

    public: // ===================================================== CONSTRUCTOR

        RealVariationalSolver();

    public: // ======================================================== MUTATORS

        void set_basis_size_conservative(std::size_t basis_size);

        void set_basis_size_destructive(std::size_t basis_size);

        double &overlap_matrix_element(std::size_t i, std::size_t j);

        double &hamiltonian_matrix_element(std::size_t i, std::size_t j);

        void compute_eigenstates();

    public: // ======================================================== ACCESORS

        bool empty() const;

        double get_eigenvalue(std::size_t i);

    private: // ========================== EIGENVALUE COMPUTATION HELPER METHODS

        double secular_objective_function(double x, double t,
                                          const Eigen::ArrayXd &beta) const;

        void find_lower_bracketing_interval(
                double &lower_bound, double &upper_bound,
                double strict_upper_bound,
                double t, const Eigen::ArrayXd &beta) const;

        double solve_secular_equation_bisection(
                double lower_bound, double upper_bound,
                double t, const Eigen::ArrayXd &beta) const;

    public: // ============================= FAST EIGENVALUE COMPUTATION METHODS

        double minimum_augmented_eigenvalue(
                const double *new_overlap_column,
                const double *new_hamiltonian_column);

    }; // class RealVariationalSolver

} // namespace zsvm

#endif // ZSVM_REAL_VARIATIONAL_SOLVER_HPP_INCLUDED

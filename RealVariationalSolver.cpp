#include "RealVariationalSolver.hpp"

// C++ standard library headers
#include <cmath> // for std::sqrt, std::pow, std::abs
#include <limits> // for std::numeric_limits<double>

// Eigen linear algebra library headers
#include <Eigen/Eigenvalues> // for Eigen::GeneralizedSelfAdjointEigenSolver


zsvm::RealVariationalSolver::RealVariationalSolver()
        : size(0), clean(false) {}


void zsvm::RealVariationalSolver::set_basis_size_conservative(
        std::size_t basis_size) {
    size = basis_size;
    overlap_matrix.conservativeResize(basis_size, basis_size);
    hamiltonian_matrix.conservativeResize(basis_size, basis_size);
    clean = false;
}


void zsvm::RealVariationalSolver::set_basis_size_destructive(
        std::size_t basis_size) {
    size = basis_size;
    overlap_matrix.resize(basis_size, basis_size);
    hamiltonian_matrix.resize(basis_size, basis_size);
    clean = false;
}


void zsvm::RealVariationalSolver::set_overlap_matrix_element(
        std::size_t i, std::size_t j, double elem) {
    overlap_matrix(i, j) = elem;
    clean = false;
}


void zsvm::RealVariationalSolver::set_hamiltonian_matrix_element(
        std::size_t i, std::size_t j, double elem) {
    hamiltonian_matrix(i, j) = elem;
    clean = false;
}


void zsvm::RealVariationalSolver::compute_eigenstates() {
    if (!clean) {
        Eigen::GeneralizedSelfAdjointEigenSolver eigen_solver(
                hamiltonian_matrix, overlap_matrix,
                Eigen::ComputeEigenvectors | Eigen::Ax_lBx);
        eigenvalues = eigen_solver.eigenvalues();
        eigenvectors = eigen_solver.eigenvectors();
        clean = true;
    }
}


double zsvm::RealVariationalSolver::get_eigenvalue(std::size_t i) {
    if (!clean) { compute_eigenstates(); }
    return eigenvalues(i);
}


double zsvm::RealVariationalSolver::secular_objective_function(
        double x, double t, const Eigen::ArrayXd &beta) const {
    return t - x - (beta / (eigenvalues.array() - x)).sum();
}


void zsvm::RealVariationalSolver::find_lower_bracketing_interval(
        double &lower_bound, double &upper_bound,
        double strict_upper_bound,
        double t, const Eigen::ArrayXd &beta) const {
    int k = 0;
    lower_bound = strict_upper_bound - std::pow(2.0, -k);
    upper_bound = strict_upper_bound - std::pow(2.0, -k - 1);
    double lower_objective = secular_objective_function(
            lower_bound, t, beta);
    double upper_objective = secular_objective_function(
            upper_bound, t, beta);
    while (true) {
        if (lower_objective > 0.0 && upper_objective > 0.0) {
            ++k;
            lower_bound = upper_bound;
            lower_objective = upper_objective;
            upper_bound = strict_upper_bound - std::pow(2.0, -k - 1);
            upper_objective = secular_objective_function(
                    upper_bound, t, beta);
        } else if (lower_objective < 0.0 && upper_objective < 0.0) {
            --k;
            upper_bound = lower_bound;
            upper_objective = lower_objective;
            lower_bound = strict_upper_bound - std::pow(2.0, -k);
            lower_objective = secular_objective_function(
                    lower_bound, t, beta);
        } else {
            break;
        }
    }
}


double zsvm::RealVariationalSolver::solve_secular_equation_bisection(
        double lower_bound, double upper_bound,
        double t, const Eigen::ArrayXd &beta) const {
    constexpr double TOLERANCE = 1.0e-14;
    while (true) {
        const double midpoint = 0.5 * (lower_bound + upper_bound);
        const double relative_difference = std::abs(
                (upper_bound - lower_bound) / midpoint);
        if (relative_difference < TOLERANCE) { return midpoint; }
        const double midpoint_objective = secular_objective_function(
                midpoint, t, beta);
        if (midpoint_objective == 0.0) {
            return midpoint;
        } else if (midpoint_objective > 0.0) {
            lower_bound = midpoint;
        } else if (midpoint_objective < 0.0) {
            upper_bound = midpoint;
        } else {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }
}


double zsvm::RealVariationalSolver::minimum_augmented_eigenvalue(
        const double *new_overlap_column,
        const double *new_hamiltonian_column) {
    if (!clean) { compute_eigenstates(); }
    Eigen::VectorXd overlap_vector(size);
    for (std::size_t i = 0; i < size; ++i) {
        overlap_vector(i) = new_overlap_column[i];
    }
    Eigen::VectorXd hamiltonian_vector(size);
    for (std::size_t i = 0; i < size; ++i) {
        hamiltonian_vector(i) = new_hamiltonian_column[i];
    }
    const Eigen::VectorXd alpha =
            eigenvectors.transpose() * overlap_vector;
    const double norm_factor = 1.0 / std::sqrt(
            new_overlap_column[size] - alpha.squaredNorm());
    const Eigen::VectorXd phi = -norm_factor * (eigenvectors * alpha);
    const Eigen::VectorXd psi = hamiltonian_matrix * phi +
                                norm_factor * hamiltonian_vector;
    const Eigen::ArrayXd beta =
            (eigenvectors.transpose() * psi).array().square();
    const double t = phi.dot(psi) + norm_factor * (
            hamiltonian_vector.dot(phi) +
            norm_factor * new_hamiltonian_column[size]);
    double lower_bound, upper_bound;
    find_lower_bracketing_interval(lower_bound, upper_bound,
                                   eigenvalues[0], t, beta);
    return solve_secular_equation_bisection(
            lower_bound, upper_bound, t, beta);
}

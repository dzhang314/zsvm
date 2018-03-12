#ifndef DZNL_AMOEBA_OPTIMIZER_HPP_INCLUDED
#define DZNL_AMOEBA_OPTIMIZER_HPP_INCLUDED

// C++ standard library headers
#include <cstddef> // for std::size_t
#include <functional> // for std::function
#include <utility> // for std::pair
#include <vector> // for std::vector

namespace dznl {

    class AmoebaOptimizer {

    private: // =============================================== MEMBER VARIABLES

        const std::size_t n;
        std::function<double(const double *)> f;
        std::vector<std::pair<double, std::vector<double>>> x;

    public: // ===================================================== CONSTRUCTOR

        AmoebaOptimizer(
                const double *initial_point,
                std::size_t num_dimensions,
                double initial_step_size,
                std::function<double(const double *)> objective_function);

    private: // ======================================================= MUTATORS

        void sort();

    public: // ======================================================= ACCESSORS

        double current_minimum(double *current_point);

    private: // ======================================= GEOMETRIC HELPER METHODS

        std::vector<double> compute_centroid() const;

        std::vector<double> compute_reflection(
                const std::vector<double> &p,
                const std::vector<double> &q) const;

        std::vector<double> compute_midpoint(
                const std::vector<double> &p,
                const std::vector<double> &q) const;

    public: // ============================================ OPTIMIZATION METHODS

        double step();

    }; // class AmoebaOptimizer

} // namespace dznl

#endif // DZNL_AMOEBA_OPTIMIZER_HPP_INCLUDED

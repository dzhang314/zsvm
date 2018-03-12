#include "AmoebaOptimizer.hpp"

// C++ standard library headers
#include <algorithm> // for std::sort


dznl::AmoebaOptimizer::AmoebaOptimizer(
        const double *initial_point,
        std::size_t num_dimensions,
        double initial_step_size,
        std::function<double(const double *)> objective_function)
        : n(num_dimensions), f(std::move(objective_function)) {
    std::vector<double> y(n);
    for (std::size_t i = 0; i < n + 1; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            y[j] = (i == j) ? initial_point[j] + initial_step_size
                            : initial_point[j];
        }
        x.emplace_back(f(y.data()), y);
    }
}


void dznl::AmoebaOptimizer::sort() {
    std::sort(x.begin(), x.end(), [](
            const std::pair<double, std::vector<double>> &p,
            const std::pair<double, std::vector<double>> &q) {
        return p.first < q.first;
    });
}


double dznl::AmoebaOptimizer::current_minimum(double *current_point) {
    sort();
    for (std::size_t i = 0; i < n; ++i) { current_point[i] = x[0].second[i]; }
    return x[0].first;
}


std::vector<double> dznl::AmoebaOptimizer::compute_centroid() const {
    std::vector<double> centroid(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) { centroid[j] += x[i].second[j]; }
    }
    for (std::size_t i = 0; i < n; ++i) { centroid[i] /= n; }
    return centroid;
}


std::vector<double> dznl::AmoebaOptimizer::compute_reflection(
        const std::vector<double> &p,
        const std::vector<double> &q) const {
    std::vector<double> reflected_point(n);
    for (std::size_t i = 0; i < n; ++i) {
        reflected_point[i] = 2 * q[i] - p[i];
    }
    return reflected_point;
}


std::vector<double> dznl::AmoebaOptimizer::compute_midpoint(
        const std::vector<double> &p,
        const std::vector<double> &q) const {
    std::vector<double> midpoint(n);
    for (std::size_t i = 0; i < n; ++i) {
        midpoint[i] = (p[i] + q[i]) / 2;
    }
    return midpoint;
}


double dznl::AmoebaOptimizer::step() {
    sort();
    const double &best_value = x[0].first;
    const std::vector<double> &best_point = x[0].second;
    const double &second_worst_value = x[n - 1].first;
    double &worst_value = x[n].first;
    std::vector<double> &worst_point = x[n].second;
    const std::vector<double> centroid = compute_centroid();
    const std::vector<double> reflected_point =
            compute_reflection(worst_point, centroid);
    const double reflected_value = f(reflected_point.data());
    if (best_value <= reflected_value &&
        reflected_value < second_worst_value) {
        worst_value = reflected_value;
        worst_point = reflected_point;
        return best_value;
    } else if (reflected_value < best_value) {
        const std::vector<double> expanded_point =
                compute_reflection(centroid, reflected_point);
        const double expanded_value = f(expanded_point.data());
        if (expanded_value < reflected_value) {
            worst_value = expanded_value;
            worst_point = expanded_point;
            return expanded_value;
        } else {
            worst_value = reflected_value;
            worst_point = reflected_point;
            return reflected_value;
        }
    } else if (reflected_value >= second_worst_value) {
        const std::vector<double> contracted_point =
                compute_midpoint(worst_point, centroid);
        const double contracted_value = f(contracted_point.data());
        if (contracted_value < worst_value) {
            worst_value = contracted_value;
            worst_point = contracted_point;
            return (contracted_value < best_value)
                   ? contracted_value : best_value;
        }
    }
    double best_shrunk_value = best_value;
    for (std::size_t i = 1; i < n + 1; ++i) {
        double &current_value = x[i].first;
        std::vector<double> &current_point = x[i].second;
        for (std::size_t j = 0; j < n; ++j) {
            current_point[j] = (best_point[j] + current_point[j]) / 2;
        }
        current_value = f(current_point.data());
        if (current_value < best_shrunk_value) {
            best_shrunk_value = current_value;
        }
    }
    return best_shrunk_value;
}

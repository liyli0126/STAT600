// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

#include "ReactionDiffusionCpp.h"
#include <stdexcept>
#include <algorithm>

ReactionDiffusionSolver1D::ReactionDiffusionSolver1D(double D_val, double rho_val,
                                                     double length, int grid_points,
                                                     double time_step)
  : D(D_val), rho(rho_val), domain_length(length),
    grid_size(grid_points), dt(time_step) {
  
  if (grid_size < 3) {
    throw std::invalid_argument("Grid size must be at least 3.");
  }
  
  dx = domain_length / (grid_size - 1);
  double cfl = D * dt / (dx * dx);
  if (cfl > 0.5) {
    Rcpp::warning("CFL condition may be violated (CFL = " + std::to_string(cfl) + ").");
  }
}

void ReactionDiffusionSolver1D::solve_tridiagonal(
    const std::vector<double>& a,
    const std::vector<double>& b,
    const std::vector<double>& c,
    const std::vector<double>& d,
    std::vector<double>& x) const {
  
  int n = b.size();
  std::vector<double> c_prime(n-1);
  std::vector<double> d_prime(n);
  
  c_prime[0] = c[0] / b[0];
  d_prime[0] = d[0] / b[0];
  
  for (int i = 1; i < n-1; ++i) {
    double denom = b[i] - a[i-1] * c_prime[i-1];
    c_prime[i] = c[i] / denom;
    d_prime[i] = (d[i] - a[i-1] * d_prime[i-1]) / denom;
  }
  d_prime[n-1] = (d[n-1] - a[n-2] * d_prime[n-2]) /
    (b[n-1] - a[n-2] * c_prime[n-2]);
  
  x[n-1] = d_prime[n-1];
  for (int i = n-2; i >= 0; --i) {
    x[i] = d_prime[i] - c_prime[i] * x[i+1];
  }
}

std::vector<double> ReactionDiffusionSolver1D::solve_crank_nicolson(
    const std::vector<double>& initial_condition, double end_time) const {
  
  if (initial_condition.size() != (size_t)grid_size) {
    throw std::invalid_argument("Initial condition size mismatch.");
  }
  
  std::vector<double> u_current = initial_condition;
  std::vector<double> u_next(grid_size);
  int num_steps = static_cast<int>(end_time / dt + 0.5);
  
  double alpha = D * dt / (2.0 * dx * dx);
  double beta = dt * rho / 2.0;
  
  std::vector<double> a(grid_size - 1, -alpha);
  std::vector<double> b(grid_size, 1.0 + 2.0 * alpha);
  std::vector<double> c(grid_size - 1, -alpha);
  b[0] = 1.0 + alpha;
  b[grid_size - 1] = 1.0 + alpha;
  
  for (int step = 0; step < num_steps; ++step) {
    std::vector<double> rhs(grid_size);
    for (int i = 1; i < grid_size - 1; ++i) {
      double diff_exp = alpha * (u_current[i-1] + u_current[i+1] - 2.0 * u_current[i]);
      double react_exp = beta * u_current[i] * (1.0 - u_current[i]);
      rhs[i] = u_current[i] + diff_exp + 2.0 * react_exp;
    }
    rhs[0] = u_current[0] + alpha * (u_current[1] - u_current[0]) +
      beta * u_current[0] * (1.0 - u_current[0]);
    rhs[grid_size-1] = u_current[grid_size-1] +
      alpha * (u_current[grid_size-2] - u_current[grid_size-1]) +
      beta * u_current[grid_size-1] * (1.0 - u_current[grid_size-1]);
    
    solve_tridiagonal(a, b, c, rhs, u_next);
    
    for (int i = 0; i < grid_size; ++i) {
      u_next[i] = std::max(0.0, std::min(1.0, u_next[i]));
    }
    u_current = u_next;
  }
  return u_current;
}
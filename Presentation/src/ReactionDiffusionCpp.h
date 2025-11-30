#ifndef REACTION_DIFFUSION_H
#define REACTION_DIFFUSION_H

#include <vector>
#include <cmath>

class ReactionDiffusionSolver1D {
private:
  double D;
  double rho;
  double domain_length;
  int grid_size;
  double dx;
  double dt;
  
  void solve_tridiagonal(const std::vector<double>& a,
                         const std::vector<double>& b,
                         const std::vector<double>& c,
                         const std::vector<double>& d,
                         std::vector<double>& x) const;
  
public:
  ReactionDiffusionSolver1D(double D_val, double rho_val,
                            double length, int grid_points, double time_step);
  
  std::vector<double> solve_crank_nicolson(const std::vector<double>& initial_condition,
                                           double end_time) const;
  
  double get_D() const { return D; }
  double get_rho() const { return rho; }
  double get_dx() const { return dx; }
  double get_dt() const { return dt; }
};

#endif

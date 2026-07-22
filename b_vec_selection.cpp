// wells_trap.cpp
#include <Rcpp.h>
#include <cmath>
#include <string>
#include <vector>

using namespace Rcpp;

inline double B_abs(double c, double x, double eps, double y, double sigma) {
  const double sigma2 = sigma * sigma;
  
  const double c2 = c * c;
  const double c3 = c2 * c;
  const double c4 = c2 * c2;
  const double c5 = c4 * c;
  const double c6 = c3 * c3;
  
  const double x2 = x * x;
  const double x3 = x2 * x;
  const double x4 = x2 * x2;
  
  const double denom = 1.0 / (24.0 * eps * eps * eps);
  
  const double top =
    (-6.0 * c3 + 3.0 * y) * x4 +
    x3 * (-36.0 * c4 - 24.0 * sigma2 * eps + 36.0 * c2 - 8.0) +
    x2 * (54.0 * (c3 - y) * (c2 - 1.0 / 3.0)) +
    x  * (-24.0 * c6 + 54.0 * c2 * eps * sigma2 + 36.0 * c3 * y - 12.0 * sigma2 * eps - 12.0 * y * y) +
    6.0 * c5 - 9.0 * c4 * y + c3 * (-18.0 * sigma2 * eps - 4.0) +
    6.0 * c2 * y + 12.0 * eps * sigma2 * y;
  
  return std::fabs(top * denom);
}

inline double Var_abs(double c, double x, double eps, double y, double sigma) {
  const double sigma2 = sigma * sigma;
  const double sigma4 = sigma2 * sigma2;
  
  const double c2 = c * c;
  const double c3 = c2 * c;
  
  const double x2 = x * x;
  const double x3 = x2 * x;
  
  const double Nx  = (-x3 + 3.0 * c2 * x - 2.0 * c3) / eps;
  const double dNx = (-3.0 * x2 + 3.0 * c2) / eps;
  const double d2Nx = -6.0 * x / eps;
  const double d3Nx = -6.0 / eps;
  
  const double A = (1.0 - 3.0 * c2) / eps;
  const double b = c - (c - c3 - y) / (1.0 - 3.0 * c2);
  
  const double Var =
    (-1.0 / 12.0) * sigma2 * d2Nx * Nx -
    (1.0 / 3.0)  * sigma2 * d2Nx * A * (x - b) +
    (1.0 / 6.0)  * sigma2 * dNx * dNx +
    (1.0 / 3.0)  * sigma2 * A * dNx +
    (-1.0 / 6.0)  * sigma4 * d3Nx;
  
  return std::fabs(Var);
}

inline double B_scalar(double c, double x, double eps, double y, double sigma, bool use_abs = true) {
  const double sigma2 = sigma * sigma;
  const double c2 = c * c;
  const double c3 = c2 * c;
  const double c4 = c2 * c2;
  const double c5 = c4 * c;
  const double c6 = c3 * c3;
  
  const double x2 = x * x;
  const double x3 = x2 * x;
  const double x4 = x2 * x2;
  
  const double denom = 1.0 / (24.0 * eps * eps * eps);
  
  const double top =
    (-6.0 * c3 + 3.0 * y) * x4 +
    x3 * (-36.0 * c4 - 24.0 * sigma2 * eps + 36.0 * c2 - 8.0) +
    x2 * (54.0 * (c3 - y) * (c2 - 1.0 / 3.0)) +
    x  * (-24.0 * c6 + 54.0 * c2 * eps * sigma2 + 36.0 * c3 * y - 12.0 * sigma2 * eps - 12.0 * y * y) +
    6.0 * c5 - 9.0 * c4 * y + c3 * (-18.0 * sigma2 * eps - 4.0) +
    6.0 * c2 * y + 12.0 * eps * sigma2 * y;
  
  const double val = top * denom;
  return use_abs ? std::fabs(val) : val;
}

inline double Var_scalar(double c, double x, double eps, double y, double sigma, bool use_abs = true) {
  const double sigma2 = sigma * sigma;
  const double sigma4 = sigma2 * sigma2;
  
  const double c2 = c * c;
  const double c3 = c2 * c;
  
  const double x2 = x * x;
  const double x3 = x2 * x;
  
  const double Nx  = (-x3 + 3.0 * c2 * x - 2.0 * c3) / eps;
  const double dNx = (-3.0 * x2 + 3.0 * c2) / eps;
  const double d2Nx = -6.0 * x / eps;
  const double d3Nx = -6.0 / eps;
  
  const double A = (1.0 - 3.0 * c2) / eps;
  const double b = c - (c - c3 - y) / (1.0 - 3.0 * c2);
  
  const double Var =
    (-1.0 / 12.0) * sigma2 * d2Nx * Nx -
    (1.0 / 3.0)  * sigma2 * d2Nx * A * (x - b) +
    (1.0 / 6.0)  * sigma2 * dNx * dNx +
    (1.0 / 3.0)  * sigma2 * A * dNx +
    (-1.0 / 6.0)  * sigma4 * d3Nx;
  
  return use_abs ? std::fabs(Var) : Var;
}

// ---------- exact grid-search replacements for the "steps" methods ---------

// [[Rcpp::export]]
NumericVector minimax_steps_cpp(NumericVector x_vals,
                                NumericVector c_grid,
                                double eps, double y, double sigma) {
  const int nx = x_vals.size();
  const int nc = c_grid.size();
  
  NumericVector out(nx);
  
  for (int i = 0; i < nx; ++i) {
    const double x = x_vals[i];
    double best = R_PosInf;
    double best_c = c_grid[0];
    
    for (int j = 0; j < nc; ++j) {
      const double c = c_grid[j];
      const double B = B_scalar(c, x, eps, y, sigma, true);
      const double V = Var_scalar(c, x, eps, y, sigma, true);
      
      const double obj = std::max(B, std::sqrt(V));
      
      if (obj < best) {
        best = obj;
        best_c = c;
      }
    }
    
    out[i] = best_c;
  }
  
  return out;
}

// [[Rcpp::export]]
NumericVector argminBargminV_steps_cpp(NumericVector x_vals,
                                       NumericVector c_grid,
                                       double eps, double y, double sigma,
                                       double tol = 0.01) {
  const int nx = x_vals.size();
  const int nc = c_grid.size();
  
  NumericVector out(nx);
  
  for (int i = 0; i < nx; ++i) {
    const double x = x_vals[i];
    
    // pass 1: minimum of V
    double Vmin = R_PosInf;
    for (int j = 0; j < nc; ++j) {
      const double V = Var_scalar(c_grid[j], x, eps, y, sigma, true);
      if (V < Vmin) Vmin = V;
    }
    
    // pass 2: among V <= Vmin + tol, choose smallest B
    double bestB = R_PosInf;
    double best_c = c_grid[0];
    bool found = false;
    
    for (int j = 0; j < nc; ++j) {
      const double c = c_grid[j];
      const double V = Var_scalar(c, x, eps, y, sigma, true);
      if (V <= Vmin + tol) {
        const double B = B_scalar(c, x, eps, y, sigma, true);
        if (!found || B < bestB) {
          bestB = B;
          best_c = c;
          found = true;
        }
      }
    }
    
    out[i] = best_c;
  }
  
  return out;
}

// [[Rcpp::export]]
NumericVector argminVargminB_steps_cpp(NumericVector x_vals,
                                       NumericVector c_grid,
                                       double eps, double y, double sigma,
                                       double tol = 0.01) {
  const int nx = x_vals.size();
  const int nc = c_grid.size();
  
  NumericVector out(nx);
  
  for (int i = 0; i < nx; ++i) {
    const double x = x_vals[i];
    
    // pass 1: minimum of B
    double Bmin = R_PosInf;
    for (int j = 0; j < nc; ++j) {
      const double B = B_scalar(c_grid[j], x, eps, y, sigma, true);
      if (B < Bmin) Bmin = B;
    }
    
    // pass 2: among B <= Bmin + tol, choose smallest V
    double bestV = R_PosInf;
    double best_c = c_grid[0];
    bool found = false;
    
    for (int j = 0; j < nc; ++j) {
      const double c = c_grid[j];
      const double B = B_scalar(c, x, eps, y, sigma, true);
      if (B <= Bmin + tol) {
        const double V = Var_scalar(c, x, eps, y, sigma, true);
        if (!found || V < bestV) {
          bestV = V;
          best_c = c;
          found = true;
        }
      }
    }
    
    out[i] = best_c;
  }
  
  return out;
}

// ---------- exact selectors for the "wells" branches ----------
// These assume you already computed the vectors in R exactly as before.

// [[Rcpp::export]]
double minimax_from_vectors_cpp(NumericVector c_grid,
                                NumericVector Bvals,
                                NumericVector Vvals) {
  const int n = c_grid.size();
  if (Bvals.size() != n || Vvals.size() != n) {
    stop("c_grid, Bvals, and Vvals must have the same length.");
  }
  
  double best = R_PosInf;
  double best_c = c_grid[0];
  
  for (int j = 0; j < n; ++j) {
    const double obj = (Bvals[j] > Vvals[j]) ? Bvals[j] : Vvals[j];
    if (obj < best) {
      best = obj;
      best_c = c_grid[j];
    }
  }
  return best_c;
}

// [[Rcpp::export]]
double argminBargminV_from_vectors_cpp(NumericVector c_grid,
                                       NumericVector Bvals,
                                       NumericVector Vvals,
                                       double tol = 0.1) {
  const int n = c_grid.size();
  if (Bvals.size() != n || Vvals.size() != n) {
    stop("c_grid, Bvals, and Vvals must have the same length.");
  }
  
  double Vmin = R_PosInf;
  for (int j = 0; j < n; ++j) {
    if (Vvals[j] < Vmin) Vmin = Vvals[j];
  }
  
  double bestB = R_PosInf;
  double best_c = c_grid[0];
  bool found = false;
  
  for (int j = 0; j < n; ++j) {
    if (Vvals[j] <= Vmin + tol) {
      if (!found || Bvals[j] < bestB) {
        bestB = Bvals[j];
        best_c = c_grid[j];
        found = true;
      }
    }
  }
  
  return best_c;
}

// [[Rcpp::export]]
double argminVargminB_from_vectors_cpp(NumericVector c_grid,
                                       NumericVector Bvals,
                                       NumericVector Vvals,
                                       double tol = 0.1) {
  const int n = c_grid.size();
  if (Bvals.size() != n || Vvals.size() != n) {
    stop("c_grid, Bvals, and Vvals must have the same length.");
  }
  
  double Bmin = R_PosInf;
  for (int j = 0; j < n; ++j) {
    if (Bvals[j] < Bmin) Bmin = Bvals[j];
  }
  
  double bestV = R_PosInf;
  double best_c = c_grid[0];
  bool found = false;
  
  for (int j = 0; j < n; ++j) {
    if (Bvals[j] <= Bmin + tol) {
      if (!found || Vvals[j] < bestV) {
        bestV = Vvals[j];
        best_c = c_grid[j];
        found = true;
      }
    }
  }
  
  return best_c;
}

inline double inv_density(double x, double eps, double y, double sigma) {
  return std::exp(2.0 * (0.5 * x * x - 0.25 * x * x * x * x - y * x) / (sigma * sigma * eps));
}

inline NumericVector uniform_grid(double a, double b, int n) {
  if (n < 2) stop("n_grid must be at least 2.");
  NumericVector g(n);
  const double h = (b - a) / (n - 1.0);
  for (int i = 0; i < n; ++i) g[i] = a + h * i;
  return g;
}

inline double trapz_uniform(const std::vector<double>& f, double h) {
  const int n = (int)f.size();
  double s = 0.5 * f.front() + 0.5 * f.back();
  for (int i = 1; i < n - 1; ++i) s += f[i];
  return s * h;
}

struct WellData {
  NumericVector x;
  std::vector<double> w; // trapezoid weight * inv_density/Z
  double p_inv;          // 1 / P_A
};

inline void get_well_interval(const std::string& split_mode,
                              const std::string& which_well,
                              double xmin, double xmax,
                              double uns_fix,
                              double& a, double& b) {
  if (split_mode == "uns_fix") {
    if (which_well == "left") {
      a = xmin;  b = uns_fix;
    } else if (which_well == "right") {
      a = uns_fix;  b = xmax;
    } else {
      stop("For split_mode='uns_fix', which_well must be 'left' or 'right'.");
    }
  } else if (split_mode == "thirds") {
    const double t = 1.0 / std::sqrt(3.0);
    if (which_well == "left") {
      a = xmin;  b = -t;
    } else if (which_well == "middle") {
      a = -t;    b = t;
    } else if (which_well == "right") {
      a = t;     b = xmax;
    } else {
      stop("For split_mode='thirds', which_well must be 'left', 'middle', or 'right'.");
    }
  } else {
    stop("Unknown split_mode. Use 'uns_fix' or 'thirds'.");
  }
  
  if (!(a < b)) stop("Invalid interval: left endpoint must be smaller than right endpoint.");
}

inline WellData make_well_data(double eps, double y, double sigma,
                               const std::string& split_mode,
                               const std::string& which_well,
                               double uns_fix,
                               double x_bound,
                               int n_grid) {
  const double xmin = -x_bound;
  const double xmax =  x_bound;
  
  if (split_mode == "uns_fix") {
    if (!(uns_fix > xmin && uns_fix < xmax)) {
      stop("uns_fix must lie inside (-x_bound, x_bound). Increase x_bound if needed.");
    }
  }
  
  // Full-grid normalization Z
  NumericVector x_full = uniform_grid(xmin, xmax, n_grid);
  const double h_full = (xmax - xmin) / (n_grid - 1.0);
  
  std::vector<double> inv_full(n_grid);
  for (int i = 0; i < n_grid; ++i) {
    inv_full[i] = inv_density(x_full[i], eps, y, sigma);
  }
  const double Z = trapz_uniform(inv_full, h_full);
  
  // Pick interval based on split_mode
  double a, b;
  get_well_interval(split_mode, which_well, xmin, xmax, uns_fix, a, b);
  
  NumericVector x_well = uniform_grid(a, b, n_grid);
  const double h_well = (b - a) / (n_grid - 1.0);
  
  std::vector<double> pvec(n_grid);
  for (int i = 0; i < n_grid; ++i) {
    pvec[i] = inv_density(x_well[i], eps, y, sigma) / Z;
  }
  
  const double P_A = trapz_uniform(pvec, h_well);
  
  std::vector<double> w(n_grid);
  for (int i = 0; i < n_grid; ++i) {
    const double tw = (i == 0 || i == n_grid - 1) ? 0.5 * h_well : h_well;
    w[i] = pvec[i] * tw;
  }
  
  WellData out;
  out.x = x_well;
  out.w = w;
  out.p_inv = 1.0 / P_A;
  return out;
}

inline void compute_BV_for_well(const NumericVector& c_grid,
                                const WellData& W,
                                double eps, double y, double sigma,
                                std::vector<double>& Bvals,
                                std::vector<double>& Vvals) {
  const int nc = c_grid.size();
  const int nx = W.x.size();
  
  Bvals.resize(nc);
  Vvals.resize(nc);
  
  for (int ic = 0; ic < nc; ++ic) {
    const double c = c_grid[ic];
    double sumB = 0.0;
    double sumV = 0.0;
    
    for (int ix = 0; ix < nx; ++ix) {
      const double x = W.x[ix];
      const double wt = W.w[ix];
      sumB += B_abs(c, x, eps, y, sigma) * wt;
      sumV += Var_abs(c, x, eps, y, sigma) * wt;
    }
    
    Bvals[ic] = std::fabs(sumB * W.p_inv);
    Vvals[ic] = std::fabs(sumV * W.p_inv);
  }
}

inline double select_minimax(const NumericVector& c_grid,
                             const std::vector<double>& B,
                             const std::vector<double>& V) {
  double best = R_PosInf;
  double best_c = c_grid[0];
  
  for (int j = 0; j < c_grid.size(); ++j) {
    const double obj = std::max(B[j], std::sqrt(V[j]));
    if (obj < best) {
      best = obj;
      best_c = c_grid[j];
    }
  }
  
  return best_c;
}

inline double select_argminBargminV(const NumericVector& c_grid,
                                    const std::vector<double>& B,
                                    const std::vector<double>& V,
                                    double tol) {
  double Vmin = R_PosInf;
  for (double v : V) if (v < Vmin) Vmin = v;
  
  double bestB = R_PosInf;
  double best_c = c_grid[0];
  bool found = false;
  
  for (int j = 0; j < c_grid.size(); ++j) {
    if (V[j] <= Vmin + tol) {
      if (!found || B[j] < bestB) {
        bestB = B[j];
        best_c = c_grid[j];
        found = true;
      }
    }
  }
  return best_c;
}

inline double select_argminVargminB(const NumericVector& c_grid,
                                    const std::vector<double>& B,
                                    const std::vector<double>& V,
                                    double tol) {
  double Bmin = R_PosInf;
  for (double b : B) if (b < Bmin) Bmin = b;
  
  double bestV = R_PosInf;
  double best_c = c_grid[0];
  bool found = false;
  
  for (int j = 0; j < c_grid.size(); ++j) {
    if (B[j] <= Bmin + tol) {
      if (!found || V[j] < bestV) {
        bestV = V[j];
        best_c = c_grid[j];
        found = true;
      }
    }
  }
  return best_c;
}

// [[Rcpp::export]]
List wells_eval_trap_cpp(NumericVector c_grid,
                         double eps, double y, double sigma,
                         std::string split_mode,
                         std::string which_well,
                         double uns_fix,
                         double x_bound = 3.0,
                         int n_grid = 501) {
  WellData W = make_well_data(eps, y, sigma,
                              split_mode, which_well,
                              uns_fix, x_bound, n_grid);
  
  std::vector<double> Bvals, Vvals;
  compute_BV_for_well(c_grid, W, eps, y, sigma, Bvals, Vvals);
  
  return List::create(
    Named("B") = wrap(Bvals),
    Named("V") = wrap(Vvals)
  );
}

// [[Rcpp::export]]
NumericVector wells_bvec_trap_cpp(std::string method,
                                  NumericVector c_grid,
                                  double eps, double y, double sigma,
                                  std::string split_mode,
                                  double uns_fix,
                                  double tol = 0.1,
                                  double x_bound = 3.0,
                                  int n_grid = 501) {
  if (method == "minimax wells") {
    if (split_mode == "uns_fix") {
      WellData WL = make_well_data(eps, y, sigma, split_mode, "left",  uns_fix, x_bound, n_grid);
      WellData WR = make_well_data(eps, y, sigma, split_mode, "right", uns_fix, x_bound, n_grid);
      
      std::vector<double> BL, VL, BR, VR;
      compute_BV_for_well(c_grid, WL, eps, y, sigma, BL, VL);
      compute_BV_for_well(c_grid, WR, eps, y, sigma, BR, VR);
      
      NumericVector out(2);
      out[0] = select_minimax(c_grid, BL, VL);
      out[1] = select_minimax(c_grid, BR, VR);
      return out;
      
    } else if (split_mode == "thirds") {
      WellData WL = make_well_data(eps, y, sigma, split_mode, "left",   uns_fix, x_bound, n_grid);
      WellData WM = make_well_data(eps, y, sigma, split_mode, "middle", uns_fix, x_bound, n_grid);
      WellData WR = make_well_data(eps, y, sigma, split_mode, "right",   uns_fix, x_bound, n_grid);
      
      std::vector<double> BL, VL, BM, VM, BR, VR;
      compute_BV_for_well(c_grid, WL, eps, y, sigma, BL, VL);
      compute_BV_for_well(c_grid, WM, eps, y, sigma, BM, VM);
      compute_BV_for_well(c_grid, WR, eps, y, sigma, BR, VR);
      
      NumericVector out(3);
      out[0] = select_minimax(c_grid, BL, VL);
      out[1] = select_minimax(c_grid, BM, VM);
      out[2] = select_minimax(c_grid, BR, VR);
      return out;
    }
    
  } else if (method == "argminBargminV wells") {
    if (split_mode == "uns_fix") {
      WellData WL = make_well_data(eps, y, sigma, split_mode, "left",  uns_fix, x_bound, n_grid);
      WellData WR = make_well_data(eps, y, sigma, split_mode, "right", uns_fix, x_bound, n_grid);
      
      std::vector<double> BL, VL, BR, VR;
      compute_BV_for_well(c_grid, WL, eps, y, sigma, BL, VL);
      compute_BV_for_well(c_grid, WR, eps, y, sigma, BR, VR);
      
      NumericVector out(2);
      out[0] = select_argminBargminV(c_grid, BL, VL, tol);
      out[1] = select_argminBargminV(c_grid, BR, VR, tol);
      return out;
      
    } else if (split_mode == "thirds") {
      WellData WL = make_well_data(eps, y, sigma, split_mode, "left",   uns_fix, x_bound, n_grid);
      WellData WM = make_well_data(eps, y, sigma, split_mode, "middle", uns_fix, x_bound, n_grid);
      WellData WR = make_well_data(eps, y, sigma, split_mode, "right",   uns_fix, x_bound, n_grid);
      
      std::vector<double> BL, VL, BM, VM, BR, VR;
      compute_BV_for_well(c_grid, WL, eps, y, sigma, BL, VL);
      compute_BV_for_well(c_grid, WM, eps, y, sigma, BM, VM);
      compute_BV_for_well(c_grid, WR, eps, y, sigma, BR, VR);
      
      NumericVector out(3);
      out[0] = select_argminBargminV(c_grid, BL, VL, tol);
      out[1] = select_argminBargminV(c_grid, BM, VM, tol);
      out[2] = select_argminBargminV(c_grid, BR, VR, tol);
      return out;
    }
    
  } else if (method == "argminVargminB wells") {
    if (split_mode == "uns_fix") {
      WellData WL = make_well_data(eps, y, sigma, split_mode, "left",  uns_fix, x_bound, n_grid);
      WellData WR = make_well_data(eps, y, sigma, split_mode, "right", uns_fix, x_bound, n_grid);
      
      std::vector<double> BL, VL, BR, VR;
      compute_BV_for_well(c_grid, WL, eps, y, sigma, BL, VL);
      compute_BV_for_well(c_grid, WR, eps, y, sigma, BR, VR);
      
      NumericVector out(2);
      out[0] = select_argminVargminB(c_grid, BL, VL, tol);
      out[1] = select_argminVargminB(c_grid, BR, VR, tol);
      return out;
      
    } else if (split_mode == "thirds") {
      WellData WL = make_well_data(eps, y, sigma, split_mode, "left",   uns_fix, x_bound, n_grid);
      WellData WM = make_well_data(eps, y, sigma, split_mode, "middle", uns_fix, x_bound, n_grid);
      WellData WR = make_well_data(eps, y, sigma, split_mode, "right",   uns_fix, x_bound, n_grid);
      
      std::vector<double> BL, VL, BM, VM, BR, VR;
      compute_BV_for_well(c_grid, WL, eps, y, sigma, BL, VL);
      compute_BV_for_well(c_grid, WM, eps, y, sigma, BM, VM);
      compute_BV_for_well(c_grid, WR, eps, y, sigma, BR, VR);
      
      NumericVector out(3);
      out[0] = select_argminVargminB(c_grid, BL, VL, tol);
      out[1] = select_argminVargminB(c_grid, BM, VM, tol);
      out[2] = select_argminVargminB(c_grid, BR, VR, tol);
      return out;
    }
  }
  
  stop("Unknown method or split_mode.");
  return NumericVector();
}
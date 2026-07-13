// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <limits>
using namespace Rcpp;
using arma::mat;
using arma::vec;

inline double cube_root(double x) {
  if (x == 0.0) return 0.0;
  return (x > 0.0 ? 1.0 : -1.0) * std::pow(std::abs(x), 1.0/3.0);
}
// [[Rcpp::export]]
NumericVector fixed_point_rcpp(NumericVector par){
  if (par.size() < 4) stop("par length < 4 in fixed_point_rcpp");
  double eps   = par[0];
  double alpha = par[1];
  double gamma = par[2];
  double beta  = par[3];
  
  double p = gamma - 1.0;
  double q = beta - alpha;
  double D = q*q / 4.0 + p*p*p / 27.0;
  double u1 = (alpha - beta) / 2.0 + std::sqrt(std::max(0.0, D));
  double u2 = (alpha - beta) / 2.0 - std::sqrt(std::max(0.0, D));
  double xstar = cube_root(u1) + cube_root(u2);
  double ystar = gamma * xstar + beta;
  
  return NumericVector::create(xstar, ystar);
}

// A implementation
// [[Rcpp::export]]
arma::mat A_rcpp(NumericVector par, NumericVector center, const std::string &method) {
  if (par.size() < 4) stop("par length < 4 in A_rcpp");
  if (center.size() < 2) stop("center length < 2 in A_rcpp");
  double eps   = par[0];
  double alpha = par[1];
  double gamma = par[2];
  double beta  = par[3];
  
  arma::mat A_mat(2,2);
  if (method == "buckwar") {
    double c1 = 0.0;
    double c2 = -1.0/eps;
    double c3 = gamma;
    double c4 = -1.0;
    A_mat(0,0) = c1; A_mat(0,1) = c2;
    A_mat(1,0) = c3; A_mat(1,1) = c4;
  } else {
    double b1 = center[0];
    double c1 = 1.0/eps - 3.0/eps * b1*b1;
    double c2 = -1.0/eps;
    double c3 = gamma;
    double c4 = -1.0;
    A_mat(0,0) = c1; A_mat(0,1) = c2;
    A_mat(1,0) = c3; A_mat(1,1) = c4;
    
  } 
  return A_mat;
}

// Computes N(x,y) where (x,y) is the current state 
arma::vec ODE(const arma::vec &state, NumericVector par, NumericVector center, const std::string &method){
  if (par.size() < 4) stop("par length < 4 in ODE");
  double eps   = par[0];
  double alpha = par[1];
  double gamma = par[2];
  double beta  = par[3];
  
  double x = state(0);
  double y = state(1);
  (void) y;
  
  arma::vec out(2);
  
  if (method == "two knees" || method == "fix" || method == "piecewise taylor"|| method == "new_piecewise taylor" || method == "test" || method == "custom"|| method == "eigen" || method == "fast eigen" || method == "mixed eigen" || method=="knee" || method == "fenichel") {
    double x_star = center[0];
    double dx = ((-std::pow(x,3)) + 3.0 * x_star*x_star * x - 2.0 * std::pow(x_star,3)) / eps;
    double dy = 0.0;
    out(0) = dx; out(1) = dy;
    return out;
  } else if (method == "fast-slow" || method == "fast" || method == "piecewise" || method == "new_piecewise" ) {
    double b1 = center[0];
    double b2 = center[1];
    double dx = ((-std::pow(x,3)) + 3.0 * b1*b1 * x - 3.0 * std::pow(b1,3) + b1 + alpha - b2) / eps;
    double dy = gamma*b1 + beta - b2;
    out(0) = dx; out(1) = dy;
    return out;
  } else if (method == "buckwar") {
    double dx = ((-std::pow(x,3)) + x + alpha) / eps;
    double dy = beta;
    out(0) = dx; out(1) = dy;
    return out;
  } 
  stop("Unknown method in ODE");
  return out;
}

// fh: single RK4 step from t=0 to t=h
// [[Rcpp::export]]
NumericVector fh_rcpp(NumericVector x, double h, NumericVector par, NumericVector center, const std::string &method) {
  if (x.size() < 2) stop("x length < 2 in fh_rcpp");
  arma::vec state(2);
  state(0) = x[0];
  state(1) = x[1];
  
  arma::vec k1 = ODE(state, par, center, method);
  arma::vec s2 = state + 0.5 * h * k1;
  arma::vec k2 = ODE(s2, par, center, method);
  arma::vec s3 = state + 0.5 * h * k2;
  arma::vec k3 = ODE(s3, par, center, method);
  arma::vec s4 = state + h * k3;
  arma::vec k4 = ODE(s4, par, center, method);
  
  arma::vec next = state + (h / 6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
  
  return NumericVector::create(next(0), next(1));
}

// mu 
// [[Rcpp::export]]
NumericVector mu_rcpp(NumericVector x, double h, NumericVector par, NumericVector center, const std::string &method) {
  if (x.size() < 2) stop("x length < 2 in mu_rcpp");
  
  double eps   = par[0];
  double alpha = par[1];
  double gamma = par[2];
  double beta = par[3];
  
  arma::mat A_mat = A_rcpp(par, center, method);
  arma::mat M = arma::expmat(A_mat * h);

  arma::vec x_vec = as<arma::vec>(x);
  arma::vec b = as<arma::vec>(center);
  if (method == "test" || method == "custom"|| method == "piecewise taylor" || method == "new_piecewise taylor"|| method == "eigen" || method == "fast eigen" || method == "mixed eigen" || method == "knee" || method=="fenichel") {  

    arma::vec Fb(2);
    Fb(0) = (b(0) - std::pow(b(0),3) + alpha - b(1)) / eps;
    Fb(1) = gamma * b(0) + beta - b(1);
    
    arma::mat A_inv = arma::inv(A_mat);

    arma::vec b_tilde = b - A_inv * Fb;
    arma::vec res = M * x_vec + (arma::eye<arma::mat>(2,2) - M) * b_tilde;
    
    return NumericVector::create(res(0), res(1));
  }else if (method == "piecewise" || method == "new_piecewise"){
    arma::vec res = M * x_vec + (arma::eye<arma::mat>(2,2) - M) * b;
    return NumericVector::create(res(0), res(1));
  }
  else{
    arma::vec res = M * x_vec + (arma::eye<arma::mat>(2,2) - M) * b;
    return NumericVector::create(res(0), res(1));
  }
}

// Omega
// [[Rcpp::export]]
arma::mat Omega_rcpp(double h, NumericVector par, NumericVector center, const std::string &method) {
  if (par.size() < 6) stop("par length < 6 in Omega_rcpp");
  double sigma12 = std::pow(par[4], 2);
  double sigma22 = std::pow(par[5], 2);
  
  arma::mat M(4,4, arma::fill::zeros);
  arma::mat A22 = A_rcpp(par, center, method);
  
  M.submat(0,0,1,1) = A22;
  M(0,2) = sigma12; M(1,3) = sigma22;
  M.submat(2,2,3,3) = -A22.t();
  
  arma::mat exph = arma::expmat(M * h);
  
  arma::mat F1 = exph.submat(0,0,1,1);
  arma::mat G1 = exph.submat(0,2,1,3);
  
  arma::mat result = G1 * F1.t();
  return result;
}

// Equivalent to numDeriv in R
arma::mat numeric_jacobian_std(std::function<NumericVector(NumericVector)> f, NumericVector x) {
  const double eps = std::numeric_limits<double>::epsilon();
  double step_base = std::pow(eps, 1.0/3.0);
  
  arma::mat J(2,2, arma::fill::zeros);
  for (int j = 0; j < 2; ++j) {
    double xj = x[j];
    double h = step_base * std::max(1.0, std::abs(xj));
    NumericVector xp = clone(x), xm = clone(x);
    xp[j] = xj + h;
    xm[j] = xj - h;
    NumericVector fp = f(xp);
    NumericVector fm = f(xm);
    J(0,j) = (fp[0] - fm[0]) / (2.0 * h);
    J(1,j) = (fp[1] - fm[1]) / (2.0 * h);
  }
  return J;
}

// [[Rcpp::export]]
NumericVector fh_rhs_2d(double t, NumericVector y, NumericVector par) {
  // y[0] is state 1 (e.g., x)
  // y[1] is state 2 (e.g., z)
  double x = y[0];
  double z = y[1];
  
  double eps = par[0];
  double b   = par[1];
  bool const_term_in_N = (par[2] == 1.0); 
  
  // Define your two differential equations
  double dx_dt; 
  if (const_term_in_N) {
    dx_dt = (1.0 / eps) * (-x*x*x - b + 3.0*b*b*x + b - 3.0*b*b*b); 
  } else {
    dx_dt = (1.0 / eps) * (-x*x*x + 3.0*b*b*x - 2.0*b*b*b);
  }
  
  double dz_dt = 0; // Example second equation
  
  // Return a vector of size 2
  NumericVector out(2);
  out[0] = dx_dt;
  out[1] = dz_dt;
  return out;
}

// [[Rcpp::export]]
NumericVector fh_sundials_2d(double x_init,
                             double z_init,
                             double h,
                             double eps,
                             double b1,
                             bool const_term_in_N,
                             double reltol = 1e-12,
                             double abstol = 1e-12) {
  
  static Environment sundialr_env =
    Environment::namespace_env("sundialr");
  
  static Function cvode_solver =
    sundialr_env["cvode"];
  
  static Environment global_env =
    Environment::global_env();
  
  static Function rcpp_rhs_ptr =
    global_env["fh_rhs_2d"];
  
  NumericVector time_vector = NumericVector::create(0.0, h);
  
  NumericVector init_values = NumericVector::create(x_init, z_init);
  
  double flag_val = const_term_in_N ? 1.0 : 0.0;
  
  NumericVector solver_params =
    NumericVector::create(eps, b1, flag_val);
  
  NumericMatrix res = cvode_solver(
    _["time_vector"]    = time_vector,
    _["IC"]             = init_values,
    _["input_function"] = rcpp_rhs_ptr,
    _["Parameters"]     = solver_params,
    _["reltolerance"]   = reltol,
    _["abstolerance"]   = abstol
  );
  
  return NumericVector::create(res(1,1), res(1,2));
}

// [[Rcpp::export]]
arma::mat jacobian_fh_2d_rcpp(double x, double z, double h, double eps, double b, bool const_term_in_N) {
  const double eps_mach = std::numeric_limits<double>::epsilon();
  const double step_base = std::pow(eps_mach, 1.0 / 3.0);
  
  NumericVector states = NumericVector::create(x, z);
  
  arma::mat J(2, 2, arma::fill::zeros);
  
  for (int j = 0; j < 2; ++j) {
    double current_val = states[j];
    double h_step = step_base * std::max(1.0, std::abs(current_val));
    
    NumericVector xp = clone(states);
    NumericVector xm = clone(states);
    
    xp[j] = current_val + h_step;
    xm[j] = current_val - h_step;
    
    NumericVector fp = fh_sundials_2d(xp[0], xp[1], h, eps, b, const_term_in_N);
    NumericVector fm = fh_sundials_2d(xm[0], xm[1], h, eps, b, const_term_in_N);
    
    // Compute central finite differences
    J(0, j) = (fp[0] - fm[0]) / (2.0 * h_step); 
    J(1, j) = (fp[1] - fm[1]) / (2.0 * h_step);
  }
  
  return J;
}


// [[Rcpp::export]]
Rcpp::List log_lik_path_rcpp(NumericMatrix df,
                         NumericVector times,
                         NumericVector theta_drift,
                         NumericVector theta_diffusion,
                         const std::string &method,
                         const std::string &fh_solver,
                         const List &b_vec) {
  int nd = theta_drift.size();
  int ns = theta_diffusion.size();
  NumericVector par(nd + ns);
  for (int i = 0; i < nd; ++i) par[i] = theta_drift[i];
  for (int j = 0; j < ns; ++j) par[nd + j] = theta_diffusion[j];
  
  const int N = df.nrow();
  if (N < 2) return 0.0;
  if (times.size() < 2) stop("times length < 2");
  double h = times[1] - times[0]; 
  double eps = par(0);
  if (par.size() < 6) stop("par must contain at least 6 elements (sigma params at indices 4,5)");
  double left_knee = -1/std::sqrt(3), right_knee = 1/std::sqrt(3);
  double loglik = 0.0;
  
  // Common variables
  NumericVector center(2);
  arma::mat Omega_mat;
  arma::mat Inv_Omega;
  double signO = 0.0, logdetO = 0.0;
  
  if (method == "fix" || method == "buckwar" || method == "knee") {
    // choose center

    center = b_vec[0];

    Omega_mat = Omega_rcpp(h, par, center, method);
    Inv_Omega = arma::inv(Omega_mat); 
    arma::log_det(logdetO, signO, Omega_mat); 
    if (signO <= 0.0) stop("Omega determinant non-positive");
    
    // accumulate sums
    double sum_log_det_Omega = 0.0;
    double sum_quad = 0.0;
    double sum_log_det_D = 0.0;
    arma::vec Z = arma::vec(2);
    for (int i = 0; i < N - 1; ++i) {
      NumericVector x_old = NumericVector::create(df(i,0), df(i,1));
      NumericVector x_new = NumericVector::create(df(i+1,0), df(i+1,1));
      
      // finv at x_new
      // NumericVector finv = fh_rcpp(x_new, -h/2.0, par, center, method);
      NumericVector finv;
      // mu_f at x_old
      NumericVector tmp;
      if (fh_solver == "rk4"){
        finv = fh_rcpp(x_new, -h/2.0, par, center, method);
        tmp = fh_rcpp(x_old,  h/2.0, par, center, method);
      }else if (fh_solver == "sundials"){
        finv = fh_sundials_2d(x_new(0),x_new(1), -h/2.0, eps, center(0), false);
        tmp = fh_sundials_2d(x_old(0),x_old(1),  h/2.0, eps, center(0), false);
      }
      NumericVector mu_f = mu_rcpp(tmp, h, par, center, method);
      
      // residual
      arma::vec z(2);
      z(0) = finv[0] - mu_f[0];
      z(1) = finv[1] - mu_f[1];
      Z += z;
      arma::vec temp = Inv_Omega * z;
      double quad = arma::dot(z, temp);
      arma::mat D;
      if (fh_solver == "rk4"){
        NumericVector cur_center = center;
        NumericVector cur_par = par;
        std::string cur_method = method;
        std::function<NumericVector(NumericVector)> fh_wrapper =
          [cur_par, cur_center, cur_method, h](NumericVector xx) -> NumericVector {
            return fh_rcpp(xx, -h/2.0, cur_par, cur_center, cur_method);
          };
        D = numeric_jacobian_std(fh_wrapper, x_new);
      }else if (fh_solver == "sundials"){
        D = jacobian_fh_2d_rcpp(x_new(0),x_new(1), -h/2,  eps, center(0), false);
      }
      
      double sigD=0.0, logdetD = 0.0;
      arma::log_det(logdetD, sigD, D);
      
      sum_log_det_Omega += logdetO;
      sum_quad += quad;
      sum_log_det_D += logdetD;
    }
    
    loglik = sum_log_det_Omega + sum_quad - 2.0 * sum_log_det_D;
    return List::create(
        Named("ll") = loglik
    );
  } 
  // Three b's
  else if (method == "piecewise" || method == "piecewise taylor") {
    double alpha = par[1];
    const List& b_list = b_vec;
    
    const NumericVector b_left   = b_list[0];
    const NumericVector b_middle = b_list[1];
    const NumericVector b_right  = b_list[2];
    
    // compute Omegas and inverses for each part
    arma::mat Omega_left  = Omega_rcpp(h, par, b_left, method);
    arma::mat Omega_mid   = Omega_rcpp(h, par, b_middle, method);
    arma::mat Omega_right = Omega_rcpp(h, par, b_right, method);
    
    arma::mat Inv_left  = arma::inv(Omega_left);
    arma::mat Inv_mid   = arma::inv(Omega_mid); 
    arma::mat Inv_right = arma::inv(Omega_right);
    
    double logdet_left=0.0, sleft=0.0;
    double logdet_mid=0.0, smid=0.0;
    double logdet_right=0.0, sright=0.0;
    arma::log_det(logdet_left, sleft, Omega_left);
    arma::log_det(logdet_mid,  smid,  Omega_mid);
    arma::log_det(logdet_right,sright,Omega_right);
    
    double sum_log_det_Omega = 0.0;
    double sum_quad = 0.0;
    double sum_log_det_D = 0.0;
    
    for (int i = 0; i < N - 1; ++i) {
      NumericVector x_old = NumericVector::create(df(i,0), df(i,1));
      NumericVector x_new = NumericVector::create(df(i+1,0), df(i+1,1));
      double xold1 = x_old[0];
      double yold1 = x_old[1];
      bool cond_left  = (xold1 < left_knee);
      bool cond_right = (xold1 >  right_knee);
      
      const NumericVector* used_b;
      arma::mat *used_Inv;
      double used_logdet;
      if (cond_left) {
        used_b  = &b_left;
        used_Inv = &Inv_left;
        used_logdet = logdet_left;
      } else if (cond_right) {
        used_b =  &b_right;
        used_Inv = &Inv_right;
        used_logdet = logdet_right;
      } else {
        used_b = &b_middle;
        used_Inv = &Inv_mid;
        used_logdet = logdet_mid;
      }
      
      NumericVector finv = fh_rcpp(x_new, -h/2.0, par, *used_b, method);
      NumericVector tmp  = fh_rcpp(x_old,  h/2.0, par, *used_b, method);
      NumericVector mu_f = mu_rcpp(tmp, h, par, *used_b, method);
      
      arma::vec z(2);
      z(0) = finv[0] - mu_f[0];
      z(1) = finv[1] - mu_f[1];
      
      arma::vec temp = (*used_Inv) * z;
      double quad = arma::dot(z, temp);
      
      NumericVector cur_center = *used_b;
      NumericVector cur_par = par;
      std::function<NumericVector(NumericVector)> fh_wrapper =
        [cur_par, cur_center, h,method](NumericVector xx) -> NumericVector {
          return fh_rcpp(xx, -h/2.0, cur_par, cur_center, method);
        };
      arma::mat D = numeric_jacobian_std(fh_wrapper, x_new);
      double sigD=0.0, logdetD=0.0;
      arma::log_det(logdetD, sigD, D);
      
      sum_log_det_Omega += used_logdet;
      sum_quad += quad;
      sum_log_det_D += logdetD;
    }
    
    
    loglik = sum_log_det_Omega + sum_quad - 2.0 * sum_log_det_D;
    return List::create(Named("ll") = loglik);
  } // Two b's
  else if (method == "two knees") {
    double alpha = par[1];
    const List& b_list = b_vec;
    
    const NumericVector b_left   = b_list[0];
    const NumericVector b_right  = b_list[1];
    
    // compute Omegas and inverses for each part
    arma::mat Omega_left  = Omega_rcpp(h, par, b_left, method);
    arma::mat Omega_right = Omega_rcpp(h, par, b_right, method);
    
    arma::mat Inv_left  = arma::inv(Omega_left);
    arma::mat Inv_right = arma::inv(Omega_right);
    
    double logdet_left=0.0, sleft=0.0;
    double logdet_right=0.0, sright=0.0;
    arma::log_det(logdet_left, sleft, Omega_left);
    arma::log_det(logdet_right,sright,Omega_right);
    
    double sum_log_det_Omega = 0.0;
    double sum_quad = 0.0;
    double sum_log_det_D = 0.0;
    
    for (int i = 0; i < N - 1; ++i) {
      NumericVector x_old = NumericVector::create(df(i,0), df(i,1));
      NumericVector x_new = NumericVector::create(df(i+1,0), df(i+1,1));
      double xold1 = x_old[0];
      double yold1 = x_old[1];
      bool cond_left  = (yold1 < -1000);

      const NumericVector* used_b;
      arma::mat *used_Inv;
      double used_logdet;
      if (cond_left) {
        used_b  = &b_left;
        used_Inv = &Inv_left;
        used_logdet = logdet_left;
      } else  {
        used_b =  &b_right;
        used_Inv = &Inv_right;
        used_logdet = logdet_right;
      } 
      NumericVector finv;
      NumericVector tmp;
      if (fh_solver == "rk4"){
        finv = fh_rcpp(x_new, -h/2.0, par, *used_b, method);
        tmp  = fh_rcpp(x_old,  h/2.0, par, *used_b, method);
      }else if (fh_solver == "sundials"){
        finv = fh_sundials_2d(x_new(0),x_new(1), -h/2.0, eps, (*used_b)[0], false);
        tmp  = fh_sundials_2d(x_old(0),x_old(1),  h/2.0, eps, (*used_b)[0], false);
      }
 
      NumericVector mu_f = mu_rcpp(tmp, h, par, *used_b, method);
      
      arma::vec z(2);
      z(0) = finv[0] - mu_f[0];
      z(1) = finv[1] - mu_f[1];
      
      arma::vec temp = (*used_Inv) * z;
      double quad = arma::dot(z, temp);
      
      arma::mat D;
      if (fh_solver == "rk4"){
        NumericVector cur_center = *used_b;
        NumericVector cur_par = par;
        std::function<NumericVector(NumericVector)> fh_wrapper =
          [cur_par, cur_center, h,method](NumericVector xx) -> NumericVector {
            return fh_rcpp(xx, -h/2.0, cur_par, cur_center, method);
          };
        arma::mat D = numeric_jacobian_std(fh_wrapper, x_new);
      }else if (fh_solver == "sundials"){
        arma::mat D = jacobian_fh_2d_rcpp(x_new(0), x_new(1), -h/2, eps, (*used_b)[0], false);
      }
      
      double sigD=0.0, logdetD=0.0;
      arma::log_det(logdetD, sigD, D);
      
      sum_log_det_Omega += used_logdet;
      sum_quad += quad;
      sum_log_det_D += logdetD;
    }
    
    
    loglik = sum_log_det_Omega + sum_quad - 2.0 * sum_log_det_D;
    return List::create(Named("ll") = loglik);
  } 
  
  // Four b's 
  else if (method == "new_piecewise" || method == "new_piecewise taylor") {
    double alpha = par[1];
    const List& b_list = b_vec;
    NumericVector b_left  = b_vec[0];
    NumericVector b_upper= b_vec[1];
    NumericVector b_lower= b_vec[2];
    NumericVector b_right = b_vec[3];
    
    // compute Omegas and inverses for each part
    arma::mat Omega_left  = Omega_rcpp(h, par, b_left, method);
    arma::mat Omega_upper   = Omega_rcpp(h, par, b_upper, method);
    arma::mat Omega_lower   = Omega_rcpp(h, par, b_lower, method);
    arma::mat Omega_right = Omega_rcpp(h, par, b_right, method);
    
    arma::mat Inv_left  = arma::inv(Omega_left);
    arma::mat Inv_upper   = arma::inv(Omega_upper); 
    arma::mat Inv_lower   = arma::inv(Omega_lower); 
    arma::mat Inv_right = arma::inv(Omega_right);
    
    double logdet_left=0.0, sleft=0.0;
    double logdet_upper=0.0, supper=0.0;
    double logdet_lower=0.0, slower=0.0;
    double logdet_right=0.0, sright=0.0;
    arma::log_det(logdet_left, sleft, Omega_left);
    arma::log_det(logdet_upper,  supper,  Omega_upper);
    arma::log_det(logdet_lower,  slower,  Omega_lower);
    arma::log_det(logdet_right,sright,Omega_right);
    
    double sum_log_det_Omega = 0.0;
    double sum_quad = 0.0;
    double sum_log_det_D = 0.0;
    
    for (int i = 0; i < N - 1; ++i) {
      NumericVector x_old = NumericVector::create(df(i,0), df(i,1));
      NumericVector x_new = NumericVector::create(df(i+1,0), df(i+1,1));
      double xold1 = x_old[0];
      double yold1 = x_old[1];
      bool cond_left  = (xold1 < left_knee);
      bool cond_right = (xold1 >  right_knee);
      bool cond_upper  = (yold1 > alpha);
      bool cond_lower = (yold1 <  alpha);
      
      const NumericVector* used_b;
      arma::mat *used_Inv;
      double used_logdet;
      if (cond_left) {
        used_b = &b_left;
        used_Inv = &Inv_left;
        used_logdet = logdet_left;
      } else if (cond_right) {
        used_b = &b_right;
        used_Inv = &Inv_right;
        used_logdet = logdet_right;
      } else {
        if (cond_upper){
          used_b = &b_upper;
          used_Inv = &Inv_upper;
          used_logdet = logdet_upper;
        }else{
          used_b = &b_lower;
          used_Inv = &Inv_lower;
          used_logdet = logdet_lower;
        }
      }
      
      NumericVector finv;
      NumericVector tmp;
      if (fh_solver == "rk4"){
        finv = fh_rcpp(x_new, -h/2.0, par, *used_b, method);
        tmp  = fh_rcpp(x_old,  h/2.0, par, *used_b, method);
      }else if (fh_solver == "sundials"){
        finv = fh_sundials_2d(x_new(0),x_new(1), -h/2.0, eps, (*used_b)[0], false);
        tmp  = fh_sundials_2d(x_old(0),x_old(1),  h/2.0, eps, (*used_b)[0], false);
      }

      NumericVector mu_f = mu_rcpp(tmp, h, par, *used_b, method);
      
      arma::vec z(2);
      z(0) = finv[0] - mu_f[0];
      z(1) = finv[1] - mu_f[1];
      
      arma::vec temp = (*used_Inv) * z;
      double quad = arma::dot(z, temp);
      
      arma::mat D;
      if (fh_solver == "rk4"){
        NumericVector cur_center = *used_b;
        NumericVector cur_par = par;
        std::function<NumericVector(NumericVector)> fh_wrapper =
          [cur_par, cur_center, h,method](NumericVector xx) -> NumericVector {
            return fh_rcpp(xx, -h/2.0, cur_par, cur_center, method);
          };
        arma::mat D = numeric_jacobian_std(fh_wrapper, x_new);
      }else if (fh_solver == "sundials"){
        D = jacobian_fh_2d_rcpp(x_new(0), x_new(1), -h/2, eps, (*used_b)[0], false);
      }

      double sigD=0.0, logdetD=0.0;
      arma::log_det(logdetD, sigD, D);
      
      sum_log_det_Omega += used_logdet;
      sum_quad += quad;
      sum_log_det_D += logdetD;
    }
    
    
    loglik = sum_log_det_Omega + sum_quad - 2.0 * sum_log_det_D;
    return List::create(Named("ll") = loglik);
  } 
  // seven b's 
  else if (method == "eigen") {
    double alpha = par[1];
    double eps = par[0];
    double gamma = par[2];
    double beta = par[3];
    const List& b_list = b_vec;
    NumericVector b_0  = b_vec[0];
    NumericVector b_1r= b_vec[1];
    NumericVector b_2r= b_vec[2];
    NumericVector b_3r = b_vec[3];
    NumericVector b_1l  = b_vec[4];
    NumericVector b_2l= b_vec[5];
    NumericVector b_3l= b_vec[6];
    
    // compute Omegas and inverses for each part
    arma::mat Omega_0  = Omega_rcpp(h, par, b_0, "eigen");
    arma::mat Omega_1r   = Omega_rcpp(h, par, b_1r, "eigen");
    arma::mat Omega_2r   = Omega_rcpp(h, par, b_2r, "eigen");
    arma::mat Omega_3r = Omega_rcpp(h, par, b_3r, "eigen");
    arma::mat Omega_1l   = Omega_rcpp(h, par, b_1l, "eigen");
    arma::mat Omega_2l   = Omega_rcpp(h, par, b_2l, "eigen");
    arma::mat Omega_3l = Omega_rcpp(h, par, b_3l, "eigen");
    
    arma::mat Inv_0  = arma::inv(Omega_0);
    arma::mat Inv_1r   = arma::inv(Omega_1r); 
    arma::mat Inv_2r   = arma::inv(Omega_2r); 
    arma::mat Inv_3r = arma::inv(Omega_3r);
    arma::mat Inv_1l   = arma::inv(Omega_1l); 
    arma::mat Inv_2l   = arma::inv(Omega_2l); 
    arma::mat Inv_3l = arma::inv(Omega_3l);
    
    double logdet_0=0.0, s0=0.0;
    double logdet_1r=0.0, s1r=0.0;
    double logdet_2r=0.0, s2r=0.0;
    double logdet_3r=0.0, s3r=0.0;
    double logdet_1l=0.0, s1l=0.0;
    double logdet_2l=0.0, s2l=0.0;
    double logdet_3l=0.0, s3l=0.0;
    arma::log_det(logdet_0, s0, Omega_0);
    arma::log_det(logdet_1r,  s1r,  Omega_1r);
    arma::log_det(logdet_2r,  s2r,  Omega_2r);
    arma::log_det(logdet_3r,s3r,Omega_3r);
    arma::log_det(logdet_1l,  s1l,  Omega_1l);
    arma::log_det(logdet_2l,  s2l,  Omega_2l);
    arma::log_det(logdet_3l,s3l,Omega_3l);
    
    double sum_log_det_Omega = 0.0;
    double sum_quad = 0.0;
    double sum_log_det_D = 0.0;
    
    for (int i = 0; i < N - 1; ++i) {
      NumericVector x_old = NumericVector::create(df(i,0), df(i,1));
      NumericVector x_new = NumericVector::create(df(i+1,0), df(i+1,1));
      double xold1 = x_old[0];
      double yold1 = x_old[1];
      
      double tau   = (1 - 3*xold1*xold1 - eps) / eps;
      double Delta = 3*xold1*xold1/eps - 1/eps + gamma/eps;
      
      double disc = tau*tau - 4*Delta;
      bool real = (disc >= 0);
      
      double root = real ? std::sqrt(disc) : 0.0;
      
      bool cond_0  = real && (tau - root >= 0);
      
      bool cond_1r = !real && (tau >= 0) && (xold1 >= 0);
      bool cond_2r = !real && (tau < 0)  && (xold1 >= 0);
      bool cond_3r = real  && (tau + root < 0) && (xold1 >= 0);
      
      bool cond_1l = !real && (tau >= 0) && (xold1 < 0);
      bool cond_2l = !real && (tau < 0)  && (xold1 < 0);
      bool cond_3l = real  && (tau + root < 0) && (xold1 < 0);
      
      const NumericVector* used_b;
      arma::mat *used_Inv;
      double used_logdet;
      if (cond_0) {
        used_b = &b_0;
        used_Inv = &Inv_0;
        used_logdet = logdet_0;
      } else if (cond_1r) {
        used_b = &b_1r;
        used_Inv = &Inv_1r;
        used_logdet = logdet_1r;
      } else if (cond_2r) {
        used_b = &b_2r;
        used_Inv = &Inv_2r;
        used_logdet = logdet_2r;
      } else if (cond_3r) {
        used_b = &b_3r;
        used_Inv = &Inv_3r;
        used_logdet = logdet_3r;
      } else if (cond_1l) {
        used_b = &b_1l;
        used_Inv = &Inv_1l;
        used_logdet = logdet_1l;
      } else if (cond_2l) {
        used_b = &b_2l;
        used_Inv = &Inv_2l;
        used_logdet = logdet_2l;
      } else if (cond_3l) {
        used_b = &b_3l;
        used_Inv = &Inv_3l;
        used_logdet = logdet_3l;
      } 
      
      NumericVector finv = fh_rcpp(x_new, -h/2.0, par, *used_b, "eigen");
      NumericVector tmp  = fh_rcpp(x_old,  h/2.0, par, *used_b, "eigen");
      NumericVector mu_f = mu_rcpp(tmp, h, par, *used_b, "eigen");
      
      arma::vec z(2);
      z(0) = finv[0] - mu_f[0];
      z(1) = finv[1] - mu_f[1];
      
      arma::vec temp = (*used_Inv) * z;
      double quad = arma::dot(z, temp);
      
      NumericVector cur_center = *used_b;
      NumericVector cur_par = par;
      std::function<NumericVector(NumericVector)> fh_wrapper =
        [cur_par, cur_center, h](NumericVector xx) -> NumericVector {
          return fh_rcpp(xx, -h/2.0, cur_par, cur_center, "eigen");
        };
        arma::mat D = numeric_jacobian_std(fh_wrapper, x_new);
        double sigD=0.0, logdetD=0.0;
        arma::log_det(logdetD, sigD, D);
        
        sum_log_det_Omega += used_logdet;
        sum_quad += quad;
        sum_log_det_D += logdetD;
    }
    
    
    loglik = sum_log_det_Omega + sum_quad - 2.0 * sum_log_det_D;
    return List::create(Named("ll") = loglik);
  } 
  // Several b's 
  else if (method == "custom" || method == "custom test" || method == "fast eigen" || method == "mixed eigen" || method == "local" || method == "fenichel") {
    
    double eps   = par[0];
    double alpha = par[1];
    double gamma = par[2];
    double beta  = par[3];
    
    double sum_log_det_Omega = 0.0;
    double sum_quad = 0.0;
    double sum_log_det_D = 0.0;
    
    for (int i = 0; i < N - 1; ++i) {
      NumericVector x_old = NumericVector::create(df(i,0), df(i,1));
      NumericVector x_new = NumericVector::create(df(i+1,0), df(i+1,1));
      double xold1 = x_old[0];
      double yold1 = x_old[1];
      
      NumericVector b = b_vec[i];
      double b1 = b[0];
      double b2 = b[1];

      arma::mat Omega = Omega_rcpp(h, par, b, method);
      
      double logdet, sign;
      arma::log_det(logdet, sign, Omega);
    
      NumericVector finv = fh_rcpp(x_new, -h/2.0, par, b, method);
      //NumericVector finv = fh_sundials_2d(x_new(0),x_new(1), -h/2.0, eps, b(0), false);
      NumericVector tmp  = fh_rcpp(x_old,  h/2.0, par, b, method);
      //NumericVector tmp  = fh_sundials_2d(x_old(0),x_old(1),  h/2.0, eps, b(0),false);
      NumericVector mu_f = mu_rcpp(tmp, h, par, b, method);
      
      arma::vec z(2);
      z(0) = finv[0] - mu_f[0];
      z(1) = finv[1] - mu_f[1];
      
      arma::vec temp = arma::solve(Omega, z);     // solves Omega * temp = z
      double quad = arma::dot(z, temp);
      
      //Jacobian D at x_new
      NumericVector cur_center = b;
      NumericVector cur_par = par;
      std::function<NumericVector(NumericVector)> fh_wrapper =
        [cur_par, cur_center, h,method](NumericVector xx) -> NumericVector {
          return fh_rcpp(xx, -h/2.0, cur_par, cur_center, method);
        };
      arma::mat D = numeric_jacobian_std(fh_wrapper, x_new);
      //arma::mat D = jacobian_fh_2d_rcpp(x_new(0), x_new(1), -h/2, eps, b(0), false);
      double sigD = 0.0, logdetD = 0.0;
      arma::log_det(logdetD, sigD, D);
      
      // accumulate
      sum_log_det_Omega += logdet;
      sum_quad += quad;
      sum_log_det_D += logdetD;
    } 
    
    loglik = sum_log_det_Omega + sum_quad - 2.0 * sum_log_det_D;
    return List::create(Named("ll") = loglik);
  }else {
    stop("Unknown method in log_lik_path_rcpp");
  }
}

// equivalent mvrnorm from MASS in R
// [[Rcpp::export]]
arma::mat rmvnorm_rcpp(int n, const arma::mat& Sigma) {
  Rcpp::Environment MASS = Rcpp::Environment::namespace_env("MASS");
  Rcpp::Function mvrnorm = MASS["mvrnorm"];
  
  int d = Sigma.n_rows;
  Rcpp::NumericVector mu(d, 0.0);
  
  Rcpp::NumericMatrix res = mvrnorm(
    Rcpp::Named("n") = n,
    Rcpp::Named("mu") = mu,
    Rcpp::Named("Sigma") = Rcpp::wrap(Sigma)
  );
  
  return Rcpp::as<arma::mat>(res);
}

// [[Rcpp::export]]
Rcpp::DataFrame one_step_pred_rcpp(const NumericMatrix x0,
                                   double h,
                                   const NumericVector par,
                                   int n,
                                   std::string method = "fix",
                                   std::string order  = "fmuf",
                                   std::string fh_solver  = "sundials",
                                   NumericVector center = NumericVector::create()) {
  
  int n_x0 = x0.nrow();
  int total = n_x0 * n;
  
  // output storage 
  Rcpp::NumericVector out_x(total);
  Rcpp::NumericVector out_y(total);
  Rcpp::NumericVector out_x0_1(total);
  Rcpp::NumericVector out_x0_2(total);
  
  int pos = 0;
  double alpha = par[1];
  double eps = par[0];
  // check user-provided center (if any)
  if (center.size() != 0 && center.size() != 2) {
    stop("If supplied, 'center' must have length 2.");
  }
  
  // main loop over initial points
  for (int j = 0; j < n_x0; ++j) {
    
    NumericVector x0j = x0(j, _); 
    
    NumericVector center_j;
    arma::mat Sigma_noise;
    arma::mat half_Sigma_noise;
    
    // If user provided a center (length 2), use it for every row
    if (center.size() == 2) {
      center_j = NumericVector::create(center[0], center[1]);
    } else {
      if (method == "fix") {
        center_j = fixed_point_rcpp(par);
      } else if (method == "buckwar") {
        center_j = NumericVector::create(0.0, 0.0);
      } else if (method == "piecewise") {
        double x = x0(j, 0);
        if (x < -1.0 / std::sqrt(3.0))
          center_j = NumericVector::create(-1.0, alpha);
        else if (x > 1.0 / std::sqrt(3.0))
          center_j = NumericVector::create( 1.0, alpha);
        else
          center_j = NumericVector::create(0.0, alpha);
      } else {
        stop("Unknown method and no center provided");
      }
    }
    
    // compute covariances using center_j
    Sigma_noise      = Omega_rcpp(h, par, center_j, method);
    half_Sigma_noise = Omega_rcpp(h / 2.0, par, center_j, method);
    
    // normal order
    if (order == "fmuf") {
      NumericVector S1;
      if (fh_solver == "sundials"){
        S1 = fh_sundials_2d(x0j(0),x0j(1), h / 2.0, eps, center_j(0), false);
      }else if (fh_solver == "rk4"){
        S1 = fh_rcpp(x0j, h / 2.0, par, center_j, method);
      }
      arma::vec muS1 = as<arma::vec>(mu_rcpp(S1, h, par, center_j, method));
      
      arma::mat xi_mat = rmvnorm_rcpp(n, Sigma_noise);
      arma::mat S2_mat = xi_mat.each_row() + muS1.t();
      
      for (int k = 0; k < n; ++k) {
        NumericVector S2k = NumericVector::create(
          S2_mat(k, 0),
          S2_mat(k, 1)
        );
        
        arma::vec S3;
        if (fh_solver == "sundials"){
          S3 = as<arma::vec>(
            fh_sundials_2d(S2k(0),S2k(1), h / 2.0, eps, center_j(0), false)
          );
        }else if (fh_solver == "rk4"){
          S3 = as<arma::vec>(
            fh_rcpp(S2k, h / 2.0, par, center_j,method)
          );
        }
          
        
        out_x[pos]    = S3(0);
        out_y[pos]    = S3(1);
        out_x0_1[pos] = x0(j, 0);
        out_x0_2[pos] = x0(j, 1);
        pos++;
      }
      
    }
    // reversed order
    else if (order == "mufmu") {
      
      arma::vec mu_half =
        as<arma::vec>(mu_rcpp(x0j, h / 2.0, par, center_j, method));
      
      arma::mat xi1 = rmvnorm_rcpp(n, half_Sigma_noise);
      arma::mat xi2 = rmvnorm_rcpp(n, half_Sigma_noise);
      
      arma::mat inner = xi1.each_row() + mu_half.t();
      
      NumericVector inner_k(2);   
      NumericVector S_inner_nv(2);
      
      for (int k = 0; k < n; ++k) {
        
        inner_k[0] = inner(k, 0);
        inner_k[1] = inner(k, 1);
        
        arma::vec S_inner =
          as<arma::vec>(
            fh_rcpp(inner_k, h, par, center_j, method)
          );
        
        S_inner_nv[0] = S_inner(0);
        S_inner_nv[1] = S_inner(1);
        
        arma::vec mu_inner =
          as<arma::vec>(
            mu_rcpp(S_inner_nv, h / 2.0, par, center_j, method)
          );
        
        arma::vec final = mu_inner + xi2.row(k).t();
        
        out_x[pos]    = final(0);
        out_y[pos]    = final(1);
        out_x0_1[pos] = x0(j, 0);
        out_x0_2[pos] = x0(j, 1);
        pos++;
      }
      
    } else {
      Rcpp::stop("Unknown order");
    }
  }
  
  // build output 
  Rcpp::CharacterVector label(total, "S3");
  
  Rcpp::CharacterVector x0_label(total);
  for (int i = 0; i < total; ++i) {
    x0_label[i] =
      std::to_string(out_x0_1[i]) + "_" +
      std::to_string(out_x0_2[i]);
  }
  
  return Rcpp::DataFrame::create(
    Rcpp::Named("x") = out_x,
    Rcpp::Named("y") = out_y,
    Rcpp::Named("x0_1") = out_x0_1,
    Rcpp::Named("x0_2") = out_x0_2,
    Rcpp::Named("label") = label,
    Rcpp::Named("x0_label") = x0_label
  );
}

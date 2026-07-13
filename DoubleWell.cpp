// doublewell_rcpp.cpp
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;



// A(par, center, part) -> scalar
static double A_cpp_scalar(NumericVector par, double b) {
  double eps = par(0);
  return (1.0/eps) * (-3.0 * b * b + 1.0);

}

static double N(double x, const arma::vec &par, double b, bool const_term_in_N) {
  double eps = par(0);
  double y   = par(1);
  double du1;
  
  if (const_term_in_N) {
    du1 = (1.0 / eps) * (-x*x*x - y + 3.0*b*b*x + b - 3.0*b*b*b);
  } else {
    du1 = (1.0 / eps) * (-x*x*x + 3.0*b*b*x - 2.0*b*b*b);
  }
  
  return du1;
}

// single-step RK4 for fh (scalar ODE)
// [[Rcpp::export]]
static double fh_scalar(double x, double h, const arma::vec &par, double b,bool const_term_in_N) {
  double k1 = N(x, par,b, const_term_in_N);
  double k2 = N(x + 0.5 * h * k1, par,b, const_term_in_N);
  double k3 = N(x + 0.5 * h * k2, par,b,const_term_in_N);
  double k4 = N(x + h * k3, par,b,const_term_in_N);
  double x_new = x + (h/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
  return x_new;
}


// [[Rcpp::export]]
NumericVector fh_rhs(double t, NumericVector y, NumericVector par) {
  double x = y[0]; // state variable
  
  double eps           = par[0];
  double b             = par[1];
  bool const_term_in_N = (par[2] == 1.0); 
  
  double du1;
  if (const_term_in_N) {
    du1 = (1.0 / eps) * (-x*x*x - b + 3.0*b*b*x + b - 3.0*b*b*b); 
  } else {
    du1 = (1.0 / eps) * (-x*x*x + 3.0*b*b*x - 2.0*b*b*b);
  }
  
  NumericVector out(1);
  out[0] = du1;
  return out;
}

// [[Rcpp::export]]
double fh_sundials(double x_init, double h, double eps, double b, bool const_term_in_N) {
  Environment sundialr_env = Environment::namespace_env("sundialr");
  Function cvode_solver = sundialr_env["cvode"];
  
  // FIX: Look up the exported C++ function pointer in R's global environment
  Environment global_env = Environment::global_env();
  Function rcpp_rhs_ptr = global_env["fh_rhs"];
  
  NumericVector time_vector = NumericVector::create(0.0, h);
  NumericVector init_values = NumericVector::create(x_init);
  
  double flag_val = const_term_in_N ? 1.0 : 0.0;
  NumericVector solver_params = NumericVector::create(eps, b, flag_val);
  
  // Call solver passing the direct Rcpp pointer object instead of a string name
  NumericMatrix res = cvode_solver(
    _["time_vector"]    = time_vector,
    _["IC"]             = init_values,       
    _["input_function"] = rcpp_rhs_ptr,      // FIX: Changed from "fh_rhs" to rcpp_rhs_ptr
    _["Parameters"]     = solver_params,     
    _["reltolerance"]   = 1e-6,              
    _["abstolerance"]   = 1e-8               
  );
  
  double x_new = res(1, 1); 
  return x_new;
}

// mu scalar
// [[Rcpp::export]]
static double mu_scalar(double x, double h, NumericVector par, double b ,bool const_term_in_N) {
  double A = A_cpp_scalar(par, b);
  double M = std::exp(A * h);
  
  if (const_term_in_N){
    return M * x + (1.0 - M) * b;  
  }else{
    double Fb = (-b*b*b+b-par(1))/par(0);
    double b_tilde = b - Fb/A;
    return M * x + (1.0 - M) * b_tilde;
  }
}
// [[Rcpp::export]]
static double Omega_scalar(double h, NumericVector par, double b) {
  double a = A_cpp_scalar(par, b);
  double sigma12 = par(2)*par(2);
  double val = (-sigma12 / (2.0 * a)) * (1.0 - std::exp(2.0 * a * h));
  return val;
}

// should be equivalent to numDeriv in R
// [[Rcpp::export]]
double jacobian_fh_scalar_rcpp(double x, double h,
                               const arma::vec &par, double b, bool const_term_in_N,std::string fh_solver) {
  double epsilon = par[0];
  // Base step
  const double eps = std::numeric_limits<double>::epsilon();
  const double step_base = std::pow(eps, 1.0/3.0);
  
  // Scale step with x
  double scale = std::max(1.0, std::abs(x));
  double h1 = step_base * scale;
  double h2 = h1 / 2.0;
  double f1p;
  double f1m;
  // ---- First central difference (h1)
  if (fh_solver == "rk4"){
    f1p = fh_scalar(x + h1, h, par, b, const_term_in_N);
    f1m = fh_scalar(x - h1, h, par, b, const_term_in_N);
  }else if (fh_solver == "sundials"){
    f1p = fh_sundials(x + h1, h, epsilon, b, const_term_in_N);
    f1m = fh_sundials(x - h1, h, epsilon, b, const_term_in_N);
  }else{
    Rcpp::stop("fh_solver_str must be 'rk4' or 'sundials'");
  }
  double D1 = (f1p - f1m) / (2.0 * h1);
  
  //---- Second central difference (h2)
  // double f2p = fh_sundials(x + h2, h, epsilon, b, const_term_in_N);
  // double f2m = fh_sundials(x - h2, h, epsilon, b, const_term_in_N);
  // double D2 = (f2p - f2m) / (2.0 * h2);

  //double D = (4.0 * D2 - D1) / 3.0;
  return D1;
}
// [[Rcpp::export]]
NumericVector sigma_MLE(NumericVector Z, NumericVector A, const NumericVector &times_r) {
  int n = Z.size();
  NumericVector result(n);
  double h = times_r[1] - times_r[0];
  
  for (int i = 0; i < n; i++) {
    result[i] = (2.0 * Z[i] * Z[i] * A[i]) /
      (std::exp(2.0 * A[i] * h) - 1.0);
  }
  
  return result;
}

// [[Rcpp::export]]
Rcpp::List log_lik_path_rcpp(const NumericMatrix &df_r,
                             const NumericVector &times_r,
                             const NumericVector &theta_drift,
                             const NumericVector &theta_diffusion,
                             Rcpp::RObject method,
                             const NumericVector &b_vec, double uns_fix,Rcpp::RObject fh_solver) {
  
  int nd = theta_drift.size();
  if (nd < 2) stop("theta_drift must have at least 2 elements.");
  
  NumericVector par(nd + 1);
  for (int i = 0; i < nd; i++) par(i) = theta_drift[i];
  par(2) = theta_diffusion(0); 
  double eps = par(0);
  double y = par(1);
  
  arma::mat data = as<arma::mat>(df_r);
  int N = data.n_rows;
  if (N < 2) stop("df must have at least 2 rows.");
  if (times_r.size() < 2) stop("times_r must have at least 2 elements.");
  
  double h = times_r[1] - times_r[0];
  int L = N - 1;
  std::string method_str = as<std::string>(method);
  std::string fh_solver_str = as<std::string>(fh_solver);
  arma::vec data_vec = data.col(0);
  arma::vec data_old = data_vec.rows(0, N - 2);
  arma::vec data_new = data_vec.rows(1, N - 1);
  arma::vec data_vec_shifted;
  arma::vec data_new_shifted;
  
  std::vector<double> Z_list(L);
  std::vector<double> b_list(L);
  std::vector<double> A_list(L);
  std::vector<double> logdensity_list(L);
  
  
  // ===================== one b =====================
  if (method_str == "uns"||  method_str == "average"|| method_str == "pi"|| method_str == "zero"|| method_str == "0.5"|| method_str == "left fix") {
    
    bool const_term_in_N = false;
    double use_b =b_vec[0];
    double Omega  = Omega_scalar(h, par, use_b);
    
    double inv_Omega  = 1.0 / std::abs(Omega);
    
    double log_Omega  = std::log(std::abs(Omega));
    
    double sum_logD = 0.0;
    double sum_quad = 0.0;
    double sum_logOmega = 0.0;
    
    for (int i = 0; i < L; i++) {
      double x_old = data_old(i);
      double x_new = data_new(i);
      
      
      // double tmp = fh_scalar(x_old, h / 2.0, par, use_b, const_term_in_N);
      // double z = fh_scalar(x_new, -h / 2.0, par, use_b, const_term_in_N) - 
      //   mu_scalar(tmp, h, par, use_b, const_term_in_N);
      // Z_list[i] = z;
      double tmp;
      double z;
      if (fh_solver_str == "rk4"){
        tmp = fh_scalar(x_old, h / 2.0, par, use_b, const_term_in_N);
        z = fh_scalar(x_new, -h / 2.0, par, use_b, const_term_in_N) -
          mu_scalar(tmp, h, par, use_b, const_term_in_N);
      }else if (fh_solver_str == "sundials"){
        tmp = fh_sundials(x_old, h / 2.0, eps, use_b, const_term_in_N);
        z = fh_sundials(x_new, -h / 2.0, eps, use_b, const_term_in_N) - 
          mu_scalar(tmp, h, par, use_b, const_term_in_N);
      }else{
        Rcpp::stop("fh_solver_str must be 'rk4' or 'sundials'");
      }
      Z_list[i] = z;
      
      double A = A_cpp_scalar(par, use_b);
      A_list[i] = A;
      
      double D = jacobian_fh_scalar_rcpp(x_new, -h / 2.0, par, use_b, const_term_in_N,fh_solver_str);
      if (std::abs(D) < 1e-12) D = (D >= 0 ? 1e-12 : -1e-12);
      
      sum_logD += std::log(std::abs(D));
      sum_quad += z * inv_Omega * z;
      sum_logOmega += log_Omega;
  
    }
    
    double ll = sum_logOmega + sum_quad - 2.0 * sum_logD;
    return List::create(
      Named("ll") = ll,
      Named("A_list") = A_list,
      Named("Z_list") = Z_list
    );
  }
  // ===================== two b's =====================
  else if (method_str == "fix" || method_str == "wrong fix" || method_str == "minimax wells 2" || method_str == "argminVargminB wells 2"|| method_str == "argminBargminV wells 2") {
    
    bool const_term_in_N;
    const_term_in_N = false;
  
    double left_b  = b_vec[0];
    double right_b = b_vec[1];
    
    
    double Omega_left  = Omega_scalar(h, par, left_b);
    double Omega_right = Omega_scalar(h, par, right_b);
    
    double inv_left  = 1.0 / std::abs(Omega_left);
    double inv_right = 1.0 / std::abs(Omega_right);
    
    double log_left  = std::log(std::abs(Omega_left));
    double log_right = std::log(std::abs(Omega_right));
    
    double sum_logD = 0.0;
    double sum_quad = 0.0;
    double sum_logOmega = 0.0;

    for (int i = 0; i < L; i++) {
      double x_old = data_old(i);
      double x_new; 
  
      x_new = data_new(i);
      
      
      bool cond = (x_old > uns_fix);
      double use_b = cond ? right_b : left_b;
      double tmp;
      double z;
      if (fh_solver_str == "rk4"){
        tmp = fh_scalar(x_old, h / 2.0, par, use_b, const_term_in_N);
        z = fh_scalar(x_new, -h / 2.0, par, use_b, const_term_in_N) -
         mu_scalar(tmp, h, par, use_b, const_term_in_N);
      }else if (fh_solver_str == "sundials"){
        tmp = fh_sundials(x_old, h / 2.0, eps, use_b, const_term_in_N);
        z = fh_sundials(x_new, -h / 2.0, eps, use_b, const_term_in_N) - 
          mu_scalar(tmp, h, par, use_b, const_term_in_N);
      }else{
        Rcpp::stop("fh_solver_str must be 'rk4' or 'sundials'");
      }

      Z_list[i] = z;
      
      // double A = A_cpp_scalar(par, use_b);
      // A_list[i] = A;
      // double denom = std::exp(2 * A * h) - 1.0;
      // if (std::abs(denom) < 1e-12) denom = 1e-12;
      // 
      // sum_sigma += (2.0 * z * z * A) / denom;
   

      double D = jacobian_fh_scalar_rcpp(x_new, -h / 2.0, par, use_b, const_term_in_N, fh_solver_str);
      if (std::abs(D) < 1e-12) D = (D >= 0 ? 1e-12 : -1e-12);
      
      sum_logD += std::log(std::abs(D));
      if (cond) {
        sum_quad += z * inv_right * z;
        sum_logOmega += log_right;
      } else {
        sum_quad += z * inv_left * z;
        sum_logOmega += log_left;
      }
    }
    
    double ll = sum_logOmega + sum_quad - 2.0 * sum_logD;
    
    return List::create(
      Named("ll") = ll,
      Named("Z_list") = Z_list
    );
  }
  // three b's
  else if (method_str == "random"|| method_str == "zero" || method_str == "all fix" || method_str == "minimax wells 3" || method_str == "argminVargminB wells 3"|| method_str == "argminBargminV wells 3") {
    

    bool const_term_in_N = false;

    double left_b  = b_vec[0];
    double middle_b = b_vec[1];
    double right_b = b_vec[2];
    
    
    double Omega_left  = Omega_scalar(h, par, left_b);
    double Omega_middle  = Omega_scalar(h, par, middle_b);
    double Omega_right = Omega_scalar(h, par, right_b);
    
    double inv_left  = 1.0 / std::abs(Omega_left);
    double inv_middle  = 1.0 / std::abs(Omega_middle);
    double inv_right = 1.0 / std::abs(Omega_right);
    
    double log_left  = std::log(std::abs(Omega_left));
    double log_middle  = std::log(std::abs(Omega_middle));
    double log_right = std::log(std::abs(Omega_right));
    
    double sum_logD = 0.0;
    double sum_quad = 0.0;
    double sum_logOmega = 0.0;
    
    for (int i = 0; i < L; i++) {
      double x_old = data_old(i);
      double x_new = data_new(i);
      
      bool cond1 = (x_old > 1/std::sqrt(3.0));
      bool cond2 = (x_old < -1/std::sqrt(3.0));
      double use_b = cond1 ? right_b : (cond2 ? left_b : middle_b);
      
      double tmp;
      double z;
      if (fh_solver_str == "rk4"){
        tmp = fh_scalar(x_old, h / 2.0, par, use_b, const_term_in_N);
        z = fh_scalar(x_new, -h / 2.0, par, use_b, const_term_in_N) -
          mu_scalar(tmp, h, par, use_b, const_term_in_N);
      }else if (fh_solver_str == "sundials"){
        tmp = fh_sundials(x_old, h / 2.0, eps, use_b, const_term_in_N);
        z = fh_sundials(x_new, -h / 2.0, eps, use_b, const_term_in_N) - 
          mu_scalar(tmp, h, par, use_b, const_term_in_N);
      }else{
        Rcpp::stop("fh_solver_str must be 'rk4' or 'sundials'");
      }
      
      Z_list[i] = z;
      
      // double A = A_cpp_scalar(par, use_b);
      // A_list[i] = A;
      // double denom = std::exp(2 * A * h) - 1.0;
      // if (std::abs(denom) < 1e-12) denom = 1e-12;
      // 
      // sum_sigma += (2.0 * z * z * A) / denom;
      
      
      double D = jacobian_fh_scalar_rcpp(x_new, -h / 2.0, par, use_b, const_term_in_N,fh_solver_str);
      if (std::abs(D) < 1e-12) D = (D >= 0 ? 1e-12 : -1e-12);
      
      sum_logD += std::log(std::abs(D));
      
      if (cond1) {
        sum_quad += z * inv_right * z;
        sum_logOmega += log_right;
      } else if (cond2) {
        sum_quad += z * inv_left * z;
        sum_logOmega += log_left;
      } else {
        sum_quad += z * inv_middle * z;
        sum_logOmega += log_middle;
      }
    }
    
    double ll = sum_logOmega + sum_quad - 2.0 * sum_logD;
    
    return List::create(
      Named("ll") = ll,
      Named("Z_list") = Z_list
    );
  }
  // ===================== multiple b's =====================
  else if (method_str == "local"  || method_str == "minimax steps" || method_str == "argminVargminB steps" || method_str == "argminBargminV steps") {
    double sum_logD = 0.0;
    double sum_quad = 0.0;
    double sum_logOmega = 0.0;
    bool const_term_in_N = false;
    for (int i = 0; i < L; i++) {
      double x_old = data_old(i);
      double x_new = data_new(i);
      
      //double use_b = closest_real_root_rcpp(x_old, par, const_term_in_N);
      double use_b = b_vec[i];
      double tmp;
      double z;
      if (fh_solver_str == "rk4"){
        tmp = fh_scalar(x_old, h / 2.0, par, use_b, const_term_in_N);
        z = fh_scalar(x_new, -h / 2.0, par, use_b, const_term_in_N) - 
          mu_scalar(tmp, h, par, use_b, const_term_in_N);
      }else if (fh_solver_str == "sundials"){
        tmp = fh_sundials(x_old, h / 2.0, eps, use_b, const_term_in_N);
        z = fh_sundials(x_new, -h / 2.0, eps, use_b, const_term_in_N) - 
          mu_scalar(tmp, h, par, use_b, const_term_in_N);
      }else{
        Rcpp::stop("fh_solver_str must be 'rk4' or 'sundials'");
      }
      
      Z_list[i] = z;
      
      double Omega = Omega_scalar(h, par, use_b);
      double absOmega = std::abs(Omega);
      if (absOmega < 1e-12) absOmega = 1e-12;
      
      double D = jacobian_fh_scalar_rcpp(x_new, -h / 2.0, par, use_b, const_term_in_N,fh_solver_str);
      if (std::abs(D) < 1e-12) D = (D >= 0 ? 1e-12 : -1e-12);
      
      sum_logD += std::log(std::abs(D));
      sum_quad += z * (1.0 / absOmega) * z;
      sum_logOmega += std::log(absOmega);
    }
    
    double ll = sum_logOmega + sum_quad - 2.0 * sum_logD;
    
    return List::create(
      Named("ll") = ll,
      Named("Z_list") = Z_list
    );
  }
  else {
    stop("Unknown method in log_lik_path_rcpp");
  }
}
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
using namespace Rcpp;


// [[Rcpp::export]]
double get_H_hat_cpp(double t, double x, double k, List valid_info) {
  NumericVector X_tx = valid_info["X_tx"];
  IntegerVector delta_c_tx = valid_info["delta_c_tx"];
  int n = X_tx.size();
  std::vector<double> u_seq;
  u_seq.reserve(n);
  for (int i = 0; i < n; ++i) {
    if (X_tx[i] <= k && delta_c_tx[i] == 1) {
      u_seq.push_back(X_tx[i]);
    }
  }
  if (u_seq.empty()) return 1.0;
  std::sort(u_seq.begin(), u_seq.end());
  u_seq.erase(std::unique(u_seq.begin(), u_seq.end()), u_seq.end());
  double H = 1.0;
  for (double u : u_seq) {
    int dN = 0;
    int risk_Y = 0;
    for (int i = 0; i < n; ++i) {
      if (X_tx[i] >= u) risk_Y++;
      if (X_tx[i] == u && delta_c_tx[i] == 1) dN++;
    }
    H *= (1.0 - double(dN) / risk_Y);
  }
  return H;
}





// [[Rcpp::export]]
double get_bivariate_S_Q_cpp(double t, double tprime,
                             double x, double xprime,
                             List valid_info_t_max,
                             List valid_info_tprime,
                             NumericVector Q_t_max,
                             NumericVector Q_tprime) {
  
  if (x==0 && xprime==0) return 1.0; 
  
  NumericVector X_tx_max = valid_info_t_max["X_tx"];
  NumericVector X_tx_prime = valid_info_tprime["X_tx"];
  IntegerVector delta_tx_max = valid_info_t_max["delta_tx"];
  IntegerVector delta_tx_prime = valid_info_tprime["delta_tx"];
  NumericVector T_tx_max = valid_info_t_max["T_tx"];
  NumericVector T_tx_prime = valid_info_tprime["T_tx"];
  
  int n = X_tx_max.size();
  double sum_inv = 0.0;
  for (int i = 0; i < n; ++i) {
    if (delta_tx_max[i] && delta_tx_prime[i]
          && Q_t_max[i] > x
          && Q_tprime[i] > xprime) {
          
          // which info to use for H-hat
          bool use_t    = (t+T_tx_max[i] >= tprime+T_tx_prime[i]);
      List infoi    = use_t ? valid_info_t_max : valid_info_tprime;
      double  tt    = use_t ? t : tprime;
      double  xx    = use_t ? x : xprime;
      double k      = use_t ? T_tx_max[i] : T_tx_prime[i];
      
      // call our C++ helper with the right info
      double H = get_H_hat_cpp(tt, xx, k, infoi);
      sum_inv += 1.0 / H;
    }
  }
  
  return sum_inv / n;
}

// [[Rcpp::export]]
double get_dMdM_cpp(double t, double tprime, double u, double uprime,
                    List valid_info_t_max, List valid_info_tprime) {
  
  IntegerVector delta_c_tx_max = valid_info_t_max["delta_c_tx"];
  IntegerVector delta_c_tx_prime = valid_info_tprime["delta_c_tx"];
  NumericVector X_tx_max   = valid_info_t_max["X_tx"];
  NumericVector X_tx_prime = valid_info_tprime["X_tx"];
  
  // Size of the t_max–at‐risk set
  int n_max = X_tx_max.size();
  // Get Y(t,t',x,x',u,u') first
  int YY = 0;
  for (int i = 0; i < n_max; ++i) {
    if (X_tx_max[i] >= u && X_tx_prime[i] >= uprime) 
      YY++;
  }
  if (YY == 0) return 0.0;
  
  // First term, needs t+u = t'+u'
  double first_term;
  if (std::abs((t+u) - (tprime+uprime)) > 1e-8) {
    first_term = 0.0;
  }
  else {
    int dNdN = 0;
    for (int i = 0; i < n_max; ++i) {
      if (X_tx_max[i] == u && X_tx_prime[i] == uprime && delta_c_tx_max[i] == 1 && delta_c_tx_prime[i] == 1) {
        dNdN++;
      }
    }
    first_term = double(dNdN)/double(n_max);
  }
  // Rcpp::Rcout << "first_term = " << first_term << "\n";
  // Compute dN^c(t,x,u), dN^c(t',x',u'), Y(t,x,u) and Y(t',x',u') first
  int dN_txu          = 0;
  int dN_txu_prime    = 0;
  int Y_txu           = 0;
  int Y_txu_prime     = 0;
  for (int i = 0; i < n_max; ++i) {
    if (X_tx_max[i]    >= u)      ++Y_txu;
    if (X_tx_prime[i]  >= uprime) ++Y_txu_prime;
    
    if (X_tx_max[i]    == u)      ++dN_txu;
    if (X_tx_prime[i]  == uprime) ++dN_txu_prime;
  }
  
  // Second term, needs t+u >= t'+u'
  double second_term;
  if ((t+u) < (tprime+uprime) || dN_txu_prime == 0 || Y_txu_prime == 0) {
    second_term = 0.0;
  }
  else {
    int dN_Yprime = 0;
    for (int i = 0; i < n_max; ++i) {
      if (X_tx_max[i] == u && delta_c_tx_max[i] == 1 && X_tx_prime[i] >= uprime) {
        dN_Yprime++;
      }
    }
    second_term = double(dN_Yprime)*double(dN_txu_prime)/(double(Y_txu_prime)*double(n_max));
  }
  //Rcpp::Rcout << "second_term = " << second_term << "\n";
  
  // Third term, needs t+u <= t'+u'
  double third_term;
  if ((t+u) > (tprime+uprime) || dN_txu== 0 || Y_txu == 0) {
    third_term = 0.0;
  }
  else {
    int dNprime_Y = 0;
    for (int i = 0; i < n_max; ++i) {
      if (X_tx_prime[i] == uprime && delta_c_tx_prime[i] == 1 && X_tx_max[i] >= u) {
        dNprime_Y++;
      }
    }
    third_term = double(dNprime_Y)*double(dN_txu)/(double(Y_txu)*double(n_max));
  }
  //Rcpp::Rcout << "third_term = " << third_term << "\n";
  
  // Fourth term
  double fourth_term;
  if (dN_txu == 0 || dN_txu_prime == 0 || Y_txu == 0 || Y_txu_prime == 0) {
    fourth_term = 0.0;
  }
  else {
    fourth_term = double(dN_txu)*double(dN_txu_prime)*double(YY)/(double(Y_txu)*double(Y_txu_prime)*double(n_max));
  }
  //Rcpp::Rcout << "fourth_term = " << fourth_term << "\n";
  
  
  return first_term-second_term-third_term+fourth_term;
}

// [[Rcpp::export]]
double get_G_hat_cpp(double t, double x, double u, List valid_info,
                     NumericVector Q_t) {
  
  NumericVector X_tx = valid_info["X_tx"];
  IntegerVector delta_tx = valid_info["delta_tx"];
  NumericVector T_tx = valid_info["T_tx"];
  int n = X_tx.size();
  
  // Quick checks
  bool anyT = false, anyQ = false, anyD = false;
  for(int i=0;i<n;i++){
    if (T_tx[i] >= u) anyT = true;
    if (Q_t[i]  >  x) anyQ = true;
    if (delta_tx[i] == 1) anyD = true;
  }
  if (!(anyT && anyQ && anyD)) return 0.0;
  
  // Identify idx where delta==1, Q>x and T>=u
  std::vector<int> idx;
  for(int i=0;i<n;i++){
    if (delta_tx[i]==1 && Q_t[i]>x && T_tx[i]>=u)
      idx.push_back(i);
  }
  if (idx.empty()) return 0.0;
  
  // Compute H_hat for each and form numerator
  double G_num = 0.0;
  for(int i : idx) {
    double k = T_tx[i];
    double H = get_H_hat_cpp(t, x, k, valid_info);
    G_num += 1.0 / H;
  }
  G_num /= n;
  
  // Compute Kaplan–Meier S(u) on {X_tx, delta_tx}
  // S_hat = product{v ≤ u} (1 − d(v)/Y(v))
  std::vector<double> events;
  events.reserve(n);
  for(int i=0;i<n;i++) if (delta_tx[i]==1) events.push_back(X_tx[i]);
  std::sort(events.begin(), events.end());
  events.erase(std::unique(events.begin(), events.end()), events.end());
  
  double S_hat = 1.0;
  for(double v : events) {
    if (v >= u) break;
    int Y = 0, d = 0;
    for(int i=0;i<n;i++){
      if (X_tx[i] >= v) Y++;
      if (X_tx[i] == v && delta_tx[i] == 1) d++;
    }
    S_hat *= (1.0 - double(d)/Y);
  }
  
  return G_num / S_hat;
}

// [[Rcpp::export]]
double get_var_inner_cpp(double t, double tprime, double x, double xprime, 
                         double u, double uprime, 
                         List valid_info_t, List valid_info_tprime, NumericVector Q_t, 
                         NumericVector Q_tprime, double bivariate_S_Q, 
                         double S_Q_max, double S_Q_prime) {
  
  double G_hat = get_G_hat_cpp(t,x,u,valid_info_t,Q_t);
  double G_hat_prime = get_G_hat_cpp(tprime,xprime,uprime,valid_info_tprime,Q_tprime);
  double H_hat = get_H_hat_cpp(t,x,u,valid_info_t);
  double H_hat_prime = get_H_hat_cpp(tprime,xprime, uprime,valid_info_tprime);
  return (bivariate_S_Q-S_Q_max*G_hat_prime-S_Q_prime*G_hat+G_hat*G_hat_prime)/(H_hat*H_hat_prime);
}


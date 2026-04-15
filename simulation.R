
seed <- 100
Rcpp::sourceCpp("helper.cpp")
source("functions.R")

check_qol_censor <- function(sample, qol) {
  perfect_qol <- which(qol == 1)
  full_follow <- which(sample$Ci == 36)
  ok_go <- sum(perfect_qol %in% full_follow) > 0
  return (ok_go)
}

repeat {
  set.seed(seed)
  
  source("Data.R")  
  
  n1 <- 200; n2 <- 200
  rho1 <- 0.5; rho2 <- 0.3; J1 <- rep(100, n1); J2 <- rep(100, n2)
  sd <- 1
  alpha <- 0.8
  rate_S1 <- 1/12; rate_S2 <- 1/12*alpha; rate_D1 <- 1/36; rate_D2 <- 1/36*alpha
  follow_full_p1 <- 0.3; follow_full_p2 <- 0.3
  followup <- 36 
  rate_C <- 1/60
  
  follow_full1 <- rbinom(n1, 1, follow_full_p1)
  follow_full2 <- rbinom(n2, 1, follow_full_p2)
  C1 <- runif(n1, 24, followup)*(1-follow_full1) + followup*follow_full1
  C2 <- runif(n2, 24, followup)*(1-follow_full2) + followup*follow_full2
  #C1 <- pmin(runif(n1, 24, followup),rexp(n1, rate_C))*(1-follow_full1) + followup*follow_full1
  #C2 <- pmin(runif(n2, 24, followup),rexp(n1, rate_C))*(1-follow_full2) + followup*follow_full2
  
  sample1_info <- gen_data(n1, J1, rho1, rho2, sd, rate_S1,rate_D1, C1, bad_qol=0.06)
  sample1 <- sample1_info[[1]]
  qol1 <- sample1_info[[2]]
  # The parameters are r_g, h_g, d_g, and b_g
  U1 <- gen_U(sample1, qol1, 0.5, 1, 0.05, 0.2)
  
  sample2_info <- gen_data(n2, J2, rho1, rho2, sd, rate_S2, rate_D2, C2, bad_qol=0.06)
  sample2 <- sample2_info[[1]]
  qol2 <- sample2_info[[2]]
  U2 <- gen_U(sample2, qol2, 0.8, 0.5, 0.025, 0.3) 
  ok1 <- check_qol_censor(sample1, qol1)
  ok2 <- check_qol_censor(sample2, qol2)
  ok_go  <- ok1 && ok2
  
  if (ok_go) {
    message("Success with seed = ", seed)
    break
  }
  
  # 4) otherwise bump seed and retry
  seed <- seed + 1000
}



# t_b+tau needs to be <= max follow-up
tau <- 3
t_vec <-seq(0, 32, by = 4)
x_seq <- seq(0, tau, 0.1)
S_estimator1_info <- get_S_estimator(t_vec, x_seq, tau, sample1, U1)
S_estimator1 <- S_estimator1_info$S_Q
info1 <- S_estimator1_info$info
S_estimator2_info <- get_S_estimator(t_vec, x_seq, tau, sample2, U2)
S_estimator2 <- S_estimator2_info$S_Q
info2 <- S_estimator2_info$info


mu_tau1 <- get_estimate_final(S_estimator1, t_vec) 
mu_tau2 <- get_estimate_final(S_estimator2, t_vec) 
start <- Sys.time()
t_integrals1 <- get_var_estimator_parallel(t_vec, x_seq, tau, sample1, U1, S_estimator1, n2/(n1+n2), info1)
t_integrals1[lower.tri(t_integrals1)] <- t(t_integrals1)[lower.tri(t_integrals1)]
row1  <- apply(t_integrals1 , 1, trapz, x = t_vec)  
var1 <- as.numeric(trapz(t_vec, row1)) 
end <- Sys.time()
t_integrals2 <- get_var_estimator_parallel(t_vec, x_seq, tau, sample2, U2, S_estimator2, n2/(n1+n2), info2)
t_integrals2[lower.tri(t_integrals2)] <- t(t_integrals2)[lower.tri(t_integrals2)]
row2  <- apply(t_integrals2 , 1, trapz, x = t_vec)  
var2 <- as.numeric(trapz(t_vec, row2)) 



result <- t(c(mu_tau1, mu_tau2, var1, var2))
t_integrals <- rbind(t_integrals1, t_integrals2)

mu_t_tau1 <- rep(NA, nrow(S_estimator1))
for (i in 1:nrow(S_estimator1)) {
  mu_t_tau1[i] <- trapz(as.numeric(colnames(S_estimator1)),as.numeric(S_estimator1[i,]))
}

mu_t_tau2 <- rep(NA, nrow(S_estimator2))
for (i in 1:nrow(S_estimator2)) {
  mu_t_tau2[i] <- trapz(as.numeric(colnames(S_estimator2)),as.numeric(S_estimator2[i,]))
}

mu_t_tau_combined <- rbind(mu_t_tau1, mu_t_tau2)

# 1000 obs: 6.153167 hours per group
# 500 obs: 31.5173 mins per group
# 200 obs: 2.555598 mins per group


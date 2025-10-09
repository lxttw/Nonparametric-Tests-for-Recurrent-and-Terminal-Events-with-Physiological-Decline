library(metamisc)
library(reshape2)
library(tidyr)
library(mhazard)
library(survival)
library(pracma)
library(foreach)
library(doParallel)
library(Rcpp)
library(RhpcBLASctl)
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)
seed <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")) + 20000
job_id <- as.integer(Sys.getenv("SLURM_ARRAY_JOB_ID"))
Rcpp::sourceCpp("helper.cpp")

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
  rho1 <- 0.3; rho2 <- 0.1; J1 <- rep(100, n1); J2 <- rep(100, n2)
  sd <- 1
  alpha <- 0.9
  rate_S1 <- 1/12; rate_S2 <- 1/12*alpha; rate_D1 <- 1/36; rate_D2 <- 1/36*alpha
  follow_full_p1 <- 0.3; follow_full_p2 <- 0.3
  followup <- 36 
  
  follow_full1 <- rbinom(n1, 1, follow_full_p1)
  follow_full2 <- rbinom(n2, 1, follow_full_p2)
  C1 <- runif(n1, 24, followup)*(1-follow_full1) + followup*follow_full1
  C2 <- runif(n2, 24, followup)*(1-follow_full2) + followup*follow_full2
  
  sample1_info <- gen_data(n1, J1, rho1, rho2, sd, rate_S1,rate_D1, C1, bad_qol=0.06)
  sample1 <- sample1_info[[1]]
  qol1 <- sample1_info[[2]]
  U1 <- gen_U(sample1, qol1, 0.7, 0.5, 0.03, 0.3)
  
  sample2_info <- gen_data(n2, J2, rho1, rho2, sd, rate_S2, rate_D2, C2, bad_qol=0.06)
  sample2 <- sample2_info[[1]]
  qol2 <- sample2_info[[2]]
  U2 <- gen_U(sample2, qol2, 0.7, 0.5, 0.03, 0.3) 
  
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


get_T_L_eta_vec <- function(t, events_list, L) {
  vapply(events_list, function(ev) {
    events_before_t <- findInterval(t, ev)
    # the very next event time (or Inf if none)
    next_ev <- if (events_before_t < length(ev)) ev[events_before_t + 1] else Inf
    # return min(next_ev, L)
    pmin(next_ev, L)
  }, numeric(1), USE.NAMES = FALSE)
}

get_Q_vec <- function(t, times, U, T_Ls) {
  n <- nrow(U)
  M <- length(times)
  
  # Find the index k so that times[k] <= t < times[k+1]
  k <- findInterval(t, times)
  if (k == 0 || k == M) stop("t out of bounds")
  # Build the truncated times and pre-compute diffs
  times_new <- times[k:M]                # length M-k+1
  times_diff <- diff(times_new)           # length M-k
  
  # Subset U
  U_sub <- U[, k:M, drop = FALSE]                # n x (M-k+1)
  U_new <- U_sub[, 1:(ncol(U_sub)-1), drop = FALSE]   # n x (M-k)
  U_last <- U_sub                               # same dim as U_sub
  
  # Compute cumulative area since t for each row
  # trunk areas: time segment * height, gives an n x (M-k) matrix
  trunk <- U_new * matrix(times_diff, nrow = n, ncol = length(times_diff), byrow = TRUE)
  # prepend zero and do row‐wise cumsum → an n x (M-k+1) matrix
  U_cum <- cbind(0, t(apply(trunk, 1, cumsum)))
  # subtract the partial at start: U_new[,1] * (t - times_new[1])
  offset <- U_new[,1] * (t - times_new[1])
  U_cum <- U_cum - offset
  
  # compute Q_i for each subject
  closest_event_before <- findInterval(T_Ls, times_new)     
  Q_vec <- U_cum[cbind(seq_len(n), closest_event_before)]+
    (T_Ls-times_new[closest_event_before])*U_last[cbind(seq_len(n), closest_event_before)]
  
  # return both U_cum and Q_vec for later use
  list(Q_t = Q_vec, U_cum = U_cum, U_new = U_new, times_new = times_new)
}

get_s_vec <- function(t, x, U_cum, U_new, times_new) {
  n <- nrow(U_cum)
  # compute the s_i (first time U_cum >= x)
  # by finding for each row the first column j with U_cum >= x
  j_s <- apply(U_cum, 1, function(u) {
    which.max(u >= x) })  # returns 1 if no such u
  
  s_vec <- numeric(n)
  if (x == 0) {
    s_vec[] <- 0
  } else {
    # those rows where U_cum never reaches x get Inf
    never <- U_cum[cbind(seq_len(n), j_s)] < x
    s_vec[never] <- Inf
    
    # for others, do linear interp
    ok <- !never
    js <- j_s[ok]
    # prev time = times_new[js-1], cum = U_cum[ok, js-1], slope = U_new[ok, js-1]
    numer <- x - U_cum[cbind(which(ok), js-1)]
    denom <- U_new[cbind(which(ok), js-1)]
    s_abs <- times_new[js-1] + numer/denom
    s_vec[ok] <- s_abs - t
  }
  return(s_vec)
}




# All these are only valid if the subject's death and censoring time > t
get_info_vec <- function(t, x, C_valid, U_cum, U_new, times_new, T_L_eta) {
  s_tx_vec <- get_s_vec(t,x,U_cum, U_new,times_new)
  T_tx <- pmin(T_L_eta-t, s_tx_vec)
  X_tx <- pmin(T_tx, C_valid-t)
  delta_tx <- as.numeric(T_tx + t <= C_valid)
  delta_c_tx <- as.numeric(C_valid <= T_tx + t)
  res <- list(T_tx, X_tx, delta_tx, delta_c_tx)
  names(res) <- c("T_tx", "X_tx", "delta_tx", "delta_c_tx")
  return (res)
}


get_S_Q_vec <- function(t, x, L, C_valid, Q_t, U_cum, U_new, times_new, T_L_eta, n_valid) {
  # probability of having cumulative utility > 0 is 1
  
  # Get info
  info <- get_info_vec(t, x, C_valid, U_cum, U_new, times_new, T_L_eta)
  if (x == 0) {
    return(list(log_n_valid = n_valid, log_info = info,
                S_Q = 1))
  }
  T_tx <- info$T_tx
  delta_tx  <- info$delta_tx
  
  # Identify the subset who both have an event (delta==1) and Q_t > x
  subset_valid  <- which(delta_tx == 1 & Q_t > x)
  if (length(subset_valid) == 0) {
    return(list(log_n_valid = n_valid, log_info = info,
                S_Q = 0))
  }
  
  # Compute H_hat for each of those
  H_vals <- vapply(subset_valid,
                   function(i) get_H_hat_cpp(t, x, T_tx[i], info),
                   numeric(1),
                   USE.NAMES = FALSE)
  
  # Sum their contributions and normalize each contributes 1 / H_hat
  S_Q <- sum(1/H_vals)/n_valid
  return(list(log_n_valid = n_valid, log_info = info,
              S_Q = S_Q))
}



get_events_list <- function(events, end_time, n, len_events) {
  result <- list()
  # exclude events after censoring/death
  for (i in 1:n) {
    events_i_temp <- events[i,2:(len_events-1)] # remove ID (1st column) and Ci (last column)
    events_i <- events_i_temp[0 <= events_i_temp & events_i_temp <= end_time[i]]
    events_i <- sort(events_i)
    result[[i]] <- events_i
  }
  return (result)
}



get_S_estimator <- function(t_vec, x_seq, tau, events, U) {
  n <- nrow(events)
  D <- events$Di
  C <- events$Ci
  end_time <- pmin(D, C)
  #end_time <- C
  len_events <- ncol(events)
  events_list <- get_events_list(events, end_time, n, len_events)
  times_U <- as.numeric(colnames(U))
  S_Q <- matrix(NA, nrow = length(t_vec), ncol = length(x_seq))
  rownames(S_Q) <- t_vec
  colnames(S_Q) <- x_seq
  info_list <- vector("list", length(t_vec))
  names(info_list) <- as.character(t_vec)
  for (i in 1:length(t_vec)) {
    t <- t_vec[i]
    L <- t+tau
    # Find everyone still at risk at time t
    valid_index <- which(end_time > t)
    n_valid <- length(valid_index)
    if (n_valid == 0) {
      # fill entire row with 0
      S_Q[i, ] <- 0
      info_list[[i]] <- rep(list(NULL), length(x_seq))
      names(info_list[[i]]) <- as.character(x_seq)
      next
    }
    events_list_valid <- events_list[valid_index]
    T_L_eta <- get_T_L_eta_vec(t, events_list_valid, L)
    U_valid <- U[valid_index,]
    C_valid <- C[valid_index]
    Q_info <- get_Q_vec(t, times_U, U_valid, T_L_eta)
    info_sub <- vector("list", length(x_seq))
    names(info_sub) <- as.character(x_seq)
    for (j in seq_along(x_seq)) {
      x  <- x_seq[j]
      info_sub[[j]] <- get_S_Q_vec(t = t, x = x, L = L, C_valid  = C_valid, 
                                   Q_t = Q_info$Q_t, U_cum = Q_info$U_cum,
                                   U_new = Q_info$U_new, times_new= Q_info$times_new,
                                   T_L_eta = T_L_eta, n_valid = n_valid)
      S_Q[i, j] <- info_sub[[j]]$S_Q
    }
    info_list[[i]] <- info_sub
  }
  # Return both
  result <- list(S_Q = S_Q, info = info_list)
  return (result)
}



get_var_outer_parallel <- function(t, tprime, x, xprime, C_valid_tprime, U_cum_t_max, 
                                   U_new_t_max, times_new_t_max,  
                                   T_L_t_max, S_estimator, S_Q_t_max,Q_t, Q_t_max,Q_tprime, 
                                   p_g_prime, p_g_t_min, valid_info_t, valid_info_tprime, tau) {
  valid_info_t_max <- get_info_vec(t, x, C_valid_tprime, U_cum_t_max, 
                                   U_new_t_max, times_new_t_max, T_L_t_max)
  bivariate_S_Q <- get_bivariate_S_Q_cpp(t,tprime, x, xprime, valid_info_t_max, valid_info_tprime,
                                         Q_t_max, Q_tprime)
  S_estimator_col <- as.numeric(colnames(S_estimator))
  S_estimator_row <- as.numeric(rownames(S_estimator))
  S_estimator_max_col <- as.numeric(colnames(S_Q_t_max))
  # Because of issue of floats
  S_Q_max <- S_Q_t_max[1, S_estimator_max_col == round(x, 5)]
  #S_Q <- S_estimator[S_estimator_row == round(t,5), S_estimator_col == round(x,5)]
  S_Q_prime <- S_estimator[S_estimator_row == round(tprime,5), S_estimator_col == round(xprime,5)]
  S_Q_S_Q_prime <- S_Q_max*S_Q_prime
  if (abs(t-tprime) > tau) return (p_g_prime/p_g_t_min*(bivariate_S_Q - S_Q_S_Q_prime))
  # u has to be one of the X(t,x) with delta^c(t,x) = 1 and u' has to be one of the X(t',x') with delta^c(t',x') = 1
  u_seq <- sort(unique(valid_info_t_max$X_tx[valid_info_t_max$delta_c_tx == 1]))
  u_seq <- u_seq[u_seq <= tau]
  uprime_seq <- sort(unique(valid_info_tprime$X_tx[valid_info_tprime$delta_c_tx == 1]))
  uprime_seq <- uprime_seq[uprime_seq <= tau]
  if(length(u_seq) == 0 | length(uprime_seq) == 0) return (p_g_prime/p_g_t_min*(bivariate_S_Q - S_Q_S_Q_prime))
  integrals <- matrix(0, nrow = length(u_seq), ncol = length(uprime_seq))
  for (i in 1:length(u_seq)) {
    u <- u_seq[i]
    for (j in 1:length(uprime_seq)) {
      uprime <- uprime_seq[j]
      dMdM_ij <- get_dMdM_cpp(t, tprime,u, uprime, valid_info_t_max, valid_info_tprime)
      if (dMdM_ij > 0) {
        integral_inner_ij <- get_var_inner_cpp(t, tprime, x, xprime, u,uprime, 
                                              valid_info_t, valid_info_tprime, Q_t, Q_tprime, 
                                              bivariate_S_Q, S_Q_max, S_Q_prime)
        integrals[i,j] <- integral_inner_ij*dMdM_ij
      }
    }
  }
  res <- p_g_prime/p_g_t_min*(bivariate_S_Q - S_Q_S_Q_prime + sum(integrals))
  return (res)
}

get_var_estimator_parallel <- function(t_vec, x_seq, tau, events, U, S_estimator, 
                                       p_g_prime, info, n_cores = 5) {
  
  # Set up cluster
  cl <- makeCluster(n_cores, type = "FORK")
  registerDoParallel(cl)
  on.exit({
    stopCluster(cl)
    registerDoSEQ()
  })
  
  
  # Precompute all the pieces that don't depend on i
  n <- nrow(events)
  D <- events$Di
  C <- events$Ci
  end_time <- pmin(D, C)
  #end_time <- C
  len_events <- ncol(events)
  events_list<- get_events_list(events, end_time, n, len_events)
  times_U <- as.numeric(colnames(U))
  # Parallelize over i = 1:length(t_vec)
  t_integrals <- foreach(i= seq_along(t_vec),
                         .combine = rbind,
                         .packages = c("pracma","survival","mhazard")) %dopar% {
                           # prepare a row of NAs
                           row_i <- rep(NA_real_, length(t_vec))
                           
                           t <- t_vec[i]
                           # skip if nobody at risk
                           if (sum(end_time > t) == 0) return(row_i)
                           # precompute everything for this t
                           valid_index_t <- which(end_time > t)
                           events_list_t <- events_list[valid_index_t]
                           L_t <- t + tau
                           T_L_t <- get_T_L_eta_vec(t, events_list_t, L_t)
                           U_t <- U[valid_index_t, , drop=FALSE]
                           C_t <- C[valid_index_t]
                           Q_t_info <- get_Q_vec(t, times_U, U_t, T_L_t)
                           n_valid_t <- length(valid_index_t)
                           p_g_t_min <- n_valid_t/n
                           # inner loop over j >= i
                           for (j in i:length(t_vec)) {
                             tprime <- t_vec[j]
                             # skip if nobody at risk at tprime
                             if (sum(end_time > tprime) == 0) next
                             valid_index_tprime <- which(end_time > tprime)
                             n_valid_tprime <- length(valid_index_tprime)
                             events_list_tprime <- events_list[valid_index_tprime]
                             L_tprime <- tprime + tau
                             T_L_t_max <- get_T_L_eta_vec(t, events_list_tprime, L_t)
                             T_L_tprime <- get_T_L_eta_vec(tprime, events_list_tprime, L_tprime)
                             U_valid_tprime <- U[valid_index_tprime,]
                             C_valid_tprime <- C[valid_index_tprime]
                             Q_info_t_max <- get_Q_vec(t, times_U, U_valid_tprime, T_L_t_max)
                             Q_info_tprime <- get_Q_vec(tprime, times_U, U_valid_tprime, T_L_tprime)
                             S_Q_t_max <- matrix(NA, nrow = 1, ncol = length(x_seq))
                             colnames(S_Q_t_max) <- x_seq
                             info_max <- vector("list", length(x_seq))
                             for (k in seq_along(x_seq)) {
                               x  <- x_seq[k]
                               info_max[[k]] <- get_S_Q_vec(t = t, x = x, L = L, C_valid = C_valid_tprime, 
                                                            Q_t = Q_info_t_max$Q_t, U_cum = Q_info_t_max$U_cum,
                                                            U_new = Q_info_t_max$U_new, times_new= Q_info_t_max$times_new,
                                                            T_L_eta = T_L_t_max, n_valid = n_valid_tprime)
                               S_Q_t_max[1,k] <- info_max[[k]]$S_Q
                             }
                             # double‐integral over x_seq × x_seq
                             inner_over_x <- vapply(x_seq, function(x) {
                               inner_xp <- vapply(x_seq, function(xprime) {
                                 # extract the correct pre‐logged info object:
                                 valid_info_t <- info[[as.character(t)]][[as.character(x)]]$log_info
                                 valid_info_tprime <- info[[as.character(tprime) ]][[as.character(xprime)]]$log_info
                                 get_var_outer_parallel(t, tprime, x, xprime,
                                                        C_valid_tprime, U_cum_t_max=Q_info_t_max$U_cum, 
                                                        U_new_t_max=Q_info_t_max$U_new, times_new_t_max=Q_info_t_max$times_new,  
                                                        T_L_t_max, S_estimator, S_Q_t_max, Q_t=Q_t_info$Q_t,
                                                        Q_t_max=Q_info_t_max$Q_t, Q_tprime=Q_info_tprime$Q_t,
                                                        p_g_prime, p_g_t_min, valid_info_t, valid_info_tprime, tau)}, 
                                 numeric(1), USE.NAMES = FALSE)
                               
                               trapz(x_seq, inner_xp)}, numeric(1), USE.NAMES = FALSE)
                             # store the outer integral
                             row_i[j] <- trapz(x_seq, inner_over_x)}
                           row_i}
  t_integrals
}

get_estimate_final <- function(S_estimator, t_vec) {
  mu_t_tau <- rep(NA, nrow(S_estimator))
  for (i in 1:nrow(S_estimator)) {
    mu_t_tau[i] <- trapz(as.numeric(colnames(S_estimator)),as.numeric(S_estimator[i,]))
  }
  return (trapz(t_vec, mu_t_tau))
}


# t_b+tau needs to be less than 36
tau <- 6
t_vec <-seq(0, 28, by = 4)
x_seq <- seq(0, tau, 0.1)
S_estimator1_info <- get_S_estimator(t_vec, x_seq, tau, sample1, U1)
S_estimator1 <- S_estimator1_info$S_Q
info1 <- S_estimator1_info$info
S_estimator2_info <- get_S_estimator(t_vec, x_seq, tau, sample2, U2)
S_estimator2 <- S_estimator2_info$S_Q
info2 <- S_estimator2_info$info


mu_tau1 <- get_estimate_final(S_estimator1, t_vec) #100.8395, 100.211
mu_tau2 <- get_estimate_final(S_estimator2, t_vec) #128.1823; 125.007

t_integrals1 <- get_var_estimator_parallel(t_vec, x_seq, tau, sample1, U1, S_estimator1, n2/(n1+n2), info1)
t_integrals1[lower.tri(t_integrals1)] <- t(t_integrals1)[lower.tri(t_integrals1)]
row1  <- apply(t_integrals1 , 1, trapz, x = t_vec)  
var1 <- as.numeric(trapz(t_vec, row1)) #1212.9027545, 665.3385
t_integrals2 <- get_var_estimator_parallel(t_vec, x_seq, tau, sample2, U2, S_estimator2, n2/(n1+n2), info2)
t_integrals2[lower.tri(t_integrals2)] <- t(t_integrals2)[lower.tri(t_integrals2)]
row2  <- apply(t_integrals2 , 1, trapz, x = t_vec)  
var2 <- as.numeric(trapz(t_vec, row2)) #1119.15332259, 442.4996



result <- t(c(mu_tau1, mu_tau2, var1, var2))

outdir    <- file.path("output", job_id)
#if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
outfile   <- file.path(outdir, paste0("result", seed, ".csv"))
write.csv(result, outfile, row.names = FALSE)
#write.csv(result, paste0("output/result", seed, ".csv"), row.names = F)

t_integrals <- rbind(t_integrals1, t_integrals2)

outdir.supp    <- file.path("output_supp", job_id)
if (!dir.exists(outdir.supp)) dir.create(outdir.supp, recursive = TRUE)
outfile.supp   <- file.path(outdir.supp, paste0("result", seed, ".csv"))
write.csv(t_integrals, outfile.supp)



mu_t_tau1 <- rep(NA, nrow(S_estimator1))
for (i in 1:nrow(S_estimator1)) {
  mu_t_tau1[i] <- trapz(as.numeric(colnames(S_estimator1)),as.numeric(S_estimator1[i,]))
}

mu_t_tau2 <- rep(NA, nrow(S_estimator2))
for (i in 1:nrow(S_estimator2)) {
  mu_t_tau2[i] <- trapz(as.numeric(colnames(S_estimator2)),as.numeric(S_estimator2[i,]))
}

mu_t_tau_combined <- rbind(mu_t_tau1, mu_t_tau2)
rownames(mu_t_tau_combined) <- c("mu_t_tau1", "mu_t_tau2")
outdir.est    <- file.path("output_estimates", job_id)
if (!dir.exists(outdir.est)) dir.create(outdir.est, recursive = TRUE)
outfile.est  <- file.path(outdir.est, paste0("result", seed, ".csv"))
write.csv(mu_t_tau_combined, outfile.est)


library(metamisc)
library(reshape2)
library(tidyr)
library(mhazard)
library(survival)
library(pracma)
library(foreach)
library(doParallel)
library(Rcpp)
Rcpp::sourceCpp("helper.cpp")
tol <- 1e-10


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
  # Use tol for floating issue
  subset_valid  <- which(delta_tx == 1 & Q_t - x > tol)
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
sum_function=function(X)
{
  temp=array(X,c(length(X),length(X)))
  lowerTriangle(temp,diag=F)=0
  apply(temp,2,sum)
}


format_data_PM=function(sample)
{
  ############format observed data for Andersen-Gill method
  #load data
  X <- pmin(sample$Di, sample$Ci)
  delta <- as.numeric(sample$Di < sample$Ci)
  events <- sample
  events$Di <- NULL
  n=length(X)
  b=length(t)
  Z <- get_events_list(events, X, n, ncol(events))
  max_len <- max(lengths(Z))
  Z_star <- t(sapply(Z, function(x) {
    length(x) <- max_len
    x
  }))
  X_Z_star=cbind(Z_star,X)
  T_start=NULL
  T_stop=NULL
  status=NULL
  ID_AG=NULL
  for(i in 1:n)
  {
    T_stop_subject=as.numeric(na.omit(X_Z_star[i,]))
    n_reccurent_subject=length(T_stop_subject)-1
    if(n_reccurent_subject>0)
    {
      T_start_subject=c(0,T_stop_subject[1:n_reccurent_subject])
    }
    if(n_reccurent_subject==0)
    {
      T_start_subject=0
    }    
    status_subject=c(rep(1,n_reccurent_subject),delta[i])
    
    T_start=c(T_start,T_start_subject)
    T_stop=c(T_stop,T_stop_subject)
    status=c(status,status_subject)  
    ID_AG=c(ID_AG,rep(i,n_reccurent_subject+1))
    
    if(length(T_start)!=length(T_stop)){print(i)}
  }
  list(ID=ID_AG,T_start=T_start,T_stop=T_stop,status=status)
}

PM_test_stat=function(data1_format,data2_format)
{
  #Andersen-Gill model based two-sample test
  n1=length(data1_format$ID)
  treatment=c(rep(0,length(data1_format$ID)),rep(1,length(data2_format$ID)))
  T_start=c(data1_format$T_start,data2_format$T_start)
  T_stop=c(data1_format$T_stop,data2_format$T_stop)
  status=c(data1_format$status,data2_format$status)
  id=c(data1_format$ID,data2_format$ID+n1)
  ag_model=coxph(Surv(T_start, T_stop, status) ~ treatment + cluster(id))
  test_stat_p=1-pchisq(ag_model$wald.test,1)
  list(test_stat_p=test_stat_p)
}

format_data_GL=function(sample, time)
{
  ############format observed data for Ghosh and Lin method
  #load data
  X <- pmin(sample$Di, sample$Ci)
  delta <- as.numeric(sample$Di < sample$Ci)
  events <- sample
  events$Di <- NULL
  n=length(X)
  Z <- get_events_list(events, X, n, ncol(events))
  max_len <- max(lengths(Z))
  Z_star <- t(sapply(Z, function(x) {
    length(x) <- max_len
    x
  }))
  n_time=length(time)
  time_array=t(array(time,c(n_time,n)))
  delta_array=array(delta,c(n,n_time))
  dN_D_i=array(as.numeric(X==time_array & delta_array==1),c(n,n_time))
  dN_D=apply(dN_D_i,2,sum)
  Y_i=array(as.numeric(X>=time_array),c(n,n_time))
  Y=apply(Y_i,2,sum)
  time_array2<<-t(array(time,c(n_time,ncol(Z_star))))
  dN_i=t(apply(Z_star,1,get_dN_for_i))
  
  dN=apply(dN_i,2,sum)
  dR_hat=dN/Y
  dch_hat=dN_D/Y
  dM_i=dN_i - Y_i*t(array(dR_hat,c(n_time,n)))
  dM_D_i=dN_D_i - Y_i*t(array(dch_hat,c(n_time,n)))
  
  #calculate ch estimate
  temp=t(array(dch_hat,c(n_time,n_time)))
  upperTriangle(temp,diag=F)=0
  ch_hat=apply(temp,1,sum)
  
  #calculate survival estimate
  s_hat=exp(-ch_hat)
  
  #calculate mu estimate
  dmu_hat=s_hat*dR_hat
  temp=t(array(dmu_hat,c(n_time,n_time)))
  upperTriangle(temp,diag=F)=0
  mu_hat=apply(temp,1,sum)  
  
  #calculate psi estimate
  temp=dM_D_i*t(array(n/Y,c(n_time,n)))
  temp2=t(apply(temp,1,sum_function))
  dpsi_i=dM_i*t(array(s_hat*n/Y,c(n_time,n)))-temp2*t(array(dmu_hat,c(n_time,n)))
  
  #calculate last observed event time
  tau=max(X*delta)
  
  list(Y=Y,tau=tau,dmu_hat=dmu_hat,dch_hat=dch_hat,dpsi_i=dpsi_i,dM_D_i=dM_D_i)
}

get_dN_for_i=function(X1)
{
  temp=(X1==time_array2)
  apply(temp,2,sum,na.rm=T)
}

GL_test_stat=function(p_GL,time,data1_format,data2_format)
{
  #Load data
  Y1=data1_format$Y
  Y2=data2_format$Y
  tau1=data1_format$tau
  tau2=data2_format$tau
  dmu_hat1=data1_format$dmu_hat
  dmu_hat2=data2_format$dmu_hat
  dch_hat1=data1_format$dch_hat
  dch_hat2=data2_format$dch_hat
  dpsi_i1=data1_format$dpsi_i
  dpsi_i2=data2_format$dpsi_i
  dM_D_i1=data1_format$dM_D_i
  dM_D_i2=data2_format$dM_D_i
  
  n_time=length(time)
  
  n1=nrow(dpsi_i1)
  n2=nrow(dpsi_i2)
  
  tau=min(tau1,tau2)
  time_used=rep(1,length(time))
  time_used[time>tau]=0
  
  #calculate test stats
  K_LR=Y1*Y2*(n1+n2)/((Y1+Y2)*n1*n2)
  Q_LR=sum((K_LR*(dmu_hat1-dmu_hat2))*time_used)
  Q_D=sum((K_LR*(dch_hat1-dch_hat2))*time_used)
  Q_WCT=p_GL*Q_LR + (1-p_GL)*Q_D
  
  #calculate variance of Q_WCT
  rho_hat1=n1/(n1+n2)
  rho_hat2=n2/(n1+n2)
  U_i1=apply((p_GL*dpsi_i1 + (1-p_GL)*dM_D_i1*t(array(n1/Y1,c(n_time,n1))))*t(array(K_LR,c(n_time,n1))),1,sum)
  U_i2=apply((p_GL*dpsi_i2 + (1-p_GL)*dM_D_i2*t(array(n2/Y2,c(n_time,n2))))*t(array(K_LR,c(n_time,n2))),1,sum)
  var_Q_WCT=rho_hat2*mean(U_i1^2) + rho_hat1*mean(U_i2^2)
  
  #calculate test statistic that has standard normal distribution under H0
  test_stat=sqrt(n1*n2/(n1+n2))*Q_WCT/sqrt(var_Q_WCT)
  test_stat_p=2*(1-pnorm(abs(test_stat)))
  list(test_stat_p=test_stat_p)
}



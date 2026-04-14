
library(metamisc)
library(reshape2)
library(tidyr)
library(mvtnorm)
tol <- 1e-10


# Simulate dataset

generate_corr <- function(n, rho1, rho2) {
  res <- diag(n)
  if (n == 1){
    return (res)
  }
  res[1:(n-1),n] <- rho2
  res[n, 1:(n-1)] <- rho2
  res[res == 0] <- rho1
  return (res)
}


gen_data <- function(n, J, rho1, rho2, sd, rate_S, rate_D,C, bad_qol) {
  qol <- rbinom(n, 1, bad_qol)  
  for (i in seq_len(n)) {
    if (qol[i]==1) {
      Ji_minus1 <- J[i]-1
      res <- c(i, seq(1001, 1000+Ji_minus1, 1), 36.6, C[i])
      names(res) <- c("ID",paste0("T", 1:Ji_minus1), "Di", "Ci")
      if (i == 1) {
        res_final <- as.data.frame(t(res))
      }
      else{
        res_final <- rbind(res_final, res)
      }
    } else {
    Ji_minus1 <- J[i]-1
    corr <- generate_corr(Ji_minus1+1, rho1, rho2)
    cov <- cor2cov(rep(sd, Ji_minus1+1), corr)
    mean <- rep(0, Ji_minus1+1)
    UV <- rmvnorm(1, mean,cov)[1,]
    UV_prime <- pnorm(UV, mean = 0, sd = 1)
    len <- length(UV_prime)
    Di <- qexp(UV_prime[len], rate = rate_D) # Check if the mean is correct
    Si <- qexp(UV_prime[1:len-1], rate = rate_S)
    Ti <- cumsum(Si)
    res <- c(i, Ti, Di, C[i]) 
    names(res) <- c("ID",paste0("T", 1:length(Ti)), "Di", "Ci")
    if (i == 1) {
      res_final <- as.data.frame(t(res))
    }
    else{
      res_final <- rbind(res_final, res)
    }
    }
  }
  return (list(res_final, qol))
}




map_single_time_U <- function(time, events_i, base_curr, drop, lag) {
  
  events_after <- events_i[events_i >= time]
  # Instead of using NA to indicate that there's no event before/after, use a super large number
  closest_event_after <- ifelse(length(events_after) == 0, Inf, 
                                min(events_after))
  events_before <- events_i[events_i <= time]
  closest_event_before <- ifelse(length(events_before) == 0, -Inf, max(events_before))
  closest_event_before <- ifelse(closest_event_before < 0, -Inf, closest_event_before)
  gap_forward <- closest_event_after - time
  gap_backward <- time - closest_event_before
  
  utility <- ifelse(gap_forward <= lag | gap_backward <= lag, base_curr*drop, base_curr)
  return (utility)
}




map_vec_time_U <- function(events_i, base_initial_i, times, drop, lag, perm_drop, base_min) {
  base_i <- rep(base_initial_i, length(times))
  events_i <- sort(events_i[events_i <= max(times) & events_i >= 0], decreasing = F)
  for (evt in events_i) {
    times_after <- which(times >= evt + lag)
    time_before <- max(which(times < evt + lag))
    base_i[times_after] <- base_i[time_before]-perm_drop
  }
  
  utilities_i <- rep(NA, length(times))
  for (j in 1:length(times)) {
    time <- times[j]
    base_curr <- max(base_i[j], base_min)
    utilities_i[j] <- map_single_time_U(time, events_i, base_curr, drop, lag)
  }

  return(utilities_i)
}



gen_U <- function(sample, qol, drop, lag, perm_drop, base_min) {
  n <- nrow(sample)
  cols <- ncol(sample)
  events <- sample[,2:(cols-1)]
  max_censor <- max(sample$Ci)
  times <- seq(0, max_censor, by = 0.05)
  utilities <- c(0.6, 0.7, 0.8, 0.9, 1)
  probs <- c(0.15, 0.15, 0.15, 0.15, 0.40)
  base_utilities <- sample(utilities, size = n, replace = TRUE, prob = probs)
  U_table <- matrix(NA_real_, nrow = n, ncol = length(times))
  for(i in seq_len(n)) {
    events_i <- as.numeric(events[i, ])   
    base_initial_i <- base_utilities[i]         
    U_table[i,] <- map_vec_time_U(events_i, base_initial_i, times, drop, lag, perm_drop,base_min)
  }
  
  colnames(U_table) <- times
  U_table[qol==1,] <- matrix(0.1, ncol = ncol(U_table), nrow = sum(qol))
  return(U_table)
}






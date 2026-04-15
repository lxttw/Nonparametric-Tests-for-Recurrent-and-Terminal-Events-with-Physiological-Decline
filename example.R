Rcpp::sourceCpp("helper.cpp")
source("functions.R")


map_utility <- function(fev1pp) {
  utility <- ifelse(fev1pp >= 80, 1, fev1pp)
  utility <- ifelse(fev1pp < 80 & fev1pp >= 50, 0.8, utility)
  utility <- ifelse(fev1pp < 50 & fev1pp >= 30, 0.5, utility)
  utility <- ifelse(fev1pp < 30, 0.3, utility)
  return (utility)
}

get_events_for_severity <- function(events, severity, select) {
  res <- matrix(data = NA, nrow = nrow(events), ncol = ncol(events)-1)
  for (i in 1:nrow(events)) {
    if (severity[i,1] != events[i,1]) {
      stop("ID Not Matching")
    }
    events_i <- as.numeric(events[i,-1])
    severity_i <- as.numeric(severity[i,-1])
    if (select == "severe") {
      events_new_i <- ifelse(severity_i %in% c(3,4), events_i, NA)
    }
    else if (select == "severe_mod") {
      events_new_i <- ifelse(severity_i %in% c(2,3,4), events_i, NA)
    }
    else {
      events_new_i <- events_i
    }
    res[i,1:ncol(res)] <- events_new_i
  }
  res <- as.data.frame(res)
  res <- cbind(events$ID, res)
  colnames(res) <- colnames(events)
  
  return (res)
}


get_utility_drop <- function(U_curr, utilities) {
  index_curr <- findInterval(U_curr, utilities)
  index_drop <- pmax(index_curr-1, 1)
  return (utilities[index_drop])
}

get_utility_for_events <- function(events_list, U, utilities) {
  U_new <- U
  for (i in 1:nrow(events)) {
    events_i <- events_list[[i]]
    times <- as.numeric(colnames(U))
    if (length(events_i) == 0) next
    starts <- findInterval(events_i, times)
    U_curr <- U[i,starts]
    U_drop <- get_utility_drop(U_curr, utilities)
    for (j in 1:length(starts)) {
      st <- starts[j]
      en <- min(st + 14, ncol(U))
      prev <- max(st - 14, 1)
      U_new[i, st:en] <- U_drop[j]
      U_new[i, prev:st] <- U_drop[j]
    }
  }
  U_new <- as.data.frame(U_new)
  return (U_new)
}

get_estimate_final_new <- function(S_estimator, t_vec) {
  mu_t_tau <- rep(NA, nrow(S_estimator))
  for (i in 1:nrow(S_estimator)) {
    mu_t_tau[i] <- trapz(as.numeric(colnames(S_estimator)),as.numeric(S_estimator[i,]))
  }
  return (list(mu_t_tau,trapz(t_vec, mu_t_tau)))
}

dat_original <- read.csv("Azith.csv")
baseline_utility <- map_utility(dat_original$fev1pp_00)
# You can change 0.3 to 0.5 or 0.8 for different baseline utility strata
dat <- dat_original[baseline_utility == 0.3,]
events <- dat[,c(2,4:14)]
Di <- dat$dayofdeath
Di <- ifelse(is.na(Di), 400, Di)
Ci <- dat$Days_In_Study
Ci <- ifelse(Ci == Di, 380, Ci)
severity <- dat[,c(2,16:26)]
group <- dat$trtgroup
fev00 <- dat$fev1pp_00
fev06 <- dat$fev1pp_06
fev12 <- dat$fev1pp_12


# wherever fev06 is missing, substitute the fev00 value
idx06_na <- is.na(fev06)
fev06[idx06_na] <- fev00[idx06_na]

# wherever fev12 is missing, substitute the MOST RECENT previous
idx12_na <- is.na(fev12)
fev12[idx12_na] <- fev06[idx12_na]

fev1pp <- matrix(NA, nrow = nrow(dat), ncol = 361)
fev1pp[,1:180]   <- fev00
fev1pp[,181:360] <- fev06
fev1pp[,361] <- fev12

U <- t(apply(fev1pp, 1, map_utility))
colnames(U) <- seq(0, 360, 1)
events_new <- get_events_for_severity(events, severity, "severe_mod")
events_new <- cbind(events_new, Di, Ci)
n <- nrow(events_new)
end_time <- pmin(events_new$Di, events_new$Ci)
len_events <- ncol(events_new)
events_list <- get_events_list(events_new, end_time, n, len_events)
utilities <- c(0.15, 0.3, 0.5, 0.8, 1)
U_new <- get_utility_for_events(events_list, U, utilities)


events_new_treatment <- events_new[group == 1, ]
events_new_control <- events_new[group == 2, ]
U_treatment <- U_new[group == 1, ]
U_control <-  U_new[group == 2, ]

events_list_treatment <- events_list[group == 1]
events_list_control <- events_list[group == 2]

tau <- 90
t_vec <-c(0, 60, 120, 180, 240)
x_seq <- seq(0, tau, 1)
S_estimator1_info <- get_S_estimator(t_vec, x_seq, tau, events_new_treatment, U_treatment)
S_estimator1 <- S_estimator1_info$S_Q
info1 <- S_estimator1_info$info
S_estimator2_info <- get_S_estimator(t_vec, x_seq, tau, events_new_control, U_control)
S_estimator2 <- S_estimator2_info$S_Q
info2 <- S_estimator2_info$info

AUC1 <- get_estimate_final_new(S_estimator1, t_vec)[[2]]
AUC2 <- get_estimate_final_new(S_estimator2, t_vec)[[2]]
mu_t_tau1 <- get_estimate_final_new(S_estimator1, t_vec)[[1]]
mu_t_tau2 <- get_estimate_final_new(S_estimator2, t_vec)[[1]]

n1 <- nrow(events_new_treatment)
n2 <- nrow(events_new_control)
t_integrals1 <- get_var_estimator_parallel(t_vec, x_seq, tau, events_new_treatment, 
                                           U_treatment, S_estimator1, n2/(n1+n2), info1)
t_integrals1[lower.tri(t_integrals1)] <- t(t_integrals1)[lower.tri(t_integrals1)]
row1  <- apply(t_integrals1 , 1, trapz, x = t_vec)  
var1 <- as.numeric(trapz(t_vec, row1)) 
t_integrals2 <- get_var_estimator_parallel(t_vec, x_seq, tau, events_new_control, 
                                           U_control, S_estimator2, n1/(n1+n2), info2)
t_integrals2[lower.tri(t_integrals2)] <- t(t_integrals2)[lower.tri(t_integrals2)]
row2  <- apply(t_integrals2 , 1, trapz, x = t_vec)  
var2 <- as.numeric(trapz(t_vec, row2)) 

diff <- (AUC1 - AUC2)/(240*tau)
sd <- (sqrt((var1+var2)/(n1*n2/(n1+n2))))/(240*tau)
paste0("(", round(diff-qnorm(0.975)*sd,5), ",", round(diff+qnorm(0.975)*sd,5), ")")
z <- diff/sd
pval <- 2*(1 - pnorm(abs(z)))
pval




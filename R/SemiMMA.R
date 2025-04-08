#' Construct Subset Data for Outcome j
#'
#' Creates a subset of data for a specific outcome.
#' @param theta A matrix with random effect
#' @param y A matrix with effect sizes
#' @param s A matrix with standard deviations
#' @param j A scalar for the outcome j
#' @return A matrix with the subset data D_j contains effect size and standard deviations
#' @export

construct_subsets <- function(y, s, theta, j){
  J <-  ncol(y)
  if(J == 2){
    temp_y = matrix(y[,-j], ncol=1)
  } else {
    temp_y = y[,-j]
  }
  temp_row = apply(temp_y, 1, function(x){sum(is.na(x)) == 0})
  if(sum(temp_row) == 1){
    result = matrix(c(y[temp_row,],s[temp_row,], theta[temp_row,]), nrow = 1)
  } else {
    result = cbind(y[temp_row,],s[temp_row,], theta[temp_row,])
  }
  return(result)
}
#' Probability of Observing Outcome j in Study i
#'
#' Calculates observation probability.
#' @param s_i A vector of within-study standard deviation for study i
#' @param y_i A vector of effect size for study i
#' @param alpha_j A scalar of parameter alpha_j in the selection model for outcome j
#' @param beta_j A vector of parameter beta_j in the selection model for outcome j
#' @param gamma_j A vector of parameter gamma_j in the selection model for outcome j
#' @return A scalar of the probability for observing the outcome j in study i
#' @export
probability_of_observation <- function(s_i, y_i, alpha_j, beta_j, gamma_j){
  p_ij <- (1+exp(alpha_j + sum(beta_j * y_i) + sum(gamma_j * s_i)))^(-1)
  return(p_ij[1])
}

#' Compute Estimating Functions for Selection Model
#'
#' Evaluates the estimating functions for a semiparametric selection model, facilitating the estimation of parameters that mitigate outcome reporting bias in multivariate meta-analysis.
#' @param param A vector of parameters in the selection model
#' @param subset_data A matrix with the subset data D_j contains effect size and standard deviations
#' @return A matrix with estimating function values
#' @export
estimating_function <- function(param, subset_data, j){
  nos <- nrow(subset_data)
  temp_y <- subset_data[,1:J]
  temp_s <- subset_data[,J+1:J]
  temp_theta <- subset_data[,2*J+1:J]
  if(is.null(dim(temp_y))){
    temp_y = matrix(temp_y, nrow = 1)
    temp_s = matrix(temp_s, nrow = 1)
    temp_theta = matrix(temp_theta, nrow = 1)
  }

  result <- numeric(3*J-1)

  alpha_j = param[1]
  beta_j = param[1+1:J]
  gamma_j = param[1+J+1:J]

  result <- matrix(0, nrow = nos, ncol = 3 * J - 1)

  for(i in 1:nos){
    temp_indiv_y <- temp_y[i,]
    temp_indiv_s <- temp_s[i,]
    temp_v <- subset_data[i,-j]

    if(sum(is.na(temp_indiv_y)) > 0){
      temp_g <- (-1) * temp_v
    } else {
      temp_g <- (1/probability_of_observation(temp_indiv_s, temp_indiv_y, alpha_j, beta_j, gamma_j) -1) * temp_v
    }

    result[i, ] <- temp_g
  }
  return(result)
}

#' Perform One-Step Parameter Update
#'
#' Updates the parameter in a single iteration
#' @param outcome_parameter A vector with previous parameters
#' @param y A matrix with effect sizes
#' @param s A matrix with standard deviations
#' @return A vector with updated parameters
#' @export
one_step_update <- function(outcome_parameter, y, s){
  nos = nrow(y)
  J = ncol(y)
  current_mu = outcome_parameter[1:J]
  current_tau = outcome_parameter[J+1:J]
  current_rho = outcome_parameter[2*J+1:(J*(J-1)/2)]
  location_index = 1
  current_sigma = diag(current_tau^2)
  for(j in 1:(J-1)){
    for(k in (j+1):J){
      current_sigma[j,k] = current_rho[location_index] * current_tau[j] * current_tau[k]
      current_sigma[k,j] = current_sigma[j,k]
      location_index = location_index + 1
    }
  }
  if(!matrixcalc::is.positive.definite(current_sigma)){
    current_sigma = Matrix::nearPD(current_sigma)$mat
  }




  init_phi <- c(runif(1,-1,1), # alpha
                runif(J,-1,1), # beta
                runif(J,-1,1)) # gamma1



  beta <- matrix(nrow = J, ncol = J)
  gamma <- matrix(nrow = J, ncol = J)
  alpha <- numeric(J)


  for(j in 1:J){
    warn_msg <- 1

    while(warn_msg > 0){
      warn_msg <- 0
      #sim_theta <- condition_expect_theta(y, s, current_mu, current_sigma)
      sim_theta <- MASS::mvrnorm(n = nos, mu = c(current_mu), Sigma=current_sigma)
      temp_subset_data <-  construct_subsets(y, s, sim_theta, j)
      temp_res <- estimating_function(init_phi, temp_subset_data, j)
      temp_res_cov <-  cov(temp_res)
      diag(temp_res_cov) <- diag(temp_res_cov) + 1e-6
      gmm_estimating <- local({
        j_local <- j
        function(theta, x) {
          estimating_function(theta, x, j_local)
        }
      })
      temp_gmm_res <- tryCatch(gmm::gmm(gmm_estimating, temp_subset_data, init_phi,
                                        vcov = "iid", weights = temp_res_cov, onlyCoefficients = T),
                               warning = function(w){warn_msg <<- 1})

    }
    temp_phi <- coef(temp_gmm_res)
    alpha[j] <- temp_phi[1]
    beta[,j] <- temp_phi[1+1:J]
    gamma[,j] <- temp_phi[1+J+1:J]
  }

  E_y_nominator = 0
  E_y2_nominator = 0
  E_yjyk_nominator = 0
  est_mu_var_nominator = 0
  denominator = 0
  for(i in 1:nos){
    temp_y = y[i,]
    temp_s = s[i,]
    temp_unobs_ind = is.na(temp_y)
    temp_obs_ind = !is.na(temp_y)
    temp_nobs = sum(temp_obs_ind)
    temp_obs_prob = 1
    if(temp_nobs == J){
      for(j in 1:J){
        temp_obs_prob = temp_obs_prob * probability_of_observation(temp_s, temp_y, alpha[j], beta[,j], gamma[,j])
      }
      E_y_nominator = temp_y / temp_obs_prob + E_y_nominator
      E_y2_nominator = temp_y^2 / temp_obs_prob + E_y2_nominator
      denominator = denominator + 1/temp_obs_prob
      temp_yjyk = numeric(J*(J-1)/2)
      count_index = 1
      for(j in 1:(J-1)){
        for(k in (j+1):J){
          temp_yjyk[count_index] = temp_y[j] * temp_y[k]
          count_index = count_index + 1
        }
      }
      E_yjyk_nominator = temp_yjyk/temp_obs_prob + E_yjyk_nominator
    }
  }
  update_mu = E_y_nominator/denominator
  update_tau = numeric(J)

  mean_s2 = colMeans(s^2)
  temp_tau = E_y2_nominator/denominator - update_mu^2 - mean_s2
  temp_tau[temp_tau <= 0] = temp_tau[temp_tau <= 0] + mean_s2[temp_tau <= 0]
  if(sum(temp_tau <= 0) > 0 ){
    temp_tau[temp_tau <= 0] = (E_y2_nominator/denominator)[temp_tau <= 0]
  }
  update_tau = sqrt(temp_tau)

  est_mu_var_nominator = 0
  denominator = 0
  for(i in 1:nos){
    temp_y = y[i,]
    temp_s = s[i,]
    temp_unobs_ind = is.na(temp_y)
    temp_obs_ind = !is.na(temp_y)
    temp_nobs = sum(temp_obs_ind)
    temp_obs_prob = 1
    if(temp_nobs == J){
      for(j in 1:J){
        temp_obs_prob = temp_obs_prob * probability_of_observation(temp_s, temp_y, alpha[j], beta[,j], gamma[,j])
      }
      est_mu_var_nominator = (temp_s^2 +  update_tau^2) / temp_obs_prob^2 + est_mu_var_nominator
      denominator = denominator + 1/temp_obs_prob
    }
  }
  est_mu_se <- sqrt(est_mu_var_nominator)/denominator

  E_yjyk = E_yjyk_nominator/denominator
  count_index = 1
  update_rho = numeric(J*(J-1)/2)
  for(j in 1:(J-1)){
    for(k in (j+1):J){
      update_rho[count_index] = (E_yjyk[count_index] - update_mu[j]*update_mu[k])/(update_tau[j]*update_tau[k])
      count_index = count_index + 1
    }
  }
  update_output_parameter = c(update_mu, update_tau, update_rho, est_mu_se)

  return(update_output_parameter)
}

#' Impute missing within-study standard deviations
#'
#'
#' @param y A matrix with effect sizes
#' @param s A matrix with within-study standard deviations having missing values
#' @return A matrix with imputed within-study standard deviations no missing values
#' @export

impute_within_study_sd <- function(s, n){
  imputed_s = s
  J <- ncol(s)

  for(k in 1:J){
    temp_missing_row <- is.na(s[,k])
    temp_k <- mean(1/s[!temp_missing_row,k]^2/n[!temp_missing_row])
    imputed_s[temp_missing_row, k] = sqrt(1/(temp_k*n[temp_missing_row]))
  }
  return(imputed_s)

}

#' Fit Semiparametric Multivariate Meta-Analysis Model
#'
#' Fits a semiparametric model  using an iterative optimization approach.
#' @param y A matrix with effect sizes
#' @param s A matrix with within-study standard deviations
#' @return A vector contains estimated mu and the standard deviations
#' @export
SemiMMA <- function(y, s){
  if (sum(is.na(s)) > 0)
    stop("Within-study sd matrix have NA. Please impute missing value first.")

  J = ncol(y)
  mean_s = colMeans(s)
  temp_var = apply(y, 2, var, na.rm=T)
  init_tau = numeric(J)
  for(j in 1:J){
    init_tau[j] = ifelse(temp_var[j] > mean_s[j]^2,
                         temp_var[j] - mean_s[j]^2,
                         sqrt(temp_var[j]))
  }
  init_mu = colMeans(y, na.rm = T)

  init_rhoB = numeric(J*(J-1)/2)
  count_index = 1
  for(j in 1:(J-1)){
    for(k in (j+1):J){
      init_rhoB[count_index] = 0.4
      count_index = count_index + 1
    }
  }
  init_para = c(init_mu, init_tau, init_rhoB, rep(runif(1),J))

  # temp_res = SQUAREM::squarem(par = init_para, fixptfn= one_step_update,  y = y, s = s,
  #                 control = list(tol = 1e-4, maxiter = 100))
  temp_res = tryCatch(SQUAREM::squarem(par = init_para, fixptfn= one_step_update,  y = y, s = s,
                                       control = list(tol = 1e-4, maxiter = 100)),
                      error=function(e){NA})

  if(is.na(temp_res[1])){
    temp_mu = NA
    temp_se = NA
  } else {
    temp_mu = temp_res$par[1:J]
    temp_se = temp_res$par[2*J+J*(J-1)/2 + 1:J]
  }
  return(c(temp_mu,temp_se))
}

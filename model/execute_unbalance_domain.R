source("assist_functions.R")

#' Execute BFL domain model under unbalanced setting
#'
#' @param test_site The test site name
#' @param sites Vector of site names
#' @sim_data_filtered_list List of data for original domains
#' @param sim_data_target_domain_list List of data for target domain
#' @param posterior_phi_full List of posterior phi
#' @param LCVA_local_model_test_obs_fit List of fitted LCVA models
#' @return A list with accuracy, confusion matrix, and CSMF accuracy
execute_unbalance_domain <- function(test_site, sites, sim_data_filtered_list, sim_data_target_domain_list,
                                      posterior_phi_full, LCVA_local_model_test_obs_fit) {
  
  library(rstan)
  library(caret)
  
  cause_list = list()
  for (i in 1:length(sites)) {
    train_site = sites[i]
    cause_list[[train_site]] <- as.numeric(sim_data_filtered_list[[train_site]]$mapping_new_to_origin)
  }
  
  stan_model <- stan_model(file = "../code/BFL_no_partial_label.stan") 
  
 #  print(paste("Processing BFL-domain model for test site:", test_site))
  
  cause_lists_train <- lapply(sites, function(site) cause_list[[site]])
  C_max <- max(unlist(cause_lists_train))
  all_causes <- sort(unique(unlist(cause_lists_train)))
  num_causes <- length(all_causes)
  
  models <- seq_along(sites)
  site_to_model <- setNames(models, sites)
  
  model_presence <- matrix(0, nrow = num_causes, ncol = length(sites), dimnames = list(all_causes, sites))
  for (site in sites) {
    for (cause in cause_list[[site]]) {
      model_presence[as.character(cause), site] <- 1
    }
  }
  
  count <- apply(model_presence, 1, function(x) sum(x != 0))
  missing_indices <- which(is.na(sim_data_target_domain_list[[test_site]]$filtered_data$data$Y.t))
  
  test_X_i <- as.matrix(sim_data_target_domain_list[[test_site]]$filtered_data$data$X[missing_indices, ])
  whole_X_i <- as.matrix(sim_data_target_domain_list[[test_site]]$filtered_data$data$X)
  num_models <- length(sites)
  
  posterior_input_partial_i <- array(0, dim = c(nrow(test_X_i), C_max, num_models))
  
  for (j in seq_along(sites)) {
    model_site <- sites[j]
    if (model_site != test_site){
      posterior_phi_j = posterior_phi_full[[model_site]]
    } else {
      out_obs <- LCVA.pred(fit = LCVA_local_model_test_obs_fit[[model_site]], X_test = whole_X_i, 
                           model = "C", Burn_in = 500, Nitr = 4000, return_likelihood = TRUE, verbose = FALSE)
      posterior_phi_j = apply(out_obs$x_given_y_prob[-c(1:2000), ,], c(2, 3), mean)
    }
    
    for (i in 1:length(cause_list[[model_site]])) {
      cause_index <- cause_list[[model_site]][i]
      posterior_input_partial_i[, cause_index, j] <- posterior_phi_j[missing_indices, i]
    }
  }
  
  stan_data_domain_and_no_partial_i <- list(
    N = nrow(test_X_i),
    C_max = C_max,
    M = num_models,
    num_causes = num_causes,
    causes = all_causes,
    model_presence = model_presence,
    count = count,
    phi = posterior_input_partial_i
  )
  
  fit <- sampling(stan_model, data = stan_data_domain_and_no_partial_i, iter = 2000, chains = 4, refresh = 0)
  posterior_samples <- rstan::extract(fit)
  
  pred_res <- get_global_post_pred_Y(
    posterior_input_partial_i,
    posterior_samples$lambda,
    posterior_samples$pi
  )
  
  Y_pred <- pred_res$posterior_pred_Y
  
  true_test_labels <- as.numeric(sim_data_target_domain_list[[test_site]]$mapping_new_to_origin[
    sim_data_target_domain_list[[test_site]]$filtered_data$data.truth$Y.t[missing_indices]
  ])
  
  acc <- get_acc(
    Y_pred, 
    true_test_labels
  )
  
  Y_pred_factor = factor(get_Y_pred(Y_pred), levels = levels(as.factor(true_test_labels)))
  
  true_test_labels <- factor(true_test_labels, levels = levels(as.factor(true_test_labels)))
  conf_matrix <- confusionMatrix(Y_pred_factor, true_test_labels)
  
  true_test_labels <- as.numeric(sim_data_target_domain_list[[test_site]]$filtered_data$data.truth$Y.t[missing_indices])
  true_label_dist <- as.numeric(table(true_test_labels)) / length(true_test_labels)
  
  aligned_prev <- align_prevalence(
    colMeans(posterior_samples$pi),
    setNames(all_causes, all_causes),
    true_label_dist,
    sim_data_target_domain_list[[test_site]]$mapping_new_to_origin[as.numeric(names(table(true_test_labels)))]
  )
  
  csmf_acc <- CSMF_acc(aligned_prev$aligned_A, aligned_prev$aligned_B)
  
  cause_level = sim_data_filtered_list[[test_site]]$filtered_data$data$cause_level
  adjust_est_prev <- aligned_prev$aligned_A
  names(adjust_est_prev) = cause_level
  
  cause_pred <- cause_level[get_Y_pred(Y_pred)]
  
  return(list(csmf = adjust_est_prev, cause_pred = cause_pred, acc = acc, conf_matrix = conf_matrix, csmf_acc = csmf_acc))

}

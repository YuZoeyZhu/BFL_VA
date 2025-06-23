source("assist_functions.R")

execute_balance_mix <- function(test_site, sites, posterior_phi_full, sim_data_filtered_list) {
  library(rstan)
  library(caret)
  
  # message(paste("Processing BFL-mix model for test site:", test_site))
  
  stan_model <- stan_model(file = "../code/BFL_partial_labeled_balanced.stan") 
  
  cause_list = list()
  for (i in 1:length(sites)) {
    train_site = sites[i]
    cause_list[[train_site]] <- as.numeric(sim_data_filtered_list[[train_site]]$mapping_new_to_origin)
  }
  
  split_res <- split_domain_and_partial(sim_data_filtered_list[[test_site]])
  domain_indices <- split_res$domain_indices
  partial_indices <- split_res$partial_indices
  domain_data <- split_res$domain_data
  
  train_X <- as.matrix(domain_data$data$X)
  train_Y <- domain_data$data.truth$Y.t
  colnames(train_X) <- NULL
  rownames(train_X) <- NULL
  
  lcva_fit <- LCVA.train(
    X = train_X, Y = train_Y,
    Domain = rep(1, length(train_Y)),
    K = 5, model = "S", Nitr = 2000, thin = 2, nchain = 3, seed = 12345, verbose = FALSE
  )
  
  cause_lists_train <- lapply(sites, function(site) cause_list[[site]])
  C_max <- max(unlist(cause_lists_train))
  all_causes <- sort(unique(unlist(cause_lists_train)))
  num_causes <- length(all_causes)
  
  model_presence <- matrix(0, nrow = num_causes, ncol = length(sites), 
                           dimnames = list(all_causes, sites))
  for (site in sites) {
    for (cause in cause_list[[site]]) {
      model_presence[as.character(cause), site] <- 1
    }
  }
  count <- apply(model_presence, 1, function(x) sum(x != 0))
  
  missing_indices <- which(is.na(sim_data_filtered_list[[test_site]]$filtered_data$data$Y.t))
  test_idx <- c(partial_indices, missing_indices)
  test_X_i <- as.matrix(sim_data_filtered_list[[test_site]]$filtered_data$data$X[test_idx, ])
  whole_X_i <- as.matrix(sim_data_filtered_list[[test_site]]$filtered_data$data$X)
  colnames(test_X_i) <- NULL
  rownames(test_X_i) <- NULL
  
  num_models <- length(sites)
  posterior_input_partial_i <- array(0, dim = c(nrow(test_X_i), C_max, num_models))
  
  for (j in seq_along(sites)) {
    model_site <- sites[j]
    if (model_site != test_site) {
      posterior_phi_j <- posterior_phi_full[[model_site]]
    } else {
      out <- LCVA.pred(fit = lcva_fit, X_test = whole_X_i, model = "C", Burn_in = 500, Nitr = 4000, return_likelihood = TRUE, verbose = FALSE)
      posterior_phi_j <- apply(out$x_given_y_prob[-c(1:2000), ,], c(2, 3), mean)
    }
    for (i in seq_along(cause_list[[model_site]])) {
      cause_index <- cause_list[[model_site]][i]
      posterior_input_partial_i[, cause_index, j] <- posterior_phi_j[test_idx, i]
    }
  }
  
  Y <- as.numeric(sim_data_filtered_list[[test_site]]$mapping_new_to_origin[
    sim_data_filtered_list[[test_site]]$filtered_data$data$Y.t[test_idx]])
  Y_known <- ifelse(is.na(Y), 0, 1)
  Y[Y_known == 0] <- 1
  
  stan_data <- list(
    N = nrow(test_X_i), C_max = C_max, M = num_models,
    num_causes = num_causes, causes = all_causes,
    model_presence = model_presence, count = count,
    phi = posterior_input_partial_i, Y_known = Y_known, Y = Y
  )
  
  fit <- sampling(stan_model, data = stan_data, iter = 2000, chains = 4, refresh = 0)
  posterior_samples <- rstan::extract(fit)
  
  pred_res <- get_global_post_pred_Y(
    posterior_input_partial_i[which(Y_known == 0), ,],
    posterior_samples$lambda, posterior_samples$pi
  )
  
  Y_pred = pred_res$posterior_pred_Y
  
  true_test_labels <- as.numeric(sim_data_filtered_list[[test_site]]$mapping_new_to_origin[
    sim_data_filtered_list[[test_site]]$filtered_data$data.truth$Y.t[missing_indices]])
  acc <- get_acc(Y_pred, true_test_labels)
  
  Y_pred_factor <- factor(get_Y_pred(Y_pred), levels = levels(as.factor(true_test_labels)))
  true_test_labels_factor <- factor(true_test_labels, levels = levels(as.factor(true_test_labels)))
  conf_matrix <- confusionMatrix(Y_pred_factor, true_test_labels_factor)
  
  true_test_labels_all <- as.numeric(sim_data_filtered_list[[test_site]]$filtered_data$data.truth$Y.t)
  true_test_label_distribution <- as.numeric(table(true_test_labels_all)) / length(true_test_labels_all)
  n0 <- length(true_test_labels_all)
  nU <- length(test_idx)
  labeled_Yt <- sim_data_filtered_list[[test_site]]$filtered_data$data.truth$Y.t[-test_idx]
  labeled_Yt <- sim_data_filtered_list[[test_site]]$mapping_new_to_origin[labeled_Yt]
  nLc <- as.numeric(table(factor(labeled_Yt, levels = all_causes)))
  
  aligned_prev <- align_prevalence(
    (nU * colMeans(posterior_samples$pi) + nLc) / n0,
    setNames(all_causes, all_causes),
    true_test_label_distribution,
    sim_data_filtered_list[[test_site]]$mapping_new_to_origin[as.numeric(names(table(true_test_labels_all)))]
  )
  
  csmf_acc <- CSMF_acc(aligned_prev$aligned_A, aligned_prev$aligned_B)
  
  cause_level = sim_data_filtered_list[[test_site]]$filtered_data$data$cause_level
  adjust_est_prev <- aligned_prev$aligned_A
  names(adjust_est_prev) = cause_level
  
  cause_pred <- cause_level[get_Y_pred(Y_pred)]
  
  return(list(csmf = adjust_est_prev, cause_pred = cause_pred, acc = acc, conf_matrix = conf_matrix, csmf_acc = csmf_acc))
  
}

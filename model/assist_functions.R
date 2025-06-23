## function: filter the sparse causes with threshold 
filter_sparse_causes <- function(sim_data, threshold = 0) {
  # Calculate the frequency of each cause of death
  cause_freq <- table(sim_data$data.truth$Y.t)
  
  # Identify causes with samples less than or equal to the threshold
  sparse_causes <- names(cause_freq[cause_freq <= threshold])
  sparse_causes <- as.numeric(sparse_causes)
  
  # Find indices of samples to keep (not in sparse causes)
  keep_indices <- which(!sim_data$data.truth$Y.t %in% sparse_causes)
  
  # Filter the data
  sim_data$data.truth$Y.t <- sim_data$data.truth$Y.t[keep_indices]
  sim_data$data$Y.t <- sim_data$data$Y.t[keep_indices]
  sim_data$data$X <- sim_data$data$X[keep_indices, ]
  sim_data$data$G <- sim_data$data$G[keep_indices]
  
  # Get unique causes after filtering and create a continuous mapping
  unique_causes <- unique(sim_data$data.truth$Y.t)
  new_labels <- seq_along(unique_causes)
  # mapping <- setNames(new_labels, unique_causes)
  mapping_origin_to_new <- setNames(new_labels, unique_causes)
  mapping_new_to_origin <- setNames(unique_causes, new_labels)
  
  # Apply the new labeling to the filtered data
  sim_data$data.truth$Y.t <- as.numeric(mapping_origin_to_new[as.character(sim_data$data.truth$Y.t)])
  sim_data$data$Y.t <- as.numeric(mapping_origin_to_new[as.character(sim_data$data$Y.t)])
  
  # Update the number of causes after filtering
  sim_data$C <- length(unique_causes)
  
  # Return a list containing the filtered data and the mapping
  list(filtered_data = sim_data, mapping_origin_to_new = mapping_origin_to_new, mapping_new_to_origin = mapping_new_to_origin)
}


# Function to check if the missing categories are a subset of observed categories
check_valid_sample <- function(full_data, missing_indices) {
  observed_data <- full_data[-missing_indices]  # Data that is still observed
  missing_data <- full_data[missing_indices]  # Data that is missing
  
  observed_categories <- unique(observed_data[!is.na(observed_data)])  # Unique categories in observed data
  missing_categories <- unique(missing_data)  # Unique categories in missing data
  
  # Check if all missing categories exist in observed categories
  return(all(missing_categories %in% observed_categories))
}

generate_missing_Yt <- function(Y.t, miss_prop) {
  
  num_samples <- round(miss_prop * length(Y.t))
  
  # Identify small sample size categories (sample size ≤ 10)
  category_counts <- table(Y.t)
  small_categories <- as.numeric(names(category_counts[category_counts <= 10]))
  
  # Initialize missing indices
  missing_indices <- integer(0)
  
  # Handle small sample size categories separately
  for (cod in small_categories) {
    cod_indices <- which(Y.t == cod)
    n_cod <- length(cod_indices)  # Number of samples in this category
    
    if (n_cod > 1) {
      num_miss <- min(round(miss_prop * n_cod), n_cod - 1)  # Ensure at least one remains
      small_miss_indices <- sample(cod_indices, num_miss, replace = FALSE)
      missing_indices <- c(missing_indices, small_miss_indices)
    }
  }
  
  
  # For other categories, ensure a valid sample
  repeat {
    general_indices <- which(!Y.t %in% small_categories)  # Exclude small sample indices
    general_missing <- sample(general_indices, num_samples - length(missing_indices), replace = FALSE)
    
    # Combine small category missing indices and general missing indices
    final_missing_indices <- c(missing_indices, general_missing)
    
    if (check_valid_sample(Y.t, final_missing_indices)) {
      break
    }
  }
  
  return(final_missing_indices)  # Return modified Y.t with missing values
}


generate_train_test_Yt_unbalanced_no_overlap_sampling <- function(Y.t, miss_prop = NULL) {
  unique_categories <- sort(unique(Y.t))
  num_categories <- length(unique_categories)
  
  total_n <- length(Y.t)
  target_unlabel_size <- round(miss_prop * total_n)
  
  is.label <- integer(0)
  no.label <- integer(0)
  
  unlabel_counts <- numeric(num_categories)
  category_sizes <- numeric(num_categories)
  
  # Step 1: Sample Beta proportion for each class and get raw unlabel counts
  for (i in seq_along(unique_categories)) {
    cod <- unique_categories[i]
    cod_indices <- which(Y.t == cod)
    category_sizes[i] <- length(cod_indices)
    
    p <- rbeta(1, 0.2, 0.2)  # Beta sampling per class
    unlabel_counts[i] <- round(length(cod_indices) * p)
  }
  
  # Step 2: Rescale to match target total unlabel size
  scaling_factor <- target_unlabel_size / sum(unlabel_counts)
  unlabel_counts_scaled <- round(unlabel_counts * scaling_factor)
  
  # Step 3: Sample unlabeled and labeled indices
  for (i in seq_along(unique_categories)) {
    cod <- unique_categories[i]
    cod_indices <- which(Y.t == cod)
    n_unlabel <- min(length(cod_indices), unlabel_counts_scaled[i])
    
    unlabel_indices <- sample(cod_indices, n_unlabel)
    label_indices <- setdiff(cod_indices, unlabel_indices)
    
    is.label <- c(is.label, label_indices)
    no.label <- c(no.label, unlabel_indices)
  }
  
  # Optional sanity check
  # cat("Proportion unlabeled:", length(no.label) / total_n, "\n")
  
  
  # Step 5: Compute prevalence
  get_full_pi <- function(indices) {
    tab <- table(factor(Y.t[indices], levels = unique_categories))
    prop <- as.numeric(tab / length(indices))
    names(prop) <- unique_categories
    return(prop)
  }
  
  actual_pi1 <- get_full_pi(is.label)
  actual_pi2 <- get_full_pi(no.label)
  
  return(list(
    training_indices = is.label,
    testing_indices = no.label,
    actual_pi1 = actual_pi1,
    actual_pi2 = actual_pi2
  ))
}




generate_missing_Yt_unbalanced <- function(Y.t, miss_prop) {
  
  # Number of unique categories
  unique_categories <- unique(Y.t)
  num_categories <- length(unique_categories)
  
  # Generate training and testing prevalence distribution
  pi1 <- as.numeric(gtools::rdirichlet(1, rep(1, num_categories)))  
  pi2 <- as.numeric(gtools::rdirichlet(1, rep(1, num_categories)))  
  
  # Obtain training samples
  total_train_samples <- round(miss_prop * length(Y.t))
  train_category_counts <- rmultinom(1, size = total_train_samples, prob = pi1)
  
  training_indices <- unlist(mapply(function(cat, n) {
    inds <- which(Y.t == cat)
    if (length(inds) == 0) return(integer(0))
    sample(inds, size = n, replace = TRUE)
  }, unique_categories, train_category_counts))
  
  # Obtain testing samples
  total_test_samples <- round((1-miss_prop) * length(Y.t))
  test_category_counts <- rmultinom(1, size = total_test_samples, prob = pi2)
  
  testing_indices <- unlist(mapply(function(cat, n) {
    inds <- which(Y.t == cat)
    if (length(inds) == 0) return(integer(0))
    sample(inds, size = n, replace = TRUE)
  }, unique_categories, test_category_counts))
  
  
  # Check which categories are missing from training set
  categories_in_train <- unique(Y.t[training_indices])
  missing_categories <- setdiff(unique_categories, categories_in_train)
  
  # For each missing category, add sample to training (extreme small causes special handle)
  if (length(missing_categories) > 0) {
    additional_indices <- unlist(lapply(missing_categories, function(cat) {
      which(Y.t == cat)
    }))
    training_indices <- c(training_indices, additional_indices)
  }
  
  pi1_actual <- table(Y.t[training_indices]) / length(training_indices)
  pi2_actual <- table(Y.t[testing_indices]) / length(testing_indices)
  
  return(list(
    training_indices = training_indices,  
    testing_indices = testing_indices,  
    actual_pi1 = pi1_actual,  
    actual_pi2 = pi2_actual   
  ))
}



get_observed_sim_data <- function(sim_data){
  obs_data = list()
  obs_data$data.truth = list()
  obs_data$data = list()
  obs_data$data.truth$Y.t = sim_data$data$Y.t[which(!is.na(sim_data$data$Y.t))]
  obs_data$data$Y.t = sim_data$data$Y.t[which(!is.na(sim_data$data$Y.t))]
  obs_data$data$X = sim_data$data$X[which(!is.na(sim_data$data$Y.t)), ]
  obs_data$data$G = rep(1, length(obs_data$data.truth$Y.t))
  obs_data$NG = 1
  
  # Get unique causes after filtering and create a continuous mapping
  unique_causes <- unique(sim_data$data.truth$Y.t)
  
  # Update the number of causes after filtering
  obs_data$C <- length(unique_causes)
  return(list(obs_data = obs_data))
}



get_global_post_pred_Y <- function(posterior_phi, posterior_lambda, posterior_pi){
  # Number of posterior samples
  n_samples <- dim(posterior_lambda)[1]
  I <- dim(posterior_phi)[1]
  C <- dim(posterior_phi)[2]
  
  posterior_pred_Y_prob <- sapply(1:n_samples, function(s) {
    apply(posterior_phi, c(1, 2), function(phi_icm) sum(phi_icm * posterior_lambda[s, , ])) * posterior_pi[s,]
  }, simplify = "array")  
  
  # Normalize across C for each (s, i)
  posterior_pred_Y_prob <- apply(posterior_pred_Y_prob, c(1, 3), function(x) x / sum(x))
    posterior_pred_Y = apply(posterior_pred_Y_prob, c(2, 3), function(x) sample(1:C, 1, prob = x))
  return(list(posterior_pred_Y_prob = posterior_pred_Y_prob, posterior_pred_Y = t(posterior_pred_Y)))
}


get_sim_data_from_phmrc_by_location <- function(site){
  phmrc = read.csv("../data/phmrc_clean.csv")
  causes = unique(phmrc$cause)
  phmrc_site = phmrc[which(phmrc$site == site), ]
  
  
  phmrc_site$cause <- factor(phmrc_site$cause, levels = causes)
  phmrc_site$Y =  as.integer(phmrc_site$cause)
  
  data = list()
  data$cause_level = causes
  data$cause = phmrc_site$cause
  data$Y.t = phmrc_site$Y
  data$X = phmrc_site[, 1:168]
  data$G = rep(1, length(phmrc_site$Y))
  
  data.truth = list()
  data.truth$Y.t = phmrc_site$Y
  
  sim_data = list()
  sim_data$data = data
  sim_data$data.truth = data.truth
  sim_data$C = 34
  sim_data$NG = 1
  
  return(sim_data)
}





get_acc <- function(pred_Yt, true_Yt){
  acc = rep(NA, ncol(pred_Yt))
  
  for(t in 1:ncol(pred_Yt)){
    # Tabulate the frequencies of each cause of death
    cause_frequencies <- table(pred_Yt[, t])
    
    # Find the cause with the highest frequency
    most_probable_cause <- which.max(cause_frequencies)
    
    most_probable_cause_number <- as.integer(names(most_probable_cause))
    
    acc[t] = ifelse(most_probable_cause_number == true_Yt[t], 1, 0)
  }
  
  return(mean(acc))
}


get_Y_pred <- function(pred_Yt){
  Y = rep(NA, ncol(pred_Yt))
  
  for(t in 1:ncol(pred_Yt)){
    # Tabulate the frequencies of each cause of death
    cause_frequencies <- table(pred_Yt[, t])
    
    # Find the cause with the highest frequency
    most_probable_cause <- which.max(cause_frequencies)
    
    most_probable_cause_number <- as.integer(names(most_probable_cause))
    
    Y[t] = most_probable_cause_number
  }
  
  return(Y)
}




align_prevalence <- function(prevalence_A, mapping_A, prevalence_B, mapping_B) {
  # Extract unique cause lists
  causes_A <- unname(mapping_A)  # Unique causes in A
  causes_B <- unname(mapping_B)  # Unique causes in B
  
  # Get the union of both cause lists
  all_causes <- sort(unique(c(causes_A, causes_B)))
  
  # Create empty prevalence vectors filled with 0s
  aligned_A <- setNames(rep(0, length(all_causes)), all_causes)
  aligned_B <- setNames(rep(0, length(all_causes)), all_causes)
  
  # Fill in existing prevalence values for A
  aligned_A[as.character(causes_A)] <- prevalence_A
  
  # Fill in existing prevalence values for B
  aligned_B[as.character(causes_B)] <- prevalence_B
  
  return(list(aligned_A = aligned_A, aligned_B = aligned_B))
}

CSMF_acc <- function(pi_hat, pi) {
  # Ensure pi_hat and pi are of the same length and that they are proper probability distributions
  if(length(pi_hat) != length(pi) || sum(pi_hat) > 1.0001 || sum(pi) > 1.0001) {
    stop("Invalid input: pi_hat and pi must be of the same length and sum to 1.")
  }
  
  # Calculate the CSMF accuracy
  C <- length(pi)
  csmf_accuracy <- 1 - (sum(abs(pi_hat - pi)) / (2 * (1 - min(pi))))
  
  return(csmf_accuracy)
}



split_domain_and_partial <- function(sim_data){
  non_missing_indices <- which(!is.na(sim_data$filtered_data$data$Y.t))
  
  # Extract cause of death labels for non-missing data
  Y_t <-sim_data$filtered_data$data$Y.t[non_missing_indices]
  
  # Count occurrences of each cause
  cause_counts <- table(Y_t)
  
  # Initialize domain and partial subsets
  domain_indices <- c()
  partial_indices <- c()
  
  # Iterate through each unique cause
  for (cause in names(cause_counts)) {
    cause_indices <- non_missing_indices[Y_t == as.numeric(cause)]  # Get indices for this cause
    
    if (length(cause_indices) == 1) {
      # If only 1 occurrence, place in domain
      domain_indices <- c(domain_indices, cause_indices)
    } else {
      # Otherwise, randomly split half and half
      split_size <- floor(length(cause_indices) / 2)
      selected_domain <- sample(cause_indices, split_size)
      
      domain_indices <- c(domain_indices, selected_domain)
      partial_indices <- c(partial_indices, setdiff(cause_indices, selected_domain))
    }
  }
  
  domain_data = list()
  domain_data$data.truth = list()
  domain_data$data = list()
  domain_data$data.truth$Y.t = sim_data$filtered_data$data.truth$Y.t[domain_indices]
  domain_data$data$Y.t = sim_data$filtered_data$data.truth$Y.t[domain_indices]
  domain_data$data$X = sim_data$filtered_data$data$X[domain_indices, ]
  domain_data$data$G = rep(1, length(domain_data$data$Y.t))
  domain_data$NG = 1
  domain_data$C = length(unique(domain_data$data$Y.t))
  
  return(list(domain_indices = domain_indices, partial_indices = partial_indices, domain_data = domain_data))
}

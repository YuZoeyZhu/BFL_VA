#' Run BFL model evaluation
#'
#' @param model_type One of "domain", "partial", "base", or "mix"
#' @param test_site Name of the test site
#' @param sites List of training site names
#' @param sim_data_filtered_list List of preprocessed simulation data
#' @param posterior_phi_full List of posterior phi from base models
#' @param LCVA_local_model_test_obs_fit (optional) Required for domain 
#' @return A list with ACC, confusion matrix, and CSMF
run_BFL_balance_model <- function(model_type, test_site, sites,
                          sim_data_filtered_list,
                          posterior_phi_full,
                          LCVA_local_model_test_obs_fit = NULL) {
  if (!(model_type %in% c("domain", "partial", "base", "mix"))) {
    stop("Invalid model_type. Choose from: domain, partial, base, mix.")
  }
  
  if (model_type == "domain") {
    source("execute_balance_domain.R")
    # Run the domain model
    res <- execute_balance_domain(
      test_site = test_site,
      sites = sites,
      posterior_phi_full = posterior_phi_full,
      LCVA_local_model_test_obs_fit = LCVA_local_model_test_obs_fit,
      sim_data_filtered_list = sim_data_filtered_list
    )
    return(res)
  }
  
  if (model_type == "partial") {
    source("execute_balance_partial.R")
    # Run the partial model
    res <- execute_balance_partial(
      test_site = test_site,
      sites = sites,
      posterior_phi_full = posterior_phi_full,
      sim_data_filtered_list = sim_data_filtered_list
    )
    return(res)
  }
  
  if (model_type == "base") {
    source("execute_balance_base.R")
    # Run the domain model
    res <- execute_balance_base(
      test_site = test_site,
      sites = sites,
      posterior_phi_full = posterior_phi_full,
      sim_data_filtered_list = sim_data_filtered_list
    )
    return(res)
  }
  
  if (model_type == "mix") {
    source("execute_balance_mix.R")
    # Run the domain model
    res <- execute_balance_mix(
      test_site = test_site,
      sites = sites,
      posterior_phi_full = posterior_phi_full,
      sim_data_filtered_list = sim_data_filtered_list
    )
    return(res)
  }
}



#' Run BFL model evaluation under unbalanced data setting
#'
#' @param model_type One of "domain", "partial", "base", or "mix"
#' @param test_site Name of the test site
#' @param sites List of training site names
#' @param sim_data_filtered_list List of preprocessed simulation data
#' @param sim_data_target_domain_list List of simulation data with unbalanced missingness for the test site
#' @param posterior_phi_full List of posterior phi from base models
#' @param LCVA_local_model_test_obs_fit (optional) Required for domain model
#' 
#' @return A list with ACC, confusion matrix, and CSMF
run_BFL_unbalance_model <- function(model_type, test_site, sites,
                                    sim_data_filtered_list,
                                    sim_data_target_domain_list,
                                    posterior_phi_full,
                                    LCVA_local_model_test_obs_fit = NULL) {
  if (!(model_type %in% c("domain", "partial", "base", "mix"))) {
    stop("Invalid model_type. Choose from: domain, partial, base, mix.")
  }
  
  if (model_type == "domain") {
    source("execute_unbalance_domain.R")
    res <- execute_unbalance_domain(
      test_site = test_site,
      sites = sites,
      sim_data_filtered_list = sim_data_filtered_list,
      sim_data_target_domain_list = sim_data_target_domain_list,
      posterior_phi_full = posterior_phi_full,
      LCVA_local_model_test_obs_fit = LCVA_local_model_test_obs_fit
    )
    return(res)
  }
  
  if (model_type == "partial") {
    source("execute_unbalance_partial.R")
    res <- execute_unbalance_partial(
      test_site = test_site,
      sites = sites,
      sim_data_filtered_list = sim_data_filtered_list,
      sim_data_target_domain_list = sim_data_target_domain_list,
      posterior_phi_full = posterior_phi_full
    )
    return(res)
  }
  
  if (model_type == "mix") {
    source("execute_unbalance_mix.R")
    res <- execute_unbalance_mix(
      test_site = test_site,
      sites = sites,
      sim_data_filtered_list = sim_data_filtered_list,
      sim_data_target_domain_list = sim_data_target_domain_list,
      posterior_phi_full = posterior_phi_full
    )
    return(res)
  }
}




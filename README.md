# README

This repository contains the reproducible materials for the paper _Bayesian Federated Cause-of-Death Classification and Quantification Under Distribution Shift_.


## Overview of the structure and code scripts 

+ preprocess-phmrc-data.R: Pre-process code for the PHMRC dataset

+ code/: Stan functions for fitting the Bayesian Federated Learning (BFL) models with or without partial label. And for the cases with partial labels, assuming with or without label shift. 
  * `BFL_no_partial_label.stan`: BFL model without partial label. Fitted for BFL_base and BFL_domain models.
  * `BFL_partial_labeled_balanced.stan`:  BFL model with partial label and assuming without label shift. Fitted for BFL_partial and BFL_mixed models.
  * `BFL_partial_labeled_unbalanced.stan`:  BFL model with partial label and assuming with label shift. Fitted for BFL_partial and BFL_mixed models.

+ model/: R functions for fitting the local LCVA models and BFL model variants with or without label shift.


+ Reproducible_examples/: R Markdown for a reproducible example with no/mild/severe label shift.

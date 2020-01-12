# Function to simulate (hypothetical) binomially distributed data for the SP13
# replication. This will then be the input to our analysis pipeline script.
# It is called/sourced from a separate script.

library("mvtnorm")
library("tidyverse")
library("boot")  # for inv.logit()

# The function generates binomially distributed data from the generative
# model we will assume for our actual data from the replication. Some of the
# parameters that go into this function can be taken from our re-analysis of
# the original data. However, the original data do not consist of trial-level
# observations. Therefore we need to make up some of the parameters, which
# is done in the script that calls this function.

# For power analyses:
# The function can also be used for simulation-based power analyses. Step 3b
# is important in that regard. Here is what it does by default:
# Rather than sampling *all* fixed effects coefficients from the covariance
# matrix, step 3b keeps the effect of the critical interaction constant
# by plugging the intended interaction effect back in (see step 3b). This makes
# sense since we are evaluating the power conditioned on a specific effect
# size, so we don't want it to vary randomly in each simuation. In fact, for
# the Type I error analysis (false positive rate), this step is *required*.
# (The relevant discussion is found in email correspondence with Florian
# (see e-mail sent by Jaeger, Florian <fjaeger@UR.Rochester.edu>; Subj: "Re:
# Pragmatic advice on power analysis to determine sample size for conceptual
# replication", sent on Sat 2018-10-13 04:46).
# To sample the critical effect from the estimated distribution, set the
# parameter 'keep_critical_effect_constant' to FALSE

simulate_binom <- function (
  Nsubj  = 2,  # Number of participants
  Nitem  = 4,  # Number of items, i.e. the words participants have to memorize
  fixef_means,        # vector of mean estimates for fixed effects
  fixef_sigma,        # covariance matrix for fixed effects (diagonals contain SE^2)
  ranef_sigma_subj,   # covariance matrix for random effects by subject
  ranef_sigma_item,   # covariance matrix for random effects by item
  keep_critical_effect_constant = TRUE,  # see comment in paragraph above
  full_output = TRUE,      # Output list with fixef/ranef dfs in addition to data?
  print_each_step = FALSE  # print output at each step to unveil inner workings
  ) {

  myprint <- function(x) {
    if (print_each_step) {
      cat("\n"); print(paste("This is:", deparse(substitute(x)))); cat("\n")
      print(x)
      }
    }

  # 1) create design matrix X for experiment
  # The factors crossed
  myfactors <- data.frame(
    movement = factor(rep(c("arms", "legs"), each = 2)),
    word_type = factor(c("arm-word", "leg-word"))
    )
  # Use contrast coding
  contrasts(myfactors$movement) <- contr.sum(2)
  contrasts(myfactors$word_type) <- contr.sum(2)
  # Generate model matrix for the 2x2 factorial design
  X <- model.matrix(~ 1 + movement * word_type, myfactors)
  myprint(X)
  # Remove attributes to avoid warnings with left_join() later on:
  attr(myfactors$movement, "contrasts") <- NULL
  attr(myfactors$word_type, "contrasts") <- NULL
  myprint(myfactors)

  # 2) Create data frame that will contain the simulated data
  # 2a) Basic structure
  n_obs <- Nsubj * Nitem * 2  # total number of observations (or #rows)
  data <- tidyr::expand(myfactors, movement, word_type, subject = 1 : Nsubj) %>%
    slice(rep(1 : n(), each = Nitem / 2)) %>%  # Multiply rows by half the item number
    select(subject, movement, word_type) %>%   # Reorder columns
    arrange(subject, movement, word_type) %>%  # Rearrange rows
    mutate(
      item  = rep(1 : Nitem, times = 2 * Nsubj),  # Ss see each item twice (once per block)
      block = rep(c(1, 2, 2, 1), each = Nitem, length.out = n_obs)  # counterbalanced movement blocks
           )
  myprint(data)
  # 2b) Randomization of items
  data <- data %>%
    # group items randomly into trials (4 word per trial)
    group_by(subject, word_type, block) %>%
    mutate(sequence = base::sample.int(Nitem / 2)) %>%  # shuffle items
    arrange(subject, movement, word_type, sequence) %>%
    ungroup() %>%
    mutate(
      trialID = rep(1 : (Nitem / 4), each = 4, length.out = n_obs),  # identify trials
      pos_in_trial = rep(1 : 4, length.out = n_obs)  # position of word in trial
      ) %>%
    # shuffle trials randomly within blocks
    nest(-subject, - block, -trialID) %>%
    group_by(subject, block) %>%
    mutate(trial_in_block = base::sample.int(Nitem / 4)) %>%
    # But because we want "trial_in_exp" (not in block), we have to add #trials in block 1 to trials in block 2
    mutate(
      block_correction = ifelse(block == 1, 0, Nitem / 4),
      trial_in_exp = trial_in_block + block_correction
      ) %>%
    unnest()
  myprint(data)
  # Remove interim columns not necessary in the output
  data <- data %>%
    select(-trialID, -trial_in_block, -block_correction, -sequence) %>%
    arrange(subject, movement, trial_in_exp, pos_in_trial)
  myprint(data)

  # 3) Fixed effects in this simulation
  # This is the input means vector for fixed effects (from original model):
  myprint(fixef_means)
  # 3a) Sample fixef for current simulation from the model's covariance matrix
  fixef <- rmvnorm(
    n = 1,  # One simulation
    mean = fixef_means,  # The mean of fixed effects
    sigma = fixef_sigma  # covariance of fixed effects
    )
  myprint(fixef)
  # 3b) Plug the intended coefficient for the interaction back in (by default).
  if (keep_critical_effect_constant) { fixef[4] <- fixef_means[4] }
  myprint(fixef)
  # 3c) Save fixed effects as df for output
  fixef_df <- tibble(coef = colnames(fixef), betas = fixef[1, ])
  myprint(fixef_df)
  # 3c) Join fixed effects with data
  fixef_cell_means <- X %*% t(fixef)  # fixed effects in each experimental cell
  myprint(fixef_cell_means)
  fixef_cell_means_df <- myfactors %>% add_column(LO_fixef = as.vector(fixef_cell_means))  # works bc reads values columnwise
  myprint(fixef_cell_means_df)
  data <- left_join(data, fixef_cell_means_df)
  myprint(data)

  # 4) By-subject adjustments (multivariate normal with mean zero)
  # 4a) Each matrix row gives subject adjustments for betas of the 4 predictors in X
  subj_adjust <- rmvnorm(
    n = Nsubj,
    mean = rep(0, ncol(ranef_sigma_subj)),
    sigma = ranef_sigma_subj
    )
  colnames(subj_adjust) <- colnames(fixef)
  myprint(subj_adjust)
  # 4b) Save as df for output
  subj_adj_df <- as_tibble(subj_adjust) %>%
    add_column(subjID = 1 : Nsubj) %>%  # add subject identifiers)
    select(last_col(), everything())
  myprint(subj_adj_df)
  # 4c) Obtain by-subject cell means by multiplying model matrix X and
  # by-subject adjustments. Each column corresponds to matrix-multiplying X
  # by a *transposed* row-vector of subj_coef, and thus gives us the expected
  # cell adjustments for each subject (4 values per subject).
  subj_cell_adj <- X %*% t(subj_adjust)
  myprint(subj_cell_adj)
  # 4d) Arrange into data frame with subject-specific log-odds per cell
  subj_cell_adj_df <- myfactors %>% slice(rep(1 : n(), times = Nsubj)) %>%
    add_column(
      subject = rep(1 : Nsubj, each = nrow(subj_cell_adj)),
      LO_subj = as.vector(subj_cell_adj)  # works bc reads values columnwise
      )
  myprint(subj_cell_adj_df)
  # 4e) Join with data:
  data <- left_join(data, subj_cell_adj_df)
  myprint(data)

  # 5) By-item adjustments (multivariate normal with mean zero)
  # 5a) Each matrix row gives item adjustments for betas of 2 predictors:
  # intercept and movement (a by-item random slope)
  item_adjust <- rmvnorm(
    n = Nitem,
    mean = rep(0, ncol(ranef_sigma_item)),
    sigma = ranef_sigma_item
    )
  colnames(item_adjust) <- colnames(fixef)[1 : 2]
  myprint(item_adjust)
  # 5b) Save as df for output
  item_adj_df <- as_tibble(item_adjust) %>%
    add_column(itemID = 1 : Nitem) %>%  # Item identifier
    select(last_col(), everything())
  myprint(item_adj_df)
  # 5c) Obtain by-item cell means as above (for subjects). Only now there are
  # only two coefficients (intercept and by-item slope for movement)
  myprint(X[c(1,3), 1:2])
  item_cell_adj <- X[c(1,3), 1:2] %*% t(item_adjust)
  myprint(item_cell_adj)
  # 5d) Arrange into data frame that shows item-specific log-odds per cell
  item_cell_adj_df <- tibble(
    movement = rep(unique(myfactors$movement), times = Nitem),
    item     = rep(1 : Nitem, each = 2),
    LO_item  = as.vector(item_cell_adj)  # works bc reads values columnwise
    )
  myprint(item_cell_adj_df)
  # 4e) Join with data:
  data <- left_join(data, item_cell_adj_df)
  myprint(data)

  # 5a) Sum up all the values in logodds space and convert to probability
  data$LO <- with(data, LO_fixef + LO_subj + LO_item)
  data$prob <- inv.logit(data$LO)
  myprint(data)
  # 5b) and generate an actual observed binary response
  data$Error <- rbinom(n = nrow(data), size = 1, prob = data$prob)

  # 6) Add nuisance variable "preced_error" (binary) which is 1 if an error
  # was made on any of the preceding words in a trial, 0 otherwise
  prec_error <- function(x) {
    cumsum <- cumsum(x)
    shift_cumsum <- c(0, cumsum[1 : (length(cumsum) - 1)])
    out <- ifelse(shift_cumsum > 0, 1, 0)
    out
  }
  data <- data %>%
    group_by(subject, trial_in_exp) %>%
    mutate(preced_error = prec_error(Error))

  # reorder columns
  data <- data %>%
    select(subject : pos_in_trial, preced_error, Error, LO_fixef : prob)

  # out
  if (! full_output) {
    data
  } else {
    list(
      data = data,
      fixef_sample = fixef_df,
      subj_adjust = subj_adj_df,
      item_adjust = item_adj_df,
      fixef_sigma = fixef_sigma,
      ranef_sigma_subj = ranef_sigma_subj,
      ranef_sigma_item = ranef_sigma_item
      )
  }
}

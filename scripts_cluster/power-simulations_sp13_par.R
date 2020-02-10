## ---- message=FALSE------------------------------------------------------
# library("knitr")
library("tidyverse")  # ggplot2, dplyr, readr, purrr, etc
library("broom")
library("lme4")
library("brms")
library("tictoc")
library("foreach")
library("doParallel")


start_time <- strftime(Sys.time())

## ------------------------------------------------------------------------
# Load original data from SP13
d <- read_csv("SP13_orig-data_total-errors_long-format.csv") %>%
  # keep critical interference conditions only
  filter(movement %in% c("arm_paradi", "leg_paradi"))
# For each experimental cell, the maximum number of errors is 48 (see
# "Appendix_B_reanalysis_original.Rmd" for justification0)
d$n <- 48


## ------------------------------------------------------------------------
# Load from disk - list containing 1) full model, 2) null model, 3) BF
bfm_max <- read_rds("sp13_bfm_max.rds")


## ------------------------------------------------------------------------
bfm_full <- bfm_max[[1]]
# summary(bfm_full)


## ------------------------------------------------------------------------
source("generate_data_fnc.R")


## ------------------------------------------------------------------------
# The info we need is in the summary of fixed effects of the model:
fixef(bfm_full)

# Vector of coefficient means
fixef_means <- fixef(bfm_full)[, 1]
fixef_means

# Covariance matrix (Sigma); we need to square the SEs to convert them to Variances
fixef_sigma <- diag(fixef(bfm_full)[, 2] ^ 2)  # we assume uncorrelated diagonal matrix
fixef_sigma

# Random effects for the model:
VarCorr(bfm_full)$subject$sd
# Specifically we extract by-subject random SDs from this matrix (they are stored
# as a "stddev" attribute of the matrix) and square them to obtain variances:
ranef_sigma_subj <- diag(VarCorr(bfm_full)$subject$sd[, 1]) ^ 2
ranef_sigma_subj


## ------------------------------------------------------------------------
# We don't have estimates for by-item random slope for movement; we will assume
# it to be the same as the by-subject random variability for the critical
# interaction, which probably is an overestimate but should make our power
# estimates more conservative
ranef_sigma_item <- diag(c(0.65, 0.06)) ^ 2  # square to obtain variances
ranef_sigma_item



## ------------------------------------------------------------------------
analyze_simulation <- function(df, print_model_summaries = FALSE) {

  # start counter
  tic()

  # Contrast coding (umeric coding allows for models without RE correlations)
  df$mv <- scale(as.numeric(df$movement), scale = FALSE) * -2
  df$wt <- scale(as.numeric(df$word_type), scale = FALSE) * -2
  df$mv_wt <- df$mv * df$wt

  # function to catch convergence issues
  converge <- function (fm) {
  # Paste in case there are several messages
  message <- paste(summary(fm)$optinfo$conv$lme4$messages, collapse = " / ")
  if( is.null(message) ) { "" } else { message }
  }

  # Function to fit a full and null model given model formula as string:
  fit_formula <- function (f_str, d = df) {
    formula <- as.formula(f_str)  # convert to formula
    full <- glmer(formula, data = d,
                  control = glmerControl(optCtrl=list(maxfun=1e6)),
                  family = "binomial"
                  )
    null <- update(full, formula = ~ . - mv_wt)

    if (print_model_summaries) {
      print(summary(full))
      print(summary(null))
    }
    list(full, null)
  }

  # Fit models without RE correlations using double-bar syntax:
  myf <- "Error ~ mv + wt + mv_wt + (1 + mv + wt + mv_wt || subject) + (1 + mv || item)"
  models <- fit_formula(myf)

  # Put fixed effect model coefficients into one data frame
  model_info <- bind_rows(
    unnest(tidy(models[[1]])) %>%
      mutate(model = "full", convergence = converge(models[[1]])),
    unnest(tidy(models[[2]])%>% filter(group == "fixed" | grepl("mv_wt", term))) %>%
      mutate(model = "null", convergence = converge(models[[2]]))
  ) %>%
    # Now add more concise info (redundantly repeated across rows, can be filtered
    # out later):
    mutate(
      BIC_full = BIC(models[[1]]),
      BIC_null = BIC(models[[2]]),
      BF_BIC = exp( (BIC(models[[2]]) - BIC(models[[1]]) ) / 2),
      conv_full = converge(models[[1]]),
      conv_null = converge(models[[2]]),
      formula_full = myf
    ) %>%
    select(-group)
  # Add the time it took
  mytoc <- toc()
  model_info$t_ellapsed <- mytoc$toc - mytoc$tic
  model_info
}


## ---- message=FALSE, warning=FALSE---------------------------------------
sim_ex <-  simulate_binom(
  Nsubj = 20,
  Nitem = 16,
  fixef_means = fixef_means,
  fixef_sigma = fixef_sigma,
  ranef_sigma_subj = ranef_sigma_subj,
  ranef_sigma_item = ranef_sigma_item,
  type1 = FALSE  # TRUE turns simulations into Type 1 error analysis
  )
# tictoc::tic()
# an_ex <- analyze_simulation(sim_ex, print_model_summaries = TRUE)
# tictoc::toc()
# an_ex


## ---- message=FALSE, warning=FALSE---------------------------------------
# analyze_simulation(sim_ex, print_model_summaries = TRUE)


## ------------------------------------------------------------------------
# Function to analyze many data sets, called with pmap
sim_many <- function(rseed, sim_id, Nsubj, sim_type, incl_col_names, ...) {
  istype1 <- if (sim_type == "type2") {
    FALSE
    } else if (sim_type == "type1") {
    TRUE
    } else{
        stop("sim_type has to be either 'type1' or 'type2'")
    }
  set.seed(rseed)
  df <- simulate_binom(Nsubj, type1 = istype1, ...)
  result <- analyze_simulation(df)
  result$sim_id   <- sim_id
  result$rseed    <- rseed
  result$sim_type <- sim_type
  result$Nsubj    <- Nsubj
  write_csv(
    result,
    path = "power_simulation_results_append.csv",
    col_names = incl_col_names,  # only needed first time it's run
    append = TRUE)
  result
}


## ------------------------------------------------------------------------
# # Number of simulations per N:
# nb_sims <- 5
# # How many different Ns?
# params <- tibble(Nsubj = c(15, 60, 96))
# params <- params[rep(1 : nrow(params), each = nb_sims), ]
# params$rseed <- sample.int(10 ^ 9, size = nrow(params))
# params$sim_type <- "type2"
# params$sim_id <- seq_len(nrow(params))
# head(params, 4)
# tail(params, 4)


## ---- message=FALSE, warning=FALSE---------------------------------------
# tictoc::tic()
# mysims <- pmap(
#   .l = params,
#   .f = sim_many,
#   # fixed simulation parameters (for clarity, comment out the ones passed through params)
#   # Nsubj =,
#   Nitem = 104,
#   fixef_means = fixef_means,
#   fixef_sigma = fixef_sigma,
#   ranef_sigma_subj = ranef_sigma_subj,
#   ranef_sigma_item = ranef_sigma_item
#   ) %>%
#   bind_rows()
# tictoc::toc()


## ------------------------------------------------------------------------
sim_till_aim <- function(
  aim,
  filter_out_nsubj = 0,  # by default this will consider all sample sizes
  filter_out_type = ""   # by default it will do type1 and type2
  ) {
  neg2zero <- function(x) ifelse(x < 0, 0, x)
  # retrieve simulations and check which ones converged (both null and full)
  if (file.exists("power_simulation_results_append.csv")) {
    sims <- read_csv("power_simulation_results_append.csv")
  } else {
    return (FALSE)
  }
  # If there are simulations, get the right info about convergence etc
  sims <- sims %>%
    select(rseed, Nsubj, sim_type, conv_full, conv_null) %>% unique() %>%
    mutate(
      converged = ( is.na(conv_full) & is.na(conv_null) ),
      Nsubj = as.numeric(Nsubj)
      ) %>%
    group_by(Nsubj, sim_type) %>%
    summarise(sims_run = n(), converged = sum(converged)) %>%
    mutate(
      aim = aim,
      needed_neg = aim - converged,
      needed = neg2zero(needed_neg)) %>%
    # Estimated convergence rate at different sample sizes
    left_join(tibble(
      Nsubj =     c(15,   60,   96,   102,  108),
      conv_rate = c(0.13, 0.40, 0.55, 0.56, 0.60))
      ) %>%
    mutate(nsims = round(needed / conv_rate)) %>%
    filter(
      (! Nsubj %in% filter_out_nsubj),
      (sim_type != filter_out_type)
      )
  sims
}


## ---- message=FALSE------------------------------------------------------
# Estimated nnumber of simulations needed to get to an aim
sim_till_aim(10000)
# Filter out simulations with sample size 15 and 60
sim_till_aim(10000, filter_out_nsubj = c(15, 60))
sim_till_aim(10000, filter_out_nsubj = c(15, 60), filter_out_type = "type2")


## ------------------------------------------------------------------------
# Function that checks how many converged simulations there are and runs new
# ones until aim is reached:
run_till_aim <- function(aim, sample_sizes = c(15, 60, 108), ...) {
  tictoc::tic()
  sims <- sim_till_aim(aim, ...)

  # 1st case: no simulations have been run, in which case sims is FALSE
  if (is.logical(sims)) {
    # so then run we'll run one of each sample size to get started
    print("Create file and run one type1 and type2 sim for each sample size")
    params <- tibble(
      Nsubj = rep(sample_sizes, each = 2),
      sim_type = rep(c("type2", "type1"), length.out = 2 * length(sample_sizes))
      )
    # include col names 1st time
    params$incl_col_names <- c(TRUE, rep(FALSE, nrow(params) - 1))
    params$rseed <- sample.int(10 ^ 9, size = nrow(params))  # random seeds
    params$sim_id <- seq_len(nrow(params))
    print("We'll now run:")
    print(params)
    # We can't parellize it yet, bc we want column names to be in row 1:
    pmap(
      .l = params,
      .f = sim_many,
      # fixed simulation parameters (for clarity, comment out the ones passed through params)
      # Nsubj =,
      Nitem = 104,
      fixef_means = fixef_means,
      fixef_sigma = fixef_sigma,
      ranef_sigma_subj = ranef_sigma_subj,
      ranef_sigma_item = ranef_sigma_item
    )

    } else {  # So there is already a file...
      # Check if we have attained the aim
      if (sum(sims$nsims) == 0) {
        return("Done, my friend!!!")
      }
      # There are simulations but we haven't reached the aim.
      # Print how many simulations will be run:
      print("We'll now run:")
      print(sims)
      # set up parameters
      params <- tibble(
        rseed = sample.int(10 ^ 9, size = sum(sims$nsims)),
        Nsubj = rep(sims$Nsubj, sims$nsims),
        sim_type = rep(sims$sim_type, sims$nsims),
        incl_col_names = FALSE
      )
      params$sim_id <- seq_len(nrow(params))

      # Now we can parallelize!
      registerDoParallel(parallel::detectCores())  # use multicore, set to the number of our cores
      foreach (i=1:nrow(params),
               .packages = c("tidyverse", "broom", "lme4", "mvtnorm", "boot", "tictoc"),
               .export = c(
                 "sim_many", "simulate_binom", "analyze_simulation",
                 "fixef_means", "fixef_sigma", "ranef_sigma_subj", "ranef_sigma_item")) %dopar% {

                   pmap(
                     .l = params[i,],
                     .f = sim_many,
                     # fixed simulation parameters (for clarity, comment out the ones passed through params)
                     # Nsubj =,
                     Nitem = 104,
                     fixef_means = fixef_means,
                     fixef_sigma = fixef_sigma,
                     ranef_sigma_subj = ranef_sigma_subj,
                     ranef_sigma_item = ranef_sigma_item
                   )
                 }
      }

  tictoc::toc()

  # recursively runs itself
  run_till_aim(aim, ...)
  }


## ---- message=FALSE------------------------------------------------------


sim_till_aim(10000)
sim_till_aim(10000, filter_out_nsubj = c(60, 108))
run_till_aim(10000, filter_out_nsubj = c(60, 108))


## ------------------------------------------------------------------------
# Save session info to file

end_time <- strftime(Sys.time())
end_time_formatted <- gsub(":", "-", gsub("[ ]", "_", end_time))
fname <- paste("run_", end_time_formatted, ".txt", sep = "")

# session info
writeLines(capture.output(sessionInfo()), fname)
# start and end time
write("\n\nSTART/END:", file = fname, append = TRUE)
write(start_time, file = fname, append = TRUE)
write(end_time, file = fname, append = TRUE)

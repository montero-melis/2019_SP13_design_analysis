## Generate simulated data for the replication.

library("brms")
library("tidyverse")
library("lme4")


# Load simulate_binom() function
source("generate_data_fnc.R")

# Load model fitted to original data so we can set simulation parameters:
fm <- readRDS("sims_etc/bayes_glmm_default.rds")

## Save model estimates in a format that serves as input for the function

# Fixed effects
fixef_means <- fixef(fm)[, 1]
fixef_means
fixef_sigma <- diag(fixef(fm)[, 2]) ^ 2  # we assume uncorrelated diagonal matrix
fixef_sigma

# For random effects by subject we only have an estimate of the random intercept:
VarCorr(fm)
# Let us assume that SDs for the other coefficients are proportional to the
# corresponding fixef SEM in the same way the intercept is:
prop <- VarCorr(fm)$subject$sd[1] / fixef(fm)[1, 2]  # random variability >3 times SEM
ranef_sigma_subj <- diag(fixef(fm)[, 2] * prop) ^ 2  # assume uncorrelated diagonal matrix
ranef_sigma_subj

# Random effects by item we need to completely make up. We let SD for the intercept
# and Movement conditions just be the same as the corresponding by-subject SDs:
ranef_sigma_item <- diag(fixef(fm)[1:2, 2] * prop) ^ 2  # assume uncorrelated diagonal matrix
ranef_sigma_item

# Simulate data
set.seed(654198461)  # make data set replicable
d_list <- simulate_binom(
  Nsubj = 60,
  Nitem = 104,
  fixef_means = fixef_means,
  fixef_sigma = fixef_sigma,
  ranef_sigma_subj = ranef_sigma_subj,
  ranef_sigma_item = ranef_sigma_item,
  full_output = TRUE,
  print_each_step = FALSE
)
# The output consists of a list of data frames:
names(d_list)
d_list
# First list contains the actual data set
d <- d_list[[1]]
head(d)


# Plot the results by subjects in probability space
# (To plot them in log-odds space use the predict function)
plot <- d %>% group_by(subject, movement, word_type) %>%
  summarise(M = mean(Error),
            Sum = sum(Error)) %>%
  mutate(sbj_wtype = paste(subject, word_type, sep = "_")) %>%
  ggplot(aes(x = movement, colour = word_type, y = M)) +
  stat_summary(fun.y = mean, geom = "point", size = 4,
               position = position_dodge(width = .25)) +
  stat_summary(fun.y = mean, geom = "line", aes(group = word_type), size = 2,
               position = position_dodge(width = .25)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2, size = 2,
               position = position_dodge(width = .25)) +
  ylim(0, 1) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position="top")
plot

plot +
  geom_jitter(height = 0, width = .05, alpha = .5) +
  geom_line(aes(group = sbj_wtype), alpha = .25)


# Coding scheme before fitting GLMM to data
contrasts(d$movement) <- contr.sum(2)
contrasts(d$word_type) <- contr.sum(2)


# Fit different lme4 models with increasingly complex random effect structure:

fm0 <- glmer(
  Error ~ 1 + movement * word_type +
    (1 | subject) +
    (1 | item),
  data = d,
  family = "binomial"
  )
summary(fm0)

fm1a <- glmer(
  Error ~ 1 + movement * word_type +
    (1 + movement | subject) +
    (1 | item),
  data = d,
  family = "binomial"
)
summary(fm1a)

fm1b <- glmer(
  Error ~ 1 + movement * word_type +
    (1 | subject) +
    (1 + movement | item),
  data = d,
  family = "binomial"
)
summary(fm1b)

fm2 <- glmer(
  Error ~ 1 + movement * word_type +
  	(1 + movement | subject) +
    (1 + movement | item),
  data = d,
  family = "binomial"
)
summary(fm2)


fm3 <- glmer(
  Error ~ 1 + movement * word_type +
  	(1 + movement + word_type | subject) +
    (1 + movement | item),
  data = d,
  family = "binomial"
)
summary(fm3)


fm4 <- glmer(
  Error ~ 1 + movement * word_type +
  	(1 + movement * word_type | subject) +
    (1 + movement | item),
  data = d,
  family = "binomial"
)
summary(fm4)

# Model summaries
summary(fm0)
summary(fm1a)
summary(fm1b)
summary(fm2)
summary(fm3)
summary(fm4)

# Compare models with anova function (LRT):
anova(fm0, fm1a, fm1b, fm2, fm3, fm4)

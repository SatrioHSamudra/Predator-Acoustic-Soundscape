library(glmmTMB)
library(readxl)
library(dplyr)
library(MASS)
library(tidyverse)
library(lme4)
library(car)
library(emmeans)
library(DHARMa)
library(MuMIn)
library(writexl)
library(scales)

setwd("~/Documents/2025_UCL/Dissertation/File")

fish <- read_excel("merged_hydromoth2.xlsx", sheet = "Sheet1", col_names = TRUE) %>% as.data.frame()

# Remove timestamp 55s
fish <- fish[fish$timestamp_s != 55, ]

# Remove rows where species_r or species_l contains 'turtle' or 'humphead snapper'
fish <- fish[!grepl("turtle|humphead snapper", fish$species_r, ignore.case = TRUE) &
             !grepl("turtle|humphead snapper", fish$species_l, ignore.case = TRUE), ]

# Convert fish metrics to scale
fish$fish_diversity_scaled <- scale(fish$`Fish Diversity (Count of Species)`)
fish$fish_biomass_scaled <- scale(fish$`Fish Total Biomass (kg)`)
fish$fish_abundance_scaled <- scale(fish$`Fish Abundance (Sum of MaxN)`)

# Filter and drop unused factor levels before modeling
fish_subset <- droplevels(subset(fish, predator_occurrence == "present"))

# MAKE A MODEL FOR PHONIC RICHNESS
fish_clean0 <- na.omit(fish_subset[, c("no_sonotypes_no_snap", "habitat",
                                       "fish_abundance_scaled", "fish_diversity_scaled",
                                       "fish_biomass_scaled", "predator_abundance",
                                       "predator_type", "site", "unit")])

model_combine1 <- glmmTMB(
  no_sonotypes_no_snap ~ habitat + 
    fish_abundance_scaled + 
    fish_diversity_scaled +
    fish_biomass_scaled + predator_abundance + predator_type +
    (1|site) + (1|unit),
  data = fish_clean0,
  family = "poisson", na.action = "na.fail"
)

dredge_results_new <- dredge(model_combine1, fixed = ~ predator_abundance + predator_type)
dredge_results_new

avg_model_new <- model.avg(dredge_results_new, subset = delta < 2)
summary(avg_model_new)

# Step 1: Generate prediction grid for predator_abundance
pred_grid <- data.frame(
  predator_abundance = seq(
    from = min(fish_clean0$predator_abundance, na.rm = TRUE),
    to   = max(fish_clean0$predator_abundance, na.rm = TRUE),
    length.out = 200
  )
)

# Step 2: Hold other predictors at 0 (or mean, if preferred)
pred_grid$fish_biomass_scaled <- 0
pred_grid$fish_diversity_scaled <- 0
pred_grid$fish_abundance_scaled <- 0

# Step 3: Extract averaged coefficients from avg_model
coefs <- summary(avg_model_new)$coefmat.full

# Step 4: Compute linear predictor using avg_model
pred_grid$lp <- with(pred_grid, 
                     coefs["cond((Int))", "Estimate"] +
                       coefs["cond(fish_abundance_scaled)", "Estimate"] * fish_abundance_scaled +
                       coefs["cond(fish_diversity_scaled)", "Estimate"] * fish_diversity_scaled +
                       coefs["cond(fish_biomass_scaled)", "Estimate"] * fish_biomass_scaled +
                       coefs["cond(predator_abundance)", "Estimate"] * predator_abundance
)

# Step 6: Convert to predicted phonic richness (Poisson scale)
pred_grid$fit <- exp(pred_grid$lp)

# Step 7: Compute CI (approximate based only on predator_abundance)
se_div <- coefs["cond(predator_abundance)", "Std. Error"]
pred_grid$se_fit <- abs(pred_grid$predator_abundance) * se_div
pred_grid$lower <- exp(pred_grid$lp - 1.96 * pred_grid$se_fit)
pred_grid$upper <- exp(pred_grid$lp + 1.96 * pred_grid$se_fit)

# Reference value: predicted phonic richness at mean fish diversity (i.e., scaled = 0)
ref_value <- pred_grid$fit[which.min(abs(pred_grid$predator_abundance - 0))]

# Calculate % change relative to reference
pred_grid$fit_pctchange   <- (pred_grid$fit - ref_value) / ref_value * 100
pred_grid$lower_pctchange <- (pred_grid$lower - ref_value) / ref_value * 100
pred_grid$upper_pctchange <- (pred_grid$upper - ref_value) / ref_value * 100


# MAKE A MODEL FOR SONOTYPE OCCURRENCE PROBABILITY
fish_clean2 <- na.omit(fish_subset[, c("sonotypes_no_snap", "habitat",
                                       "fish_abundance_scaled", "fish_diversity_scaled",
                                       "fish_biomass_scaled", "predator_abundance",
                                       "predator_type", "site", "unit")])

model_combine2 <- glmmTMB(
  sonotypes_no_snap ~ habitat + 
    fish_abundance_scaled + 
    fish_diversity_scaled +
    fish_biomass_scaled + predator_abundance + predator_type +
    (1|site) + (1|unit),
  data = fish_clean2,
  family = "binomial", na.action = "na.fail"
)

dredge_results_prob <- dredge(model_combine2, fixed = ~ predator_abundance + predator_type)
dredge_results_prob

avg_model_prob <- model.avg(dredge_results_prob, subset = delta < 2)
summary(avg_model_prob)

# Step 1: Generate prediction grid for predator_abundance
pred_grid1 <- data.frame(
  predator_abundance = seq(
    from = min(fish_clean2$predator_abundance, na.rm = TRUE),
    to   = max(fish_clean2$predator_abundance, na.rm = TRUE),
    length.out = 200
  )
)

# Step 2: Hold other predictors at 0 (or mean, if preferred)
pred_grid1$fish_biomass_scaled <- 0
pred_grid1$fish_diversity_scaled <- 0
pred_grid1$fish_abundance_scaled <- 0

# Step 3: Extract averaged coefficients from avg_model
coefs <- summary(avg_model_prob)$coefmat.full

# Step 4: Compute linear predictor using avg_model
pred_grid1$lp <- with(pred_grid1, 
                     coefs["cond((Int))", "Estimate"] +
                       coefs["cond(fish_abundance_scaled)", "Estimate"] * fish_abundance_scaled +
                       coefs["cond(fish_diversity_scaled)", "Estimate"] * fish_diversity_scaled +
                       coefs["cond(fish_biomass_scaled)", "Estimate"] * fish_biomass_scaled +
                       coefs["cond(predator_abundance)", "Estimate"] * predator_abundance
)

# Step 6: Convert to predicted phonic richness (Poisson scale)
pred_grid1$fit <- exp(pred_grid1$lp)

# Step 7: Compute CI (approximate based only on predator_abundance)
se_div <- coefs["cond(predator_abundance)", "Std. Error"]
pred_grid1$se_fit <- abs(pred_grid1$predator_abundance) * se_div
pred_grid1$lower <- exp(pred_grid1$lp - 1.96 * pred_grid1$se_fit)
pred_grid1$upper <- exp(pred_grid1$lp + 1.96 * pred_grid1$se_fit)

# Reference value: predicted phonic richness at mean fish diversity (i.e., scaled = 0)
ref_value <- pred_grid1$fit[which.min(abs(pred_grid1$predator_abundance - 0))]

# Calculate % change relative to reference
pred_grid1$fit_pctchange   <- (pred_grid1$fit - ref_value) / ref_value * 100
pred_grid1$lower_pctchange <- (pred_grid1$lower - ref_value) / ref_value * 100
pred_grid1$upper_pctchange <- (pred_grid1$upper - ref_value) / ref_value * 100


# MAKE A MODEL FOR SNAPPING SHRIMP
fish_clean3 <- na.omit(fish_subset[, c("snap", "habitat",
                                       "fish_abundance_scaled", "fish_diversity_scaled",
                                       "fish_biomass_scaled", "predator_abundance",
                                       "predator_type", "site", "unit")])

model_combine3 <- glmmTMB(
  snap ~ habitat + 
    fish_abundance_scaled + 
    fish_diversity_scaled +
    fish_biomass_scaled + predator_abundance + predator_type +
    (1|site) + (1|unit),
  data = fish_clean3,
  family = "binomial", na.action = "na.fail"
)

dredge_results_snap <- dredge(model_combine3, fixed = ~ predator_abundance + predator_type)
dredge_results_snap

avg_model_snap <- model.avg(dredge_results_snap, subset = delta < 2)
summary(avg_model_snap) # Neither predator abundance nor type featured in the model-averaged estimates
                        # Hence not visualise on the final graph


# Combine all acoustic metrics
pred_grid$metric <- "Phonic Richness"
pred_grid1$metric <- "Sonotype Occurrence"

combined_df <- bind_rows(
  dplyr::select(pred_grid1, predator_abundance, fit_pctchange, lower_pctchange, upper_pctchange, metric),
  dplyr::select(pred_grid,  predator_abundance, fit_pctchange, lower_pctchange, upper_pctchange, metric)
)

combined_df$metric <- factor(combined_df$metric, levels = c(
  "Sonotype Occurrence",
  "Phonic Richness"
))

ggplot(combined_df, aes(x = predator_abundance, y = fit_pctchange)) +
  geom_ribbon(aes(ymin = lower_pctchange, ymax = upper_pctchange), fill = "#1f78b4", alpha = 0.4) +
  geom_line(color = "black", size = 1.2) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2,
             labeller = as_labeller(c(
               "Sonotype Occurrence" = "A) Sonotype Occurrence",
               "Phonic Richness" = "B) Phonic Richness"))) +
  annotate("text", x = 6, y = 0, label = "*", size = 10, hjust = 1, vjust = 0.5) +
  labs(
    x = "Predator Abundance",
    y = "Percent Change (%)"
  ) +
  scale_y_continuous(
    limits = c(-60, 10),
    breaks = seq(-60, 10, 20)
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12.5, hjust = -0.045, face = "bold"),
    text = element_text(size = 12.5),
    axis.text = element_text(size = 12.5, colour = "black"),
    axis.ticks.length = unit(0.1, "cm"),
    plot.caption = element_text(hjust = 0.5, size = 7, face = "italic"),
    panel.background = element_rect(color = "black", size = 1)
  )

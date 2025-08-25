#############################################################################
# Predator abundance and type effects on acoustic metrics
# Author: PTSH0
# Project: MSc Biodiversity and Global Change dissertation
# Date: 2025-08-01
# 
# This script reproduces GLMM analyses of predator effects (abundance, type)
# on fish- and shrimp-generated acoustic metrics.
#############################################################################

# ---- Load  packages ----
library(readxl)
library(dplyr)
library(MASS)
library(tidyverse)
library(glmmTMB) 
library(emmeans)
library(DHARMa)
library(MuMIn)
library(scales)
library(ggplot2)

# ---- Set working directory ----
# Use project root or relative path for reproducibility
# setwd("~/Documents/2025_UCL/Dissertation/File")
# Prefer relative path for GitHub projects:
data_location <- "data/merged_hydromoth2.xlsx"
output_location <- "outputs/figures/predator_abundance.png"

# ---- Import and clean data ----
fish <- read_excel(data_location, sheet = "Sheet1", col_names = TRUE) %>% 
  as.data.frame()

fish <- fish[fish$timestamp_s != 55, ]  # remove 55s
fish <- fish[
  !grepl("turtle|humphead snapper", fish$species_r, ignore.case = TRUE) &
  !grepl("turtle|humphead snapper", fish$species_l, ignore.case = TRUE),
]

# scale fish metric variables
fish$fish_diversity_scaled  <- scale(fish$`Fish Diversity (Count of Species)`)
fish$fish_biomass_scaled    <- scale(fish$`Fish Total Biomass (kg)`)
fish$fish_abundance_scaled  <- scale(fish$`Fish Abundance (Sum of MaxN)`)

# keep only predator presence
fish_subset <- droplevels(subset(fish, predator_occurrence == "present"))

# ---- Run model dredging for all acoustic metrics ----
run_model <- function(data, response, family, label) {
  clean <- na.omit(data[, c(response, "habitat",
                            "fish_abundance_scaled", "fish_diversity_scaled",
                            "fish_biomass_scaled", "predator_abundance",
                            "predator_type", "site", "unit")])
  
  model <- glmmTMB(
    as.formula(paste(response, "~ habitat + fish_abundance_scaled + fish_diversity_scaled + fish_biomass_scaled + predator_abundance + predator_type + (1|site) + (1|unit)")),
    data = clean,
    family = family,
    na.action = "na.fail"
  )
  
  dredge_res <- dredge(model, fixed = ~ predator_abundance + predator_type)
  avg_model  <- model.avg(dredge_res, subset = delta < 2)
  
  # prediction grid
  pred_grid <- data.frame(
    predator_abundance = seq(min(clean$predator_abundance, na.rm = TRUE),
                             max(clean$predator_abundance, na.rm = TRUE),
                             length.out = 200)
  )
  pred_grid$fish_biomass_scaled   <- 0
  pred_grid$fish_diversity_scaled <- 0
  pred_grid$fish_abundance_scaled <- 0
  
  coefs <- summary(avg_model)$coefmat.full
  
  # linear predictor for explanatory variables
  pred_grid$lp <- with(pred_grid, 
                       coefs["cond((Int))", "Estimate"] +
                         coefs["cond(fish_abundance_scaled)", "Estimate"] * fish_abundance_scaled +
                         coefs["cond(fish_diversity_scaled)", "Estimate"] * fish_diversity_scaled +
                         coefs["cond(fish_biomass_scaled)", "Estimate"] * fish_biomass_scaled +
                         coefs["cond(predator_abundance)", "Estimate"] * predator_abundance
  )
  
  pred_grid$fit   <- exp(pred_grid$lp)
  se_div          <- coefs["cond(predator_abundance)", "Std. Error"]
  pred_grid$se_fit <- abs(pred_grid$predator_abundance) * se_div
  pred_grid$lower <- exp(pred_grid$lp - 1.96 * pred_grid$se_fit)
  pred_grid$upper <- exp(pred_grid$lp + 1.96 * pred_grid$se_fit)
  
  # % change relative to reference
  ref_value <- pred_grid$fit[which.min(abs(pred_grid$predator_abundance - 0))]
  pred_grid$fit_pctchange   <- (pred_grid$fit   - ref_value) / ref_value * 100
  pred_grid$lower_pctchange <- (pred_grid$lower - ref_value) / ref_value * 100
  pred_grid$upper_pctchange <- (pred_grid$upper - ref_value) / ref_value * 100
  
  pred_grid$metric <- label
  return(pred_grid)
}

# ---- Run final models ----
pred_sonotype <- run_model(fish_subset, "sonotypes_no_snap", "binomial", "Sonotype Occurrence")
pred_phonic   <- run_model(fish_subset, "no_sonotypes_no_snap", "poisson",  "Phonic Richness")
# snapping shrimp model not plotted since predator effects not retained in the best-performing models

combined_df <- bind_rows(pred_sonotype, pred_phonic)
combined_df$metric <- factor(combined_df$metric, levels = c("Sonotype Occurrence", "Phonic Richness"))

# ---- Plot final graph ----
pred_abundance <- ggplot(combined_df, aes(x = predator_abundance, y = fit_pctchange)) +
  geom_ribbon(aes(ymin = lower_pctchange, ymax = upper_pctchange),
              fill = "#1f78b4", alpha = 0.4) +
  geom_line(color = "black", size = 1.2) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2,
             labeller = as_labeller(c(
               "Sonotype Occurrence" = "A) Sonotype Occurrence",
               "Phonic Richness"     = "B) Phonic Richness"
             ))) +
  annotate("text", x = 6, y = 0, label = "*", size = 10,
           hjust = 1, vjust = 0.5) +
  labs(x = "Predator Abundance", y = "Percent Change (%)") +
  scale_y_continuous(limits = c(-60, 10), breaks = seq(-60, 10, 20)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12.5, hjust = -0.045, face = "bold"),
    text = element_text(size = 12.5),
    axis.text = element_text(size = 12.5, colour = "black"),
    axis.ticks.length = unit(0.1, "cm"),
    panel.background = element_rect(color = "black", size = 1)
  )
pred_abundance

# save plot
ggsave(output_location, pred_abundance, width = 8, height = 5, dpi = 300)

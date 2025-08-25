########################################################
# Predator occurrence effects on acoustic metrics
# Author: PTSH0
# Project: MSc Biodiversity and Global Change thesis
# Date: 2025-08-01
#
# This script reproduces GLMM analyses comparing three
# acoustic metrics between predator present vs. absent.
########################################################

# ---- Load packages ----
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
library(patchwork)

# ---- Set working directory and paths ----
# setwd("~/Documents/2025_UCL/Dissertation/File") <- my project directory
# Alternative path for GitHub projects
data_location <- "data/merged_hydromoth2.xlsx"
output_location <- "outputs/figures/predator_occurrence.png"

# ---- Import and clean data ----
fish <- read_excel("merged_hydromoth2.xlsx", sheet = "Sheet1", col_names = TRUE) %>% as.data.frame()

# Remove timestamp = 55s
fish <- fish[fish$timestamp_s != 55, ]

# Remove rows where either species col has turtle or humphead snapper
fish <- fish[
  !grepl("turtle|humphead snapper", fish$species_r, ignore.case = TRUE) &
  !grepl("turtle|humphead snapper", fish$species_l, ignore.case = TRUE),
]

# Convert "na" string to proper NA
fish <- fish %>%
  mutate(
    predator_type_1 = na_if(predator_type_1, "na"),
    predator_type_2 = na_if(predator_type_2, "na")
  )

# Harmonise predator type columns into a new variable
fish <- fish %>%
  mutate(
    predator_type_1 = tolower(str_trim(as.character(predator_type_1))),
    predator_type_2 = tolower(str_trim(as.character(predator_type_2))),
    predator_type_1 = na_if(predator_type_1, "na"),
    predator_type_2 = na_if(predator_type_2, "na"),
    predator_type = case_when(
      predator_type_1 == "shark" & (is.na(predator_type_2) | predator_type_2 == "shark") ~ "shark",
      predator_type_2 == "shark" & (is.na(predator_type_1) | predator_type_1 == "shark") ~ "shark",
      predator_type_1 == "non-shark" & (is.na(predator_type_2) | predator_type_2 == "non-shark") ~ "non-shark",
      predator_type_2 == "non-shark" & (is.na(predator_type_1) | predator_type_1 == "non-shark") ~ "non-shark",
      ((predator_type_1 == "shark" & predator_type_2 == "non-shark") |
         (predator_type_1 == "non-shark" & predator_type_2 == "shark")) ~ "mix",
      TRUE ~ NA_character_
    )
  )

# ---- Scale fish metrics variables ----
fish$fish_diversity_scaled  <- scale(fish$`Fish Diversity (Count of Species)`)
fish$fish_biomass_scaled    <- scale(fish$`Fish Total Biomass (kg)`)
fish$fish_abundance_scaled  <- scale(fish$`Fish Abundance (Sum of MaxN)`)

# Keep only predator present/absent
fish_clean <- droplevels(subset(fish, predator_occurrence %in% c("present", "absent")))
fish_clean$predator_occurrence_num <- ifelse(fish_clean$predator_occurrence == "present", 1, 0)

# ---- Incorporate log10 transformation for plotting ----
signed_log10_trans <- function() {
  trans_new(
    name = "signed_log10",
    transform = function(x) ifelse(x > 0, log10(x), ifelse(x < 0, -log10(abs(x)), 0)),
    inverse  = function(x) ifelse(x > 0, 10^x, ifelse(x < 0, -10^(-x), 0)),
    breaks   = function(x) {
      brks <- c(-1000, -100, -10, -1, 0, 1, 10, 100, 1000)
      brks[brks >= min(x, na.rm = TRUE) & brks <= max(x, na.rm = TRUE)]
    },
    format   = label_percent(scale = 1)
  )
}

# ---- Run model dredging for all acoustic metrics ----
run_predator_models <- function(data, response, family, label) {
  global_model <- glmmTMB(
    as.formula(paste(response, "~ habitat + fish_diversity_scaled + fish_abundance_scaled + 
                     fish_biomass_scaled + predator_occurrence + (1|site) + (1|unit)")),
    data = data,
    family = family
  )
  options(na.action = "na.fail")
  dredge_results <- dredge(global_model, fixed = ~ predator_occurrence)
  avg_model <- model.avg(dredge_results, subset = delta < 2)
  
  # extract top models + weights
  top_models <- get.models(dredge_results, subset = delta < 2)
  model_table <- as.data.frame(dredge_results) %>% filter(delta < 2)
  weights <- model_table$weight
  
  # compute emmeans for predator_occurrence
  emm_list <- lapply(top_models, function(m) {
    tryCatch(as.data.frame(emmeans(m, ~ predator_occurrence, type = "response")),
             error = function(e) NULL)
  })
  valid_idx <- which(sapply(emm_list, function(x) !is.null(x) && nrow(x) > 0))
  emm_dfs <- emm_list[valid_idx]
  weights_valid <- weights[valid_idx]
  
  # standardise col name: always rename response column to "prob"
  emm_dfs <- lapply(emm_dfs, function(df) {
    names(df)[1] <- "predator_occurrence"
    if ("prob" %in% names(df)) {
      df <- df
    } else if ("rate" %in% names(df)) {
      names(df)[names(df) == "rate"] <- "prob"
    }
    df
  })
  
  emm_combined <- bind_rows(Map(function(df, w){df$weight <- w; df}, emm_dfs, weights_valid))
  
  # summarise across models
  result <- emm_combined %>% 
    group_by(predator_occurrence) %>% 
    summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop")
  
  baseline <- result %>% filter(predator_occurrence == "absent") %>% pull(prob)
  
  # clean dataframe with 95% + 75% CIs
  outcome <- result %>%
    filter(predator_occurrence %in% c("absent","present")) %>%
    mutate(
      log_fold_change = log10(prob / baseline),
      se_log = SE / prob / log(10),
      
      # 95% CI
      lower_log_95 = log_fold_change - qnorm(0.975) * se_log,
      upper_log_95 = log_fold_change + qnorm(0.975) * se_log,
      
      # 75% CI
      lower_log_75 = log_fold_change - qnorm(0.875) * se_log,
      upper_log_75 = log_fold_change + qnorm(0.875) * se_log,
      
      # back-transform
      pct_change = (10^log_fold_change - 1) * 100,
      lower_95   = (10^lower_log_95 - 1) * 100,
      upper_95   = (10^upper_log_95 - 1) * 100,
      lower_75   = (10^lower_log_75 - 1) * 100,
      upper_75   = (10^upper_log_75 - 1) * 100,
      
      group_color  = case_when(
        predator_occurrence == "absent"  ~ "#9E9AC8",
        predator_occurrence == "present" ~ "#74C476"
      ),
      Indicator = label
    )
  return(outcome)
}

# ---- Run final models for each acoustic metric ----
df1 <- run_predator_models(fish_clean, "no_sonotypes_no_snap", "poisson",  "Phonic Richness")
df2 <- run_predator_models(fish_clean, "sonotypes_no_snap",   "binomial", "Sonotype Occurrence")
df3 <- run_predator_models(fish_clean, "snap",                "binomial", "Snapping Shrimp Activity")

combined_df <- bind_rows(df1, df2, df3) %>%
  mutate(
    Indicator = factor(Indicator, levels = c("Sonotype Occurrence", "Phonic Richness", "Snapping Shrimp Activity")),
    predator_occurrence = factor(predator_occurrence, levels = c("absent", "present"))
  )

# ---- Plot the graph ----
pred_occurrence <- ggplot(combined_df, aes(x = pct_change, y = predator_occurrence, color = predator_occurrence)) +
  # 95% CI
  geom_errorbarh(aes(xmin = lower_95, xmax = upper_95), height = 0, size = 1.3) +
  # 75% CI
  geom_errorbarh(aes(xmin = lower_75, xmax = upper_75), height = 0, size = 2) +
  # Points
  geom_point(size = 4) +
  facet_wrap(~Indicator, ncol = 1, scales = "free_y") +
  scale_x_continuous(trans = signed_log10_trans(), limits = c(-100, 1000)) +
  scale_y_discrete(limits = rev) +
  labs(x = "Decreased over absent (%)", y = "Predator") +
  scale_color_manual(values = c("absent" = "#9E9AC8", "present" = "#74C476"),
                     name = "Predator") +
  guides(color = guide_legend(reverse = TRUE)) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", hjust = 0),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
pred_occurrence

ggsave(output_path, pred_occurrence, width = 6, height = 8, dpi = 300)

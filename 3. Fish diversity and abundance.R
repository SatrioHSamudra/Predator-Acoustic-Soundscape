#############################################################
# Fish diversity and abundance effects on acoustic metrics
# Author: PTSH0
# Project: MSc Biodiversity and Global Change thesis
# Date: 2025-08-01
#
# This script reproduces GLMM analyses comparing three
# acoustic metrics on different fish diversity and abundance
#############################################################

# ---- Load packages ----
library(glmmTMB)
library(readxl)
library(dplyr)
library(tidyr)
library(MuMIn)
library(ggplot2)
library(readr)
library(grid)  

# ---- Set paths ----
data_location <- file.path("data", "dataset_final.xlsx")
output_location   <- file.path("results", "figures")

# ---- Read data and filter data----
dat <- read_excel(data_path, sheet = "Sheet1", col_names = TRUE) %>% as.data.frame()

# Keep only rows with predator present/absent
fish_clean <- dat %>%
  filter(predator_occurrence %in% c("present", "absent")) %>%
  droplevels()

# Factors for random effects + habitat
if (!is.factor(fish_clean$site))   fish_clean$site   <- factor(fish_clean$site)
if (!is.factor(fish_clean$unit))   fish_clean$unit   <- factor(fish_clean$unit)
if (!is.factor(fish_clean$habitat)) fish_clean$habitat <- factor(fish_clean$habitat)

# Required for MuMIn::dredge() to work
op <- options(na.action = "na.fail")
on.exit(options(op), add = TRUE)

# -------------------------------------------------------------------
# PART A — Fish Diversity (scaled) as focal predictor
# -------------------------------------------------------------------

# 1) Fit global models (PR = Phonic Richness, SO = Sonotype Occurrence, SS = Snapping Shrimp)
m_PR_div <- glmmTMB(
  no_sonotypes_no_snap ~ habitat + fish_diversity_scaled + fish_abundance_scaled + fish_biomass_scaled + predator_occurrence +
    (1 | site) + (1 | unit),
  data = fish_clean,
  family = poisson
)

m_SO_div <- glmmTMB(
  sonotypes_no_snap ~ habitat + fish_diversity_scaled + fish_abundance_scaled + fish_biomass_scaled + predator_occurrence +
    (1 | site) + (1 | unit),
  data = fish_clean,
  family = binomial
)

m_SS_div <- glmmTMB(
  snap ~ habitat + fish_diversity_scaled + fish_abundance_scaled + fish_biomass_scaled + predator_occurrence +
    (1 | site) + (1 | unit),
  data = fish_clean,
  family = binomial
)

# 2) Dredge and average top models (delta AIC < 2)
dr_PR_div <- dredge(m_PR_div, fixed = ~ predator_occurrence)
dr_SO_div <- dredge(m_SO_div, fixed = ~ predator_occurrence)
dr_SS_div <- dredge(m_SS_div, fixed = ~ predator_occurrence)

avg_PR_div <- model.avg(dr_PR_div, subset = delta < 2)
avg_SO_div <- model.avg(dr_SO_div, subset = delta < 2)
avg_SS_div <- model.avg(dr_SS_div, subset = delta < 2)

c_PR_div <- summary(avg_PR_div)$coefmat.full
c_SO_div <- summary(avg_SO_div)$coefmat.full
c_SS_div <- summary(avg_SS_div)$coefmat.full

# 3) Build prediction grid for fish_diversity_scaled
x_div <- seq(
  from = min(fish_clean$fish_diversity_scaled, na.rm = TRUE),
  to   = max(fish_clean$fish_diversity_scaled, na.rm = TRUE),
  length.out = 200
)

p_present <- mean(fish_clean$predator_occurrence == "present", na.rm = TRUE)

# Pick out intercepts, slopes, and SEs from the model results
get_b0   <- function(cm) cm[grep("Intercept", rownames(cm)), "Estimate"]
get_b1   <- function(cm, varname) cm[grep(varname, rownames(cm)), "Estimate"]
get_seb1 <- function(cm, varname) cm[grep(varname, rownames(cm)), "Std. Error"]
get_pred <- function(cm) { i <- grep("predator_occurrencepresent", rownames(cm)); if(length(i)==0) 0 else cm[i, "Estimate"] }

# phonic richness ~ Poisson distribution
b0  <- get_b0(c_PR_div)
b1  <- get_b1(c_PR_div, "fish_diversity_scaled")
se1 <- get_seb1(c_PR_div, "fish_diversity_scaled")
bp  <- get_pred(c_PR_div)

lp_PR_div <- b0 + b1 * x_div + p_present * bp
fit_PR_div <- exp(lp_PR_div)
se_lp_PR_div <- abs(x_div) * se1
lower_PR_div <- exp(lp_PR_div - 1.96 * se_lp_PR_div)
upper_PR_div <- exp(lp_PR_div + 1.96 * se_lp_PR_div)

ref_idx <- which.min(abs(x_div - 0))
ref_val <- fit_PR_div[ref_idx]
df_PR_div <- tibble(
  x = x_div,
  fit = fit_PR_div,
  lower = lower_PR_div,
  upper = upper_PR_div,
  metric = "Phonic Richness"
) %>%
  mutate(
    fit_pctchange   = (fit - ref_val)/ref_val*100,
    lower_pctchange = (lower - ref_val)/ref_val*100,
    upper_pctchange = (upper - ref_val)/ref_val*100
  )

# sonotype occurrence ~ Binomial distribution
b0  <- get_b0(c_SO_div)
b1  <- get_b1(c_SO_div, "fish_diversity_scaled")
se1 <- get_seb1(c_SO_div, "fish_diversity_scaled")
bp  <- get_pred(c_SO_div)

lp_SO_div <- b0 + b1 * x_div + p_present * bp
fit_SO_div <- plogis(lp_SO_div)
se_lp_SO_div <- abs(x_div) * se1
lower_SO_div <- plogis(lp_SO_div - 1.96 * se_lp_SO_div)
upper_SO_div <- plogis(lp_SO_div + 1.96 * se_lp_SO_div)

ref_val <- fit_SO_div[ref_idx]
df_SO_div <- tibble(
  x = x_div,
  fit = fit_SO_div,
  lower = lower_SO_div,
  upper = upper_SO_div,
  metric = "Sonotype Occurrence"
) %>%
  mutate(
    fit_pctchange   = (fit - ref_val)/ref_val*100,
    lower_pctchange = (lower - ref_val)/ref_val*100,
    upper_pctchange = (upper - ref_val)/ref_val*100
  )

# snapping shrimp ~ Binomial distribution
b0  <- get_b0(c_SS_div)
b1  <- get_b1(c_SS_div, "fish_diversity_scaled")
se1 <- get_seb1(c_SS_div, "fish_diversity_scaled")
bp  <- get_pred(c_SS_div)

lp_SS_div <- b0 + b1 * x_div + p_present * bp
fit_SS_div <- plogis(lp_SS_div)
se_lp_SS_div <- abs(x_div) * se1
lower_SS_div <- plogis(lp_SS_div - 1.96 * se_lp_SS_div)
upper_SS_div <- plogis(lp_SS_div + 1.96 * se_lp_SS_div)

ref_val <- fit_SS_div[ref_idx]
df_SS_div <- tibble(
  x = x_div,
  fit = fit_SS_div,
  lower = lower_SS_div,
  upper = upper_SS_div,
  metric = "Snapping Shrimp Activity"
) %>%
  mutate(
    fit_pctchange   = (fit - ref_val)/ref_val*100,
    lower_pctchange = (lower - ref_val)/ref_val*100,
    upper_pctchange = (upper - ref_val)/ref_val*100
  )

div_all <- bind_rows(df_SO_div, df_PR_div, df_SS_div) %>%
  mutate(metric = factor(metric,
                         levels = c("Sonotype Occurrence","Phonic Richness","Snapping Shrimp Activity")),
         focal = "Fish Diversity (scaled)")

# -------------------------------------------------------------------
# Part B — Fish Abundance (scaled) as focal predictor
# -------------------------------------------------------------------

m_PR_abn <- glmmTMB(
  no_sonotypes_no_snap ~ habitat + fish_diversity_scaled + fish_abundance_scaled + fish_biomass_scaled + predator_occurrence +
    (1 | site) + (1 | unit),
  data = fish_clean,
  family = poisson
)

m_SO_abn <- glmmTMB(
  sonotypes_no_snap ~ habitat + fish_diversity_scaled + fish_abundance_scaled + fish_biomass_scaled + predator_occurrence +
    (1 | site) + (1 | unit),
  data = fish_clean,
  family = binomial
)

m_SS_abn <- glmmTMB(
  snap ~ habitat + fish_diversity_scaled + fish_abundance_scaled + fish_biomass_scaled + predator_occurrence +
    (1 | site) + (1 | unit),
  data = fish_clean,
  family = binomial
)

dr_PR_abn <- dredge(m_PR_abn, fixed = ~ predator_occurrence)
dr_SO_abn <- dredge(m_SO_abn, fixed = ~ predator_occurrence)
dr_SS_abn <- dredge(m_SS_abn, fixed = ~ predator_occurrence)

avg_PR_abn <- model.avg(dr_PR_abn, subset = delta < 2)
avg_SO_abn <- model.avg(dr_SO_abn, subset = delta < 2)
avg_SS_abn <- model.avg(dr_SS_abn, subset = delta < 2)

c_PR_abn <- summary(avg_PR_abn)$coefmat.full
c_SO_abn <- summary(avg_SO_abn)$coefmat.full
c_SS_abn <- summary(avg_SS_abn)$coefmat.full

x_abn <- seq(
  from = min(fish_clean$fish_abundance_scaled, na.rm = TRUE),
  to   = max(fish_clean$fish_abundance_scaled, na.rm = TRUE),
  length.out = 200
)

# PR ~ Poisson distribution
b0  <- get_b0(c_PR_abn)
b1  <- get_b1(c_PR_abn, "fish_abundance_scaled")
se1 <- get_seb1(c_PR_abn, "fish_abundance_scaled")
bp  <- get_pred(c_PR_abn)

lp_PR_abn <- b0 + b1 * x_abn + p_present * bp
fit_PR_abn <- exp(lp_PR_abn)
se_lp_PR_abn <- abs(x_abn) * se1
lower_PR_abn <- exp(lp_PR_abn - 1.96 * se_lp_PR_abn)
upper_PR_abn <- exp(lp_PR_abn + 1.96 * se_lp_PR_abn)

ref_idx <- which.min(abs(x_abn - 0))
ref_val <- fit_PR_abn[ref_idx]
df_PR_abn <- tibble(
  x = x_abn, fit = fit_PR_abn, lower = lower_PR_abn, upper = upper_PR_abn,
  metric = "Phonic Richness"
) %>%
  mutate(
    fit_pctchange   = (fit - ref_val)/ref_val*100,
    lower_pctchange = (lower - ref_val)/ref_val*100,
    upper_pctchange = (upper - ref_val)/ref_val*100
  )

# SO ~ Binomial distribution
b0  <- get_b0(c_SO_abn)
b1  <- get_b1(c_SO_abn, "fish_abundance_scaled")
se1 <- get_seb1(c_SO_abn, "fish_abundance_scaled")
bp  <- get_pred(c_SO_abn)

lp_SO_abn <- b0 + b1 * x_abn + p_present * bp
fit_SO_abn <- plogis(lp_SO_abn)
se_lp_SO_abn <- abs(x_abn) * se1
lower_SO_abn <- plogis(lp_SO_abn - 1.96 * se_lp_SO_abn)
upper_SO_abn <- plogis(lp_SO_abn + 1.96 * se_lp_SO_abn)

ref_val <- fit_SO_abn[ref_idx]
df_SO_abn <- tibble(
  x = x_abn, fit = fit_SO_abn, lower = lower_SO_abn, upper = upper_SO_abn,
  metric = "Sonotype Occurrence"
) %>%
  mutate(
    fit_pctchange   = (fit - ref_val)/ref_val*100,
    lower_pctchange = (lower - ref_val)/ref_val*100,
    upper_pctchange = (upper - ref_val)/ref_val*100
  )

# SS ~ Binomial distribution
b0  <- get_b0(c_SS_abn)
b1  <- get_b1(c_SS_abn, "fish_abundance_scaled")
se1 <- get_seb1(c_SS_abn, "fish_abundance_scaled")
bp  <- get_pred(c_SS_abn)

lp_SS_abn <- b0 + b1 * x_abn + p_present * bp
fit_SS_abn <- plogis(lp_SS_abn)
se_lp_SS_abn <- abs(x_abn) * se1
lower_SS_abn <- plogis(lp_SS_abn - 1.96 * se_lp_SS_abn)
upper_SS_abn <- plogis(lp_SS_abn + 1.96 * se_lp_SS_abn)

ref_val <- fit_SS_abn[ref_idx]
df_SS_abn <- tibble(
  x = x_abn, fit = fit_SS_abn, lower = lower_SS_abn, upper = upper_SS_abn,
  metric = "Snapping Shrimp Activity"
) %>%
  mutate(
    fit_pctchange   = (fit - ref_val)/ref_val*100,
    lower_pctchange = (lower - ref_val)/ref_val*100,
    upper_pctchange = (upper - ref_val)/ref_val*100
  )

abn_all <- bind_rows(df_SO_abn, df_PR_abn, df_SS_abn) %>%
  mutate(metric = factor(metric,
                         levels = c("Sonotype Occurrence","Phonic Richness","Snapping Shrimp Activity")),
         focal = "Fish Abundance (scaled)")

# -------------------------------------------------------------------
# Combine fish diversity and abundance, plot the graph
# -------------------------------------------------------------------
combined_all <- bind_rows(
  div_all %>% mutate(panel = "Fish Diversity (scaled)"),
  abn_all %>% mutate(panel = "Fish Abundance (scaled)")
)

facet_labs <- c(
  "Sonotype Occurrence"      = "A) Sonotype Occurrence",
  "Phonic Richness"          = "B) Phonic Richness",
  "Snapping Shrimp Activity" = "C) Snapping Shrimp Activity"
)

# Diversity figure
fish_div <- ggplot(subset(combined_all, panel == "Fish Diversity (scaled)"),
                aes(x = x, y = fit_pctchange)) +
  geom_ribbon(aes(ymin = lower_pctchange, ymax = upper_pctchange), fill = "#1f78b4", alpha = 0.35) +
  geom_line(color = "black", linewidth = 1) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1,
             labeller = as_labeller(facet_labs)) +
  labs(x = "Fish Diversity (scaled)", y = "Percent Change (%)") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12.5, hjust = -0.045, face = "bold"),
    text = element_text(size = 12.5),
    axis.text = element_text(size = 12.5, colour = "black"),
    axis.ticks.length = unit(0.1, "cm"),
    panel.background = element_rect(color = "black", size = 1)
  )
fish_div

# Abundance figure
fish_abn <- ggplot(subset(combined_all, panel == "Fish Abundance (scaled)"),
                aes(x = x, y = fit_pctchange)) +
  geom_ribbon(aes(ymin = lower_pctchange, ymax = upper_pctchange), fill = "#1f78b4", alpha = 0.35) +
  geom_line(color = "black", linewidth = 1) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1,
             labeller = as_labeller(facet_labs)) +
  labs(x = "Fish Abundance (scaled)", y = "Percent Change (%)") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12.5, hjust = -0.045, face = "bold"),
    text = element_text(size = 12.5),
    axis.text = element_text(size = 12.5, colour = "black"),
    axis.ticks.length = unit(0.1, "cm"),
    panel.background = element_rect(color = "black", size = 1)
  )
fish_abn

# Save
ggsave(file.path(output_location, "acoustic_metrics_vs_fish_diversity.png"),
       plot = fish_div, width = 6, height = 9, dpi = 300)

ggsave(file.path(output_location, "acoustic_metrics_vs_fish_abundance.png"),
       plot = fish_abn, width = 6, height = 9, dpi = 300)
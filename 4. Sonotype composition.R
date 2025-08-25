######################################################
# Sonograms composition
# Author: PTSH0
# Project: MSc Biodiversity and Global Change thesis
# Date: 2025-08-01
#
# This script reproduces bar plots for each sonotype 
# composition on different habitat types
######################################################

# ---- Packages ----
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(grid) 

# --- Set paths and output files ---
data_location <- file.path("data", "merged_hydromoth2.xlsx")
output_location   <- file.path("results", "figures")

# --- Data ---
fish <- read_excel(data_path)

# Mutate data
fish <- fish %>%
  mutate(
    Fish_Abundance_Sum_of_MaxN      = suppressWarnings(as.numeric(`Fish Abundance (Sum of MaxN)`)),
    Fish_Diversity_Count_of_Species = suppressWarnings(as.numeric(`Fish Diversity (Count of Species)`)),
    habitat = as.factor(habitat)
  )

# ---- Plot theme ----
plot_theme <- theme_bw() +
  theme(
    strip.background = element_rect(fill = "black", color = "black", size = 2.5),
    strip.text       = element_text(size = 9, colour = "white", face = "bold"),
    legend.position  = "bottom",
    legend.text      = element_text(size = 10),
    legend.title     = element_text(size = 6.5),
    legend.key.size  = unit(1, "cm"),
    legend.spacing   = unit(0.5, "cm"),
    panel.background = element_rect(colour = "black", size = 0.5),
    panel.border     = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(colour = "black", linewidth = 0.7),
    plot.title       = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.x     = element_text(size = 10, face = "bold", hjust = 0.475, margin = margin(t = 5)),
    axis.title.y     = element_text(size = 10, face = "bold"),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.text.y      = element_text(colour = "black", size = 10)
  )

# ---- 1) Diversity summary (mean ± SE by habitat) ----
diversity_summary <- fish %>%
  group_by(habitat) %>%
  summarise(
    mean_div = mean(Fish_Diversity_Count_of_Species, na.rm = TRUE),
    se_div   = sd(Fish_Diversity_Count_of_Species, na.rm = TRUE) /
      sqrt(sum(!is.na(Fish_Diversity_Count_of_Species))),
    .groups = "drop"
  )

print(diversity_summary)

fish_diversity <- ggplot(diversity_summary, aes(x = habitat, y = mean_div)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(ymin = mean_div - se_div, ymax = mean_div + se_div), width = 0.2) +
  labs(y = "Fish Diversity (Mean ± SE)", x = "Habitat", title = "Fish Diversity by Habitat") +
  theme_minimal()

print(fish_diversity)
ggsave(file.path(output_location, "fish_diversity_by_habitat.png"),
       plot = fish_diversity, width = 6, height = 4, dpi = 300)

# ---- 2) Fish Abundance (mean ± SE by habitat) ----
abundance_summary <- fish %>%
  group_by(habitat) %>%
  summarise(
    mean_abun = mean(Fish_Abundance_Sum_of_MaxN, na.rm = TRUE),
    se_abun   = sd(Fish_Abundance_Sum_of_MaxN, na.rm = TRUE) /
      sqrt(sum(!is.na(Fish_Abundance_Sum_of_MaxN))),
    .groups = "drop"
  )

fish_abundance <- ggplot(abundance_summary, aes(x = habitat, y = mean_abun)) +
  geom_col(fill = "darkgreen") +
  geom_errorbar(aes(ymin = mean_abun - se_abun, ymax = mean_abun + se_abun), width = 0.2) +
  labs(y = "Fish Abundance (Mean ± SE)", x = "Habitat", title = "Fish Abundance by Habitat") +
  theme_minimal()

print(fish_abundance)
ggsave(file.path(output_location, "fish_abundance_by_habitat.png"),
       plot = fish_abundance, width = 6, height = 4, dpi = 300)

# ---- 3) Sum count for each sonotype per habitat ----
sonotype_cols <- c("grunt","pop","oink","pulse_purr","pulse_knock","snap")
fish[sonotype_cols] <- lapply(fish[sonotype_cols], function(x) suppressWarnings(as.numeric(x)))

sonotype_summary <- fish %>%
  group_by(habitat) %>%
  summarise(across(all_of(sonotype_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(all_of(sonotype_cols), names_to = "sonotype", values_to = "count")

# Cstom fill palette for habitat types
habitat_cols <- c("#7ac67b", "#a7a4c6", "firebrick", "#1f78b4")

sonograms <- sonotype_summary %>%
  group_split(sonotype) %>%
  setNames(unique(sonotype_summary$sonotype)) %>%
  map(~ ggplot(.x, aes(x = habitat, y = count, fill = habitat)) +
        geom_col() +
        scale_fill_manual(values = habitat_cols) +
        labs(y = NULL, x = NULL) +
        plot_theme +
        theme_classic() +
        theme(
          legend.position = "bottom",
          legend.text  = element_text(face = "bold", size = 18),
          legend.title = element_text(face = "bold", size = 20),
          axis.text.x  = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title   = element_blank(),
          panel.border = element_blank(),
          panel.grid   = element_blank(),
          axis.line.x  = element_line(),
          axis.line.y.right = element_line(),
          axis.line.y.left  = element_blank(),
          axis.text.y  = element_text(size = 25)
        ) +
        scale_y_continuous(position = "right", expand = c(0, 0))
  )
sonograms

# Print each plot
invisible(lapply(sonograms, print))

# Save each sonotype plot
purrr::iwalk(
  plots,
  ~ ggsave(filename = file.path(output_location, paste0(.y, "_sumcount_by_habitat.png")),
    plot = .x, width = 4.0, height = 5.5, dpi = 300))

cat("Save all figure to:", normalizePath(output_location), "\n")

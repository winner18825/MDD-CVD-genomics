# =============================================================================
# Figure | Bidirectional MR forest plot (IVW)
# =============================================================================
# Faceted forest plot of the inverse-variance-weighted (IVW) MR estimates,
# split by direction (MDD -> CVD vs CVD -> MDD) and by MDD
# subtype. Reads MR_IVW_results_summary.txt produced by 08_bidirectional_MR.R.
#
# Public packages only:  install.packages(c("data.table", "ggplot2"))
# =============================================================================

source("config.R")
library(data.table)
library(ggplot2)

IVW_FILE <- file.path(RESULTS_DIR, "mr", "MR_IVW_results_summary.txt")
OUT_DIR  <- file.path(RESULTS_DIR, "figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

df <- fread(IVW_FILE, sep = "\t")
cat("Loaded", nrow(df), "IVW results.\n")

plot_data <- copy(df)

# ---- Trait-label abbreviations ----------------------------------------------
# key = raw trait name (upper case / underscores), value = display abbreviation.
# Longer keys are matched first to avoid partial replacement.
abbr_map <- c(
  "MDD_ALL"        = "MDD",
  "eoMDD"          = "eoMDD",
  "loMDD"          = "loMDD",
  "CAD"            = "CAD",
  "MI"             = "MI",
  "AF"             = "AF",
  "HF"             = "HF",
  "AIS"            = "AIS",
  "Hypertension"   = "HTN",
  "PAD"            = "PAD"
)

prettify_label <- function(x) {
  x <- as.character(x)
  keys <- names(abbr_map)[order(nchar(names(abbr_map)), decreasing = TRUE)]
  for (k in keys)
    x <- gsub(paste0("\\b", k, "\\b"), abbr_map[[k]], x, ignore.case = TRUE)
  # fallback: any remaining underscore name -> initials
  needs_fb <- grepl("_", x) | grepl("[A-Z]{4,}", x)
  if (any(needs_fb)) {
    x[needs_fb] <- vapply(x[needs_fb], function(s) {
      w <- strsplit(s, "[_ ]+")[[1]]; w <- w[nzchar(w)]
      paste0(toupper(substr(w, 1, 1)), collapse = "")
    }, character(1))
  }
  gsub("_", " ", x)
}

plot_data[, exposure_label := prettify_label(exposure)]
plot_data[, outcome_label  := prettify_label(outcome)]
plot_data[, pair_label     := paste0(exposure_label, " \u2192 ", outcome_label)]

# ---- Faceting factors -------------------------------------------------------
plot_data[, direction_type := ifelse(
  exposure %in% MDD_TRAITS, "MDD \u2192 CVD", "CVD \u2192 MDD")]
plot_data[, direction_type := factor(direction_type,
                                     levels = c("CVD \u2192 MDD", "MDD \u2192 CVD"))]

plot_data[, MDD_subtype := fcase(
  grepl("eoMDD", exposure) | grepl("eoMDD", outcome), "eoMDD",
  grepl("loMDD", exposure) | grepl("loMDD", outcome), "loMDD",
  default = "Overall")]
plot_data[, MDD_subtype := factor(MDD_subtype,
                                       levels = c("Overall", "eoMDD", "loMDD"))]

# ---- Significance grading ---------------------------------------------------
bonf_threshold <- 0.05 / nrow(plot_data)
plot_data[, sig_color := fcase(
  pval < bonf_threshold, "Bonferroni",
  pval < 0.05,           "Nominal",
  default              = "NS")]

plot_data[, sig_color := factor(sig_color,
                                levels = c("Bonferroni", "Nominal", "NS"))]

# order rows within each facet by p-value
plot_data <- plot_data[order(direction_type, MDD_subtype, pval)]
plot_data[, plot_y := factor(paste0(pair_label, "  "),
                             levels = rev(unique(paste0(pair_label, "  "))))]

# ---- Plot -------------------------------------------------------------------
sig_colors <- c(Bonferroni = "#C0392B", Nominal = "#E67E22", NS = "#7F8C8D")
or_max <- min(max(plot_data$OR_uci95, na.rm = TRUE) * 1.05, 5)
or_min <- max(min(plot_data$OR_lci95, na.rm = TRUE) * 0.95, 0.2)

p <- ggplot(plot_data, aes(x = OR, y = plot_y)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = OR_lci95, xmax = OR_uci95, color = sig_color),
                 height = 0.20, linewidth = 0.7) +
  geom_point(aes(color = sig_color, fill = sig_color),
             shape = 23, size = 3.2, stroke = 0.5) +
  scale_color_manual(values = sig_colors, name = "Significance", drop = FALSE) +
  scale_fill_manual(values = sig_colors, name = "Significance", drop = FALSE) +
  scale_x_continuous(trans = "log10",
                     breaks = c(0.3, 0.5, 0.7, 1, 1.5, 2, 3),
                     limits = c(or_min, or_max),
                     expand = expansion(mult = c(0.02, 0.02))) +
  facet_grid(direction_type ~ MDD_subtype,
             scales = "free_y", space = "free_y", switch = "y") +
  labs(x = "Odds Ratio (95% CI), log scale", y = NULL,
       title = paste("Bidirectional Mendelian Randomization:",
                     "MDD \u2194 Cardiovascular Diseases (IVW)")) +
  theme_bw(base_size = 12) +
  theme(
    plot.title         = element_text(hjust = 0.5, face = "bold", size = 14,
                                      margin = margin(b = 10)),
    strip.text         = element_text(face = "bold", size = 12.5, color = "white"),
    strip.text.y       = element_text(angle = 0),
    strip.background.x = element_rect(fill = "#2C3E50", color = NA),
    strip.background.y = element_rect(fill = "#34495E", color = NA),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
    panel.spacing      = unit(0.6, "lines"),
    legend.position    = "bottom",
    legend.title       = element_text(face = "bold", size = 12),
    axis.text          = element_text(size = 11, color = "black"),
    axis.title.x       = element_text(face = "bold", size = 12,
                                      margin = margin(t = 8)),
    plot.margin        = margin(15, 15, 10, 15)
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 1.5, size = 3.2)),
         fill  = guide_legend(override.aes = list(size = 3.2)))

# ---- Save -------------------------------------------------------------------
h <- max(8, nrow(plot_data) * 0.30 + 3)
ggsave(file.path(OUT_DIR, "Figure_MR_forest.pdf"), p,
       width = 11, height = h, dpi = 300)
ggsave(file.path(OUT_DIR, "Figure_MR_forest.png"), p,
       width = 11, height = h, dpi = 300, bg = "white")

cat("Saved: Figure_MR_forest.pdf / .png in", OUT_DIR, "\n")
cat(sprintf("Bonferroni threshold: P < %.3e\n", bonf_threshold))
print(plot_data[, .N, by = sig_color])

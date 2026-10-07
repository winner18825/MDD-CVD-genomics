# ============================================================
# Genetic correlation heatmap: LDSC vs HDL
# Final style: two side-by-side heatmaps
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(grid)

source("config.R")
OUT_DIR <- "results/figures"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# 1. Define the original input order of MDD and CVD traits
mdd_order_input <- c("MDD", "eoMDD", "loMDD")
cvd_order_input <- c("AF", "AIS", "CAD", "HF", "HTN", "MI", "PAD")

mdd_order_plot <- c("MDD", "eoMDD", "loMDD")
cvd_order_plot <- c("AF", "AIS", "CAD", "HF", "HTN", "MI", "PAD")

# 2. Raw genetic correlation data (kept in the original input order)
raw_dat <- data.frame(
  MDD = rep(mdd_order_input, each = 7),
  cvd = rep(cvd_order_input, times = 3),
  ldsc_rg = c(-0.1049, 0.1158, 0.2068, 0.2014, 0.2285, 0.1457, 0.4241,
              -0.1166, 0.211, 0.2149, 0.2438, 0.2345, 0.1723, 0.4663,
              -0.0368, 0.028, 0.13, 0.0907, 0.1773, 0.0734, 0.3006),
  ldsc_p  = c(0.0001,0.001,1.2078e-32,9.0811e-16,2.1808e-32,5.2904e-11,6.4339e-10,
              0.0011,7.4314e-06,5.4268e-17,1.2558e-16,1.481e-20,4.4229e-08,5.3745e-08,
              0.3037,0.5998,5.1591e-07,0.0123,1.1143e-09,0.0231,0.0002),
  hdl_rg  = c(-0.1440, 0.1904, 0.2425, 0.3063, 0.2651, 0.2604, 0.4096, 
              -0.1636, 0.1848, 0.2099, 0.2736, 0.2453, 0.2320, 0.4431, 
              -0.0863, 0.0616, 0.1506, 0.1137, 0.2130, 0.1219, 0.2961),
  hdl_p   = c(1.71E-07,0.2540,0.0127,0.2330,3.33E-88,0.1878,7.84E-17,
              1.51E-06,6.97E-06,2.69E-24,9.26E-28,2.58E-30,8.98E-08,2.86E-09,
              0.0321,0.2450,2.14E-11,0.0001,1.74E-18,0.0079,1.57E-05),
  stringsAsFactors = FALSE
)

# 3. Significance labels
bonf <- 0.05 / nrow(raw_dat)

# significance annotation: ** Bonferroni, * nominal
add_sig <- function(p) ifelse(p < bonf, "**", ifelse(p < 0.05, "*", ""))

# 4. Convert data into long format for ggplot visualization
plot_dat <- bind_rows(
    raw_dat %>%
      transmute(MDD, cvd, method = "LDSC", rg = ldsc_rg, p  = ldsc_p),
    raw_dat %>%
      transmute(MDD, cvd, method = "HDL", rg = hdl_rg, p  = hdl_p)) %>%
    mutate( label = paste0(sprintf("%.2f", rg), add_sig(p)),
    
    cvd = factor(cvd, levels = cvd_order_plot),
    MDD = factor(MDD, levels = rev(mdd_order_plot)),
    
    method = factor(method, levels = c("LDSC", "HDL")),
    text_col = ifelse(abs(rg) >= 0.22, "white", "black")
  )

# 5. Define symmetric color scale limits based on the maximum absolute rg value
lim <- ceiling(max(abs(plot_dat$rg)) * 10) / 10

# 6. Plot
p <- ggplot(plot_dat, aes(x = cvd, y = MDD, fill = rg)) +
     geom_tile(color = "white", linewidth = 0.9, width = 0.95, height = 0.95) +
     geom_text(aes(label = label, color = text_col), fontface = "bold",
               size = 2,show.legend = FALSE) +
     facet_grid(. ~ method) +
     scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#C62828", midpoint = 0,
                          limits = c(-0.6, 0.6), breaks = seq(-0.6, 0.6, by = 0.2),
                          labels = c("-0.6","-0.4", "-0.2", "0","0.2", "0.4", "0.6"), name = expression(italic(r)[g])) +
     scale_color_identity() +
     labs(x = NULL, y = NULL, caption = "*P < 0.05    **P < 0.05/21 (Bonferroni)") +
     coord_fixed() +
     theme_minimal(base_size = 13) +
     theme(
       panel.grid = element_blank(),
       strip.text = element_text(size = 16, face = "bold"),
       strip.background = element_blank(),
       axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8, face = "bold"),
       axis.text.y = element_text(size = 8, face = "bold"),
       panel.spacing.x = unit(1.0, "lines"),
       legend.position = "right",
       legend.title = element_text(size = 13),
       legend.text  = element_text(size = 11),
       plot.caption = element_text(hjust = 0.5, size = 11, 
                                   color = "grey35",margin = margin(t = 12)),
       plot.margin = margin(10, 15, 10, 15)
  )
# 7. Save
ggsave(filename = file.path(OUT_DIR, "Figure_rg_heatmap.pdf"), plot = p,
       width = 8.8, height = 4.6, dpi = 300)

ggsave(filename = file.path(OUT_DIR, "Figure_rg_heatmap.png"), plot = p,
       width = 8.8, height = 4.6, dpi = 300, bg = "white")

cat("Saved:\n",
    file.path(OUT_DIR, "Figure_rg_heatmap_two_panel.pdf"), "\n",
    file.path(OUT_DIR, "Figure_rg_heatmap_two_panel.png"), "\n")

# =============================================================================
# Figure | SMR pleiotropic gene network
# =============================================================================
# Network diagram linking MDD subtypes (left) to pleiotropic genes
# (centre) to cardiovascular diseases (right). Edge thickness encodes the
# number of tissues in which the SMR association was significant; edge colour
# encodes the direction of the SMR effect.
#
# This script is self-contained: the edges below are the published SMR
# results. Only public packages are required.
#   install.packages(c("ggplot2", "dplyr"))
# =============================================================================

library(ggplot2)
library(dplyr)

OUT_DIR <- "results/figures"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Edge data (gene -- trait associations) ------------------------------
edges <- data.frame(
  gene = c(
    "IP6K2","IP6K2",
    "MST1","MST1","MST1","MST1",
    "RNF123","RNF123",
    "TMEM106B","TMEM106B","TMEM106B","TMEM106B",
    "ZNF652","ZNF652"
  ),
  trait = c(
    "HF","MDD",
    "CAD","HF","MDD","loMDD",
    "CAD","loMDD",
    "CAD","HF","MDD","loMDD",
    "CAD","MDD"
  ),
  n_tissues = c(
    3,3,
    2,2,1,6,
    23,1,
    3,7,9,6,
    3,2
  ),
  direction = c(
    "-","+",
    "+","+","+","-",
    "+","+",
    "+","+","+","+",
    "-","-"
  ),
  stringsAsFactors = FALSE
)
# ---- 2. Node coordinates (symmetric about y = 0) ----------------------------
MDD_nodes <- data.frame(
  name = c("MDD", "loMDD"),
  x = rep(0, 2),
  y = seq(2, -2, length.out = 2),
  type = rep("MDD", 2),
  stringsAsFactors = FALSE
)

gene_list <- c("IP6K2", "MST1", "RNF123", "TMEM106B", "ZNF652")
gene_nodes <- data.frame(
  name = gene_list, x = 5,
  y = seq(3, -3, length.out = length(gene_list)),
  type = "gene", stringsAsFactors = FALSE)

cvd_list <- c("CAD",  "HF")
cvd_nodes <- data.frame(
  name = cvd_list, x = 10,
  y = seq(2.5, -2.5, length.out = length(cvd_list)),
  type = "cvd", stringsAsFactors = FALSE)

all_nodes <- bind_rows(MDD_nodes, gene_nodes, cvd_nodes)

# ---- 3. Edge coordinates, with label positions spread to avoid overlap ------
edge_coords <- edges %>%
  left_join(gene_nodes %>% select(name, x, y), by = c("gene" = "name")) %>%
  rename(x_gene = x, y_gene = y) %>%
  left_join(all_nodes %>% select(name, x, y), by = c("trait" = "name")) %>%
  rename(x_trait = x, y_trait = y)

set.seed(42)   # fixed seed -> reproducible label placement
edge_coords <- edge_coords %>%
  group_by(trait) %>%
  mutate(base_ratio = if (n() == 1) 0.5 else seq(0.35, 0.65, length.out = n())) %>%
  ungroup() %>%
  mutate(ratio = base_ratio + runif(n(), min = -0.03, max = 0.03))

edge_coords$lbl_x <- with(edge_coords, x_gene * (1 - ratio) + x_trait * ratio)
edge_coords$lbl_y <- with(edge_coords, y_gene * (1 - ratio) + y_trait * ratio)

# ---- 4. Colours -------------------------------------------------------------
dir_colors <- c("+" = "#D73027", "-" = "#4575B4")

# ---- 5. Plot ----------------------------------------------------------------
p <- ggplot() +
  # connecting lines
  geom_segment(data = edge_coords,
               aes(x = x_gene, y = y_gene, xend = x_trait, yend = y_trait,
                   colour = direction, linewidth = n_tissues),
               alpha = 0.40, lineend = "round") +
  scale_colour_manual(values = dir_colors, name = "SMR effect direction",
                      labels = c("+" = "Positive (risk \u2191)",
                                 "-" = "Negative (risk \u2193)")) +
  scale_linewidth_continuous(range = c(0.5, 5.5),
                             name = "No. significant tissues",
                             breaks = c(1, 3, 8, 17)) +
  # tissue-count badge on each edge
  geom_point(data = edge_coords, aes(x = lbl_x, y = lbl_y),
             shape = 21, size = 9.5, fill = "white",
             colour = "grey50", stroke = 0.4) +
  geom_text(data = edge_coords,
            aes(x = lbl_x, y = lbl_y, label = n_tissues, colour = direction),
            size = 4.6, fontface = "bold", show.legend = FALSE) +
  # MDD nodes (circles)
  geom_point(data = MDD_nodes, aes(x = x, y = y),
             shape = 21, size = 24, fill = "#E64B35",
             colour = "white", stroke = 1.5) +
  geom_text(data = MDD_nodes, aes(x = x, y = y, label = name),
            colour = "white", size = 4.8, fontface = "bold") +
  # gene nodes (rectangles)
  geom_rect(data = gene_nodes,
            aes(xmin = x - 0.95, xmax = x + 0.95,
                ymin = y - 0.38, ymax = y + 0.38),
            fill = "#3C5488", colour = "white", linewidth = 0.8) +
  geom_text(data = gene_nodes, aes(x = x, y = y, label = name),
            colour = "white", size = 5.0, fontface = "bold.italic") +
  # cardiovascular disease nodes (circles)
  geom_point(data = cvd_nodes, aes(x = x, y = y),
             shape = 21, size = 24, fill = "#00A087",
             colour = "white", stroke = 1.5) +
  geom_text(data = cvd_nodes, aes(x = x, y = y, label = name),
            colour = "white", size = 4.8, fontface = "bold") +
  # column headers
  annotate("text", x = 0, y = max(all_nodes$y) + 1.2,
           label = "MDD\nsubtypes", size = 5.5, fontface = "bold",
           colour = "#E64B35", lineheight = 0.9) +
  annotate("text", x = 5, y = max(all_nodes$y) + 1.2,
           label = "Pleiotropic\ngenes", size = 5.5, fontface = "bold",
           colour = "#3C5488", lineheight = 0.9) +
  annotate("text", x = 10, y = max(all_nodes$y) + 1.2,
           label = "Cardiovascular\ndiseases", size = 5.5, fontface = "bold",
           colour = "#00A087", lineheight = 0.9) +
  coord_cartesian(xlim = c(-2, 12),
                  ylim = c(min(all_nodes$y) - 1.0, max(all_nodes$y) + 2.0)) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title       = element_text(size = 15, face = "bold"),
    legend.text        = element_text(size = 13),
    legend.key.size    = unit(0.9, "cm"),
    legend.key.width   = unit(1.2, "cm"),
    legend.spacing.x   = unit(0.4, "cm"),
    legend.box.spacing = unit(0.5, "cm"),
    plot.margin = margin(15, 15, 15, 15)
  ) +
  guides(
    colour    = guide_legend(order = 1,
                             override.aes = list(linewidth = 4, alpha = 1)),
    linewidth = guide_legend(order = 2)
  )

# ---- 6. Save ----------------------------------------------------------------
ggsave(file.path(OUT_DIR, "Figure_SMR_network.pdf"), p,
       width = 13, height = 9.5, dpi = 300)
ggsave(file.path(OUT_DIR, "Figure_SMR_network.png"), p,
       width = 13, height = 9.5, dpi = 300, bg = "white")

cat("Saved: Figure_SMR_network.pdf / .png in", OUT_DIR, "\n")
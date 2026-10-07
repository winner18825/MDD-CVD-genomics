# ============================================================
# Colocalization circular tree plot: visualize shared genetic signals (PP.H4)
# across traits using a circular dendrogram layout

# Colocalization circular tree plot
# Load packages for data import, graph construction, and visualization
library(readxl)
library(dplyr)
library(igraph)
library(ggraph)
library(ggplot2)

COLOC_FILE <- file.path(RESULTS_DIR, "coloc", "coloc_summary.txt")
OUT_DIR  <- file.path(RESULTS_DIR, "figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

df <- fread(COLOC_FILE, sep = "\t")
# Define threshold for strong colocalization evidence
pph4_cutoff <- 0.7

# Colors represent the three primary trait branches
branch_cols <- c("MDD" = "#F05A4A", "loMDD" = "#45B8D8", "eoMDD" = "#27A98B")
trait1_order <- c("MDD", "loMDD", "eoMDD")
trait2_order <- c("AF", "AIS", "CAD", "HF", "HTN", "MI", "PAD")

# Standardize variables, remove incomplete records, and classify colocalization status
dat <- dat %>%
  mutate(
    trait1 = trimws(as.character(trait1)),trait2 = trimws(as.character(trait2)),
    lead_SNP = trimws(as.character(lead_SNP)),PP.H4 = as.numeric(PP.H4)) %>%
  filter(!is.na(trait1),!is.na(trait2),
         !is.na(lead_SNP), !is.na(PP.H4))

trait1_order2 <- c(intersect(trait1_order, unique(dat$trait1)),
                   setdiff(unique(dat$trait1), trait1_order))
trait2_order2 <- c(intersect(trait2_order, unique(dat$trait2)),
                   setdiff(unique(dat$trait2), trait2_order))
dat <- dat %>%
  mutate(trait1 = factor(trait1, levels = trait1_order2),
          trait2 = factor(trait2, levels = trait2_order2)) %>%
  arrange(trait1,trait2,desc(PP.H4)) %>%
  mutate(leaf_id = sprintf("SNP_%03d", row_number()),
         coloc = PP.H4 > pph4_cutoff,
         coloc_class = ifelse(coloc, "Colocalized (PP.H4 > 0.7)", "Not colocalized"))
cat("Total results:", nrow(dat), "\n")
cat("Colocalized:", sum(dat$coloc), "\n")
cat( "Not colocalized:", sum(!dat$coloc), "\n")

node_root <- data.frame(name = "ROOT", node_class = "root",
                        label = "MDD", branch = NA_character_, stringsAsFactors = FALSE)
node_trait1 <- dat %>%
  distinct(trait1) %>%
  mutate(trait1_chr = as.character(trait1),name = paste0("T1_", trait1_chr),
         node_class = "trait1",label = case_when(trait1_chr == "MDD" ~ "MDD",trait1_chr == "loMDD" ~ "loMDD",
         trait1_chr == "eoMDD" ~ "eoMDD",TRUE ~ trait1_chr), branch = trait1_chr) %>%
  select(name, node_class, label, branch)

node_trait2 <- dat %>%distinct(trait1,trait2) %>%
  mutate(trait1_chr = as.character(trait1),trait2_chr = as.character(trait2),
         name = paste0("T2_",trait1_chr,"_",trait2_chr),node_class = "trait2",
         label = ifelse(trait2_chr == "Hypertension", "HTN", trait2_chr), branch = trait1_chr) %>%
  select(name, node_class, label, branch)

node_leaf <- dat %>%
  transmute(name = leaf_id, node_class = "snp", label = as.character(lead_SNP),
            branch = as.character(trait1))

nodes <- bind_rows(node_root, node_trait1,
                   node_trait2, node_leaf)

# Construct edges connecting root -> trait1 -> trait2 -> SNP nodes
edge1 <- dat %>% distinct(trait1) %>%
  transmute(from = "ROOT", to = paste0("T1_", as.character(trait1)),
            branch = as.character(trait1))

edge2 <- dat %>% distinct(trait1, trait2) %>%
  transmute(from = paste0("T1_", as.character(trait1)),
            to = paste0("T2_", as.character(trait1), "_", as.character(trait2)),
            branch = as.character(trait1))

edge3 <- dat %>%
  transmute(from = paste0("T2_", as.character(trait1), "_", as.character(trait2)),
            to = leaf_id, branch = as.character(trait1))

edges <- bind_rows(edge1, edge2, edge3)
# Build graph object for circular tree visualization
g <- graph_from_data_frame(d = edges, directed = TRUE, vertices = nodes)

# Generate circular dendrogram coordinates
lay <- create_layout(g,layout = "dendrogram",circular = TRUE)
leaf_meta <- dat %>% transmute(name = leaf_id, SNP = as.character(lead_SNP),
                               PP.H4 = PP.H4, coloc = coloc, coloc_class = coloc_class)
leaf_df   <- as.data.frame(lay) %>% filter(leaf == TRUE) %>%
             select(name, x, y) %>% left_join(leaf_meta,by = "name")
# Highlight colocalized SNPs with larger points
leaf_df   <- leaf_df %>% mutate(point_size = ifelse(coloc, 3.2, 1.0))
label_df  <- leaf_df %>% 
  filter(coloc == TRUE) %>%
  mutate(label_x = x * 1.10, label_y = y * 1.10,
         label_angle = ggraph::node_angle(x, y))

# Plot branches, trait labels, SNP points, and colocalization annotations
p <- ggraph(lay) +
geom_edge_diagonal(
  aes(edge_colour = branch), linewidth = 0.32,
  alpha = 0.32, show.legend = TRUE) +
  scale_edge_colour_manual(
    name = NULL, values = c("MDD"   = "#F05A4A", "loMDD" = "#45B8D8", "eoMDD" = "#27A98B"),
    breaks = c("MDD", "loMDD", "eoMDD"),
    labels = c("MDD (overall)", "loMDD", "eoMDD")) +
geom_node_text(
  aes(filter = node_class == "trait1", label = label, colour = branch),
  fontface = "bold", size = 3.5) +
  scale_colour_manual(
    values = c("MDD"   = "#F05A4A", "loMDD" = "#45B8D8", "eoMDD" = "#27A98B"),
    guide = "none") +
geom_node_text(
  aes(filter = node_class == "trait2", label = label), colour = "#333333",
  fontface = "bold", size = 2.4) +
geom_point(data = leaf_df,
  aes(x = x, y = y, fill = coloc_class, size = point_size),
  inherit.aes = FALSE, shape = 21, colour = "white", stroke = 0.2) +
scale_fill_manual(
  name = NULL,values = c("Colocalized (PP.H4 > 0.7)" = "#D89B1D","Not colocalized" = "#B9DDF4"),
  breaks = c("Colocalized (PP.H4 > 0.7)", "Not colocalized"),drop = FALSE,
  guide = guide_legend(order = 1, override.aes = list(shape = 21, fill = c("#D89B1D", "#B9DDF4"),
                                                      colour = "black", size = c(5, 3), stroke = 0.3))) +

scale_size_identity(guide = "none") +
geom_text(data = label_df,
  aes(x = label_x, y = label_y, label = SNP, angle = label_angle),
  inherit.aes = FALSE, hjust = "outward", colour = "#222222",
  fontface = "bold", size = 2.25, check_overlap = FALSE) +
coord_fixed(clip = "off") +

labs(caption = paste0("Total colocalization pairs: ",nrow(dat),
                      "   |   Colocalized: ", sum(dat$coloc),
                      "   |   Not colocalized: ",sum(!dat$coloc))) +
theme_void() +
  theme(
    legend.position = "bottom", legend.box = "vertical",
    legend.direction = "vertical", legend.justification = "left",
    legend.text = element_text(size = 9),
    plot.caption = element_text(size = 7, colour = "#777777", hjust = 0),
    plot.margin = margin(t = 80, r = 130, b = 60, l = 130)) +
guides(
  fill = guide_legend(order = 1,
                      override.aes = list(shape = 21,size = c(4, 2),colour = NA)),
  edge_colour = guide_legend(order = 2,
                            override.aes = list(linewidth = 1.5, alpha = 1)))

# Export high-resolution figures
ggsave(file.path(OUT_DIR, "coloc_circular_tree.png"), plot = p,
       width = 12, height = 12, units = "in", dpi = 600, bg = "white")
ggsave(file.path(OUT_DIR, "coloc_circular_tree.pdf"), plot = p,
       width = 12, height = 12, units = "in", device = cairo_pdf, bg = "white")
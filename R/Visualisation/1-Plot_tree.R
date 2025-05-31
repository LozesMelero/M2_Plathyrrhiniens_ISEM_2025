################################################################################
# Name: 1-Plot_tree.R
# Authors: Lucas Buffan & Frédéric Lozes-Méléro
# Contact: lucas.l.buffan@gmail.com, frederic.lozes-melero@orange.fr
# Aim: Script for plotting tree
################################################################################

library(ape)
library(treeio)
library(ggtree)
library(tidyverse)

# Load consensus tree
platy_tree <- read.mega("Results/Phylogeny/MrBayes/TED/10_rogueless/Platy_TED_FBD_datation_3clock_10rogueless_40M.con.tre")
# Fix nomenclature issue with 'Aotus' dindensis (not an Aotus per se)
platy_tree@phylo$tip.label[which(platy_tree@phylo$tip.label == "Aotus_dindensis")] <- "'Aotus'_dindensis"

# Some plotting parameters
tree_height <- 59.60 # from FigTree
gsc_ymin <- -10 # Height of the geoscale
n_tips <- length(platy_tree@phylo$tip.label)

# Plot (flat)
flat_tree <- platy_tree %>% ggtree() +
  geom_tiplab(size = 1.5) +
  # Add 95% HPD on node age
  geom_range("age_0.95HPD", size = 1, alpha = 0.3, color = "purple") +
  # Grey ribbons delimiting geological subdivisions
  annotate(geom = "rect", xmin = tree_height-56, xmax = tree_height-47.8, ymin = 0, ymax = n_tips, fill = "grey", alpha = 0.2) +
  annotate(geom = "rect", xmin = tree_height-37.71, xmax = tree_height-33.9, ymin = 0, ymax = n_tips, fill = "grey", alpha = 0.2) +
  annotate(geom = "rect", xmin = tree_height-27.82, xmax = tree_height-23.04, ymin = 0, ymax = n_tips, fill = "grey", alpha = 0.2) +
  annotate(geom = "rect", xmin = tree_height-15.97, xmax = tree_height-11.63, ymin = 0, ymax = n_tips, fill = "grey", alpha = 0.2) +
  annotate(geom = "rect", xmin = tree_height-5.33, xmax = tree_height-2.58, ymin = 0, ymax = n_tips, fill = "grey", alpha = 0.2) +
  # MECO + EOT indications
  annotate(geom = "rect", xmin = tree_height-41.5, xmax = tree_height-41, ymin = 0, ymax = n_tips, fill = "red", alpha = 0.2) + # MECO
  annotate(geom = "text", x = tree_height-42, y = n_tips-20, label = "MECO", size = 5, colour = "red", angle = "90") +
  geom_segment(x = tree_height-33.9, xend = tree_height-33.9, y = 0, yend = 240, linetype="dashed", color = "red", linewidth = 0.2) + # EOT
  annotate(geom = "text", x = tree_height-34.5, y = n_tips-20, label = "EOT", size = 5, colour = "red", angle = "90") +
  # Manual GTS
  annotate(geom = "rect", xmin = 0, xmax = tree_height-56, ymin = gsc_ymin, ymax = 0, fill = "#FDA75F", colour = "black") +
  annotate(geom = "text", label = "Palaeocene", x = 1.8, y = gsc_ymin/2) +
  annotate(geom = "rect", xmin = tree_height-56, xmax = tree_height-47.8, ymin = gsc_ymin, ymax = 0, fill = "#FCA773", colour = "black") +
  annotate(geom = "text", label = "Early Eocene", x = 7.7, y = gsc_ymin/2) +
  annotate(geom = "rect", xmin = tree_height-47.8, xmax = tree_height-37.71, ymin = gsc_ymin, ymax = 0, fill = "#FDC091", colour = "black") +
  annotate(geom = "text", label = "Mid Eocene", x = 17, y = gsc_ymin/2) +
  annotate(geom = "rect", xmin = tree_height-37.71, xmax = tree_height-33.9, ymin = gsc_ymin, ymax = 0, fill = "#FDCDA1", colour = "black") +
  annotate(geom = "text", label = "Late Eocene", x = 23.8, y = gsc_ymin/2) +
  annotate(geom = "rect", xmin = tree_height-33.9, xmax = tree_height-27.8, ymin = gsc_ymin, ymax = 0, fill = "#FED99A", colour = "black") +
  annotate(geom = "text", label = "Early Oligocene", x = 28.5, y = gsc_ymin/2) +
  annotate(geom = "rect", xmin = tree_height-27.8, xmax = tree_height-23.04, ymin = gsc_ymin, ymax = 0, fill = "#FEE6AA", colour = "black") +
  annotate(geom = "text", label = "Late Oligocene", x = 34.2, y = gsc_ymin/2) +
  annotate(geom = "rect", xmin = tree_height-23.04, xmax = tree_height-15.97, ymin = gsc_ymin, ymax = 0, fill = "#FFFF41", colour = "black") +
  annotate(geom = "text", label = "Early Miocene", x = 40, y = gsc_ymin/2) +
  annotate(geom = "rect", xmin = tree_height-15.97, xmax = tree_height-11.63, ymin = gsc_ymin, ymax = 0, fill = "#FFFF59", colour = "black") +
  annotate(geom = "text", label = "Mid Miocene", x = 46, y = gsc_ymin/2) +
  annotate(geom = "rect", xmin = tree_height-11.63, xmax = tree_height-5.33, ymin = gsc_ymin, ymax = 0, fill = "#FFFF73", colour = "black") +
  annotate(geom = "text", label = "Late Miocene", x = 51, y = gsc_ymin/2) +
  annotate(geom = "rect", xmin = tree_height-5.33, xmax = tree_height-2.58, ymin = gsc_ymin, ymax = 0, fill = "#FFFFBF", colour = "black") +
  annotate(geom = "text", label = "Pliocene", x = 55.7, y = gsc_ymin/2) +
  annotate(geom = "rect", xmin = tree_height-2.58, xmax = tree_height, ymin = gsc_ymin, ymax = 0, fill = "#FEF2E0", colour = "black") +
  annotate(geom = "text", label = "Q", x = 58, y = gsc_ymin/2) +
  # Time axis
  scale_x_continuous(breaks = c(seq(0, 50, 10), tree_height),
                     labels = c(tree_height, seq(50, 0, -10))) +
  labs(x = "Time (Ma)") +
  theme(axis.ticks.x = element_line(),
        axis.text.x = element_text())
  
ggsave("./Figures/Trees/Platy_con_tree_10_rogueless_flat.pdf", plot = flat_tree, height = 300, width = 300, units = "mm")


# Draft plot (rounded)
platy_tree %>% ggtree(layout = "fan", open.angle = 40) + # leave an open angle for GTS
  geom_tiplab(size = 1.5) +
  geom_range("age_0.95HPD", size = 1, alpha = 0.3, color = "purple") + # add 95% HPD on node age
  # Grey ribbons
  annotate(geom = "rect", xmin = tree_height-56, xmax = tree_height-47.8, ymin = 0, ymax = n_tips, fill = "grey", alpha = 0.2) +
  annotate(geom = "rect", xmin = tree_height-37.71, xmax = tree_height-33.9, ymin = 0, ymax = n_tips, fill = "grey", alpha = 0.2) +
  annotate(geom = "rect", xmin = tree_height-27.82, xmax = tree_height-23.04, ymin = 0, ymax = n_tips, fill = "grey", alpha = 0.2) +
  annotate(geom = "rect", xmin = tree_height-15.97, xmax = tree_height-11.63, ymin = 0, ymax = n_tips, fill = "grey", alpha = 0.2) +
  annotate(geom = "rect", xmin = tree_height-5.33, xmax = tree_height-2.58, ymin = 0, ymax = n_tips, fill = "grey", alpha = 0.2) +
  annotate(geom = "rect", xmin = tree_height-41.5, xmax = tree_height-41, ymin = 0, ymax = n_tips, fill = "red", alpha = 0.2) + # MECO
  geom_segment(x = tree_height-33.9, xend = tree_height-33.9, y = 0, yend = 240, linetype="dashed", color = "red", linewidth = 0.05) + # EOT
  


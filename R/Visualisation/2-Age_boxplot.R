################################################################################
# Name: 2-Age_boxplot.R
# Authors: Lucas Buffan & Frédéric Lozes-Méléro
# Contact: lucas.l.buffan@gmail.com, frederic.lozes-melero@orange.fr
# Aim: Script for representing comparative divergence ages for key nodes
################################################################################

library(tidyverse)
library(readxl)

age_tbl <- read_xlsx("./Results/Phylogeny/comparative_age_table.xlsx")
# Assign studies a letter for clarity
studies <- colnames(age_tbl)[-which(colnames(age_tbl) == "Node")]
studies[length(studies)] <- "Marivaux et al. (2023)"
colnames(age_tbl) <- c("Node", LETTERS[1:7])

## Set up plotting dataset -----------------------------------------------------
  # Function to extract data from `age_tbl`

extract_dat <- function(age, what){ 
  # ---
  # `age` in the "median [lower_95hpd − upper_95hpd]" format
  # `what` must be in c("median", "lower", "upper")
  # ---
  if(is.na(age)){
    return(NA)
  }
  else{
    spl <- strsplit(age, split = " ")[[1]]
    if(what == "median"){
      return(as.numeric(spl[1]))
    }
    else if(what == "lower"){
      low <- spl[2]
      low <- gsub("\\[", "", low) # remove initial opening bracket
      return(as.numeric(low))
    }
    else if(what == "upper"){
      up <- spl[4]
      up <- gsub("\\]", "", up)
      return(as.numeric(up))
    }
  }
}
  # Median ages
median_table <- as.data.frame(apply(X = age_tbl[,2:ncol(age_tbl)], FUN = extract_dat, MARGIN = c(1,2), what = "median"))
median_table$node <- age_tbl$Node
med_table <- median_table  %>% 
  pivot_longer(
    cols = -node,
    names_to = "study",
    values_to = "median_age")
  # Min HPD
lower_table <- apply(X = age_tbl[,2:ncol(age_tbl)], FUN = extract_dat, MARGIN = c(1,2), what = "lower")
low_table <- lower_table %>% 
  as_data_frame() %>% 
  pivot_longer(
    cols = everything(),
    names_to = "study",
    values_to = "min_HPD")
  # Max HPD
upper_table <- apply(X = age_tbl[,2:ncol(age_tbl)], FUN = extract_dat, MARGIN = c(1,2), what = "upper")
up_table <- upper_table %>% 
  as_data_frame() %>% 
  pivot_longer(
    cols = everything(),
    names_to = "study",
    values_to = "max_HPD")
  # Proper plotting df
plot_df <- med_table %>% 
  mutate(min_HPD = low_table$min_HPD, max_HPD = up_table$max_HPD)

## Plot ------------------------------------------------------------------------

div_plot <- plot_df %>% 
  ggplot(aes(x = study, y = median_age)) +
  geom_point(size = 1) +
  scale_x_discrete(labels = studies) +
  labs(x = NULL, y = "Age (Ma)") +
  geom_errorbar(aes(ymin = min_HPD, ymax = max_HPD), size = .3) +
  facet_wrap(.~node, nrow = 3, ncol = 4, scales = "free_y") +
  theme(axis.text.x = element_text(size = 5, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 5),
        axis.ticks = element_line(size = .3),
        strip.text = element_text(size = 8))

ggsave("./Figures/Trees/divergence_ages.pdf", plot = div_plot, height = 150, width = 170, units = "mm")

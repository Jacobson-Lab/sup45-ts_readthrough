# ---------------------------------------
# Figure 1E-F
# ---------------------------------------
# Table key:
## Sample = sample
## distance = distance of footprint P-site to either start or stop codon
## reg = region the distance is relative to
## count = raw footprint count at a particular distance
## norm = count converted to RPKM and then normalized to mRNA abundance (RPKM from RNA-Seq RSEM) of a gene in a sample

# 1E: ADE2; 5'-UTR length = 71, CDS length = 1713, 3'-UTR length = 57, PTC = 190-1 (minus 1 because "distance" is nt position - AUG position, which is 1)
# 1F: CAN1: 5'-UTR length = NA, CDS length = 1770, 3'-UTR length = NA, PTC = 139-1 (minus 1 because "distance" is nt position - AUG position, which is 1)
#           (CAN1 lacks UTR annotations. Use 50 for UTR lengths.)

library(ggplot2)
library(dplyr)

df <- read.table("../Data/Data_Figure1E.txt", header = TRUE)  # read in 1E or 1F
df$Sample <- recode_factor(df$Sample, sup45_wt_25C = "SUP45, 25 °C", sup45_ts_25C = "sup45-ts, 25 °C", sup45_wt_37C = "SUP45, 37 °C", sup45_ts_37C = "sup45-ts, 37 °C")

# Split dataframe into smaller dataframes based on sample
xgene <- split(x = df, f = df$Sample)

# Set length + PTC variables (for plotting purpose)
gene_name <- "ADE2"
l5 <- 71      # 5'-UTR length
lc <- 1713    # CDS length
l3 <- 57      # 3'-UTR length
ptc <- 190-1  # position of the first nt of PTC

# Main + inset (zooom in after PTC) plots
scalefactor <- max(df$norm) # for setting y-axis and calculate box coordinates
pcol <- c("#ed1c24", "#f7941d", "#00a651", "#8dc63f") # set colors
pmain <- list()  # for storing main plots
pinset <- list() # for storing inset plots
for (i in 1:4) {
  pmain[[i]] <- ggplot() +
    geom_line(data = xgene[[i]], aes(x = distance + 1, y = norm), size = 0.5, color = pcol[i]) + # Line plot for normalized count data. Add 1 to distance to convert to nt position
    geom_rect(aes(xmin = ptc, xmax = lc, ymin = -0.05*scalefactor, ymax = 0.05*scalefactor), # Draw box around region after PTC for zoom-in
              fill = NA, color = "black", size = 0.2) +
    geom_segment(aes(x = ptc, y = 0.05*scalefactor, xend = ptc+30, yend = 0.25*scalefactor), color = "black", size = 0.2) + # left line segment connecting box to inset
    geom_segment(aes(x = lc, y = 0.05*scalefactor, xend = lc-30, yend = 0.25*scalefactor), color = "black", size = 0.2) + # right line segment connecting box to inset
    scale_x_continuous(limits = c(-l5, lc+l3), expand = expansion(mult = 0.01), labels = scales::comma) + 
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)), labels = scales::comma) +
    coord_cartesian(ylim = c(0, scalefactor), clip = "off") +
    xlab("") + ylab("") +
    theme_bw(base_size = 8) + 
    theme(panel.grid = element_blank(), 
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          plot.margin = unit(c(0.1, 0.1, -0.1, 0.1), "cm"))
  
  pinset[[i]] <- ggplot() +
    geom_line(data = xgene[[i]], aes(x = distance + 1, y = norm), size = 0.5, color = pcol[i]) +
    coord_cartesian(xlim = c(ptc, lc), ylim = c(0, 0.05*scalefactor)) +
    scale_x_continuous(expand = expansion(mult = 0)) +
    scale_y_continuous(expand = expansion(mult = c(0.01))) +
    annotate("text", x = (lc-ptc)/2, y = 0.05*scalefactor, label = names(xgene)[i], fontface = "italic", 
             hjust = 0.35, vjust = 1.2, color = pcol[i], size = 8 / .pt) +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_blank(), legend.position = "none",
          panel.grid = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
          plot.margin = grid::unit(c(0, 0, -1, -1), "mm"), panel.spacing = unit(0, "cm"))
}

# Combine main and inset plots
pcomb <- list()
for (i in 1:4) {
  pcomb[[i]] <- pmain[[i]] + 
    annotation_custom(grob = ggplotGrob(pinset[[i]]), 
                      ymin = 0.25*scalefactor, ymax = 0.95*scalefactor, xmin = ptc+30, xmax = lc-30)
}

# Create boxes representing ORF and UTRs with AUG, PTC, NTC marks
gbox <- ggplot() +
  geom_rect(aes(xmin = -l5, xmax = lc+l3, ymin = 2, ymax = 4), fill = "#2e3192", color = NA) + # UTRs
  geom_rect(aes(xmin = 1, xmax = lc, ymin = 1, ymax = 5), fill = "#2e3192", color = NA) +  # CDS
  geom_rect(aes(xmin = 1-3, xmax = 3+3, ymin = 1, ymax = 5), fill = "#7cc576", color = NA) + # AUG (make box bigger than just 3 nt for visualization)
  geom_rect(aes(xmin = ptc-1-3, xmax = ptc-1+3+3, ymin = 1, ymax = 5), fill = "#fff568", color = NA) + # PTC (make box bigger than just 3 nt for visualization)
  geom_rect(aes(xmin = lc-3, xmax = lc+3+3, ymin = 1, ymax = 5), fill = "#f26c4f", color = NA) + # NTC (make box bigger than just 3 nt for visualization)
  scale_x_continuous(limits = c(-l5, lc+l3), expand = expansion(mult = 0.01)) +
  annotate("text", x = lc/2, y = 3, label = gene_name, fontface = "italic", 
           hjust = 0.5, vjust = 0.5, color = "white", size = 8 / .pt) +
  theme_bw() + theme(axis.text.x = element_text(size = 6), axis.ticks.x = element_line(), 
                     axis.text.y = element_text(color = NA), axis.ticks.y = element_line(color = NA), axis.title = element_blank(),
                     panel.grid = element_blank(), panel.border = element_blank(),
                     plot.margin = unit(c(-0.1, 0.1, 0.1, 0.1), "cm"))

# Arrange plots
f1ef <- cowplot::plot_grid(pcomb[[1]], pcomb[[2]], pcomb[[3]], pcomb[[4]], gbox, 
                           ncol = 1, rel_heights = c(rep(0.92/4, times = 4), 0.08), align = "v", axis = "lr")

# Save plot
ggsave(filename = "Figure1E_ade2.pdf", plot = f1ef, path = "../Plots", width = 3.5, height = 3, units = "in", dpi = 500)
ggsave(filename = "Figure1F_can1.pdf", plot = f1ef, path = "../Plots", width = 3.5, height = 3, units = "in", dpi = 500)


# ---------------------------------------
# Figure 5A
# ---------------------------------------
# Table key:
## Sample = sample
## Var = variable (in this case, nucleotide identities)
## samp_median = median readthrough efficiency of all genes in the sample
## median = median readthrough efficiency of genes containing a particular variable
## p_val = p-value from Wilcoxon's rank sum test comparing median and samp_median
## median_diff = median - samp_median
## p_hw = a value used for plotting; corresponds to the size of the heatmap tile; related to p_sig.
##        If p-value < 0.05 (significant), tile's height and width = 0.9. Otherwise, tile's height and width = 0.6.
## p_sig = string used to label significance of p-value in plot; related to p_hw.
##         The p-value < 0.05 (significant) is labeled as such, and p-value >= 0.05 is labeled "ns"
## position = feature of interest from the analysis (in this case, nt positions relative to the stop codon)

library(ggplot2)
library(dplyr)

df <- read.table("../Data/Data_Figure5A.txt", header = TRUE)
df$Sample <- factor(df$Sample, levels = unique(df$Sample)[c(1, 2, 5, 7, 3, 4, 6, 8)])
df$position <- recode_factor(df$position, random_factor = "Random", nt_m03 = "-3", nt_m02 = "-2", nt_m01 = "-1", stop_codon = "Stop", 
                             nt_p04 = "+4", nt_p05 = "+5", nt_p06 = "+6", nt_p07 = "+7", nt_p08 = "+8", nt_p09 = "+9")

# Plot
f5a <- ggplot(df) +
  geom_tile(aes(x = Var, y = reorder(Sample, desc(Sample)), fill = median_diff, height = p_hw, width = p_hw, size = p_sig), color = NA) +
  facet_grid(.~position, scales = "free", space = "free_x") +
  scale_fill_gradient2(low = "blue" , mid = "white", high = "red", midpoint = 0, na.value = "grey80",
                       limits = c(-1, 1), oob = scales::squish,
                       name = "variable's median RE - overall's median RE  ") +
  scale_size_manual(values = c(0, 3), name = "Significance  ") +
  xlab("") + ylab("") +
  theme_bw(base_size = 8) + 
  theme(panel.grid = element_blank(), panel.spacing = unit(0, "cm"), 
        legend.position = "top", legend.box = "vertical", legend.text = element_text(margin = ggplot2::margin(r = 2, unit = "pt")),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.text.y = element_text(face = "italic"),
        strip.text.y = element_text(angle = 0), strip.background = element_rect(fill = "white"), aspect.ratio = 8) +
  guides(fill = guide_colourbar(order = 1, barwidth = unit(3, "cm"), barheight = unit(0.25, "cm")), 
         size = guide_legend(order = 2, override.aes = list(color = "white"), keywidth = unit(0.3, "cm"), keyheight = unit(0.3, "cm")))

# Save plot
ggsave(filename = "Figure5A.pdf", plot = f5a, path = "../Plots", width = 6.5, height = 3, units = "in", dpi = 500)

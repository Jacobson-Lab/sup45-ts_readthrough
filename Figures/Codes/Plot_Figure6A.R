# ---------------------------------------
# Figure 6A
# ---------------------------------------
# Table key:
## Sample = sample
## Var = variable (in this case, P-site codon identities)
## samp_median = median readthrough efficiency of all genes in the sample
## median = median readthrough efficiency of genes containing particular variable
## p_val = p-value from Wilcoxon's rank sum test comparing median and samp_median
## median_diff = median - samp_median
## p_hw = a value used for plotting; corresponds to the size of the heatmap tile. 
##        If p-value < 0.05 (significant), tile's height and width = 0.9. Otherwise, tile's height and width = 0.6.
## p_sig = string used to label significance of p-value in plot; related to p_hw.
##         The p-value < 0.05 (significant) is labeled as such, and p-value >= 0.05 is labeled "ns"
## position = feature of interest from the analysis (in this case, P-site codon)
## aa = amino acid encoded by Var

library(ggplot2)

df <- read.table("../Data/Data_Figure6A.txt", header = TRUE)
df$Sample <- factor(df$Sample, levels = unique(df$Sample)[c(1, 2, 5, 7, 3, 4, 6, 8)])
df$aa <- factor(df$aa, levels = c("F", "S", "Y", "*", "C", "W", "L", "P", "H", "Q", "R", "I", "M", "T", "N", "K", "V", "A", "D", "E", "G"))

# Plot
f6a <- ggplot(df) +
  geom_tile(aes(x = Var, y = reorder(Sample, desc(Sample)), fill = median_diff, height = p_hw, width = p_hw, size = p_sig), color = NA) +
  facet_grid(.~aa, scales = "free", space = "free_x") +
  scale_fill_gradient2(low = "blue" , mid = "white", high = "red", midpoint = 0, na.value = "grey80",
                       limits = c(-3, 3), oob = scales::squish,
                       name = "variable's median RE - overall's median RE  ") +
  scale_size_manual(values = c(0, 3), name = "Significance ") +
  xlab("P-site codon") + ylab("") +
  theme_bw(base_size = 8) + 
  theme(panel.grid = element_blank(), panel.spacing = unit(0, "cm"),
        legend.position = "top", legend.box = "vertical", legend.text = element_text(margin = ggplot2::margin(r = 2, unit = "pt")),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(face = "italic"),
        strip.text.y = element_text(angle = 0), strip.background = element_rect(fill = "white"), aspect.ratio = 8) +
  guides(fill = guide_colourbar(order = 1, barwidth = unit(2.5, "cm"), barheight = unit(0.25, "cm")), 
         size = guide_legend(order = 2, override.aes = list(color = "white"), keywidth = unit(0.3, "cm"), keyheight = unit(0.3, "cm")))

# Save plot
ggsave(filename = "Figure6A.pdf", plot = f6a, path = "../Plots", width = 6.5, height = 3, units = "in", dpi = 500)

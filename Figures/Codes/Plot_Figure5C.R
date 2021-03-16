# ---------------------------------------
# Figure 5C
# ---------------------------------------
# Table key:
## Sample = sample
## Var = variable (in this case, stop codon and nucleotide +4 identities)
## samp_median = median readthrough efficiency of all genes in the sample
## median = median readthrough efficiency of genes containing particular variable
## p_val = p-value from Wilcoxon's rank sum test comparing median and samp_median
## median_diff = median - samp_median
## p_hw = a value used for plotting; corresponds to the size of the heatmap tile. 
##        If p-value < 0.05 (significant), tile's height and width = 0.9. Otherwise, tile's height and width = 0.6.
## p_sig = string used to label significance of p-value in plot; related to p_hw.
##         The p-value < 0.05 (significant) is labeled as such, and p-value >= 0.05 is labeled "ns"
## position = feature of interest from the analysis (in this case, stop codon and nucleotide +4)

library(ggplot2)
library(dplyr)

df <- read.table("../Data/Data_Figure5C.txt", header = TRUE)
df$Sample <- factor(df$Sample, levels = unique(df$Sample)[c(1, 2, 5, 7, 3, 4, 6, 8)])

# Create columns of stop codon and nt +4 identities (for plot facet)
df$stop <- gsub('.{1}$', '', df$Var)  # replace the last character in Var with blank (keep first 3 characters) and store in a new column named "stop"
df$ntp04 <- gsub('^.{3}', '', df$Var) # replace the first 3 characters in Var with blank (keep alst character) and store in a new column named "nt_p04"

# Plot
f5c <- ggplot(df) +
  geom_tile(aes(x = ntp04, y = reorder(Sample, desc(Sample)), fill = median_diff, height = p_hw, width = p_hw, size = p_sig), color = NA) +
  facet_grid(.~stop, scales = "free", space = "free_x") +
  scale_fill_gradient2(low = "blue" , mid = "white", high = "red", midpoint = 0, na.value = "grey80",
                       limits = c(-2, 2), oob = scales::squish,
                       name = "variable's median RE - overall's median RE  ") +
  scale_size_manual(values = c(0, 3), name = "Significance  ") +
  xlab("+4") + ylab("") +
  theme_bw(base_size = 8) + 
  theme(panel.grid = element_blank(), panel.spacing = unit(0, "cm"), 
        legend.position = "top", legend.box = "vertical", legend.text = element_text(margin = ggplot2::margin(r = 2, unit = "pt")),
        axis.text.y = element_text(face = "italic"),
        strip.text.y = element_text(angle = 0), strip.background = element_rect(fill = "white"), aspect.ratio = 8) +
  guides(fill = guide_colourbar(order = 1, barwidth = unit(2, "cm"), barheight = unit(0.25, "cm"), title.hjust = 0.5),
         size = guide_legend(order = 2, override.aes = list(color = "white"), keywidth = unit(0.3, "cm"), keyheight = unit(0.3, "cm")))

# Save plot
ggsave(filename = "Figure5C.pdf", plot = f5c, path = "../Plots", width = 4, height = 3, units = "in", dpi = 500)

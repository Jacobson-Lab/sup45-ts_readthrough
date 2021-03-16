# ---------------------------------------
# Figure 5B
# ---------------------------------------
# Table key:
## nt_position = nt position relative to stop codon
## V1 = variable (in this case, nucleotide identities)
## N = raw count of genes in the high or low readthrough "group" containing a particular "variable" at a particular "position" in a "sample"
## fraction = fraction of genes in the high or low readthrough "group" containing a particular "variable" at a particular "position" in a "sample"
## ratio = fraction / fraction from the reference gene set
## Group = readthrough group
## sample = sample
## log2_ratio = log2 of ratio
## RefvsG_padj = adjusted p-value for Reference vs Group (High or Low) comparison obtained from Chi-square test of independence
## HvsL_padj = adjusted p-value for High vs Low comparison obtained from Chi-square test of independence
## overall_p = overall p-value from Chi-square test of independence
## [RefvsG or HvsL or overall]_sig = string used to label significance of p-value in plot; related to [RefvsG or HvsL or overall]_hw.
##                                   The padj < 0.05 (significant) is labeled as such, and padj >= 0.05 is labeled "ns"
## [RefvsG or HvsL]_hw = a value used for plotting; corresponds to the size of the heatmap tile; related to [RefvsG or HvsL or overall]_sig.
##                       If padj < 0.05 (significant), tile's height and width = 0.9. Otherwise, tile's height and width = 0.6.

library(ggplot2)
library(dplyr)

df <- read.table("../Data/Data_Figure5B.txt", header = TRUE)
df$sample <- factor(df$sample, levels = unique(df$sample)[c(1, 2, 5, 7, 3, 4, 6, 8)])
df$nt_position <- recode_factor(df$nt_position, random_factor = "Random", nt_m03 = "-3", nt_m02 = "-2", nt_m01 = "-1", stop_codon = "Stop", 
                                nt_p04 = "+4", nt_p05 = "+5", nt_p06 = "+6", nt_p07 = "+7", nt_p08 = "+8", nt_p09 = "+9")

# Plot
f5b <- ggplot(df) +
  geom_tile(aes(x = V1, y = reorder(sample, desc(sample)), fill = log2_ratio, height = RefvsG_hw, width = RefvsG_hw, size = RefvsG_sig), color = NA) +
  facet_grid(Group~nt_position, scales = "free", space = "free_x") +
  scale_fill_gradient2(low = "blue" , mid = "white", high = "red", midpoint = 0, na.value = "grey70",
                       limits = c(-1.5, 1.5), oob = scales::squish,
                       name = expression(paste(log[2], "(", frac(Group, Reference), ")  "))) +
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
ggsave(filename = "Figure5B.pdf", plot = f5b, path = "../Plots", width = 6.5, height = 4.5, units = "in", dpi = 500)

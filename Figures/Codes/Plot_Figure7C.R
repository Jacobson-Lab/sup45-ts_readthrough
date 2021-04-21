# ---------------------------------------
# Figure 7C
# ---------------------------------------
# Table key:
## sample = sample
## Var1 = variable 1 (in this case, P-site codon)
## Var2 = variable 2 (in this case, stop codon)
## value = Chi-square residuals from Chi-square test of association between Var1 and Var2. Positive = attraction, negative = repulsion.
## pval = p-value obtained from Chi-square test of association between Var1 and Var2
## psig = string used to label significance of p-value in plot; related to p_hw.
##        The p-value < 0.05 (significant) is labeled as such, and p-value >= 0.05 is labeled "ns"
## p_hw = a value used for plotting; corresponds to the size of the heatmap tile; related to p_sig.
##        If p-value < 0.05 (significant), tile's height and width = 0.9. Otherwise, tile's height and width = 0.6.
## aa = amino acid encoded by the codon in Var1

library(ggplot2)
library(ggpubr)

df_all <- read.table("../Data/Data_Figure7C_All.txt", header = TRUE)
df_all$sample <- factor(df_all$sample, levels = unique(df_all$sample)[c(1, 2, 5, 7, 3, 4, 6, 8)])
df_all$aa <- factor(df_all$aa, levels = c("F", "S", "Y", "*", "C", "W", "L", "P", "H", "Q", "R", "I", "M", "T", "N", "K", "V", "A", "D", "E", "G"))

# Plot
f7c <- ggplot(df_all) +
  geom_tile(aes(x = Var1, y = reorder(sample, desc(sample)), fill = value, height = p_hw, width = p_hw)) +
  facet_grid(Var2~aa, scales = "free_x", space = "free") +
  scale_fill_gradient2(low = "blue" , mid = "white", high = "red", midpoint = 0, na.value = "grey70",
                       limits = c(-4, 4), oob = scales::squish,
                       name = expression(paste(chi^2, " residuals "))) +
  xlab("\nP-site codon") + ylab("") +
  theme_bw(base_size = 8) + 
  theme(panel.grid = element_blank(), panel.spacing = unit(0, "cm"), axis.text.y = element_text(face = "italic"),
        legend.position = "top", legend.key.height = unit(0.25, "cm"), strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Save plot
ggsave(filename = "Figure7C.pdf", plot = f7c, path = "../Plots", height = 3.5, width = 7, units = "in", dpi = 500)

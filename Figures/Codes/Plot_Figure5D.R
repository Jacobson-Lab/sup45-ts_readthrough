# ---------------------------------------
# Figure 5D
# ---------------------------------------
# Table key:
## sample = sample
## Var1 = variable 1 (in this case, nt +4)
## Var2 = variable 2 (in this case, stop codon)
## value = Chi-square residuals from Chi-square test of association between Var1 and Var2. Positive = attraction, negative = repulsion.
## pval = p-value obtained from Chi-square test of association between Var1 and Var2
## psig = string used to label significance of p-value in plot; related to p_hw.
##        The p-value < 0.05 (significant) is labeled as such, and p-value >= 0.05 is labeled "ns"
## p_hw = a value used for plotting; corresponds to the size of the heatmap tile; related to p_sig.
##        If p-value < 0.05 (significant), tile's height and width = 0.9. Otherwise, tile's height and width = 0.6.

library(ggplot2)
library(ggpubr)

df_all <- read.table("../Data/Data_Figure5D_All.txt", header = TRUE)
df_all$sample <- factor(df_all$sample, levels = unique(df_all$sample)[c(1, 2, 5, 7, 3, 4, 6, 8)])

df_High <- read.table("../Data/Data_Figure5D_High.txt", header = TRUE)
df_High$sample <- factor(df_High$sample, levels = unique(df_High$sample)[c(1, 2, 5, 7, 3, 4, 6, 8)])

df_Low <- read.table("../Data/Data_Figure5D_Low.txt", header = TRUE)
df_Low$sample <- factor(df_Low$sample, levels = unique(df_Low$sample)[c(1, 2, 5, 7, 3, 4, 6, 8)])

# Plot
f5d_plot <- function(df, x_axis_label = "") {
  p <- ggplot(df) +
    geom_tile(aes(x = Var1, y = reorder(sample, desc(sample)), fill = value, height = p_hw, width = p_hw)) +
    facet_grid(.~Var2) +
    scale_fill_gradient2(low = "blue" , mid = "white", high = "red", midpoint = 0, na.value = "grey70",
                         limits = c(-4, 4), oob = scales::squish,
                         name = expression(paste(chi^2, " residuals "))) +
    xlab(x_axis_label) + ylab("") +
    theme_bw(base_size = 8) + 
    theme(panel.grid = element_blank(), panel.spacing = unit(0, "cm"), axis.text.y = element_text(face = "italic"),
          legend.position = "top", legend.key.height = unit(0.25, "cm"), strip.background = element_rect(fill = "white"))
  return(p)
}
  
p_all <- f5d_plot(df = df_all, x_axis_label = "nt +4")
p_High <- f5d_plot(df = df_High, x_axis_label = "nt +4")
p_Low <- f5d_plot(df = df_Low, x_axis_label = "nt +4")

# Arrange plots
f5d <- ggarrange(plotlist = lapply(list(p_all, p_High, p_Low), "+", coord_equal()), 
                 common.legend = TRUE, nrow = 1, 
                 labels = c("All", "High", "Low"), label.x = 0.66, hjust = 0.5, font.label = list(size = 8, face = "plain"))

# Save plot
ggsave(filename = "Figure5D.pdf", plot = f5d, path = "../Plots", height = 2.5, width = 7, units = "in", dpi = 500)

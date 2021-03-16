# ---------------------------------------
# Figure S2A
# ---------------------------------------
# Table key:
## gene = gene
## sample = sample
## ptc_re = efficiency of PTC readthrough
## rep = replicate number

library(ggplot2)
#library(ggpubr)
library(dplyr)
library(reshape2)

df <- read.table("../Data/Data_FigureS2A.txt", header = TRUE)
df$sample <- gsub(", ", "\n", df$sample)
df$sample <- factor(df$sample, levels = unique(df$sample))

# Calculate average
df2 <- dcast(df, gene + sample ~ rep, value.var = "ptc_re")
df2$mean <- rowMeans(df2[, c("rep1", "rep2")], na.rm = TRUE)

# Plot
fs2a <- ggplot(data = df, aes(x = sample, y = ptc_re)) +
  #geom_bar(data = df2, aes(x = sample, y = mean, color = sample), stat = "identity", fill = "white") + # Option 1: Add mean bar using df2
  stat_summary(fun = mean, geom = "bar", aes(fill = sample), color = NA) + # Option 2: Add mean bar calculated directly from df in ggplot
  geom_point(position = position_dodge2(width = 0.75), size = 0.75, color = "black") +
  geom_text(data = df[c(1, 16), ], aes(label = gene, x = levels(sample)[2], y = max(ptc_re)), # subset df because geom_text over-plots (it plots every line of data)
            fontface = "italic", size = 8 / .pt, hjust = 0) +
  facet_grid(.~gene) +
  scale_fill_manual(values = c("#ed1c24", "#f7941d", "#00a651", "#8dc63f")) +
  xlab("Sample") + ylab("# of in-frame footprints\nDownstream of PTC / Upstream of PTC") +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(), strip.background = element_blank(), strip.text = element_blank(),
        legend.position = "none", axis.text.x = element_text(face = "italic", angle = 90, hjust = 1, vjust = 0.5))

# Save plot
ggsave(filename = "FigureS2A.pdf", plot = fs2a, path = "../Plots/", height = 2.5, width = 3.5, unit = "in", dpi = 500)


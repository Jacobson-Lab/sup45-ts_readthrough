# ---------------------------------------
# Figure 4B
# ---------------------------------------

library(ggplot2)
library(dplyr)

df <- read.table("../Data/Data_Figure4B_AUROC_average.txt", header = TRUE)
df$sample <- recode_factor(df$sample, sup45_wt_25C = "SUP45\n25 °C", sup45_ts_25C = "sup45-ts\n25 °C",
                           sup45_wt_37C = "SUP45\n37 °C",  sup45_ts_37C = "sup45-ts\n37 °C", 
                           sup45_wt = "SUP45-D", sup45_d = "sup45-d", 
                           rli1_wt = "RLI1-D", rli1_d = "rli1-d")

# Set color for each sample (order is by the recode_factor)
p_col <- c("#ed1c24", "#f7941d", "#00a651", "#8dc63f", "#2e3192", "#00aeef", "#662d91", "#ec008c")

# Plot
f4b <- ggplot(data = df, aes(x = sample, y = ROC, fill = sample)) +
  geom_bar(stat = "identity", position = "dodge", color = NA) +
  geom_errorbar(aes(ymin = ROC-sd, ymax = ROC+sd), width = 0.2, position = position_dodge(0.9)) +
  geom_hline(yintercept = 0.5, color = "black", linetype = 2) +
  ylab("AUROC") + xlab("Sample") +
  coord_cartesian(ylim = c(0, 1)) +
  scale_fill_manual(values = p_col) +
  theme_bw(base_size = 8) + 
  theme(panel.grid = element_blank(), axis.text.x = element_text(face = "italic", angle = 90, hjust = 1, vjust = 0.5), 
        legend.position = "none")

# Save plot
ggsave(filename = "Figure4B_AUROC.pdf", plot = f4b, path = "../Plots", width = 3, height = 1.5, units = "in", dpi = 300)

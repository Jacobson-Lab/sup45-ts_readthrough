# ---------------------------------------
# Figure 7
# ---------------------------------------
# Table key:
## Sample = sample
## log_rte = readthrough efficiency (log2 transformed) as calculated by the formula in Figure 3A
## value = value for particular feature
## feature = feature (X variable) tested (in this case, l_utr3 = 3'-UTR length)

library(ggplot2)
library(ggpubr)
library(dplyr)

df <- read.table("../Data/Data_Figure7.txt", header = TRUE)
df$Sample <- recode_factor(df$Sample, sup45_wt_25C = "SUP45, 25 °C", sup45_wt_37C = "SUP45, 37 °C", sup45_wt = "SUP45-D", rli1_wt = "RLI1-D", 
                           sup45_ts_25C = "sup45-ts, 25 °C", sup45_ts_37C = "sup45-ts, 37 °C", sup45_d = "sup45-d", rli1_d = "rli1-d")

# Plot
f7 <- ggplot(df, aes(x = value, y = log_rte)) +
  geom_point(alpha = 1, color = "black", size = 0.25) +
  geom_smooth(method = "lm", formula = y~x) +
  stat_cor(method = "spearman", geom = "label", label.x.npc = "middle", label.y = 5.5, hjust = 0.5, vjust = 0.5, size = 2, cor.coef.name = "rho") +
  facet_wrap(~Sample, nrow = 2) +
  scale_x_log10() +
  xlab("3'-UTR length (nt)") + ylab("Readthrough Efficiency") + 
  theme_bw(base_size = 8) + 
  theme(panel.grid = element_blank(), panel.spacing.x = unit(0.2,"cm"), panel.spacing.y = unit(0.2, "cm"),
        strip.text = element_text(face = "italic"), strip.background = element_rect(fill = "white"), aspect.ratio = 1)

# Save plot
ggsave(filename = "Figure7.pdf", plot = f7, path = "../Plots/", height = 3, width = 6, units = "in", dpi = 500)
ggsave(filename = "Figure7_eps.eps", plot = f7, path = "../Plots/", height = 3, width = 6, units = "in", dpi = 500)


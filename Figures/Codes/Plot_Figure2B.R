# ---------------------------------------
# Figure 2B
# ---------------------------------------
# Table key:
## sample = sample
## psite_region = 5'-UTR, CDS, Extension, or Distal 3'-UTR
## frame = reading frame; 0, 1, or 2; 0 = in-frame translation in the CDS
## count = Number of reads in a particular reading frame in a region of a sample
## sum = Sum of read counnts in the three reading frame in each region of a sample
## fraction = count / sum
## dataset = data set the sample belongs to (sup45-ts, sup45-d, or rli1)

library(ggplot2)
#library(ggpubr)
library(dplyr)
library(reshape2)

df <- read.table("../Data/Data_Figure2B_raw_count.txt", header = TRUE)
df$Sample <- gsub("_rep[0-9]", "", df$sample)
df$rep <- gsub(".*_r", "r", df$sample)
df$Sample <- recode_factor(df$Sample, sup45_wt_25C = "SUP45\n25 °C", sup45_ts_25C = "sup45-ts\n25 °C",
                           sup45_wt_37C = "SUP45\n37 °C",  sup45_ts_37C = "sup45-ts\n37 °C", 
                           sup45_wt = "SUP45-D", sup45_d = "sup45-d", 
                           rli1_wt = "RLI1-D", rli1_d = "rli1-d")
df$psite_region <- recode_factor(df$psite_region, `5utr` = "5'-UTR", cds = "CDS", `3utr_extension` = "Extension", `3utr_distal` = "Distal 3'-UTR")
df$frame <- recode_factor(df$frame, `0` = "Frame 0", `1` = "Frame 1", `2` = "Frame 2")

# Calculate average
df2 <- dcast(df, psite_region + frame + Sample ~ rep, value.var = "fraction")
df2$mean <- rowMeans(df2[, c("rep1", "rep2", "rep3")], na.rm = TRUE)

# Plot
p_col <- c("#ed1c24", "#f7941d", "#00a651", "#8dc63f", "#2e3192", "#00aeef", "#662d91", "#ec008c")
f2b <- ggplot(data = df, aes(x = Sample, y = fraction, fill = Sample)) +
  #geom_bar(data = df2, aes(x = Sample, y = mean, color = Sample), stat = "identity", fill = "white") + # Option 1: Add mean bar using df2
  stat_summary(fun = mean, geom = "bar", aes(fill = Sample), color = NA) + # Option 2: Add mean bar calculated directly from df in ggplot
  geom_point(position = position_dodge2(width = 0.75), size = 0.75, color = "black") +
  geom_hline(yintercept = 1/3, color = "gray30", linetype = 2, size = 0.5) +
  facet_grid(psite_region~frame, scales = "fixed") +
  scale_fill_manual(values = p_col) +
  xlab("Sample") + ylab("Fraction of ribosome footprints") +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(), strip.background = element_rect(fill = "white"), strip.text.y = element_text(angle = 270),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"), legend.position = "none")

# Save plot
ggsave(filename = "Figure2B.pdf", plot = f2b, path = "../Plots/", height = 5.5, width = 5.5, unit = "in", dpi = 500)

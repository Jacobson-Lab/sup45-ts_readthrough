# ---------------------------------------
# Figure S1
# ---------------------------------------

library(ggplot2)
library(scales)
library(ggpubr)
library(dplyr)

# Function to combine replicate RPKM into the same table
rep_rpkm <- function(path_rep1, path_rep2, sample_name) {
  rep1 <- read.table(path_rep1, header = TRUE)
  rep2 <- read.table(path_rep2, header = TRUE)
  df_FPKM <- data.frame(r1 = rep1$FPKM, r2 = rep2$FPKM)
  df_FPKM$sample <- sample_name
  return(df_FPKM)
}

# Ribo-seq
folder_path <- "../../RSEM_results/Ribo-Seq_RSEM/"
files <- list.files(path = folder_path, pattern = "RPF-.*.isoforms.results", full.names = TRUE)
df <- list()
df$w25 <- rep_rpkm(path_rep1 = files[grepl(pattern = ".*WTA25.*", files)], path_rep2 = files[grepl(pattern = ".*WTB25.*", files)], sample_name = "SUP45, 25 °C")
df$w37 <- rep_rpkm(path_rep1 = files[grepl(pattern = ".*WTA37.*", files)], path_rep2 = files[grepl(pattern = ".*WTB37.*", files)], sample_name = "SUP45, 37 °C")
df$m25 <- rep_rpkm(path_rep1 = files[grepl(pattern = ".*mutA25.*", files)], path_rep2 = files[grepl(pattern = ".*mutB25.*", files)], sample_name = "sup45-ts, 25 °C")
df$m37 <- rep_rpkm(path_rep1 = files[grepl(pattern = ".*mutA37.*", files)], path_rep2 = files[grepl(pattern = ".*mutB37.*", files)], sample_name = "sup45-ts, 37 °C")
df <- bind_rows(df)
df$sample <- factor(df$sample, levels = c("SUP45, 25 °C", "SUP45, 37 °C", "sup45-ts, 25 °C", "sup45-ts, 37 °C"))
df$seqtype <- "Ribosome profiling"

p_ribo <- ggplot(df, aes(x = r1, y = r2)) +
  geom_point(size = 0.5) + 
  geom_abline(slope = 1, color = "red") +
  facet_grid(seqtype~sample) +
  stat_cor(method = "pearson", geom = "label", label.x.npc = "left", label.y = 4.7, 
           hjust = 0, vjust = 0.5, size = 2, cor.coef.name = "r", aes(label = ..r.label..)) +
  scale_x_continuous(trans = "log10", labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(trans = "log10", labels = trans_format("log10", math_format(10^.x))) +
  coord_cartesian(xlim = c(10^-2, 10^5), ylim = c(10^-2, 10^5)) +
  xlab("Replicate 1 (RPKM)") + ylab("Replicate 2 (RPKM)") +
  theme_bw(base_size = 8) + 
  theme(strip.placement = "outside", strip.background = element_rect(fill = "white"), strip.text.x = element_text(hjust = 0.5, face = "italic"),
        panel.grid = element_blank(), panel.spacing.x = unit(0.2,"cm"), aspect.ratio = 1)

# RNA-seq
folder_path <- "../../RSEM_results/RNA-Seq_RSEM/"
files <- list.files(path = folder_path, pattern = "TotalRNA.*.isoforms.results", full.names = TRUE)
df <- list()
df$w25 <- rep_rpkm(path_rep1 = files[grepl(pattern = ".*WTA25.*", files)], path_rep2 = files[grepl(pattern = ".*WTB25.*", files)], sample_name = "SUP45, 25 °C")
df$w37 <- rep_rpkm(path_rep1 = files[grepl(pattern = ".*WTA37.*", files)], path_rep2 = files[grepl(pattern = ".*WTB37.*", files)], sample_name = "SUP45, 37 °C")
df$m25 <- rep_rpkm(path_rep1 = files[grepl(pattern = ".*mutA25.*", files)], path_rep2 = files[grepl(pattern = ".*mutB25.*", files)], sample_name = "sup45-ts, 25 °C")
df$m37 <- rep_rpkm(path_rep1 = files[grepl(pattern = ".*mutA37.*", files)], path_rep2 = files[grepl(pattern = ".*mutB37.*", files)], sample_name = "sup45-ts, 37 °C")
df <- bind_rows(df)
df$sample <- factor(df$sample, levels = c("SUP45, 25 °C", "SUP45, 37 °C", "sup45-ts, 25 °C", "sup45-ts, 37 °C"))
df$seqtype <- "RNA-Seq"

p_rna <- ggplot(df, aes(x = r1, y = r2)) +
  geom_point(size = 0.5) + 
  geom_abline(slope = 1, color = "red") +
  facet_grid(seqtype~sample) +
  stat_cor(method = "pearson", geom = "label", label.x.npc = "left", label.y = 4.7, 
           hjust = 0, vjust = 0.5, size = 2, cor.coef.name = "r", aes(label = ..r.label..)) +
  scale_x_continuous(trans = "log10", labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(trans = "log10", labels = trans_format("log10", math_format(10^.x))) +
  coord_cartesian(xlim = c(10^-2, 10^5), ylim = c(10^-2, 10^5)) +
  xlab("Replicate 1 (RPKM)") + ylab("Replicate 2 (RPKM)") +
  theme_bw(base_size = 8) + 
  theme(strip.placement = "outside", strip.background = element_rect(fill = "white"), strip.text.x = element_text(hjust = 0.5, face = "italic"),
        panel.grid = element_blank(), panel.spacing.x = unit(0.2,"cm"), aspect.ratio = 1)

# Arrange plots
fs1 <- cowplot::plot_grid(p_ribo, p_rna, nrow = 2)

# Save plot
ggsave(filename = "FigureS1.pdf", plot = fs1, path = "../Plots", width = 7, height = 4, units = "in", dpi = 300)
#ggsave(filename = "FigureS1.eps", plot = fs1, path = "../Plots", width = 7, height = 4, units = "in", dpi = 300)


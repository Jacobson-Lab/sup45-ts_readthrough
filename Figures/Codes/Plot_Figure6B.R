# ---------------------------------------
# Figure 6B
# ---------------------------------------
# Table key:
## codon_m01 = P-site codon (amino acid -1 relative to stop codon)
## all_count = count of each codon out of all 6575 genes
## oa_count = count of each codon out of genes with UTR annotations
## all_fraction = all_count converted to fraction (each value / sum of all values)
## oa_fraction = oa_count converted to fraction (each value / sum of all values)
## aa = amino acid
## tAI = tRNA adaptation index (from Pechmann S, Frydman J, Nat Struct Mol Biol ,2013; Reis M dos, Savva R, Wernisch L, Nucleic Acids Res, 2004)
## anti_codon = specify wobble interactions (from Johansson MJO, Esberg A, Huang B, Björk GR, Byström AS, Mol Cell Biol, 2008)
## column 9-14: anti-codon modifications (from Johansson MJO, Esberg A, Huang B, Björk GR, Byström AS, Mol Cell Biol, 2008)

library(ggplot2)

pc2 <- read.table("../Data/Data_Figure6B.txt", header = TRUE)
pc2$aa <- factor(pc2$aa, levels = c("F", "S", "Y", "*", "C", "W", "L", "P", "H", "Q", "R", "I", "M", "T", "N", "K", "V", "A", "D", "E", "G"))

# Plot
pos <- seq(-0.008, -0.056, -0.008)  # plot position of anti-codon modifications
ss <- 4 # "x" shape (ggplot's shape #4) that indicates anti-codon modification
sz <- 1 # size of "x" shape that indicates anti-codon modification

pp <- ggplot(pc2) + 
  # frequency (bar height) and tAI (bar fill color)
  geom_bar(aes(x = codon_m01, y = oa_fraction, fill = tAI), stat = "identity", color = "black", size = 0.2) +
  scale_y_continuous(breaks = c(0, 0.025, 0.050, 0.075, 0.100)) + # specify y-axis number labels
  scale_fill_gradient2(low = "firebrick2", mid = "white", high = "forestgreen", midpoint = 0.47) + # specify color for tAI
  # anti-codon wobble (shape)
  geom_point(aes(x = codon_m01, y = pos[1], shape = anti_codon), size = sz) + 
  scale_shape_manual(name = "   Wobble", values = c(2, 7, 0, 15, 17), na.translate = FALSE) + # specify shapes for wobble pairs
  # anti-codon modification information
  geom_point(aes(x = codon_m01, y = pos[2], color = pseudoU), shape = ss, na.rm = TRUE, size = sz) +
  geom_point(aes(x = codon_m01, y = pos[3], color = m), shape = ss, na.rm = TRUE, size = sz) +
  geom_point(aes(x = codon_m01, y = pos[4], color = m5), shape = ss, na.rm = TRUE, size = sz) +
  geom_point(aes(x = codon_m01, y = pos[5], color = mcm5), shape = ss, na.rm = TRUE, size = sz) +
  geom_point(aes(x = codon_m01, y = pos[6], color = ncm5), shape = ss, na.rm = TRUE, size = sz) +
  geom_point(aes(x = codon_m01, y = pos[7], color = s2), shape = ss, na.rm = TRUE, size = sz) +
  scale_color_manual(values = rep("black", times = 6), na.value = NA) +
  # group codons by amino acid identities
  facet_grid(.~aa, scales = "free_x", space = "free_x") + 
  xlab("\nP-site codon") + ylab("Anti-codon      Frequency in the\nmodification     reference genome\n") +
  theme_bw(base_size = 8) + 
  theme(panel.grid = element_blank(), panel.spacing = unit(0,"cm"), strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), aspect.ratio = 20,
        legend.position = "bottom", legend.box = "horizontal", legend.text = element_text(margin = ggplot2::margin(r = 2, unit = "pt"))) +
  guides(fill = guide_colourbar(order = 1, barwidth = grid::unit(2.5, "cm")), 
         shape = guide_legend(order = 2), color = guide_none())

# Labels for anti-codon modifications
mod_df <- data.frame(aa = rep(c("F", "S", "Y", "C", "W", "L", "P", "H", "Q", "R", "I", "M", "T", "N", "K", "V", "A", "D", "E", "G"), each = 6),
                     x = c(rep(pos[1], times = 6), rep(pos[1], times = (20*6)-6)),
                     y = c(pos[2:7], rep(NA, times = (20*6)-6)), 
                     mod = c("psi * '   '", "m * '   '", "m^5 * '   '", "mcm^5 * '   '", "ncm^5 * '   '", "s^2 * '   '", rep(NA, times = (20*6)-6)))
mod_df$aa <- factor(mod_df$aa, levels = c("F", "S", "Y", "*", "C", "W", "L", "P", "H", "Q", "R", "I", "M", "T", "N", "K", "V", "A", "D", "E", "G"))

# Add labels for anti-codon modifications to plot
f6b <- pp + geom_text(data = mod_df, aes(x = x, y = y, label = mod), hjust = 1, parse = TRUE, size = 2) + 
  coord_cartesian(clip = "off")

# Save plot
ggsave(filename = "Figure6B.pdf", plot = f6b, height = 3, width = 6.5, units = "in", dpi = 500, path = "../Plots/")  

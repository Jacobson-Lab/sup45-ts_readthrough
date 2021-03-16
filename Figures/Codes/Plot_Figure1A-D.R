# ---------------------------------------
# Figure 1A-D
# ---------------------------------------

library(ggplot2)
library(dplyr)

df <- read.table("../Data/Data_Figure1A-D.txt", header = TRUE)
df$sample <- recode_factor(df$sample, sup45_wt_25C = "SUP45, 25 °C", sup45_ts_25C = "sup45-ts, 25 °C",
                           sup45_wt_37C = "SUP45, 37 °C",  sup45_ts_37C = "sup45-ts, 37 °C", 
                           sup45_wt = "SUP45-D", sup45_d = "sup45-d", 
                           rli1_wt = "RLI1-D", rli1_d = "rli1-d")

# Set color for each sample (order is by the recode_factor)
f1_col <- c("#ed1c24", "#f7941d", "#00a651", "#8dc63f", "#2e3192", "#00aeef", "#662d91", "#ec008c")

# Define vertical line for every in-frame nt
vline <- data.frame(reg = c(rep("Distance from start (nt)", times = 21), rep("Distance from stop (nt)", times = 43)), 
                    xin = c(seq(-21, 39, 3), seq(-41, 85, 3)))

# Define vertical line for start and stop
rline <- data.frame(reg = c("Distance from start (nt)", "Distance from stop (nt)"), xin = c(0, 1))

# Box + line coordinates for zooming
boxx <- data.frame(reg = c("Distance from start (nt)", "Distance from stop (nt)"),
                   xmin = c(NA, 1), xmax = c(NA, 80), ymin = c(NA, -0.01), ymax = c(NA, 0.01))
bll <- data.frame(reg = c("Distance from start (nt)", "Distance from stop (nt)"),
                  x = c(NA, 1), y = c(NA, 0.01), xend = c(NA, 5), yend = c(NA, 0.03))
blr <- data.frame(reg = c("Distance from start (nt)", "Distance from stop (nt)"),
                  x = c(NA, 80), y = c(NA, 0.01), xend = c(NA, 78), yend = c(NA, 0.03))

# Plot
f1_set <- unique(df$dataset)
f1_main <- list()
f1_inset <- list()
for (i in 1:length(f1_set)) {
  f1_main[[f1_set[i]]] <- ggplot() + 
    geom_vline(data = vline, aes(xintercept = xin), linetype = "dotted", color = "grey", size = 0.5) +
    geom_line(data = df[df$dataset == f1_set[i], ], aes(x = distance, y = fraction2, color = sample), size = 0.7) +
    geom_vline(data = rline, aes(xintercept = xin), color = "red", size = 0.5) +
      facet_grid(.~reg, scales = "free_x", space = "free_x", switch = "x") +
    geom_rect(data = boxx, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = NA, color = "black", size = 0.2) +
    geom_segment(data = bll, aes(x = x, y = y, xend = xend, yend = yend), color = "black", size = 0.2) +
    geom_segment(data = blr, aes(x = x, y = y, xend = xend, yend = yend), color = "black", size = 0.2) +
    scale_x_continuous(expand = expansion(mult = 0.01)) +
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.1))) +
    scale_color_manual(values = c(f1_col[(2*i-1):(2*i)])) +
    ylab("Normalized reads") + xlab("") + labs(color = "") +
    coord_cartesian(ylim = c(0, 0.21), clip = "off") +
    theme_bw(base_size = 8, base_family = "sans") + 
    theme(strip.text.y = element_text(angle = 0), strip.placement = "outside", strip.background = element_rect(fill = "white", color = NA),
          panel.grid = element_blank(), panel.spacing.y = unit(0, "cm"), panel.spacing.x = unit(0.2, "cm"),
          plot.margin = ggplot2::margin(r = c(0.1, 0.1, 0.1, 0.1), unit = "cm"), 
          legend.position = c(0.45, 0.8), legend.background = element_blank(), legend.key = element_rect(fill = NA),
          legend.spacing.y = unit(0.5, "cm"), legend.text = element_text(face = "italic")) +
    guides(color = guide_legend(override.aes = list(size = 2)))
  
  f1_inset[[f1_set[i]]] <- ggplot() + 
    geom_vline(xintercept = c(seq(-44, 85, 3)), linetype = "dotted", color = "grey", size = 0.5) +
    geom_line(data = df[df$dataset == f1_set[i] & df$reg == "Distance from stop (nt)", ],
              aes(x = distance, y = fraction2, color = sample), size = 0.7) +
    scale_x_continuous(expand = expansion(mult = 0)) +
    scale_y_continuous(expand = expansion(mult = c(0.01))) +
    scale_color_manual(values = c(f1_col[(2*i-1):(2*i)])) +
    coord_cartesian(xlim = c(1, 80), ylim = c(0, 0.005)) +
    theme_void() + 
    theme(legend.position = "none", plot.margin = grid::unit(c(0, 0, 0, 0), "mm"), plot.background = element_rect(color = "black"))
}

# Combine main + inset plots
f1_comb <- list()
for (i in 1:length(f1_set)) {
  f1_comb[[i]] <- f1_main[[i]] +
    bobfunctions2::gg_inset(grob = ggplotGrob(f1_inset[[i]]), data = data.frame(reg = "Distance from stop (nt)"), ymin = 0.03, ymax = 0.21, xmin = 5, xmax = 78)
}

# Arrange different data sets vertically
cowplot::plot_grid(f1_comb[[1]], f1_comb[[2]], f1_comb[[3]], f1_comb[[4]], ncol = 1)

# Save plot
ggsave(filename = "Figure1A-D.pdf", plot = last_plot(), path = "../Plots", width = 7, height = 7, units = "in", dpi = 300)

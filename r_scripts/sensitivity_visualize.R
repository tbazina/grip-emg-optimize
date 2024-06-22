# Visalize sensitivity analysis results
# Install ggmagnify package from r-universe
# install.packages(
#   "ggmagnify",
#   repos = c("https://hughjonesd.r-universe.dev",  "https://cloud.r-project.org")
#   )
library(tidyverse)
library(janitor)
library(magrittr)
library(ggsci)
library(ggmagnify)
library(ggfx)

# Load data
# Position 3 or 4
position = 3
s_res <- read_csv(
  paste0('results/sensitivity_preliminary_lh_pos', position, '.csv')
  )
s_res

# Preliminary plots from lhs sensitivity analysis to obtain proper parameter ranges
prel_sens <- s_res %>% 
  # Pivot window and decay to longer format
  pivot_longer(
    cols = c(window, decay), names_to = 'parameter', values_to = 'value'
    ) %>% 
  ggplot(aes(x = value, y = corrs)) +
  facet_wrap(~ parameter, scales = 'free_x', ncol = 1) +
  geom_point(aes(fill = parameter), shape = 21, stroke = 0, size = 0.3) +
  scale_x_continuous(
    name = 'Exponential decay/window size',
    # Reduce spacing to the left and right side
    expand = expansion(mult = c(0.01, 0.01), add = c(0., 0.)),
  ) +
  scale_y_continuous(
    name = 'Maximum lagged cross-correlation',
    breaks = seq(0, 1, 0.05),
  ) +
  scale_fill_nejm() +
  theme_bw() + theme(
    legend.position = 'none',
    legend.title = element_text(
      colour="black", size = 5, margin = margin(0, 2, 0, 0, 'mm')
      ),
    legend.box = 'horizontal',
    legend.direction = 'horizontal',
    legend.box.spacing = unit(0.5, 'mm'),
    legend.spacing.y = unit(0, 'mm'),
    legend.spacing.x = unit(0, 'mm'),
    legend.margin = margin(0, 0, 0, 0, 'mm'),
    legend.key.spacing = unit(0, 'mm'),
    legend.key.size = unit(2, 'mm'),
    legend.text = element_text(
      colour="black", size = 5, margin = margin(0, 0.3, 0, 0, 'mm')
      ),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, 'mm'),
    panel.background = element_blank(),
    panel.spacing.y = unit(0, 'mm'),
    panel.spacing.x = unit(0, 'mm'),
    axis.title = element_text(face="bold", size = 5),
    axis.text = element_text(
      color="black", size = 5, margin = margin(0.0, 0.0, 0.0, 0.0, 'mm')
      ),
    # Remove x axis text and title
    axis.line = element_line(linewidth = 0.1, colour = "black"),
    plot.title = element_text(hjust = 0.5, size = 5),
    panel.grid = element_line(colour = 'grey', size = 0.1),
    panel.border = element_rect(linewidth = 0.1),
    strip.background = element_rect(linewidth = 0.01),
    strip.text = element_text(
      colour = 'black', size = 5.0, margin = margin(b = 0.3, t = 0.3, unit='mm')
    ),
    axis.ticks = element_line(linewidth = 0.1),
    axis.ticks.y = element_blank(),
    axis.ticks.length = unit(0.1, 'lines')
  ) +
  # Add at the end to inherit theme settings
  geom_magnify(
    aes(from = (
      parameter == 'window' & corrs > 0.9 & value <= max(value)
      )),
    to = list(xmin = 50., ymin = 0.6, xmax = max(s_res$window), ymax = 0.75),
    axes = 'xy',
    shadow = T,
    expand = c(0.01, 0.01),
    corners = 0.05,
    proj = 'corresponding',
    recompute = T,
    scale.inset = 1,
    # Borders and projection lines
    linetype = 1,
    linewidth = 0.2,
  ) +
  geom_magnify(
    aes(from = (
      parameter == 'decay' & corrs > 0.9 & value <= max(value)
      )),
    to = list(xmin = 0.00, ymin = 0.55, xmax = 0.01, ymax = 0.75),
    axes = 'xy',
    shadow = T,
    expand = c(0.01, 0.01),
    corners = 0.05,
    proj = 'corresponding',
    recompute = T,
    scale.inset = 1,
    # Borders and projection lines
    linetype = 1,
    linewidth = 0.2,
  ) 
prel_sens

# Save plot
ggsave(
  filename = paste0('plots/preliminary_sensitivity_pos', position, '.png'),
  plot = prel_sens, width = 18, height = 10, units = 'cm', dpi = 320
  )


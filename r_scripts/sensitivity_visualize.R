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
# Position 3 and 4, preliminary latin and sobol
position = 3
s_res <- read_csv(
  paste0('results/sensitivity_preliminary_lh_pos', position, '.csv')
  ) %>%
  mutate(position = position) 
sobol_res <- read_csv(
  paste0(
    'results/sensitivity_groupsTrue_sobol_samples65536_resamples65536_pos',
    position, '.csv'
    )
  ) %>%
  mutate(position = position)
position = 4
s_res <- s_res %>% 
  bind_rows(
    read_csv(
      paste0('results/sensitivity_preliminary_lh_pos', position, '.csv')
      ) %>%
      mutate(position = position)
    )
sobol_res <- sobol_res %>% 
  bind_rows(
    read_csv(
      paste0(
        'results/sensitivity_groupsTrue_sobol_samples65536_resamples65536_pos',
        position, '.csv'
        )
      ) %>%
      mutate(position = position)
    )

# Custome palette using nejm and starting from second color
scale_fill_nejm_skip <- function(..., aesthetics = 'fill') {
  discrete_scale(
    aesthetics = aesthetics,
    # Create palette and skip first
    palette = \(x) pal_nejm()(8)[seq(2, x+1)]
    )
}

# Preliminary plots from lhs sensitivity analysis to obtain proper parameter ranges
prel_sens <- s_res %>% 
  mutate(
    # Rename position to Position: x
    position = case_when(
      position == 3 ~ 'Position: 1',
      position == 4 ~ 'Position: 2'
      )
  ) %>% 
  # Pivot window and decay to longer format
  pivot_longer(
    cols = c(window, decay), names_to = 'parameter', values_to = 'value'
    ) %>% 
  mutate(
    # Rename window to Window size and decay to Decay factor
    parameter = case_when(
      parameter == 'window' ~ 'Window size',
      parameter == 'decay' ~ 'Decay factor'
      ),
    # Wrap parameter
    parameter = str_wrap(parameter, width = 10)
  ) %>% 
  ggplot(aes(x = value, y = corrs)) +
  facet_wrap(position ~ parameter, scales = 'free_x', ncol = 1, strip.position = 'right') +
  geom_point(aes(fill = parameter), shape = 21, stroke = 0, size = 0.3) +
  scale_x_continuous(
    name = 'Exponential decay/window size',
    # Reduce spacing to the left and right side
    expand = expansion(mult = c(0.01, 0.01), add = c(0., 0.)),
  ) +
  scale_y_continuous(
    name = 'Mean peak cross-correlation',
    breaks = seq(0, 1, 0.1),
  ) +
  scale_fill_nejm_skip() +
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
      parameter == 'Window\nsize' & corrs > 0.86 & value <= 495 & value >= 150
      )),
    to = list(xmin = 150., ymin = 0.5, xmax = 495, ymax = 0.75),
    axes = '',
    shadow = T,
    expand = c(0.0, 0.0),
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
      parameter == 'Decay\nfactor' & value <= 0.02 & corrs > 0.86
      )),
    to = list(xmin = 0.00, ymin = 0.5, xmax = 0.02, ymax = 0.80),
    axes = '',
    shadow = T,
    expand = c(0.0, 0.0),
    corners = 0.05,
    proj = 'facing',
    recompute = T,
    scale.inset = 1,
    # Borders and projection lines
    linetype = 1,
    linewidth = 0.2,
  ) 
prel_sens

# Save plot
ggsave(
  filename = paste0('plots/preliminary_sensitivity_decay_window.png'),
  plot = prel_sens, width = 12, height = 6, units = 'cm', dpi = 360
  )

# Plot results from sobol sensitivity analysis - bar + CI
# Plot total ST and first order S1 indices
sobol_res_long <- sobol_res %>% 
  pivot_longer(
    cols = c(ST, S1), names_to = 'indice', values_to = 'sensitivity'
    ) %>%
  mutate(
    sensitivity_conf = case_when(
      indice == 'ST' ~ ST_conf,
      indice == 'S1' ~ S1_conf
    )
  ) %>% 
  select(-ST_conf, -S1_conf) %>% 
  group_by(position, indice) %>%
  mutate(
    indice_sum = sum(sensitivity)
  ) %>% 
  ungroup() %>% 
  mutate(
    # Rename variables to full names
    vars = case_when(
      vars == 'freqs' ~ 'Spectral mask',
      vars == 'decay' ~ 'Decay factor',
      vars == 'window' ~ 'Window size',
      ),
    # Rename sensitivity indices to full names
    indice = case_when(
      indice == 'ST' ~ 'Total Sobol indices',
      indice == 'S1' ~ 'First order Sobol indices'
    ),
    position_vars = factor(paste0(position, '_', vars)),
    position = case_when(
      position == 3 ~ 'Position: 1',
      position == 4 ~ 'Position: 2'
    ),
    position = ordered(position, levels = c('Position: 1', 'Position: 2')),
    # Rename position levels to Position: x
    vars = ordered(vars, levels = c('Spectral mask', 'Decay factor', 'Window size')),
  )

sobol_plt <- sobol_res_long %>% 
  ggplot(aes(y = vars, x = sensitivity, fill = vars)) +
  facet_grid(position ~ indice, scales = 'fixed') +
  geom_bar(stat = 'identity', width = 0.5) +
  geom_errorbar(
    aes(
      xmin = sensitivity - sensitivity_conf,
      xmax = sensitivity + sensitivity_conf
      ),
    width = 0.2
  ) +
  # Add text values for each sensitivity index value
  geom_label(
    aes(label = round(sensitivity, 2)),
    fill = 'white', color = 'black',
    nudge_x = -0.1,
    vjust = 0.5, size = 1.5
  ) +
  # Add text value for the sum of indices
  geom_label(
    aes(x = 0.55, y = 1, label = str_wrap(paste0('Sum: ', round(indice_sum, 2)), 5)),
    fill = 'white', color = 'black',
    vjust = 0.5, size = 1.5
  ) +
  scale_y_discrete(
    name = 'Variable',
    # Wrap too long labels
    labels = function(x) str_wrap(x, 10)
  ) +
  scale_x_continuous(
    name = 'Variance explained'
    # breaks = seq(0, 1, 0.05),
  ) +
  scale_fill_nejm() +
  theme_bw() + theme(
    legend.position = 'best',
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
    axis.title = element_text(face="bold", size = 4),
    axis.text = element_text(
      color="black", size = 4, margin = margin(0.0, 0.0, 0.0, 0.0, 'mm')
      ),
    # Remove x axis text and title
    axis.line = element_line(linewidth = 0.1, colour = "black"),
    plot.title = element_text(hjust = 0.5, size = 5),
    panel.grid = element_line(colour = 'grey', size = 0.1),
    panel.border = element_rect(linewidth = 0.1),
    strip.background = element_rect(linewidth = 0.01),
    strip.text = element_text(
      colour = 'black', size = 4, margin = margin(b = 0.3, t = 0.3, unit='mm')
    ),
    axis.ticks = element_line(linewidth = 0.1),
    axis.ticks.y = element_blank(),
    axis.ticks.length = unit(0.1, 'lines')
    )

sobol_plt
# Save plot
ggsave(
  filename = paste0('plots/preliminary_sensitivity_sobol_group.png'),
  plot = sobol_plt, width = 6, height = 6, units = 'cm', dpi = 360
  )

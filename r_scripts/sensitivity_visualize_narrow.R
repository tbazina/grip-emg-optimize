# Visalize sensitivity analysis results on narrowed ranges
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
library(ggh4x)
library(datapasta)

# Load data
# Position 3 and 4, preliminary latin and sobol
position = 3
lh_res <- read_csv(
  paste0('results/sensitivity_preliminary_lh_pos', position, '_narrow22.csv')
  ) %>%
  mutate(position = position) 
fast_res <- read_csv(
  paste0(
    'results/sensitivity_groupsFalse_rbd_fast_samples65536_resamples16384_pos',
    position, '_narrow_bounds.csv'
    )
  ) %>%
  mutate(position = position)
position = 4
lh_res <- lh_res %>% 
  bind_rows(
    read_csv(
      paste0('results/sensitivity_preliminary_lh_pos', position, '_narrow22.csv')
      ) %>%
      mutate(position = position)
    )
fast_res <- fast_res %>% 
  bind_rows(
    read_csv(
      paste0(
        'results/sensitivity_groupsFalse_rbd_fast_samples65536_resamples16384_pos',
        position, '_narrow_bounds.csv'
        )
      ) %>%
      mutate(position = position)
    )

# Get names of all variables
fast_res %>% pull(vars) %>% unique()

# Check sum of sensitivity index values for 4 variables with maximum values
fast_res %>% 
  group_by(position) %>%
  filter(S1 > 0.0152)
  # slice_max(S1, n = 50) %>%
  arrange((S1), .by_group = T) %>%
  summarise(
    sum_sensitivity = sum(S1)
  )

# Plot results from RBD-FAST sensitivity analysis - bar + CI
# Plot total first order S1 indices
fast_res_long <- fast_res %>% 
  group_by(position) %>%
  # filter(S1 > 0.0152) %>%
  # slice_max(S1, n = 50) %>%
  arrange(desc(S1), .by_group = T) %>%
  rename(
    S1_index = S1
  ) %>% 
  group_by(position) %>%
  mutate(
    indice_sum = sum(S1_index)
  ) %>% 
  ungroup() %>% 
  mutate(
    # Rename sensitivity indices to full names
    # indice = case_when(
    #   indice == 'ST' ~ 'Total Sobol indices',
    #   indice == 'S1' ~ 'First order Sobol indices'
    # ),
    position_vars = factor(paste0(position, '_', vars)),
    # Rename position levels to Position: x
    position = case_when(
      position == 3 ~ 'Position: 3',
      position == 4 ~ 'Position: 4'
    ),
    position = ordered(position, levels = c('Position: 3', 'Position: 4')),
    # vars = ordered(vars, levels = c('Frequencies', 'Decay factor', 'Window size')),
  )

# Get unique names for top 8 variables
vars_select <- fast_res_long %>% pull(vars) %>% unique()

fast_plt <- fast_res_long %>%
  ggplot(aes(y = vars, x = S1_index)) +#, fill = vars)) +
  facet_wrap(position ~ ., scales = 'fixed') +
  geom_bar(stat = 'identity', width = 0.5) +
  geom_errorbar(
    aes(
      xmin = S1_index - S1_conf,
      xmax = S1_index + S1_conf
      ),
    width = 0.4
  ) +
  # Add text values for each sensitivity index value
  # geom_label(
  #   aes(label = round(sensitivity, 2)),
  #   fill = 'white', color = 'black',
  #   nudge_x = -0.1,
  #   vjust = 0.5, size = 1.5
  # ) +
  # # Add text value for the sum of indices
  # geom_label(
  #   aes(x = 0.55, y = 1, label = str_wrap(paste0('Sum: ', round(indice_sum, 2)), 5)),
  #   fill = 'white', color = 'black',
  #   vjust = 0.5, size = 1.5
  # ) +
  scale_y_discrete(
    name = 'Variable',
    # Wrap too long labels
    # labels = function(x) str_wrap(x, 10)
  ) +
  scale_x_continuous(
    name = 'Variance contributed',
    # trans = scales::pseudo_log_trans(sigma = 1e-2),
    breaks = seq(0, 1, 0.02),
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
    panel.grid = element_line(colour = 'grey', linewidth = 0.1),
    panel.border = element_rect(linewidth = 0.1),
    strip.background = element_rect(linewidth = 0.01),
    strip.text = element_text(
      colour = 'black', size = 4, margin = margin(b = 0.3, t = 0.3, unit='mm')
    ),
    axis.ticks = element_line(linewidth = 0.1),
    axis.ticks.y = element_blank(),
    axis.ticks.length = unit(0.1, 'lines')
    )

fast_plt
# Save plot
ggsave(
  filename = paste0('plots/preliminary_sensitivity_rbd_fast_freqs_narrow20.png'),
  plot = fast_plt, width = 8, height = 12, units = 'cm', dpi = 360
  )


# Preliminary plots from lhs sensitivity analysis to obtain proper parameter ranges
prel_sens <- lh_res %>% 
  group_by(position) %>%
  mutate(
    mean_corr = mean(corrs)
  ) %>% 
  mutate(
    # Rename position to Position: x
    position = case_when(
      position == 3 ~ 'Pos: 3',
      position == 4 ~ 'Pos: 4'
      )
  ) %>% 
  # rename('2 Hz' = 'f2_0') %>%
  # Pivot window and decay to longer format
  pivot_longer(
    # cols = all_of(vars_select[!vars_select %in% c('decay', 'window')]),
    cols = all_of(vars_select[!vars_select %in% c(
      'f10_01', 'f22_02', 'f42_04', 'f44_04', 'f52_05', 'f54_05', 'f60_06',
      'f62_06', 'f64_06', 'f66_06', 'f68_07', 'f82_08', 'f90_09', 'f100_1',
      'f112_22', 'f86_08', 'f94_09')
      ]),
    # cols = all_of(vars_select),
    # cols = c('decay', 'window'),
    # cols = 1:31,
    # cols = 32:62,
    # cols = 63:93,
    # cols = 94:124,
    # cols = 125:155,
    # cols = 156:186,
    # cols = 187:217,
    # cols = 218:248,
    names_to = 'parameter', values_to = 'value'
    ) %>% 
  # mutate(
  #   # Rename window to Window size and decay to Decay factor
  #   parameter = case_when(
  #     parameter == 'window' ~ 'Window size',
  #     parameter == 'decay' ~ 'Decay factor'
  #     ),
  #   # Wrap parameter
  #   parameter = str_wrap(parameter, width = 10)
  # ) %>% 
  ggplot(aes(x = value, y = corrs)) +
  facet_wrap(parameter ~ position, scales = 'free', ncol = 1, strip.position = 'right') +
  # facet_grid(parameter ~ position, scales = 'fixed') +
  geom_point(aes(fill = parameter), shape = 21, stroke = 0, size = 0.2) +
  geom_line(
    stat = 'smooth', method = 'gam', se = F,
    color = 'black', alpha = 0.6, linewidth = 0.5
    ) +
  geom_line(aes(y = mean_corr), color = 'black', linetype = 'dashed', linewidth = 0.4) +
  scale_x_continuous(
    name = 'Amplification/attenuation',
    # breaks = seq(0, 5, 0.5),
    # Reduce spacing to the left and right side
    expand = expansion(mult = c(0.01, 0.01), add = c(0., 0.)),
  ) +
  scale_y_continuous(
    name = 'Peak cross-correlation',
    breaks = seq(0, 1, 0.0025),
    # breaks = seq(0, 1, 0.5),
  ) +
  scale_fill_nejm() +
  scale_color_nejm() +
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
    panel.grid = element_line(colour = 'grey', linewidth = 0.1),
    panel.border = element_rect(linewidth = 0.1),
    strip.background = element_rect(linewidth = 0.01),
    strip.text = element_text(
      colour = 'black', size = 5.0, margin = margin(b = 0.3, t = 0.3, unit='mm')
    ),
    axis.ticks = element_line(linewidth = 0.1),
    axis.ticks.y = element_blank(),
    axis.ticks.length = unit(0.1, 'lines')
  ) 
  # Add at the end to inherit theme settings
  # geom_magnify(
  #   aes(from = (
  #     parameter == 'Window\nsize' & corrs > 0.86 & value <= 495 & value >= 150
  #     )),
  #   to = list(xmin = 150., ymin = 0.5, xmax = 495, ymax = 0.75),
  #   axes = '',
  #   shadow = T,
  #   expand = c(0.0, 0.0),
  #   corners = 0.05,
  #   proj = 'corresponding',
  #   recompute = T,
  #   scale.inset = 1,
  #   # Borders and projection lines
  #   linetype = 1,
  #   linewidth = 0.2,
  # ) +
  # geom_magnify(
  #   aes(from = (
  #     parameter == 'Decay\nfactor' & value <= 0.02 & corrs > 0.86
  #     )),
  #   to = list(xmin = 0.00, ymin = 0.5, xmax = 0.02, ymax = 0.80),
  #   axes = '',
  #   shadow = T,
  #   expand = c(0.0, 0.0),
  #   corners = 0.05,
  #   proj = 'facing',
  #   recompute = T,
  #   scale.inset = 1,
  #   # Borders and projection lines
  #   linetype = 1,
  #   linewidth = 0.2,
  # ) 
prel_sens

# TODO: run narrow3 again - accidental overwrite
ggsave(
  filename = paste0('plots/preliminary_sensitivity_lh_narrow20.png'),
  plot = prel_sens, width = 8.9, height = 20, units = 'cm', dpi = 360
  )

# Plot FINAL results from RBD-FAST sensitivity analysis - bar + CI
# Plot ALL first order S1 indices
fast_res_long <- fast_res %>% 
  # Filter out decay and window from vars
  filter(vars != 'decay' & vars != 'window') %>%
  group_by(position) %>%
  arrange((S1), .by_group = T) %>%
  rename(
    S1_index = S1
  ) %>% 
  group_by(position) %>%
  mutate(
    indice_sum = sum(S1_index)
  ) %>% 
  ungroup() %>% 
  mutate(
    # Rename position levels to Position: x
    position = case_when(
      position == 3 ~ 'Position: 3',
      position == 4 ~ 'Position: 4'
    ),
    # Remove decimal part of frequency from name
    vars = str_split(vars, '_') %>% map_chr(\(x) pluck(x, 1)),
    # Remove starting f from vars name
    vars = str_remove(vars, 'f'),
    # Convert vars to ordered factor in descending order
    vars = as.integer(vars),
    # vars = ordered(vars),
    # vars = ordered(vars),
    position = ordered(position, levels = c('Position: 3', 'Position: 4')),
  )

# Statistics on final sensitivity indices
fast_res_long %>% 
  group_by(position) %>%
  arrange(desc(S1_index), .by_group = T) %>%
  slice_max(order_by = S1_index, n = 27) %>%
  # Compute cumsum of S1_index per group
  mutate(
    cum_var_contributed = cumsum(S1_index),
  ) %>%
  # filter(vars <= 208) %>%
  # filter(S1_index > 0.0152) %>%
  summarise(
    var_contributed = sum(S1_index),
  )

# Plot graph for selecting variables by cumulative total variance contributed
fast_elbow_plt <- fast_res_long %>%
  group_by(position) %>%
  arrange(desc(S1_index), .by_group = T) %>%
  mutate(
    cum_var_contributed = cumsum(S1_index),
    num_vars = row_number(),
  ) %>%
  ungroup() %>% 
  mutate(
    # Rename position levels to Position: x
    position = case_when(
      position == 'Position: 3' ~ 3,
      position == 'Position: 4' ~ 4
    ),
    position = ordered(position, levels = c(3, 4)),
    ) %>% 
  ggplot(aes(x = num_vars, y = cum_var_contributed)) +
  # facet_wrap(position ~ ., scales = 'free_y', ncol = 1, strip.position = 'right') +
  # Sensitivity index
  geom_point(aes(fill = position), shape = 21, stroke = 0, size = 0.5) +
  scale_x_continuous(
    name = 'Num. spectral components',
    breaks = seq(0, 400, 5),
    expand = expansion(mult = 0., add = 2),
    # Wrap too long labels
    # labels = function(x) str_wrap(x, 10)
  ) +
  scale_y_continuous(
    name = 'Cumulative SI',
    # trans = scales::pseudo_log_trans(sigma = 1e-2),
    breaks = seq(0, 1, 0.1),
    # limits = c(-0.001, 0.21),
    # expand = expansion(mult = 0.0, add = 0.03)
  ) +
  scale_fill_nejm() +
  guides(
    fill = guide_legend('Position: ')
    ) +
  theme_bw() + theme(
    legend.position = 'inside',
    legend.position.inside = c(0.8, 0.7),
    legend.title = element_text(
      colour="black", size = 5, margin = margin(0, 0, 0, 0, 'mm')
      ),
    legend.box = 'horizontal',
    legend.direction = 'horizontal',
    legend.box.spacing = unit(0.5, 'mm'),
    legend.box.background = element_rect(colour = 'black', fill = 'white'),
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
      color="black", size = 3, margin = margin(0.0, 0.0, 0.0, 0.0, 'mm')
      ),
    axis.text.x = element_text(
      color="black", size = 3, margin = margin(0.0, 0.0, 0.0, 0.0, 'mm'),
      angle = 45, hjust = 0.5, vjust = 0.5
      ),
    # Remove x axis text and title
    axis.line = element_line(linewidth = 0.1, colour = "black"),
    plot.title = element_text(hjust = 0.5, size = 5),
    panel.grid = element_line(colour = 'grey', linewidth = 0.1),
    panel.border = element_rect(linewidth = 0.1),
    strip.background = element_rect(linewidth = 0.01),
    strip.text = element_text(
      colour = 'black', size = 4, margin = margin(b = 0.3, t = 0.3, unit='mm')
    ),
    axis.ticks = element_line(linewidth = 0.1),
    axis.ticks.y = element_blank(),
    axis.ticks.length = unit(0.1, 'lines')
    )

fast_elbow_plt
# Save plot
ggsave(
  filename = paste0('plots/final_sensitivity_rbd_fast_elbow_plt.png'),
  plot = fast_elbow_plt, width = 7, height = 1.7, units = 'cm', dpi = 360
  )
  
  
# Plot final sensitivity indices
fast_plt <- fast_res_long %>%
  # Contribution of frequencies > 202 is < 2 %
  group_by(position) %>%
  arrange(desc(S1_index), .by_group = T) %>%
  slice_max(order_by = S1_index, n = 27) %>%
  mutate(
    # Sum all S1 indices per position
    indice_sum = sum(S1_index),
  ) %>% 
  ungroup() %>% 
  mutate(
    # Convert vars to ordered factor
    vars = factor(vars, levels = unique(sort(vars))),
  ) %>% 
  # filter(vars <= 250) %>%
  ggplot(aes(x = vars, y = S1_index)) +
  facet_wrap(position ~ ., scales = 'free_y', ncol = 1, strip.position = 'right') +
  # Sensitivity index
  geom_bar(aes(fill = position), stat = 'identity', width = 0.8) +
  # CI on sensitivity index
  geom_errorbar(
    aes(
      ymin = S1_index - S1_conf,
      ymax = S1_index + S1_conf,
      ),
    width = 0.7,
    linewidth = 0.1,
    color = 'black'
  ) +
  # Add text value for the sum of indices
  geom_label(
    aes(x = 4 , y = 0.09, label = str_wrap(paste0('Sum: ', round(indice_sum, 2)), 5)),
    fill = 'white', color = 'black',
    vjust = 0.5, size = 1.5
  ) +
  scale_x_discrete(
    name = 'Spectral component [Hz]',
    # breaks = seq(0, 500, 4),
    # expand = expansion(mult = 0., add = 2),
    # Wrap too long labels
    # labels = function(x) str_wrap(x, 10)
  ) +
  scale_y_continuous(
    name = 'Variance contributed',
    trans = scales::pseudo_log_trans(sigma = 1e-2),
    breaks = seq(0, 200, 20) * 1e-3,
    # limits = c(-0.001, 0.21),
    expand = expansion(mult = 0.0, add = 0.03)
  ) +
  scale_fill_nejm() +
  scale_color_nejm() +
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
    panel.spacing.y = unit(0.5, 'mm'),
    panel.spacing.x = unit(0, 'mm'),
    axis.title = element_text(face="bold", size = 4),
    axis.text = element_text(
      color="black", size = 3, margin = margin(0.0, 0.0, 0.0, 0.0, 'mm')
      ),
    axis.text.x = element_text(
      color="black", size = 3, margin = margin(0.0, 0.0, 0.0, 0.0, 'mm'),
      angle = 45, hjust = 0.5, vjust = 0.5
      ),
    # Remove x axis text and title
    axis.line = element_line(linewidth = 0.1, colour = "black"),
    plot.title = element_text(hjust = 0.5, size = 5),
    panel.grid = element_line(colour = 'grey', linewidth = 0.1),
    panel.border = element_rect(linewidth = 0.1),
    strip.background = element_rect(linewidth = 0.01),
    strip.text = element_text(
      colour = 'black', size = 4, margin = margin(b = 0.3, t = 0.3, unit='mm')
    ),
    axis.ticks = element_line(linewidth = 0.1),
    axis.ticks.length = unit(0.1, 'lines')
    )

fast_plt
# Save plot
ggsave(
  filename = paste0('plots/final_sensitivity_rbd_fast_freqs.png'),
  plot = fast_plt, width = 7, height = 4, units = 'cm', dpi = 360
  )

# Manual bounds on frequency mask
bounds <- tibble::tribble(
~Variable, ~Trend,      ~Bounds,
  "decay",    "↘", "0 - 0.0005",
 "window",   "↗↘",  "275 - 330",
   "2 Hz",   "↗↘",    "0 - 0.1",
   "4 Hz",   "↗↘",    "0 - 0.2",
   "6 Hz",   "↗↘",    "0 - 0.4",
   "8 Hz",   "↗↘",    "0 - 0.5",
  "10 Hz",   "↗↘",      "0 - 1",
  "12 Hz",   "↗↘",      "0 - 1",
  "14 Hz",   "↗↘",      "0 - 1",
  "16 Hz",   "↗↘",    "0 - 1.5",
  "18 Hz",   "↗↘",      "0 - 2",
  "20 Hz",   "↗↘",    "0 - 2.5",
  "22 Hz",   "↗↘",    "0.5 - 3",
  "24 Hz",   "↗↘",  "0.5 - 2.5",
  "26 Hz",   "↗↘",    "0.5 - 3",
  "28 Hz",   "↗↘",      "1 - 3",
  "30 Hz",   "↗↘",      "1 - 3",
  "32 Hz",   "↗↘",  "1.5 - 3.5",
  "34 Hz",   "↗↘",  "1.5 - 3.5",
  "36 Hz",   "↗↘",  "1.5 - 3.5",
  "38 Hz",   "↗↘",  "1.5 - 3.5",
  "40 Hz",   "↗↘",  "1.5 - 3.5",
  "42 Hz",   "↗↘",    "1.5 - 3",
  "44 Hz",   "↗↘",    "1.5 - 3",
  "46 Hz",   "↗↘",      "1 - 2",
  "48 Hz",   "↗↘",    "0.5 - 2",
  "50 Hz",   "↗↘",   "0 - 0.75",
  "52 Hz",   "↗↘",  "0.5 - 2.5",
  "54 Hz",   "↗↘",    "1 - 3.5",
  "56 Hz",   "↗↘",    "1 - 2.5",
  "58 Hz",   "↗↘",    "1.5 - 3",
  "60 Hz",   "↗↘",    "1.5 - 3",
  "62 Hz",   "↗↘",  "2.5 - 4.5",
  "64 Hz",   "↗↘",    "1.5 - 4",
  "66 Hz",   "↗↘",    "2 - 4.5",
  "68 Hz",   "↗↘",    "2 - 4.5",
  "70 Hz",   "↗↘",  "1.5 - 3.5",
  "72 Hz",   "↗↘",    "2.5 - 4",
  "74 Hz",   "↗↘",    "3.5 - 5",
  "76 Hz",   "↗↘",    "2 - 4.5",
  "78 Hz",   "↗↘",    "2.5 - 4",
  "80 Hz",   "↗↘",    "1.5 - 4",
  "82 Hz",   "↗↘",    "2.5 - 5",
  "84 Hz",   "↗↘",  "2.5 - 4.5",
  "86 Hz",   "↗↘",      "2 - 5",
  "88 Hz",   "↗↘",  "2.5 - 4.5",
  "90 Hz",   "↗↘",    "3.5 - 6",
  "92 Hz",   "↗↘",      "3 - 5",
  "94 Hz",   "↗↘",  "3.5 - 6.5",
  "96 Hz",   "↗↘",  "3.5 - 5.5",
  "98 Hz",   "↗↘",  "2.5 - 4.5",
 "100 Hz",   "↗↘",    "2.5 - 6",
 "102 Hz",   "↗↘",    "3.5 - 5",
 "104 Hz",   "↗↘",  "3.5 - 5.5",
 "106 Hz",   "↗↘",      "3 - 6",
 "108 Hz",   "↗↘",    "3 - 5.5",
 "110 Hz",   "↗↘",  "4.5 - 6.5",
 "112 Hz",   "↗↘",  "4.5 - 6.5",
 "114 Hz",   "↗↘",  "4.5 - 6.5",
 "116 Hz",   "↗↘",  "4.5 - 6.5",
 "118 Hz",   "↗↘",  "4.5 - 6.5",
 "120 Hz",   "↗↘",  "4.5 - 6.5",
 "122 Hz",   "↗↘",      "3 - 6",
 "124 Hz",   "↗↘",    "4 - 6.5",
 "126 Hz",   "↗↘",    "4 - 6.5",
 "128 Hz",   "↗↘",    "4 - 6.5",
 "130 Hz",   "↗↘",    "4 - 6.5",
 "132 Hz",   "↗↘",    "4 - 6.5",
 "134 Hz",   "↗↘",    "4 - 6.5",
 "136 Hz",   "↗↘",    "4 - 6.5",
 "138 Hz",   "↗↘",    "4 - 6.5",
 "140 Hz",   "↗↘",    "4 - 6.5",
 "142 Hz",   "↗↘",    "4 - 6.5",
 "144 Hz",   "↗↘",    "4 - 6.5",
 "146 Hz",   "↗↘",    "4 - 6.5",
 "148 Hz",   "↗↘",    "4 - 6.5",
 "150 Hz",   "↗↘",    "4 - 6.5",
 "152 Hz",   "↗↘",    "4 - 6.5",
 "154 Hz",   "↗↘",    "4 - 6.5",
 "156 Hz",   "↗↘",    "4 - 6.5",
 "158 Hz",   "↗↘",    "4 - 6.5",
 "160 Hz",   "↗↘",    "4 - 6.5",
 "162 Hz",   "↗↘",    "4 - 6.5",
 "164 Hz",   "↗↘",    "4 - 6.5",
 "166 Hz",   "↗↘",    "4 - 6.5",
 "168 Hz",   "↗↘",    "4 - 6.5",
 "170 Hz",   "↗↘",    "4 - 6.5",
 "172 Hz",   "↗↘",    "4 - 6.5",
 "174 Hz",   "↗↘",    "4 - 6.5",
 "176 Hz",   "↗↘",    "4 - 6.5",
 "178 Hz",   "↗↘",    "4 - 6.5",
 "180 Hz",   "↗↘",    "4 - 6.5",
 "182 Hz",   "↗↘",    "4 - 6.5",
 "184 Hz",   "↗↘",    "4 - 6.5",
 "186 Hz",   "↗↘",    "4 - 6.5",
 "188 Hz",   "↗↘",    "4 - 6.5",
 "190 Hz",   "↗↘",    "4 - 6.5",
 "192 Hz",   "↗↘",    "4 - 6.5",
 "194 Hz",   "↗↘",    "4 - 6.5",
 "196 Hz",   "↗↘",    "4 - 6.5",
 "198 Hz",   "↗↘",    "4 - 6.5",
 "200 Hz",   "↗↘",    "4 - 6.5",
 "202 Hz",   "↗↘",    "4 - 6.5"
)

bounds_freq <- bounds %>% 
  # Filter decay and window
  filter(Variable != 'decay' & Variable != 'window') %>% 
  rename(
    vars = Variable, trend = Trend, bounds = Bounds
  ) %>% 
  mutate(
    # Remove Hz at the end of Variable
    vars = gsub(' Hz', '', vars) %>% as.integer(),
  ) %>% 
  # Split bounds to upper and lower using - as separator
  separate_wider_delim(
    bounds, delim = ' - ', names = c('lower_bound', 'upper_bound')
  ) %>% 
  # Convert upper_bound and lower_bound to numeric
  mutate(
    lower_bound = as.numeric(lower_bound),
    upper_bound = as.numeric(upper_bound)
  )

# Plot bounds on variables, sensitivity indices are computed for
fast_plt <- bounds_freq %>%
  filter(vars <= 202) %>% 
  # Add mean bounds
  mutate(
    mean_bound = (upper_bound - lower_bound) / 2 + lower_bound,
    # Add both position vars since bounds are common for both
    position = 'Position: 3, 4'
  ) %>%
  ggplot(aes(x = vars, y = mean_bound)) +
  facet_wrap(position ~ ., scales = 'free_y', ncol = 1, strip.position = 'right') +
  geom_line(
    stat = 'smooth',
    # method = 'gam',
    method = 'loess',
    formula = y ~ x,
    span = 0.1,
    # formula = y ~ s(x, bs = "cs"),
    color = pal_nejm()(8)[3],
    linewidth = 0.5,
  ) +
  geom_errorbar(
    aes(
      ymin = lower_bound,
      ymax = upper_bound
      ),
    width = 1.5,
    linewidth = 0.1,
    color = 'black'
  ) +
  # Add text value for the sum of indices
  scale_x_continuous(
    name = 'Spectral component [Hz]',
    breaks = seq(0, 500, 4),
    expand = expansion(mult = 0., add = 1),
    # Wrap too long labels
    # labels = function(x) str_wrap(x, 10)
  ) +
  scale_y_continuous(
    name = 'Bounds',
    breaks = seq(0, 8, 1),
    limits = c(0, 6.5),
    expand = expansion(mult = 0.0, add = 0.3)
  ) +
  scale_fill_nejm() +
  scale_color_nejm() +
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
      color="black", size = 3, margin = margin(0.0, 0.0, 0.0, 0.0, 'mm')
      ),
    axis.text.x = element_text(
      color="black", size = 3, margin = margin(0.0, 0.0, 0.0, 0.0, 'mm'),
      angle = 45, hjust = 0.5, vjust = 0.5
      ),
    # Remove x axis text and title
    axis.line = element_line(linewidth = 0.1, colour = "black"),
    plot.title = element_text(hjust = 0.5, size = 5),
    panel.grid = element_line(colour = 'grey', linewidth = 0.1),
    panel.border = element_rect(linewidth = 0.1),
    strip.background = element_rect(linewidth = 0.01),
    strip.text = element_text(
      colour = 'black', size = 4, margin = margin(b = 0.3, t = 0.3, unit='mm')
    ),
    axis.ticks = element_line(linewidth = 0.1),
    axis.ticks.length = unit(0.1, 'lines')
    )

fast_plt
# Save plot
ggsave(
  filename = paste0('plots/final_sensitivity_rbd_fast_freqs_bounds.png'),
  plot = fast_plt, width = 7, height = 1.7, units = 'cm', dpi = 360
  )

lh_res_final <- lh_res %>% 
  # Filter out decay and window from vars
  select(-c(window, decay)) %>%
  pivot_longer(
    cols = starts_with('f'),
    names_to = 'vars',
    values_to = 'mask'
  ) %>%
  mutate(
    # Rename position levels to Position: x
    position = case_when(
      position == 3 ~ 'Pos: 3',
      position == 4 ~ 'Pos: 4'
    ),
    # Remove decimal part of frequency from name
    vars = str_split(vars, '_') %>% map_chr(\(x) pluck(x, 1)),
    # Remove starting f from vars name
    vars = str_remove(vars, 'f'),
    # Convert vars to ordered factor in descending order
    vars = as.integer(vars),
    # vars = ordered(vars),
    # vars = ordered(vars),
    position = ordered(position, levels = c('Pos: 3', 'Pos: 4')),
  )

# Final sensitivity plots from lhs sensitivity analysis to obtain proper parameter ranges
final_sens <- lh_res_final %>% 
  group_by(position, vars) %>%
  mutate(
    mean_corr = mean(corrs)
  ) %>% 
  filter(vars <= 96) %>%
  # filter(vars > 96 & vars <= 202) %>%
  # Pivot window and decay to longer format
  ggplot(aes(x = mask, y = corrs, group = position)) +
  facet_wrap(
    vars ~ ., scales = 'free_x', strip.position = 'top', ncol = 8,
    labeller = as_labeller(\(x) paste0(x, ' Hz'))
    ) +
  # facet_grid(vars ~ position, scales = 'free') +
  geom_point(aes(fill = position), shape = 21, stroke = 0, size = 0.2) +
  geom_line(
    stat = 'smooth', method = 'gam', se = F,
    color = 'black', alpha = 0.6, linewidth = 0.5
    ) +
  geom_line(aes(y = mean_corr), color = 'black', linetype = 'dashed', linewidth = 0.4) +
  scale_x_continuous(
    name = 'Amplification/attenuation',
    # breaks = seq(0, 5, 0.5),
    # Reduce spacing to the left and right side
    expand = expansion(mult = c(0.01, 0.01), add = c(0., 0.)),
  ) +
  scale_y_continuous(
    name = 'Peak cross-correlation',
    breaks = seq(0, 1, 0.0025),
    # breaks = seq(0, 1, 0.5),
  ) +
  scale_fill_nejm() +
  scale_color_nejm() +
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
    panel.grid = element_line(colour = 'grey', linewidth = 0.1),
    panel.border = element_rect(linewidth = 0.1),
    strip.background = element_rect(linewidth = 0.01),
    strip.text = element_text(
      colour = 'black', size = 5.0, margin = margin(b = 0.3, t = 0.3, unit='mm')
    ),
    axis.ticks = element_line(linewidth = 0.1),
    axis.ticks.y = element_blank(),
    axis.ticks.length = unit(0.1, 'lines')
  ) 

final_sens

ggsave(
  filename = paste0('plots/final_sensitivity_lh_narrow96_2.png'),
  plot = final_sens, width = 29.7, height = 17.5, units = 'cm', dpi = 360
  )
ggsave(
  filename = paste0('plots/final_sensitivity_lh_narrow202_2.png'),
  plot = final_sens, width = 29.7, height = 17.5, units = 'cm', dpi = 360
  )

# Descriptive statistics per postion on lh_res
lh_res %>% 
  group_by(position) %>%
  summarise(
    count = n(),
    mean_corr = mean(corrs),
    sd_corr = sd(corrs),
    min_corr = min(corrs),
    max_corr = max(corrs)
  )


# Evolution of bounds across all sensitivity analysis steps
bounds_steps <- tibble::tribble(
  ~Variable, ~Trend,     ~Bounds, ~Trend,     ~Bounds, ~Trend,     ~Bounds, ~Trend,     ~Bounds, ~Trend,     ~Bounds, ~Trend,     ~Bounds, ~Trend,     ~Bounds, ~Trend,     ~Bounds, ~Trend,     ~Bounds, ~Trend,      ~Bounds, ~Trend,      ~Bounds, ~Trend,     ~Bounds, ~Trend,     ~Bounds, ~Trend,   ~Bounds, ~Trend,     ~Bounds, ~Trend,      ~Bounds, ~Trend,     ~Bounds, ~Trend,     ~Bounds, ~Trend,      ~Bounds, ~Trend,      ~Bounds,
    "decay",    "↘", "0 - 0.005",    "↘", "0 - 0.003",    "↘", "0 - 0.001",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,    "↘", "0 - 0.0005",     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↘", "0 - 0.0005",    "↘", "0 - 0.0005",
   "window",   "↗↘", "200 - 400",   "↗↘", "250 - 350",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,   "↗↘",  "275 - 330",     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,   "↗↘",  "275 - 330",   "↗↘",  "275 - 330",
     "2 Hz",     NA,          NA,    "↘",   "0 - 0.3",     NA,          NA,     NA,          NA,    "↘",   "0 - 0.2",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,    "↘",   "0 - 0.1",     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↘",    "0 - 0.1",   "↗↘",    "0 - 0.1",
     "4 Hz",     NA,          NA,     NA,          NA,    "↘",   "0 - 0.5",     NA,          NA,     NA,          NA,     NA,          NA,    "↘",   "0 - 0.2",     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↘",    "0 - 0.2",   "↗↘",    "0 - 0.2",
     "6 Hz",    "↘",  "0 - 1.25",     NA,          NA,     NA,          NA,    "↘",  "0 - 0.75",     NA,          NA,     NA,          NA,     NA,          NA,    "↘",   "0 - 0.4",     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↘",    "0 - 0.4",   "↗↘",    "0 - 0.4",
     "8 Hz",    "↘",   "0 - 2.5",     NA,          NA,    "↘",     "0 - 1",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,    "↘",    "0 - 0.5",     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↘",    "0 - 0.5",   "↗↘",    "0 - 0.5",
    "10 Hz",    "↘",   "0 - 2.5",     NA,          NA,     NA,          NA,    "↘",     "0 - 1",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,    "↘",     "0 - 1",    "↘",      "0 - 1",   "↗↘",      "0 - 1",
    "12 Hz",     NA,          NA,    "↘",     "0 - 3",     NA,          NA,    "↘",     "0 - 1",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↘",      "0 - 1",   "↗↘",      "0 - 1",
    "14 Hz",     NA,          NA,   "↗↘",   "0 - 3.5",     NA,          NA,    "↘",     "0 - 1",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,   "↗↘",      "0 - 1",   "↗↘",      "0 - 1",
    "16 Hz",     NA,          NA,     NA,          NA,    "↘",     "0 - 3",     NA,          NA,     NA,          NA,    "↘",     "0 - 2",     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,   "↗↘",    "0 - 1.5",     NA,          NA,     NA,          NA,   "↗↘",    "0 - 1.5",   "↗↘",    "0 - 1.5",
    "18 Hz",     NA,          NA,     NA,          NA,     NA,          NA,    "↘",     "0 - 2",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,   "↗↘",      "0 - 2",   "↗↘",      "0 - 2",
    "20 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,   "↗↘",     "0 - 3",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,   "↗↘", "0 - 2.5",     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,   "↗↘",    "0 - 2.5",   "↗↘",    "0 - 2.5",
    "22 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,   "↗↘", "0.5 - 3.5",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,   "↗↘",    "0.5 - 3",     NA,          NA,     NA,          NA,     NA,        NA,   "↗↘",   "0.5 - 3",     NA,           NA,     NA,          NA,     NA,          NA,   "↗↘",    "0.5 - 3",   "↗↘",    "0.5 - 3",
    "24 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,   "↗↘",   "0.5 - 3",     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,   "↗↘", "0.5 - 2.5",     NA,          NA,   "↗↘",  "0.5 - 2.5",   "↗↘",  "0.5 - 2.5",
    "26 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,   "↗↘", "0.5 - 3.5",     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,   "↗↘",    "0.5 - 3",   "↗↘",    "0.5 - 3",
    "28 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,   "↗↘",  "0.5 - 3.5",     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,   "↗↘",      "1 - 3",   "↗↘",      "1 - 3",
    "30 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,   "↗↘", "0.5 - 3.5",     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,   "↗↘",      "1 - 3",   "↗↘",      "1 - 3",
    "32 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,   "↗↘", "0.5 - 3.5",     NA,           NA,     NA,           NA,    "↗",   "1.5 - 4",     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,   "↗↘",  "1.5 - 3.5",   "↗↘",  "1.5 - 3.5",
    "34 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,   "↗↘", "0.5 - 3.5",     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,   "↗↘",    "1 - 3.5",     NA,          NA,     NA,          NA,   "↗↘",  "1.5 - 3.5",   "↗↘",  "1.5 - 3.5",
    "36 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,   "↗↘", "0.5 - 3.5",     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,   "↗↘",     "1 - 3",   "↗↘",  "1.5 - 3.5",   "↗↘",  "1.5 - 3.5",
    "38 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,   "↗↘", "1.5 - 4.5",     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,   "↗↘", "1.5 - 4",     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,   "↗↘",  "1.5 - 3.5",   "↗↘",  "1.5 - 3.5",
    "40 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,   "↗↘", "0.5 - 3.5",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,   "↗↘",   "1 - 3",     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,   "↗↘",  "1.5 - 3.5",   "↗↘",  "1.5 - 3.5",
    "42 Hz",     NA,          NA,     NA,          NA,   "↗↘",     "0 - 5",   "↗↘", "0.5 - 3.5",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,   "↗↘", "0.5 - 3.5",     NA,           NA,     NA,          NA,     NA,          NA,   "↗↘",    "1 - 3.5",   "↗↘",    "1.5 - 3",
    "44 Hz",     NA,          NA,     NA,          NA,     NA,          NA,   "↗↘", "0.5 - 3.5",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,   "↗↘", "0.5 - 3.5",     NA,           NA,     NA,          NA,     NA,          NA,   "↗↘",    "1 - 3.5",   "↗↘",    "1.5 - 3",
    "46 Hz",    "↘",   "0 - 4.5",     NA,          NA,     NA,          NA,     NA,          NA,   "↗↘", "0.5 - 3.5",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,   "↗↘", "0.75 - 2.5",     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,   "↗↘", "0.75 - 2.5",   "↗↘",      "1 - 2",
    "48 Hz",     NA,          NA,    "↘",     "0 - 3",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,   "↗↘",  "0.5 - 2.5",     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,   "↗↘", "0.5 - 2.25",     NA,          NA,     NA,          NA,   "↗↘", "0.5 - 2.25",   "↗↘",    "0.5 - 2",
    "50 Hz",    "↘",   "0 - 2.5",     NA,          NA,    "↘",   "0 - 1.5",     NA,          NA,     NA,          NA,     NA,          NA,    "↘",  "0 - 0.75",     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,   "↗↘",   "0 - 0.75",   "↗↘",   "0 - 0.75",
    "52 Hz",     NA,          NA,     NA,          NA,   "↗↘", "0.5 - 3.5",     NA,          NA,     NA,          NA,   "↗↘", "0.5 - 2.5",     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,   "↗↘", "0.5 - 2.5",     NA,           NA,     NA,          NA,     NA,          NA,   "↗↘",  "0.5 - 2.5",   "↗↘",  "0.5 - 2.5",
    "54 Hz",     NA,          NA,     NA,          NA,     NA,          NA,   "↗↘", "0.5 - 3.5",     NA,          NA,     NA,          NA,   "↗↘",     "1 - 4",   "↗↘",     "1 - 4",   "↗↘",     "1 - 4",     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,   "↗↘",      "1 - 4",   "↗↘",    "1 - 3.5",
    "56 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,   "↗↘", "0.5 - 3.5",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,   "↗↘",   "1 - 3",     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,   "↗↘",   "1 - 2.75",   "↗↘",    "1 - 2.5",
    "58 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,   "↗↘",     "1 - 4",     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,   "↗↘", "1.5 - 3.5",     NA,           NA,     NA,          NA,     NA,          NA,   "↗↘",  "1.5 - 3.5",   "↗↘",    "1.5 - 3",
    "60 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,   "↗↘",     "1 - 4",     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,   "↗↘",     "1 - 4",     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,   "↗↘",  "1.5 - 3.5",   "↗↘",    "1.5 - 3",
    "62 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,    "↗",   "2 - 5.5",   "↗↘", "2.5 - 4.5",     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,   "↗↘", "2.5 - 4.5",     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,   "↗↘",  "2.5 - 4.5",   "↗↘",  "2.5 - 4.5",
    "64 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,   "↗↘",   "1 - 4.5",     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,   "↗↘",   "1 - 4.5",     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,   "↗↘",  "1.5 - 4.5",   "↗↘",    "1.5 - 4",
    "66 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,   "↗↘",   "1 - 4.5",     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,   "↗↘",   "1 - 4.5",     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,   "↗↘",  "1.5 - 4.5",   "↗↘",    "2 - 4.5",
    "68 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,    "↗",   "2 - 5.5",     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,   "↗↘",     "2 - 5",     NA,          NA,     NA,        NA,   "↗↘",     "2 - 5",     NA,           NA,     NA,          NA,     NA,          NA,   "↗↘",    "2 - 4.5",   "↗↘",    "2 - 4.5",
    "70 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,   "↗↘",      "1 - 4",     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,   "↗↘", "1.5 - 3.5",     NA,          NA,   "↗↘",  "1.5 - 3.5",   "↗↘",  "1.5 - 3.5",
    "72 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,   "↗↘", "1.5 - 4.5",     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,   "↗↘",      "2 - 4",   "↗↘",    "2.5 - 4",
    "74 Hz",     NA,          NA,    "↗",     "1 - 5",     NA,          NA,     NA,          NA,    "↗",     "3 - 7",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,   "↗↘",  "3.5 - 5.5",     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,   "↗↘",    "3.5 - 5",   "↗↘",    "3.5 - 5",
    "76 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,    "↗", "1.5 - 5.5",     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,   "↗↘",   "2 - 5",     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,   "↗↘",    "2 - 4.5",   "↗↘",    "2 - 4.5",
    "78 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,   "↗↘", "1.5 - 4.5",     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,   "↗↘", "2.5 - 4.5",   "↗↘",  "2.5 - 4.5",   "↗↘",    "2.5 - 4",
    "80 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,   "↗↘", "1.5 - 4.5",     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,   "↗↘",    "1.5 - 4",   "↗↘",    "1.5 - 4",
    "82 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,    "↗",   "2 - 5.5",     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,   "↗↘",    "2 - 5.5",     NA,          NA,     NA,          NA,   "↗↘",  "2.5 - 5.5",   "↗↘",    "2.5 - 5",
    "84 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,   "↗↘",      "2 - 5",     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,   "↗↘",  "2.5 - 4.5",   "↗↘",  "2.5 - 4.5",
    "86 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,    "↗",   "2 - 5.5",     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,   "↗↘",     "2 - 5",   "↗↘",     "2 - 5",   "↗↘",      "2 - 5",   "↗↘",      "2 - 5",
    "88 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,   "↗↘",      "2 - 5",     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,   "↗↘",    "2 - 4.5",   "↗↘",  "2.5 - 4.5",
    "90 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,    "↗",     "3 - 7",     NA,          NA,     NA,          NA,     NA,          NA,   "↗↘",      "3 - 6",     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,   "↗↘",     "3 - 6",     NA,           NA,     NA,          NA,     NA,          NA,   "↗↘",      "3 - 6",   "↗↘",    "3.5 - 6",
    "92 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,    "↗",     "3 - 7",    "↘",   "3 - 5.5",     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,   "↗↘",   "3 - 5",     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,   "↗↘",      "3 - 5",   "↗↘",      "3 - 5",
    "94 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,    "↗",     "3 - 7",     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,   "↗↘",   "3 - 6.5",   "↗↘",   "3 - 6.5",   "↗↘",  "3.5 - 6.5",   "↗↘",  "3.5 - 6.5",
    "96 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,    "↗",   "2 - 5.5",     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",      "3 - 6",   "↗↘",  "3.5 - 5.5",
    "98 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,    "↗", "2.5 - 5.5",     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,   "↗↘",  "2.5 - 4.5",   "↗↘",  "2.5 - 4.5",
   "100 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,    "↗",    "2 - 5.5",     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,    "↗",   "2 - 5.5",   "↗↘",      "2 - 6",   "↗↘",    "2.5 - 6",
   "102 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,    "↗",   "3 - 5.5",     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,   "↗↘",    "3 - 5.5",   "↗↘",    "3.5 - 5",
   "104 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,    "↗", "2.5 - 5.5",     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",      "3 - 6",   "↗↘",  "3.5 - 5.5",
   "106 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,    "↗", "2.5 - 5.5",     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",      "3 - 6",   "↗↘",      "3 - 6",
   "108 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,    "↗", "2.5 - 5.5",     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",      "3 - 6",   "↗↘",    "3 - 5.5",
   "110 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,    "↗",  "2.5 - 5.5",     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,    "↗",   "3.5 - 6",     NA,          NA,    "↗",  "4.5 - 6.5",   "↗↘",  "4.5 - 6.5",
   "112 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,    "↗",    "3 - 5.5",     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,    "↗",   "3.5 - 6",    "↗",   "3.5 - 6",    "↗",  "4.5 - 6.5",   "↗↘",  "4.5 - 6.5",
   "114 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,    "↗",    "3 - 5.5",     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,    "↗",      "4 - 6",     NA,          NA,     NA,          NA,    "↗",  "4.5 - 6.5",   "↗↘",  "4.5 - 6.5",
   "116 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,    "↗", "2.5 - 5.5",     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",    "4 - 6.5",   "↗↘",  "4.5 - 6.5",
   "118 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,    "↗",   "3 - 5.5",     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",    "4 - 6.5",   "↗↘",  "4.5 - 6.5",
   "120 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗", "3 - 5.5",     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",    "4 - 6.5",   "↗↘",  "4.5 - 6.5",
   "122 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",    "3 - 6.5",   "↗↘",      "3 - 6",
   "124 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,    "↗", "2.5 - 5.5",     NA,        NA,     NA,          NA,     NA,           NA,    "↗",   "3.5 - 6",     NA,          NA,    "↗",    "4 - 6.5",   "↗↘",    "4 - 6.5",
   "126 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,    "↗",   "3.5 - 6",     NA,          NA,    "↗",    "4 - 6.5",   "↗↘",    "4 - 6.5",
   "128 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,    "↗",   "3 - 5.5",     NA,           NA,     NA,          NA,     NA,          NA,    "↗",    "4 - 6.5",   "↗↘",    "4 - 6.5",
   "130 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,    "↗",      "3 - 6",     NA,          NA,     NA,          NA,    "↗",    "4 - 6.5",   "↗↘",    "4 - 6.5",
   "132 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗", "3 - 5.5",     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",    "4 - 6.5",   "↗↘",    "4 - 6.5",
   "134 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,    "↗",      "3 - 6",     NA,          NA,     NA,          NA,    "↗",    "4 - 6.5",   "↗↘",    "4 - 6.5",
   "136 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,    "↗",      "3 - 6",     NA,          NA,     NA,          NA,    "↗",    "4 - 6.5",   "↗↘",    "4 - 6.5",
   "138 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",    "3 - 6.5",   "↗↘",    "4 - 6.5",
   "140 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,    "↗",     "3 - 6",    "↗",    "4 - 6.5",   "↗↘",    "4 - 6.5",
   "142 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",    "4 - 6.5",   "↗↘",    "4 - 6.5",
   "144 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",    "4 - 6.5",   "↗↘",    "4 - 6.5",
   "146 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",    "4 - 6.5",   "↗↘",    "4 - 6.5",
   "148 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",    "4 - 6.5",   "↗↘",    "4 - 6.5",
   "150 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",    "4 - 6.5",   "↗↘",    "4 - 6.5",
   "152 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",    "4 - 6.5",   "↗↘",    "4 - 6.5",
   "154 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",    "4 - 6.5",   "↗↘",    "4 - 6.5",
   "156 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",    "4 - 6.5",   "↗↘",    "4 - 6.5",
   "158 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",    "4 - 6.5",   "↗↘",    "4 - 6.5",
   "160 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",    "4 - 6.5",   "↗↘",    "4 - 6.5",
   "162 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",    "4 - 6.5",   "↗↘",    "4 - 6.5",
   "164 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",    "4 - 6.5",   "↗↘",    "4 - 6.5",
   "166 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",    "4 - 6.5",   "↗↘",    "4 - 6.5",
   "168 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",    "4 - 6.5",   "↗↘",    "4 - 6.5",
   "170 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",    "4 - 6.5",   "↗↘",    "4 - 6.5",
   "172 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",    "4 - 6.5",   "↗↘",    "4 - 6.5",
   "174 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",    "4 - 6.5",   "↗↘",    "4 - 6.5",
   "176 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",    "3 - 6.5",   "↗↘",    "4 - 6.5",
   "178 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",    "4 - 6.5",   "↗↘",    "4 - 6.5",
   "180 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",    "4 - 6.5",   "↗↘",    "4 - 6.5",
   "182 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",    "4 - 6.5",   "↗↘",    "4 - 6.5",
   "184 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",    "3 - 6.5",   "↗↘",    "4 - 6.5",
   "186 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",    "4 - 6.5",   "↗↘",    "4 - 6.5",
   "188 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",    "4 - 6.5",   "↗↘",    "4 - 6.5",
   "190 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",    "4 - 6.5",   "↗↘",    "4 - 6.5",
   "192 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",    "4 - 6.5",   "↗↘",    "4 - 6.5",
   "194 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",    "4 - 6.5",   "↗↘",    "4 - 6.5",
   "196 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",    "4 - 6.5",   "↗↘",    "4 - 6.5",
   "198 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",    "4 - 6.5",   "↗↘",    "4 - 6.5",
   "200 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",    "4 - 6.5",   "↗↘",    "4 - 6.5",
   "202 Hz",     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,          NA,     NA,           NA,     NA,           NA,     NA,          NA,     NA,          NA,     NA,        NA,     NA,          NA,     NA,           NA,     NA,          NA,     NA,          NA,    "↗",    "4 - 6.5",   "↗↘",    "4 - 6.5"
  ) %>% 
  as_tibble(.name_repair = 'unique') %>% 
  select(-starts_with("Trend")) %>% 
  # Rename columns starting with "Bounds..." to Step
  rename_with(~str_replace(., "Bounds...", "bounds_")) %>% 
  # Reduce number in column names using formula num//2 + 1
  rename_with(
    .fn = function(name) {
      num <- str_split(name, "_", simplify = T) %>% 
        # Extract second column
        extract(, 2) %>% 
        as.numeric() %>%
        {paste0("bounds_", . %/% 2 + 1)}
      print(num)
      return(num)
    },
    starts_with("bounds")) %>%
  rename(vars = Variable) %>%
  # Pivot longer on bounds
  pivot_longer(
    cols = starts_with("bounds"),
    names_to = "step",
    values_to = "bounds"
  ) %>% 
  drop_na() %>% 
  # Rename bounds_x to just x and convert to number
  mutate(
    step = str_remove(step, "bounds_") %>% as.integer(),
  ) %>% 
  # Split bounds to upper and lower using - as separator
  separate_wider_delim(
    bounds, delim = ' - ',
    names = c('low', 'up'),
    names_sep = '_'
  ) %>% 
  # Convert upper_bound and lower_bound to numeric
  mutate(across(starts_with('bounds'), as.numeric)) %>% 
  # Add mean bounds
  mutate(
    mean_bound = (bounds_up - bounds_low) / 2 + bounds_low,
    # Add both position vars since bounds are common for both
    position = 'Position: 3, 4'
  ) 

# Panel plot of bounds across all sensitivity analysis steps
# Counter for plot breaks - only change last plot and only while saving plot
# Run code from here since it changes a global variable
plt_breaks_count <- 1
breaks_fun <- function(x) {
  plt_breaks_count <<- plt_breaks_count + 1L
  if (plt_breaks_count < 80) {
    return(seq(0, 8, 2))
    }
  else {return(seq(0, 8, 0.5))}
}

bounds_steps_plt <- bounds_steps %>% 
  # Filter decay and window
  filter(vars != 'decay' & vars != 'window') %>% 
  mutate(
    # Remove Hz at the end of Variable
    vars = gsub(' Hz', '', vars) %>% as.integer(),
  ) %>% 
  ggplot(aes(x = vars, y = mean_bound)) +
  facet_grid(. ~ step, scales = 'fixed') +
  force_panelsizes(cols = c(rep(5, 19), 30)) +
  # Bounds
  geom_errorbar(
    aes(
      ymin = bounds_low,
      ymax = bounds_up
      ),
    width = 1.5,
    linewidth = 0.2,
    color = 'black'
  ) +
  # Smooth line
  geom_line(
    data = \(x) x %>% filter(step == 21),
    stat = 'smooth',
    method = 'loess',
    formula = y ~ x,
    span = 0.1,
    color = pal_nejm()(8)[1],
    linewidth = 0.5,
    alpha = 0.8
  ) +
  coord_flip() +
  scale_x_continuous(
    name = 'Spectral component [Hz]',
    breaks = seq(0, 500, 8),
    expand = expansion(mult = 0., add = 1),
    transform = 'reverse',
    sec.axis = sec_axis(
      ~ ., 
      # name = 'Spectral component [Hz]',
      breaks = seq(0, 500, 8)
      )
    # Wrap too long labels
    # labels = function(x) str_wrap(x, 10)
  ) +
  scale_y_continuous(
    name = 'Bounds',
    # breaks = seq(0, 8, 1),
    breaks = breaks_fun,
    limits = c(0, NA),
    expand = expansion(mult = 0.0, add = 0.3)
  ) +
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
    axis.title = element_text(face="bold", size = 6),
    axis.text = element_text(
      color="black", size = 5, margin = margin(0.0, 0.0, 0.0, 0.0, 'mm')
      ),
    axis.text.x = element_text(
      color="black", size = 5, margin = margin(0.0, 0.0, 0.0, 0.0, 'mm'),
      angle = 30, hjust = 0.5, vjust = 0.5
      ),
    # Remove x axis text and title
    axis.line = element_line(linewidth = 0.1, colour = "black"),
    plot.title = element_text(hjust = 0.5, size = 5),
    panel.grid = element_line(colour = 'grey', linewidth = 0.1),
    panel.border = element_rect(linewidth = 0.1),
    strip.background = element_rect(linewidth = 0.01),
    strip.text = element_text(
      colour = 'black', size = 6, margin = margin(b = 0.3, t = 0.3, unit='mm')
    ),
    axis.ticks = element_line(linewidth = 0.1),
    axis.ticks.length = unit(0.1, 'lines')
  )

bounds_steps_plt
# Save plot
ggsave(
  filename = paste0('plots/final_sensitivity_rbd_fast_freqs_bounds_steps.png'),
  plot = bounds_steps_plt, width = 18, height = 5, units = 'cm', dpi = 360
  )

# Visualize the spectral mask for the last step
bounds_steps %>%
  # Select last step
  filter(step == 21) %>%
  View()


# Plot final spectral mask
mask_plt <- bounds_steps %>%
  # Filter last step and select vars and mean_bound
  select(vars, mean_bound, step, position) %>%
  rename(mask = mean_bound) %>%
  filter(step == 21) %>%
  mutate(
    # Rename window to Window size and decay to Decay rate
    vars = ifelse(vars == 'window', 'Window size', vars),
    vars = ifelse(vars == 'decay', 'Decay rate', vars),
    # Remove Hz at the end of vars
    vars_int = gsub(' Hz', '', vars) %>% as.numeric(),
  ) %>% 
  ggplot(aes(x = vars_int, y = mask)) +
  # Mask values
  geom_step(
    color = pal_nejm()(8)[2], na.rm = T, stat = 'identity', linewidth = 0.2,
    # position = position_nudge(x = -1),
    direction = 'mid',
    ) +
  geom_errorbar(
    aes(ymin = mask, ymax = mask), color = pal_nejm()(8)[1], linewidth = 0.3, width = 1.9
    ) +
  geom_hline(yintercept = 1, color = 'black', linetype = 'dashed', linewidth = 0.3) +
  # Add labels for decay and window
  geom_label(
    data = . %>% filter(vars == 'Window size'),
    aes(x = 4, y = 5, label = paste(vars, round(mask), sep = ': ')),
    fill = 'white', color = 'black',
    vjust = 0.5, size = 2.2, hjust = 0.0, label.size = 0.1
  ) +
  geom_label(
    data = . %>% filter(vars == 'Decay rate'),
    aes(x = 4, y = 4, label = paste(vars, mask, sep = ': ')),
    fill = 'white', color = 'black',
    vjust = 0.5, size = 2.2, hjust = 0.0, label.size = 0.1
  ) +
  # Add text Amplify and Attenuate above/below horizontal line yintercept = 1
  annotate(
    geom = 'text',
    x = 200, y = 1.4, label = 'Amplification',
    color = 'black', size = 2.2, hjust = 1.0
  ) +
  annotate(
    geom = 'text',
    x = 200, y = 0.6, label = 'Attenuation',
    color = 'black', size = 2.2, hjust = 1.0
  ) +
  scale_x_continuous(
    name = 'Spectral variables [Hz]',
    breaks = seq(0, 202, 10),
    expand = expansion(mult = 0., add = 2),
    # Wrap too long labels
    # labels = function(x) str_wrap(x, 10)
  ) +
  scale_y_continuous(
    name = 'Spectral mask [/]',
    # trans = scales::pseudo_log_trans(sigma = 1e-1),
    breaks = seq(0, 6.5, 0.5),
    limits = c(0, 5.5),
    expand = expansion(mult = 0.0, add = 0.25),
    # Add another y axis
    sec.axis = sec_axis(~ ., breaks = seq(0, 6.5, 0.5))
  ) +
  scale_fill_nejm() +
  scale_color_nejm() +
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
    panel.spacing.y = unit(0.5, 'mm'),
    panel.spacing.x = unit(0, 'mm'),
    axis.title = element_text(face="bold", size = 6),
    axis.text = element_text(
      color="black", size = 5, margin = margin(0.0, 0.0, 0.0, 0.0, 'mm')
      ),
    axis.text.x = element_text(
      color="black", size = 5, margin = margin(0.0, 0.0, 0.0, 0.0, 'mm'),
      # angle = 45, 
      hjust = 0.5, vjust = 0.5
      ),
    # Remove x axis text and title
    axis.line = element_line(linewidth = 0.1, colour = "black"),
    plot.title = element_text(hjust = 0.5, size = 5),
    panel.grid = element_line(colour = 'grey', linewidth = 0.1),
    panel.border = element_rect(linewidth = 0.1),
    strip.background = element_rect(linewidth = 0.01),
    strip.text = element_text(
      colour = 'black', size = 4, margin = margin(b = 0.3, t = 0.3, unit='mm')
    ),
    axis.ticks = element_line(linewidth = 0.1),
    axis.ticks.length = unit(0.1, 'lines')
    )

mask_plt
# Save plot
ggsave(
  filename = paste0('plots/window_spectral_mask.png'),
  plot = mask_plt, width = 3.5, height = 1.2, units = 'in', dpi = 420
  )

# Export mask data
bounds_steps %>%
  # Filter last step and select vars and mean_bound
  select(vars, mean_bound, step) %>%
  rename(mask = mean_bound) %>%
  filter(step == 21) %>%
  # Remove step column
  select(-step) %>%
  mutate(
    # Remove Hz at the end of vars
    vars = gsub(' Hz', '', vars),
    # Round mask for vars == window
    mask = ifelse(vars == 'window', round(mask), mask),
    # Put vars decay and window at the end of vector
    vars = factor(vars, levels = c(seq(2, 202, 2), 'decay', 'window'))
  ) %>%
  arrange(vars) %>%
  # Export to csv
  write_csv(file = 'data/spectral_mask_optimal.csv')


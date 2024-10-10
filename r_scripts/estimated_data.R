# Plot and analyze batch estimated EMG and grip force data with error metrics
library(tidyverse)
library(magrittr)
library(ggsci)
library(xtable)
library(clipr)

# Load data from position 3 and 4
est_dat <- read_csv('data/processed_approximated_data_dwns8_pos3_td_60_grid_obs.csv')
# Bind_rows with position 4
est_dat <- est_dat %>% 
  bind_rows(
    read_csv('data/processed_approximated_data_dwns8_pos4_td_60_grid_obs.csv')
  )
# Expand file_names to subject, position, replicate, age and gender
est_dat <- est_dat %>% 
  rename(file_name = file_names) %>%
  mutate(
    # Get initials from file_name
    subject = str_sub(file_name, 1, 2),
    # Get position from file_name
    position = str_sub(file_name, 4, 4),
    # Get replicate from file_name
    replicate = str_sub(file_name, 6, 6),
    # Get age from file_name
    age = str_sub(file_name, 8, 9),
    # Get gender from file_name
    gender = str_sub(file_name, 11, 11)
  )

# Select only error metrics, subject, position and replicate for RBD
rbd_design <- est_dat %>% 
  select(subject, position, replicate, wmape, rmse, mae, mae_mvc) %>% 
  distinct() %>% 
  # Arrange by initials, position and replicate
  arrange(subject, position, replicate) %>%
  # Make subjects and position factors
  mutate(
    subject = factor(subject),
    position = factor(position)
  ) %>% 
  mutate(
    wmape_mean = mean(wmape),
  ) %>% 
  # Estimate the effect of position
  group_by(position) %>%
  mutate(
    wmape_position = mean(wmape),
    effect_position = wmape_position - wmape_mean
  ) %>% 
  ungroup() %>% 
  # Estimate the effect of subject
  group_by(subject) %>%
  mutate(
    wmape_subject = mean(wmape),
    effect_subject = wmape_subject -wmape_mean
  ) %>% 
  ungroup() %>% 
  # Estimate the mean of replicate in each position/subject block
  group_by(subject, position) %>%
  mutate(
    wmape_replicate = mean(wmape),
  ) %>% 
  ungroup()
  
rbd_design %>% 
  select(
    subject, position, wmape, wmape_mean, wmape_position, effect_position,
    wmape_subject, effect_subject, wmape_replicate
    ) %>%
  distinct() %>% 
  View()

# Xtable formatting
sanitize_titles_xt <- function(x){
  paste0('\\bfseries{\\small ', x, '}')
}
sanitize_rownames_xt <- function(x){
  paste0('\\footnotesize ', x)
}
# Analyze RBD with analysis o variance and generate latex table
aov(wmape ~ Position + Subject, 
    data = rbd_design %>% rename(Subject = subject, Position = position)) %>% 
  summary() %>% xtable(
    caption = 'RBD results - non-significant effect of measuring position on wMAPE.',
    label = 'tab:rbd-aov',
  ) %>% 
  print() %>% 
  write_clip(object_type = 'character')

# Export rbd_design to latex table
pos_3_blk <- rbd_design %>% 
  select(
    subject, position, replicate, wmape, wmape_mean, wmape_position, effect_position,
    wmape_subject, effect_subject
    ) %>%
  distinct() %>% 
  filter(position == '3') %>%
  select(subject, replicate, wmape, wmape_position, effect_position) %>% 
  pivot_wider(
    names_from = subject, id_cols = replicate, values_from = wmape, unused_fn = unique
    ) %>% 
  mutate(replicate = paste0('R', replicate)) %>%
  rename('M' = wmape_position, 'E' = effect_position) %>% 
  column_to_rownames('replicate')

pos_4_blk <- rbd_design %>% 
  select(
    subject, position, replicate, wmape, wmape_mean, wmape_position, effect_position,
    wmape_subject, effect_subject
    ) %>%
  distinct() %>% 
  filter(position == '4') %>%
  select(subject, replicate, wmape, wmape_position, effect_position) %>% 
  pivot_wider(
    names_from = subject, id_cols = replicate, values_from = wmape, unused_fn = unique
    ) %>% 
  mutate(replicate = paste0('R', replicate)) %>%
  rename('M' = wmape_position, 'E' = effect_position) %>% 
  column_to_rownames('replicate')

subj_eff_blk <- rbd_design %>% 
  select(
    subject, wmape_subject, effect_subject
  ) %>%
  distinct() %>% 
  rename('M' = wmape_subject, 'E' = effect_subject) %>%
  pivot_longer(
    cols = c('M', 'E'), names_to = 'Subject', values_to = 'Value'
  ) %>%
  pivot_wider(
    names_from = subject, id_cols = 'Subject', values_from = 'Value'
  ) %>%
  column_to_rownames('Subject') %>% 
  # Add 0 values to Position mean and Position effect columns to match the other blocks
  mutate('M' = 0, 'E' = 0)
  
# Convert to list of dataframes
rbd_xtable_list = list(
  pos_3_blk, pos_4_blk, subj_eff_blk
)
# Add subheadings
attr(rbd_xtable_list, 'subheadings') <- c('Position 3', 'Position 4', 'Subject')
# Add overall mean as message
attr(rbd_xtable_list, 'message') <- paste0(
  'Overall mean wMAPE: ', round(unique(rbd_design$wmape_mean), 2), ' \\%')

str(rbd_xtable_list)

rbd_xtable_list %>% xtableList(
  align = "LRRRRRRRRRRRRR||RR",
  # align = paste(
  #   c(
  #     rep('p{0.04\\columnwidth}', 1),
  #     rep('p{0.05\\columnwidth}', 13), 
  #     '||', 
  #     rep('p{0.05\\columnwidth}', 2)
  #     ),
  #   collapse = ''
  #   ),
  digits = rep(1, 16),
  ) %>%
  print(
    tabular.environment = 'tabulary',
    width="\\columnwidth",
    size = 'scriptsize',
    floating.environment = 'table',
    floating = T,
    booktabs = T,
    sanitize.text.function = sanitize_rownames_xt,
    sanitize.colnames.function = sanitize_titles_xt,
    sanitize.subheadings.function = sanitize_titles_xt,
  ) %>%
  write_clip(object_type = 'character')

rbd_design %>% summary()

est_dat %>% 
  select(subject, position, replicate, wmape, rmse, mae, mae_mvc) %>% 
  distinct() %>% 
  # Arrange by initials, position and replicate
  arrange(subject, position, replicate) %>%
  # Make subjects and position factors
  mutate(
    subject = factor(subject),
    position = factor(position)
  ) %>% 
  mutate(
    wmape_mean = mean(wmape),
    mae_mvc_mean = mean(mae_mvc)
  ) %>% View()

# Plot approximations and processed EMG
sec_axis_coeff <- 200
est_dat %>% 
  filter(
    file_name %in% c('lmp3m1_23_m', 'dsp3m2_22_m', 'pbp4m1_24_m')
    ) %>% 
  mutate(
    file_name = case_when(
      file_name == 'lmp3m1_23_m' ~ 'lm - P3 - R1',
      file_name == 'dsp3m2_22_m' ~ 'ds - P3 - R2',
      file_name == 'pbp4m1_24_m' ~ 'pb - P4 - R1'
    )
  ) %>% 
  group_by(file_name) %>%
  mutate(
    time_t = time_t - min(time_t)
  ) %>% 
  ungroup() %>% 
  mutate(
    emg_processed = emg_processed * sec_axis_coeff
  ) %>% 
  rename(
    'Measured Grip Force' = grip, 
    'Estimated Grip Force' = grip_approx,
    'Processed EMG' = emg_processed
    ) %>%
  pivot_longer(
    cols = c(
      'Measured Grip Force', 'Estimated Grip Force', 'Processed EMG'
      ),
    names_to = "measure",
    values_to = "value"
  ) %>%
  # Ordering of plots using factors
  mutate(measure = factor(measure, levels = c(
    'Processed EMG', 'Measured Grip Force', 'Estimated Grip Force'
    ))) %>%
  ggplot(aes(x = time_t, y = value)) +
  facet_wrap(~file_name, scales = 'free_x') +
  geom_line(aes(color =  measure), linewidth = 0.3) +
  scale_y_continuous(
    breaks = seq(0, 350, 50),
    sec.axis = sec_axis(
      ~ . / sec_axis_coeff,
      name = "Processed EMG [mV]",
      breaks = seq(0, 2, 0.25)
      )) +
  scale_x_continuous(
    breaks = seq(-5, 35, 5),
    minor_breaks = seq(-5, 35, 0.5),
    expand = c(0, 0.2)
      ) +
  scale_color_nejm() +
  labs(
    x = "Time [s]",
    y = "Grip force [N]"
  ) +
  theme_bw() + theme(
    legend.position = 'inside',
    legend.position.inside = c(0.91, 0.75),
    legend.title = element_blank(),
    legend.box = 'vertical',
    legend.direction = 'vertical',
    legend.box.spacing = unit(0.5, 'mm'),
    legend.spacing.y = unit(0.5, 'mm'),
    legend.spacing.x = unit(0.5, 'mm'),
    legend.margin = margin(0, 0, 0, 0, 'mm'),
    legend.key.spacing = unit(0.5, 'mm'),
    legend.key.size = unit(3, 'mm'),
    legend.text = element_text(
      colour="black", size = 6, margin = margin(0, 0., 1, 0, 'mm')
      ),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, 'mm'),
    panel.background = element_blank(),
    panel.spacing.y = unit(0, 'mm'),
    panel.spacing.x = unit(1.5, 'mm'),
    panel.grid.minor.x = element_line(
      color = 'azure4', linewidth = 0.1, linetype = 'dotdash'
      ),
    axis.title = element_text(face="bold", size = 6),
    axis.text = element_text(
      color="black", size = 6, margin = margin(0.0, 0.0, 0.0, 0.0, 'mm')
      ),
    # Remove x axis text and title
    axis.line = element_line(linewidth = 0.2, colour = "black"),
    plot.title = element_text(hjust = 0.5, size = 5),
    panel.grid = element_line(colour = 'grey', linewidth = 0.1),
    panel.border = element_rect(linewidth = 0.2),
    strip.background = element_rect(linewidth = 0.01),
    strip.text = element_text(
      colour = 'black', size = 6.0, margin = margin(b = 0.3, t = 0.3, unit='mm')
    ),
    axis.ticks = element_line(linewidth = 0.2),
    # axis.ticks.y = element_blank(),
    axis.ticks.length = unit(0.1, 'lines')
  )
ggsave(
  'plots/proc_emg_grip_approx_example.png',
  width = 7, height = 1.5, units = 'in', dpi = 360
  )

est_dat %>% str()

# Plot raw EMG signal
est_dat %>% 
  filter(
    file_name %in% c('lmp3m1_23_m', 'dsp3m2_22_m', 'pbp4m1_24_m')
    ) %>% 
  mutate(
    file_name = case_when(
      file_name == 'lmp3m1_23_m' ~ 'lm - P3 - R1',
      file_name == 'dsp3m2_22_m' ~ 'ds - P3 - R2',
      file_name == 'pbp4m1_24_m' ~ 'pb - P4 - R1'
    )
  ) %>% 
  group_by(file_name) %>%
  mutate(
    time_t = time_t - min(time_t)
  ) %>% 
  ungroup() %>% 
  rename(
    'Raw EMG' = emg, 
    ) %>%
  pivot_longer(
    cols = c(
      'Raw EMG'
      ),
    names_to = "measure",
    values_to = "value"
  ) %>%
  # Ordering of plots using factors
  mutate(measure = factor(measure, levels = c('Raw EMG'))) %>%
  ggplot(aes(x = time_t, y = value)) +
  facet_wrap(~file_name, scales = 'free_x') +
  geom_line(aes(color =  measure), linewidth = 0.3) +
  scale_y_continuous(
    breaks = seq(-5, 5, 0.5)
    ) +
  scale_x_continuous(
    breaks = seq(-5, 35, 5),
    # minor_breaks = seq(-5, 35, 0.5),
    expand = c(0, 0.2)
      ) +
  scale_color_nejm() +
  labs(
    x = "Time [s]",
    y = "Raw EMG [mV]"
  ) +
  theme_bw() + theme(
    legend.position = 'inside',
    legend.position.inside = c(0.95, 0.5),
    legend.title = element_blank(),
    legend.box = 'vertical',
    legend.direction = 'vertical',
    legend.box.spacing = unit(0.5, 'mm'),
    legend.spacing.y = unit(0.5, 'mm'),
    legend.spacing.x = unit(0.5, 'mm'),
    legend.margin = margin(0, 0, 0, 0, 'mm'),
    legend.key.spacing = unit(0.5, 'mm'),
    legend.key.size = unit(3, 'mm'),
    legend.text = element_text(
      colour="black", size = 6, margin = margin(0, 0., 1, 0, 'mm')
      ),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, 'mm'),
    panel.background = element_blank(),
    panel.spacing.y = unit(0, 'mm'),
    panel.spacing.x = unit(1.5, 'mm'),
    # panel.grid.minor.x = element_line(
    #   color = 'azure4', linewidth = 0.1, linetype = 'dotdash'
    #   ),
    axis.title = element_text(face="bold", size = 6),
    axis.text = element_text(
      color="black", size = 6, margin = margin(0.0, 0.0, 0.0, 0.0, 'mm')
      ),
    # Remove x axis text and title
    axis.line = element_line(linewidth = 0.2, colour = "black"),
    plot.title = element_text(hjust = 0.5, size = 5),
    panel.grid = element_line(colour = 'grey', linewidth = 0.1),
    panel.border = element_rect(linewidth = 0.2),
    strip.background = element_rect(linewidth = 0.01),
    strip.text = element_text(
      colour = 'black', size = 6.0, margin = margin(b = 0.3, t = 0.3, unit='mm')
    ),
    axis.ticks = element_line(linewidth = 0.2),
    # axis.ticks.y = element_blank(),
    axis.ticks.length = unit(0.1, 'lines')
  )
ggsave(
  'plots/raw_emg_example.png',
  width = 6.65, height = 1., units = 'in', dpi = 360
  )

# Get unique values for correlation lag and summary statistics
est_dat %>% str()
est_dat %>% 
  select(
    file_name, subject, position, replicate, age, gender, max_corr, max_corr_lag,
  ) %>% 
  # Transfer lag to ms if data is sampled with 993 Hz
  mutate(
    max_corr_lag = max_corr_lag / 993 * 1000
  ) %>%
  group_by(file_name) %>%
  distinct() %>% 
  ungroup() %>%
  group_by(position) %>%
  select(
    file_name, position, max_corr, max_corr_lag
  ) %>%
  split(.$position) %>%
  map(~summary(.x, digits = 3))
  summarise(
    n = n(),
    mean_max_corr = mean(max_corr),
    median_max_corr = median(max_corr),
    min_max_corr = min(max_corr),
    max_max_corr = max(max_corr),
    sd_max_corr = sd(max_corr),
    mean_max_corr_lag = mean(max_corr_lag),
    median_max_corr_lag = median(max_corr_lag),
    sd_max_corr_lag = sd(max_corr_lag)
  )

  

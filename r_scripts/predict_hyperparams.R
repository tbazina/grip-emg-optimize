# Plot and analyze hyperparameter tuning results.
library(tidyverse)
library(magrittr)
library(ggsci)
library(xtable)
library(clipr)

# Load data from both positions and all hyperparameters.
hyp_dat <- read_csv('results/optim_grid_pred_wnd1.2_pos3.csv') %>%
  bind_rows(
    read_csv('results/optim_grid_pred_wnd1.3_pos3.csv'),
    read_csv('results/optim_grid_pred_wnd1.4_pos3.csv'),
    read_csv('results/optim_grid_pred_wnd1.5_pos3.csv'),
    read_csv('results/optim_grid_pred_wnd1.6_pos3.csv'),
    read_csv('results/optim_grid_pred_wnd1.7_pos3.csv'),
    # read_csv('results/optim_grid_pred_wnd1.8_pos3.csv'),
    read_csv('results/optim_grid_pred_wnd1.2_pos4.csv'),
    read_csv('results/optim_grid_pred_wnd1.3_pos4.csv'),
    read_csv('results/optim_grid_pred_wnd1.4_pos4.csv'),
    read_csv('results/optim_grid_pred_wnd1.5_pos4.csv'),
    read_csv('results/optim_grid_pred_wnd1.6_pos4.csv'),
    read_csv('results/optim_grid_pred_wnd1.7_pos4.csv'),
    # read_csv('results/optim_grid_pred_wnd1.8_pos4.csv')
    read_csv('results/optim_grid_pred_smooth1.1_pos3.csv'),
    read_csv('results/optim_grid_pred_smooth1.1_pos4.csv'),
  ) %>% select(-index) %>% 
  mutate(
    position = factor(position, levels = c(3, 4)),
  )

# remove p* from filename
hyp_dat <- hyp_dat %>% 
  mutate(
    filename = str_remove(filename, 'p[0-9]')
  )

# Pivot wider using filename and wmape_pred
hyp_dat <- hyp_dat %>% 
  pivot_wider(names_from = filename, values_from = wmape_pred, names_prefix = 'wmape_')

# Compute mean across all positions for single hyperparameter combination
hyp_dat <- hyp_dat %>% 
  group_by(
    batch_wnd_coeff, batch_smooth_coeff, thin_step, num_delays_predict, force_rank
  ) %>%
  mutate(
    mean_wmape = mean(c_across(starts_with('wmape_'))),
  ) %>%
  ungroup() %>% 
  arrange(batch_wnd_coeff, batch_smooth_coeff, thin_step, num_delays_predict, force_rank)

# Compute minimum values per each hyperparameter for plotting
hyp_dat <- hyp_dat %>%
  drop_na() %>%
  group_by(batch_wnd_coeff) %>%
  mutate(min_wmape_wnd = min(mean_wmape)) %>%
  ungroup() %>%
  group_by(batch_smooth_coeff) %>%
  mutate(min_wmape_smooth = min(mean_wmape)) %>%
  ungroup() %>%
  group_by(thin_step) %>%
  mutate(min_wmape_thin = min(mean_wmape)) %>%
  ungroup() %>%
  group_by(num_delays_predict) %>%
  mutate(min_wmape_delays = min(mean_wmape)) %>%
  ungroup() %>%
  group_by(force_rank) %>%
  mutate(min_wmape_rank = min(mean_wmape)) %>%
  ungroup()

# Panel plot mean_wmape vs all hyperparameters
hyp_plot <- hyp_dat %>% 
  # Rename hyperparameters to readable names
  rename(
    'Window coefficient' = batch_wnd_coeff,
    'Smooth coefficient' = batch_smooth_coeff,
    'Thinning step' = thin_step,
    'No. time delays' = num_delays_predict,
    'No. modes' = force_rank
  ) %>% 
  pivot_longer(
    cols = c('Window coefficient', 'Smooth coefficient', 'Thinning step',
             'No. time delays', 'No. modes'),
    names_to = 'hyperparameter', values_to = 'hyp_val') %>% 
  mutate(
    min_wmape = case_when(
      hyperparameter == 'Window coefficient' ~ min_wmape_wnd,
      hyperparameter == 'Smooth coefficient' ~ min_wmape_smooth,
      hyperparameter == 'Thinning step' ~ min_wmape_thin,
      hyperparameter == 'No. time delays' ~ min_wmape_delays,
      hyperparameter == 'No. modes' ~ min_wmape_rank
    )
  ) %>% 
  ggplot(aes(x = hyp_val, y = mean_wmape)) +
  facet_grid(~hyperparameter, scales = 'free_x') +
  geom_point(
    shape = 21,
    fill = pal_nejm()(8)[4],
    color = 'black',
    size = 0.3,
    stroke = 0.1
  ) +
  scale_x_continuous(
    breaks = c(
      1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 
      3, 4, 5, 6, 7, 8, 9, 10),
  ) +
  scale_y_continuous(
    limits = c(NA, 20)) +
  geom_line(
    data = . %>% select(hyperparameter, hyp_val, min_wmape) %>% distinct(),
    aes(y = min_wmape, alpha = 0.2),
    stat = 'smooth',
    method = 'loess',
    span = 0.65,
    se = F,
    color = pal_nejm()(8)[1],
    linewidth = 0.3,
    alpha = 0.8
    ) +
  labs(
    x = 'Hyperparameter value',
    y = 'Mean wMAPE (%)',
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
    panel.spacing.y = unit(0.5, 'mm'),
    panel.spacing.x = unit(1.0, 'mm'),
    axis.title = element_text(face="bold", size = 5),
    axis.text = element_text(
      color="black", size = 4, margin = margin(0.0, 0.0, 0.0, 0.0, 'mm')
      ),
    axis.text.x = element_text(
      color="black", size = 4, margin = margin(0.0, 0.0, 0.0, 0.0, 'mm'),
      angle = 45, hjust = 0.5, vjust = 0.5
      ),
    # Remove x axis text and title
    axis.line = element_line(linewidth = 0.1, colour = "black"),
    plot.title = element_text(hjust = 0.5, size = 5),
    panel.grid = element_line(colour = 'grey', linewidth = 0.1),
    panel.border = element_rect(linewidth = 0.1),
    strip.background = element_rect(linewidth = 0.01),
    strip.text = element_text(
      colour = 'black', size = 4.5, margin = margin(b = 0.3, t = 0.3, unit='mm')
    ),
    axis.ticks = element_line(linewidth = 0.1),
    axis.ticks.length = unit(0.1, 'lines')
    )

hyp_plot
# Save plot
ggsave(
  filename = paste0('plots/hyperparameter_tuning.png'),
  plot = hyp_plot, width = 8.8, height = 2.5, units = 'cm', dpi = 360
  )

# Get hyperparameters with wmape < 18.5
hyp_dat %>% 
  filter(mean_wmape < 18.5) %>% 
  select(
    position, batch_wnd_coeff, batch_smooth_coeff, thin_step, num_delays_predict, force_rank,
    mean_wmape, starts_with('wmape_')
  ) %>% 
  group_by(batch_wnd_coeff, batch_smooth_coeff, thin_step, num_delays_predict,
           force_rank) %>%
  mutate(
    median_wmape = median(c_across(starts_with('wmape_')), na.rm = T),
    mean_median_wmape = median_wmape + mean_wmape
  ) %>% 
  ungroup() %>%
  relocate(mean_median_wmape, median_wmape, mean_wmape) %>%
  arrange(mean_median_wmape) %>%
  View()

# Optimal hyperparameters
# batch_wnd_coeff - 1.3
# batch_smooth_coeff - 1.1
# thin_step - 7
# num_delays_predict - 8
# force_rank - 4

# Extract wMAPE for optimal hyperparameters
opt_dat <- hyp_dat %>% 
  filter(
    batch_wnd_coeff == 1.3,
    batch_smooth_coeff == 1.1,
    thin_step == 7,
    num_delays_predict == 8,
    force_rank == 4
  ) %>% 
  select(
    position, starts_with('wmape_')
  ) %>%
  pivot_longer(
    cols = starts_with('wmape_'),
    names_to = 'file_name',
    values_to = 'wmape'
  ) %>%
  mutate(
    file_name = str_remove(file_name, 'wmape_')
  ) %>%
  mutate(
    # Get initials from file_name
    subject = str_sub(file_name, 1, 2),
    # Get replicate from file_name
    replicate = str_sub(file_name, 4, 4),
    # Get age from file_name
    age = str_sub(file_name, 6, 7),
    # Get gender from file_name
    gender = str_sub(file_name, 9, 9)
  )

# Select only error metrics, subject, position and replicate for RBD
rbd_design <- opt_dat %>% 
  select(subject, position, replicate, wmape) %>% 
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
    caption = 'RBD results for predicting grip force - non-significant effect of measuring position on wMAPE.',
    label = 'tab:rbd-aov-pred',
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
  'Overall prediction mean wMAPE: ', round(unique(rbd_design$wmape_mean), 2), ' \\%')

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

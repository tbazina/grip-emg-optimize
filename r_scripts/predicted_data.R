# Plot and analyze batch predicted grip force data
library(tidyverse)
library(magrittr)
library(ggsci)
library(xtable)
library(clipr)

# Load data for 3 selected measurements
pred_dat <- read_csv('data/predicted_data_optimal_pbp4m1_24_m.csv') %>% 
  bind_rows(
    read_csv('data/predicted_data_optimal_lmp3m1_23_m.csv'),
    read_csv('data/predicted_data_optimal_dsp3m2_22_m.csv')
  )

# Expand file_names to subject, position, replicate, age and gender
pred_dat <- pred_dat %>% 
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

pred_dat %>% str()

# Plot smoothed grip approximations, grip and predicted grip
pred_dat %>% 
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
    'Measured Grip Force' = grip, 
    'Smoothed Grip Force Estimation' = grip_approx_smooth,
    'Predicted Grip Force' = grip_approx_predict
    ) %>%
  pivot_longer(
    cols = c(
      'Measured Grip Force', 'Smoothed Grip Force Estimation'
      ),
    names_to = "measure",
    values_to = "value"
  ) %>%
  # Ordering of plots using factors
  mutate(measure = factor(measure, levels = c(
    'Smoothed Grip Force Estimation', 'Measured Grip Force' 
    ))) %>%
  ggplot(aes(x = time_t, y = value)) +
  facet_wrap(~file_name, scales = 'free_x') +
  geom_line(aes(color =  measure), linewidth = 0.3) +
  geom_point(
    aes(y = `Predicted Grip Force`, fill = 'Predicted Grip Force'), 
    size = 0.5,
    stroke = 0.05,
    shape = 21,
    alpha = 0.7,
    color = 'black',
    # fill = pal_nejm()(8)[1],
    show.legend = T
    ) +
  scale_y_continuous(
    breaks = seq(0, 350, 50),
      ) +
  scale_x_continuous(
    breaks = seq(-5, 35, 5),
    minor_breaks = seq(-5, 35, 0.5),
    expand = c(0, 0.2)
      ) +
  scale_color_discrete(
    type = pal_nejm()(8)[c(3, 2, 1)],
  ) +
  scale_fill_manual(
    values = c('Predicted Grip Force' = pal_nejm()(8)[1])
  ) +
  labs(
    x = "Time [s]",
    y = "Grip force [N]"
  ) +
  # guides(fill = guide_legend(override.aes = list(shape = 21))) +
  theme_bw() + theme(
    legend.position = 'inside',
    legend.position.inside = c(0.87, 0.8),
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
  'plots/grip_approx_smooth_predict_example.png',
  width = 6.65, height = 1.5, units = 'in', dpi = 360
  )

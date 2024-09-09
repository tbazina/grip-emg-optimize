# Plot and analyze processed EMG and grip force data with correlations
library(tidyverse)
library(magrittr)
library(ggsci)

# Load data from position 3/4
proc_dat <- read_csv('data/processed_data_pos3.csv')
grip_approx <- read_csv('data/grip_approximation.csv', col_names = F)

# Panel plot all the data
proc_dat %>%
  # pivot longer using emg and grip_force
  rename(proc = emg_processed, file_name = file_names) %>%
  mutate(
    # Get initials from file_name
    initials = str_sub(file_name, 1, 2),
    # Get position from file_name
    position = str_sub(file_name, 4, 4),
    # Get repetition from file_name
    repetition = str_sub(file_name, 6, 6),
    # Get age from file_name
    age = str_sub(file_name, 8, 9),
    # Get gender from file_name
    gender = str_sub(file_name, 11, 11)
  ) %>% 
  # Only filter md-3-1
  filter(initials == 'md', position == 3, repetition == 1) %>%
  # Remove minimu time from time ti fix x axis
  mutate(time_t = time_t - min(time_t)) %>%
  pivot_longer(
    cols = c(emg, proc, grip),
    names_to = "measure",
    values_to = "value"
    ) %>%
  # Set each id by pasting initials, position and repetition
  mutate(id = paste(initials, position, repetition, measure, sep = '-')) %>%
  ggplot(aes(x = time_t, y = value, color = initials)) +
  facet_wrap(~ id, scales = 'free') +
  geom_line(linewidth = 0.2) +
  scale_color_d3(
    palette = "category20b",
  ) +
  labs(
    x = "Timestamp [s]",
    y = "EMG [mV] / Grip force [N] / Processed EMG [/]",
    title = "Position 4 - EMG, Grip Force and Processed EMG Data"
  ) +
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
    panel.spacing.x = unit(2, 'mm'),
    axis.title = element_text(face="bold", size = 5),
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
      colour = 'black', size = 5.0, margin = margin(b = 0.3, t = 0.3, unit='mm')
    ),
    axis.ticks = element_line(linewidth = 0.1),
    # axis.ticks.y = element_blank(),
    axis.ticks.length = unit(0.1, 'lines')
  )
# Save only md-3-1 processed and raw emg
ggsave(
  'plots/emg_grip_proc_position_3_md-3-1.png',
  width = 10, height = 2, units = 'cm', dpi = 320
  )

ggsave(
  'plots/emg_grip_proc_position_4.png',
  width = 18, height = 10, units = 'cm', dpi = 320
  )

# Join grip approx and processed data
proc_dat %>% 
  filter(file_names == 'mdp3m1_24_m') %>% 
  mutate(
    grip_approx = grip_approx$X1,
    time_t = time_t - min(time_t)
  ) %>% 
  rename('Measured Grip Force' = grip, 'Estimated Grip Force' = grip_approx) %>%
  pivot_longer(
    cols = c('Measured Grip Force', 'Estimated Grip Force'),
    names_to = "measure",
    values_to = "value"
  ) %>%
  ggplot(aes(x = time_t, y = value)) +
  geom_line(aes(color =  measure), linewidth = 0.5) +
  scale_y_continuous(
    breaks = seq(0, 350, 50)
      ) +
  scale_x_continuous(
    breaks = seq(0, 35, 5)
      ) +
  scale_color_nejm() +
  labs(
    x = "Timestamp [s]",
    y = "Grip force / estimate [N]"
  ) +
  theme_bw() + theme(
    legend.position = 'inside',
    legend.position.inside = c(0.85, 0.9),
    legend.title = element_blank(),
    legend.box = 'vertical',
    legend.direction = 'vertical',
    legend.box.spacing = unit(5.5, 'mm'),
    legend.spacing.y = unit(5.5, 'mm'),
    legend.spacing.x = unit(5.5, 'mm'),
    legend.margin = margin(0, 0, 0, 0, 'mm'),
    legend.key.spacing = unit(5.5, 'mm'),
    legend.key.size = unit(5, 'mm'),
    legend.text = element_text(
      colour="black", size = 9, margin = margin(0, 0., 1, 0, 'mm')
      ),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, 'mm'),
    panel.background = element_blank(),
    panel.spacing.y = unit(0, 'mm'),
    panel.spacing.x = unit(0, 'mm'),
    axis.title = element_text(face="bold", size = 8.5),
    axis.text = element_text(
      color="black", size = 9, margin = margin(0.0, 0.0, 0.0, 0.0, 'mm')
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
ggsave(
  'plots/grip_approx_mdp3m1_24_m.png',
  width = 7, height = 2, units = 'in', dpi = 420
  )

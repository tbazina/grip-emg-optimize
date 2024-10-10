# Actuators paper vizualizations
library(tidyverse)
library(magrittr)
library(ggsci)
library(xtable)
library(clipr)
library(datapasta)
library(scales)

robot_angl <- tibble::tribble(
  ~pressure, ~j1, ~j2, ~j3,
  0L,  3L,  5L,  3L,
  1L,  7L,  8L, 11L,
  2L, 12L, 14L, 17L,
  3L, 18L, 20L, 23L,
  4L, 24L, 25L, 36L,
  5L, 30L, 31L, 40L,
  6L, 37L, 37L, 45L,
  7L, 43L, 40L, 53L
) %>%  mutate(rep = 1, rob = 'left') %>% 
  bind_rows(
    tibble::tribble(
      ~pressure, ~j1, ~j2, ~j3,
      0L,  2L,  3L,  3L,
      1L,  5L,  8L,  9L,
      2L, 10L, 12L, 18L,
      3L, 16L, 18L, 24L,
      4L, 22L, 23L, 34L,
      5L, 28L, 30L, 40L,
      6L, 35L, 35L, 48L,
      7L, 42L, 40L, 53L
) %>% mutate(rep = 2, rob = 'left'),
tibble::tribble(
 ~pressure, ~j1, ~j2, ~j3,
  0L,  4L,  4L,  4L,
  1L, 12L, 12L, 18L,
  2L, 19L, 21L, 26L,
  3L, 28L, 25L, 40L,
  4L, 36L, 35L, 46L,
  5L, 45L, 41L, 55L,
  6L, 53L, 47L, 63L,
  7L, 61L, 52L, 67L
 ) %>% mutate(rep = 1, rob = 'right'),
tibble::tribble(
~pressure, ~j1, ~j2, ~j3,
 0L,  5L,  5L,  2L,
 1L, 11L, 12L, 13L,
 2L, 19L, 18L, 24L,
 3L, 27L, 26L, 35L,
 4L, 37L, 32L, 46L,
 5L, 44L, 40L, 54L,
 6L, 54L, 46L, 62L,
 7L, 61L, 53L, 69L
) %>% mutate(rep = 2, rob = 'right')
) %>% 
  full_join(
    bind_rows(
    tibble::tribble(
    ~pressure,  ~L1,  ~L2,  ~L3,  ~L4,
     0L, 16.4, 48.5, 34.5, 23.6,
     1L, 16.9, 48.3, 35.2, 25.1,
     2L, 17.2, 48.9, 35.6, 24.5,
     3L, 17.5, 49.1, 35.7, 25.1,
     4L, 17.6, 49.6, 36.1, 26.1,
     5L, 17.5, 50.2, 36.4, 26.6,
     6L, 17.9, 50.6,   37, 27.1,
     7L,  18, 50.5, 38.5, 27.2
     ) %>% mutate(rep = 1, rob = 'left'),
    tibble::tribble(
    ~pressure, ~L1,  ~L2,  ~L3,  ~L4,
     0L,    16.3, 52.7, 39.3, 29.2,
     1L,    17.6, 52.3, 41.3, 29.4,
     2L,    17.9, 53.1, 42.5, 29.4,
     3L,      18, 53.7,   43,   31,
     4L,    18.4, 54.5,   44, 31.5,
     5L,      19, 54.9, 45.2, 32.2,
     6L,      19, 56.9, 45.5, 32.6,
     7L,      20, 57.5, 45.9, 33.3
       ) %>% mutate(rep = 1, rob = 'right')
    ),
    by = c('pressure', 'rep', 'rob')
    ) %>% 
# replication as character
mutate(rep = as.character(rep))

# Mean not necessary, instead use fit
# Bind rows with mean across each rep
# robot_angl <- robot_angl %>%
#   bind_rows(
#     robot_angl %>%
#     # Get average across reps for same robot
#       group_by(rob, pressure) %>%
#       summarise(
#         j1 = mean(j1),
#         j2 = mean(j2),
#         j3 = mean(j3),
#         rep = 'mean'
#       ) %>%
#       ungroup()
#   )

# Linear regression for all revolute joints and robots across repetitions
robot_angl %>% 
  filter(rob == 'left') %>%
  # filter(rob == 'right') %>%
  # lm(j1 ~ pressure, data = .) %>%
  # lm(j2 ~ pressure, data = .) %>%
  lm(j3 ~ pressure, data = .) %>%
  summary()

# Plot revolut joints and lm with CI for all repetitions and robots
plt_rev <- robot_angl %>%
  # Select only necessary parameters
  select(-L1, -L2, -L3, -L4) %>%
  # Change left/right to index/little
  mutate(rob = ifelse(rob == 'left', 'little', 'index')) %>%
  # rename j1, j2, j3 to j3, j5, j7
  rename(j2 = j1, j4 = j2, j6 = j3) %>%
  # pivot longer all joints starting with j
  pivot_longer(
    cols = starts_with('j'),
    names_to = 'joint',
    values_to = 'value'
  ) %>%
  # Rename joint names to theta_1, theta_2, theta_3
  mutate(joint = str_replace(joint, 'j', 'theta[')) %>%
  # Add ) at end of joint names
  mutate(joint = paste0(joint, ']')) %>%
  # plot the data
  ggplot(aes(x = pressure, y = value, color=rob)) +
  facet_grid(~joint, scales = 'fixed',
             labeller = label_parsed
             ) +
  # Smooth line with CI
  geom_smooth(
    aes(color = interaction(rob, joint, sep=':')),
    method = 'lm',
    formula = y ~ x,
    se = T,
    # color = 'black',
    linewidth = 0.3,
    alpha = 0.4,
    show.legend = F
  ) +
  geom_point(
    aes(
      fill = interaction(rob, joint, sep=':'),
      # color = rep
      ),
    size = 1.0,
    stroke = 0.2,
    shape = 21,
    alpha = 0.9,
    color = 'black',
    # fill = pal_nejm()(8)[1],
    show.legend = T
    ) +
  scale_color_nejm() +
  # scale_color_discrete(
  #   type = pal_nejm()(8)[c(8, 7)],
  # ) +
  scale_fill_nejm(
    labels = label_parse(),
  ) +
  scale_x_continuous(
    breaks = seq(0, 7, 1),
    # minor_breaks = seq(0, 1, length.out = 22)^1.8,
    # expand = c(0, 0.04)
      ) +
  scale_y_continuous(
    breaks = seq(0, 70, 5),
    # minor_breaks = seq(0, 1, length.out = 22)^1.8,
    expand = c(0, 0)
      ) +
  labs(
    x = expression("Input pressure "*italic(p)*" [bar]"),
    y = expression("Joint Angle "*italic(theta)*" [°]")
  ) +
  # Customize legend
  guides(
    # Remove legend for color
    color = 'none',
    fill = guide_legend(
      title = 'Actuator : Revolute Joint',
      byrow = F,
      reverse = F
    )
  ) +
  theme_bw() + theme(
    legend.position = 'top',
    legend.title = element_text(
      size = 7, margin = margin(0, 0, 1, 0, 'mm')),
    legend.title.position = 'left',
    legend.box = 'horizontal',
    legend.direction = 'horizontal',
    legend.box.spacing = unit(0.5, 'mm'),
    legend.spacing.y = unit(0.5, 'mm'),
    legend.spacing.x = unit(0.5, 'mm'),
    legend.margin = margin(0, 0, 0, 0, 'mm'),
    legend.key.spacing = unit(0.5, 'mm'),
    legend.key.size = unit(3, 'mm'),
    legend.text = element_text(
      colour="black", size = 7, margin = margin(0, 0., 1, 0, 'mm')
      ),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, 'mm'),
    panel.background = element_blank(),
    panel.spacing.y = unit(0, 'mm'),
    panel.spacing.x = unit(1.5, 'mm'),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_line(
    #   color = 'azure4', linewidth = 0.15, linetype = 'solid'
    #   ),
    axis.title = element_text(face="bold", size = 8),
    axis.text = element_text(
      color="black", size = 7, margin = margin(0.0, 0.0, 0.0, 0.0, 'mm'),
      ),
    # axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1.),
    # Remove x axis text and title
    axis.line = element_line(linewidth = 0.2, colour = "black"),
    plot.title = element_text(hjust = 0.5, size = 5),
    panel.grid = element_line(colour = 'grey', linewidth = 0.1),
    panel.border = element_rect(linewidth = 0.2),
    strip.background = element_rect(linewidth = 0.01),
    strip.text = element_text(
      color = 'black', size = 7.0, margin = margin(b = 0.3, t = 0.3, unit='mm')
    ),
    axis.ticks = element_line(linewidth = 0.2),
    # axis.ticks.y = element_blank(),
    axis.ticks.length = unit(0.1, 'lines')
  )
plt_rev
ggsave(
  filename = 'plots/actuators_revolute_joints.png',
  plot = plt_rev,
  device = grDevices::png,
  width = 16.5, height = 7, units = 'cm', dpi = 360
  )

# Linear regression for all prismatic joints and robots (no repetitions)
robot_angl %>% 
  drop_na() %>%
  # d is link offset, subtract minimum of L to compute d range
  group_by(rob) %>%
  mutate(
    d1 = L1 - min(L1),
    d3 = L2 - min(L2),
    d5 = L3 - min(L3),
    d7 = L4 - min(L4)
    ) %>%
  ungroup() %>%
  # filter(rob == 'left') %>%
  filter(rob == 'right') %>%
  # lm(d1 ~ pressure-1, data = .) %>%
  # lm(d3 ~ pressure-1, data = .) %>%
  # lm(d5 ~ pressure-1, data = .) %>%
  lm(d7 ~ pressure-1, data = .) %>%
  summary()

# Print minimum link lengths
robot_angl %>% 
  drop_na() %>%
  # d is link offset, subtract minimum of L to compute d range
  group_by(rob) %>%
  summarise(
    L1_min = min(L1),
    L2_min = min(L2),
    L3_min = min(L3),
    L4_min = min(L4)
    ) %>%
  ungroup()

# Plot primatic joints and lm with CI for all repetitions and robots
plt_prism <- robot_angl %>%
  # Select only necessary parameters
  select(-j1, -j2, -j3) %>%
  drop_na() %>%
  group_by(rob) %>%
  mutate(
    d1 = L1 - min(L1),
    d3 = L2 - min(L2),
    d5 = L3 - min(L3),
    d7 = L4 - min(L4)
    ) %>%
  ungroup() %>%
  # Change left/right to index/little
  mutate(rob = ifelse(rob == 'left', 'little', 'index')) %>%
  # rename to proper numbering
  # rename(d1 = L1, d2 = L2, d4 = L3, d6 = L4) %>%
  # pivot longer all joints starting with j
  pivot_longer(
    cols = starts_with('d'),
    names_to = 'joint',
    values_to = 'value'
  ) %>%
  # Rename joint names to d_1, d_2, d_4, d_6
  mutate(joint = str_replace(joint, 'd', 'Delta*italic(d)[')) %>%
  # Add ] at end of joint names
  mutate(joint = paste0(joint, ']')) %>%
  # plot the data
  ggplot(aes(x = pressure, y = value, color=rob)) +
  facet_grid(~joint, scales = 'fixed',
             labeller = label_parsed
             ) +
  # Smooth line with CI
  geom_smooth(
    aes(color = interaction(rob, joint, sep=':')),
    method = 'lm',
    formula = y ~ x-1,
    se = T,
    # color = 'black',
    linewidth = 0.3,
    alpha = 0.4,
    show.legend = F
  ) +
  geom_point(
    aes(
      fill = interaction(rob, joint, sep=':'),
      # color = rep
      ),
    size = 1.0,
    stroke = 0.2,
    shape = 21,
    alpha = 0.9,
    color = 'black',
    # fill = pal_nejm()(8)[1],
    show.legend = T
    ) +
  scale_color_nejm() +
  # scale_color_discrete(
  #   type = pal_nejm()(8)[c(8, 7)],
  # ) +
  scale_fill_nejm(
    labels = label_parse(),
  ) +
  scale_x_continuous(
    breaks = seq(0, 7, 1),
    # minor_breaks = seq(0, 1, length.out = 22)^1.8,
    # expand = c(0, 0.04)
      ) +
  scale_y_continuous(
    # breaks = seq(0, 70, 5),
    # minor_breaks = seq(0, 1, length.out = 22)^1.8,
    expand = c(0, 0)
      ) +
  labs(
    x = expression("Input pressure "*italic(p)*" [bar]"),
    y = expression("Joint elongation "*Delta*italic(d)*" [mm]")
  ) +
  # Customize legend
  guides(
    # Remove legend for color
    color = 'none',
    fill = guide_legend(
      title = 'Actuator : Prismatic Joint',
      byrow = F,
      reverse = F
    )
  ) +
  theme_bw() + theme(
    legend.position = 'top',
    legend.title = element_text(
      size = 7, margin = margin(0, 0, 1, 0, 'mm')),
    legend.title.position = 'left',
    legend.box = 'horizontal',
    legend.direction = 'horizontal',
    legend.box.spacing = unit(0.5, 'mm'),
    legend.spacing.y = unit(0.5, 'mm'),
    legend.spacing.x = unit(0.5, 'mm'),
    legend.margin = margin(0, 0, 0, 0, 'mm'),
    legend.key.spacing = unit(0.5, 'mm'),
    legend.key.size = unit(3, 'mm'),
    legend.text = element_text(
      colour="black", size = 7, margin = margin(0, 0., 1, 0, 'mm')
      ),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, 'mm'),
    panel.background = element_blank(),
    panel.spacing.y = unit(0, 'mm'),
    panel.spacing.x = unit(1.5, 'mm'),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_line(
    #   color = 'azure4', linewidth = 0.15, linetype = 'solid'
    #   ),
    axis.title = element_text(face="bold", size = 8),
    axis.text = element_text(
      color="black", size = 7, margin = margin(0.0, 0.0, 0.0, 0.0, 'mm'),
      ),
    # axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1.),
    # Remove x axis text and title
    axis.line = element_line(linewidth = 0.2, colour = "black"),
    plot.title = element_text(hjust = 0.5, size = 5),
    panel.grid = element_line(colour = 'grey', linewidth = 0.1),
    panel.border = element_rect(linewidth = 0.2),
    strip.background = element_rect(linewidth = 0.01),
    strip.text = element_text(
      color = 'black', size = 7.0, margin = margin(b = 0.3, t = 0.3, unit='mm')
    ),
    axis.ticks = element_line(linewidth = 0.2),
    # axis.ticks.y = element_blank(),
    axis.ticks.length = unit(0.1, 'lines')
  )
plt_prism
ggsave(
  filename = 'plots/actuators_prismatic_joints.png',
  plot = plt_prism,
  device = grDevices::png,
  width = 16.5, height = 7, units = 'cm', dpi = 360
  )

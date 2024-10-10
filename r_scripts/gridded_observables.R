# Plot and gridded identity observables
library(tidyverse)
library(magrittr)
library(ggsci)
library(xtable)
library(clipr)

# Load data for gridded identity observables
grid_obs <- read_csv('data/gridded_identity_observables_example.csv')

# Rename columns to full names
# grid_obs <- grid_obs %>% 
#   rename(
#     'Time delay 1' = time_delay_1,
#     'Time delay 30' = time_delay_30,
#     'Time delay 60' = time_delay_60,
#   )

# Combine time delays in pairs of 2
names_comb <- combn(names(grid_obs), 2, simplify = T)
# Concatenate combinations of time delays from names_comb using grid_obs dataframe
grid_obs_comb <- bind_rows(
  grid_obs %>% select(names_comb[,1]) %>% rename(
    x_td = !!names_comb[1,1],
    y_td = !!names_comb[2,1]
    ) %>% mutate(
      time_delays = paste(str_sub(names_comb[,1], 12), collapse = "-")
      ),
  grid_obs %>% select(names_comb[,2]) %>% rename(
    x_td = !!names_comb[1,2],
    y_td = !!names_comb[2,2]
    ) %>% mutate(
      time_delays = paste(str_sub(names_comb[,2], 12), collapse = "-")),
  grid_obs %>% select(names_comb[,3]) %>% rename(
    x_td = !!names_comb[1,3],
    y_td = !!names_comb[2,3]
    ) %>% mutate(
      time_delays = paste(str_sub(names_comb[,3], 12), collapse = "-"))
)

grid_obs_comb %>% 
  mutate(time_delays = factor(str_c('Time delays ', time_delays))) %>%
  ggplot(aes(x = x_td, y = y_td)) +
  facet_wrap(~time_delays, scales = 'fixed') +
  coord_fixed() +
  # geom_line(aes(color =  measure), linewidth = 0.3) +
  geom_point(
    size = 0.4,
    stroke = 0.05,
    shape = 21,
    alpha = 0.8,
    color = 'black',
    fill = pal_nejm()(8)[1],
    show.legend = T
    ) +
  scale_x_continuous(
    breaks = seq(0, 1, 0.2),
    minor_breaks = seq(0, 1, length.out = 22)^1.8,
    expand = c(0, 0.04)
      ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.2),
    minor_breaks = seq(0, 1, length.out = 22)^1.8,
    expand = c(0, 0.04)
      ) +
  scale_color_discrete(
    type = pal_nejm()(8)[c(3, 2, 1)],
  ) +
  # scale_fill_manual(
  #   values = c('Predicted Grip Force' = pal_nejm()(8)[1])
  # ) +
  labs(
    x = "Time delay",
    y = "Time delay"
  ) +
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
    panel.grid.major = element_blank(),
    panel.grid.minor = element_line(
      color = 'azure4', linewidth = 0.15, linetype = 'solid'
      ),
    axis.title = element_text(face="bold", size = 6),
    axis.text = element_text(
      color="black", size = 6, margin = margin(0.0, 0.0, 0.0, 0.0, 'mm'),
      ),
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1.),
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
  'plots/gridded_identity_observables.png',
  width = 3.5, height = 1.5, units = 'in', dpi = 360
  )

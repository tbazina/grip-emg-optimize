# Plot and gridded identity observables
library(tidyverse)
library(magrittr)
library(ggsci)
library(xtable)
library(clipr)
library(patchwork)
library(plotly)

# Load data for gridded identity observables
grid_obs <- read_csv('data/gridded_identity_observables_example.csv')

# XY projection
plot_xy <- ggplot(grid_obs, aes(x = time_delay_1, y = time_delay_30)) +
  coord_fixed() +
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
    limits = c(-0.02, 1.02),
    breaks = seq(0, 1, 0.2),
    minor_breaks = seq(0, 1, length.out = 22)^1.8,
    expand = c(0, 0.0),
      ) +
  scale_y_continuous(
    limits = c(-0.02, 1.02),
    breaks = seq(0, 1, 0.2),
    minor_breaks = seq(0, 1, length.out = 22)^1.8,
    expand = c(0, 0.0)
      ) +
  labs(
    title = "XY projection",
    x = expression(bold(e)[td(1)]),
    subtitle = expression(bold(e)[td(30)]),
  ) +
  theme_bw()

# XZ projection
plot_xz <- ggplot(grid_obs, aes(x = time_delay_1, y = time_delay_60)) +
  coord_fixed() +
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
    limits = c(-0.02, 1.02),
    breaks = seq(0, 1, 0.2),
    minor_breaks = seq(0, 1, length.out = 22)^1.8,
    expand = c(0, 0.0),
      ) +
  scale_y_continuous(
    limits = c(-0.02, 1.02),
    breaks = seq(0, 1, 0.2),
    minor_breaks = seq(0, 1, length.out = 22)^1.8,
    expand = c(0, 0.0)
      ) +
  labs(
    title = "XZ projection",
    x = expression(bold(e)[td(1)]),
    subtitle = expression(bold(e)[td(60)]),
  ) +
  theme_bw() 

# YZ projection
plot_yz <- ggplot(grid_obs, aes(x = time_delay_30, y = time_delay_60)) +
  coord_fixed() +
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
    limits = c(-0.02, 1.02),
    breaks = seq(0, 1, 0.2),
    minor_breaks = seq(0, 1, length.out = 22)^1.8,
    expand = c(0, 0.0),
      ) +
  scale_y_continuous(
    limits = c(-0.02, 1.02),
    breaks = seq(0, 1, 0.2),
    minor_breaks = seq(0, 1, length.out = 22)^1.8,
    expand = c(0, 0.0)
      ) +
  labs(
    title = "YZ projection",
    x = expression(bold(e)[td(30)]),
    subtitle = expression(bold(e)[td(60)]),
  ) +
  theme_bw()  

# Arrange the plots in a "3D-like" layout and remove margins
combined_plot <- plot_xy + plot_xz + plot_yz + plot_layout(ncol = 3)
combined_plot <- combined_plot & theme(
  plot.margin = margin(0.0, 0.0, 0.0, 0.0, 'mm'),
  plot.title = element_text(
    size = 6, margin = margin(0.0, 0.0, 0.0, 0.0, 'mm'), 
    hjust = 0.55, vjust = -1.5
    ),
  plot.title.position = 'panel',
  plot.subtitle = element_text(
    size = 6, margin = margin(0.0, 0.0, 0.0, 0.0, 'mm'), hjust = -0.1
    ),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_line(
    color = 'azure4', linewidth = 0.15, linetype = 'solid'
    ),
  panel.grid = element_line(colour = 'grey', linewidth = 0.1),
  panel.border = element_rect(linewidth = 0.2),
  axis.title.x = element_text(size = 6),
  axis.title.y = element_blank(),
  axis.text = element_text(
    color="black", size = 6, margin = margin(0.0, 0.0, 0.0, 0.0, 'mm'),
    ),
  axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1.),
  # Remove x axis text and title
  axis.line = element_line(
    linewidth = 0.2, colour = "black", arrow = arrow(
      angle = 10, type = "open", length = unit(0.05, "inches"))
    ),
  axis.ticks = element_line(linewidth = 0.2),
  axis.ticks.length = unit(0.1, 'lines')
  )
# Add left margin to the second and third plots, remove y axis ticks and labels
combined_plot[[2]] <- combined_plot[[2]] + theme(
  plot.margin = margin(0.0, 0.0, 0.0, 1.5, 'mm'),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank()
  )
combined_plot[[3]] <- combined_plot[[3]] + theme(
  plot.margin = margin(0.0, 0.0, 0.0, 1.5, 'mm'),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank()
  )
ggsave(
  'plots/gridded_indicator_observables_new.png',
  plot = combined_plot,
  width = 3.4, height = 1.5, units = 'in', dpi = 360
  )


fig <- plot_ly(grid_obs) %>%
  # Add markers on each 2d projection
  # add_markers(
  #   x = ~time_delay_1,
  #   y = ~time_delay_30,
  #   z = 0,
  #   mode = 'markers',
  #   marker = list(
  #     size = 2,
  #     color = 'red',
  #     opacity = 0.8
  #   )
  # ) %>%
  # add_markers(
  #   x = ~time_delay_1,
  #   y = 0,
  #   z = ~time_delay_60,
  #   mode = 'markers',
  #   marker = list(
  #     size = 2,
  #     color = 'black',
  #     opacity = 0.8
  #   )
  # ) %>%
  # add_markers(
  #   x = 0,
  #   y = ~time_delay_30,
  #   z = ~time_delay_60,
  #   mode = 'markers',
  #   marker = list(
  #     size = 2,
  #     color = 'blue',
  #     opacity = 0.8
  #   )
  # ) %>%
  add_markers(
    x = ~time_delay_1,
    y = ~time_delay_30,
    z = ~time_delay_60,
    type = 'scatter3d',
    mode = 'markers',
    marker = list(
      size = 2,
      color = 'black',
      colorscale = 'Viridis',
      opacity = 0.8
    )
  ) %>%
  layout(
    scene = list(
      xaxis = list(title = 'Time delay 1'),
      yaxis = list(title = 'Time delay 30'),
      zaxis = list(title = 'Time delay 60')
    )
  )
  
fig

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

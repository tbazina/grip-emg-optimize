# Clean input EMG ch1 and grip force timeseries for numpy.loadtxt

library(ggplot2)
library(ggsci)
library(tidyverse)
library(zoo)

for (meas in list.files('measurements/')) {
  print(paste('File: ', meas))
  print(paste('File: ', meas))
  print(paste('File: ', meas))
  print(paste('File: ', meas))
  print(paste('File: ', meas))
  print(paste('File: ', meas))
  
  emg_grip_dat_raw <- read_csv(paste0('measurements/', meas))
  
  emg_grip_dat_raw <- emg_grip_dat_raw %>% transmute(
    grip_force = `/gdx/grip_force_stream/grip_force`,
    grip_force_t = `/gdx/grip_force_stream/header/stamp`,
    emg_ch_1 = `/shimmer0/emg_stream/emg_ch1`,
    emg_ch_1_t = `/shimmer0/emg_stream/header/stamp`
    ) 
  grip_dat <- emg_grip_dat_raw %>% 
    select(grip_force, grip_force_t) %>%
    rename(measure_time = grip_force_t) %>% 
    drop_na()
  emg_dat <- emg_grip_dat_raw %>% 
    select(emg_ch_1, emg_ch_1_t) %>%
    rename(measure_time = emg_ch_1_t) %>% 
    drop_na()
  
  emg_grip_dat <- full_join(
    grip_dat, emg_dat, by = 'measure_time'
      ) %>% arrange(measure_time)
  
  emg_grip_dat$measure_time[emg_grip_dat$measure_time %in% unique(emg_grip_dat$measure_time[duplicated(emg_grip_dat$measure_time)])]
  duplicated_time = unique(emg_grip_dat$measure_time[duplicated(emg_grip_dat$measure_time)])
  
  test1 = zoo(emg_grip_dat$grip_force, emg_grip_dat$measure_time)
  test1 <- na.approx(test1, na.rm = F)
  emg_grip_dat$grip_force_inter <- coredata(test1)
  
  # Export data
  export_file_name <- str_sub(meas, end = -5)
  emg_grip_dat %>% select(-grip_force) %>% drop_na() %>% 
    write_csv(paste0('new_data/', export_file_name, '_sl.csv'))
}

meas <- list.files('measurements/')
str_sub(meas[1], end = -5)
paste0('new_data/', str_sub(meas[1], end = -5), '_sl.csv')

emg_grip_dat_raw %>% str()

emg_grip_dat %>% 
  mutate(
    test = (measure_time == time(test1))
  ) %>%
  filter(test == F)

plot(test1)

emg_grip_dat %>% 
  select(-grip_force) %>% 
  filter(measure_time >= 1670925211 & measure_time <= 1670925211.05) %>%
  pivot_longer(
    c(emg_ch_1, grip_force_inter), 
    names_to = 'instrument', 
    values_to = 'measured_value'
    ) %>%
  drop_na() %>%
  arrange(measure_time) %>% 
  ggplot(aes(y = measured_value, x = measure_time)) +
  facet_wrap(
    . ~ instrument,
    scales = 'free_y',
    strip.position = 'left',
    ncol = 1
    ) +
  geom_line(
    aes(color = instrument),
    # color = 'black',
    # stat = 'smooth',
    stat = 'identity',
    # formula = y ~ s(x, bs = "cs"),
    # method = 'gam',
    # se = T,
    alpha = 0.2,
    size = 0.8,
    # method = 'loess',
    # formula = y ~ x, span = 0.5
  ) +
  geom_point(
    aes(fill = instrument),
    size = 0.6,
    shape = 21,
    stroke = 0.0
  ) +
  scale_y_continuous(
    # name = 'Euclidean distance (error) [cm]',
    # breaks = seq(0, 1, 0.1),
    # limits = c(0, NA)
  ) +
  scale_x_continuous(
    name = expression(time*' [s]'),
    breaks = seq(1657796000, 1657796090, 0.01),
    expand = c(0.02, 0.02)
    # limits = c(2.5, -2.5),
    ) +
  scale_fill_nejm(
    guide = 'none'
    # name = expression('Joint')
  ) +
  scale_color_nejm(
    name = expression('Joint')
  ) +
  # ggtitle(expression('Model comparisons - real vs zero '*italic(beta))) +
  theme_bw() + theme(
    panel.border = element_rect(),
    panel.background = element_blank(),
    panel.spacing = unit(2, 'mm'),
    legend.position = 'top',
    legend.box.spacing = unit(2, 'mm'),
    legend.spacing = unit(2, 'mm'),
    legend.margin = margin(0.5, 0.5, 1, 0.5, 'mm'),
    legend.title = element_text(face = 'bold', size = 10),
    axis.title = element_text(
      size = 10,
      margin = margin(0, 0, 0, 0, 'mm')
      ),
    axis.text.x = element_text(
      colour="black", size = 8,
      margin = margin(0, 0, 0, 0, 'mm')
      ),
    axis.text.y = element_text(colour="black", size = 8),
    axis.title.y = element_blank(),
    axis.line = element_line(size=0.5, colour = "black"),
    plot.title = element_text(
      hjust = 0.5, vjust = 0, size = 12,
      margin = margin(0, 0, 0.5, 0, 'mm')
      ),
    plot.margin = margin(0, 1, 0, 1, 'mm'),
  )

ccf(
  emg_grip_dat %>% select(-grip_force) %>% drop_na() %>% pull(grip_force_inter),
  emg_grip_dat %>% select(-grip_force) %>% drop_na() %>% pull(emg_ch_1),
  lag.max = 1000
  )


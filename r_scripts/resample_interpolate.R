# Resample raw dynamometer data (200 Hz to 1000 Hz) and interpolate to match EMG data timestamps
library(tidyverse)
library(janitor)
library(magrittr)
library(ggsci)
library(zoo)

# Change option digits.secs to show 6 decimals for seconds (microseconds)
options(digits.secs = 6)
getOption("digits.secs")

# Get initials, measuring position, repetition, age and gender from filename
emg_grip_dat_raw <- list.files(
  path = "measurements_june_2024/csv", pattern = "*_m.csv", full.names = T
  ) %>% 
  as_tibble_col(column_name = "path") %>% 
  mutate(
    filename = str_sub(path, start = 28, end = -5),
    initials = str_sub(filename, start = 1, end = 2),
    position = str_sub(filename, start = 4, end = 4) %>% as.factor(),
    repetition = str_sub(filename, start = 6, end = 6) %>% as.factor(),
    age = str_sub(filename, start = 8, end = 9) %>% as.integer(),
    gender = str_sub(filename, start = 11, end = 11)
  ) %>% 
  mutate(
    # Reads correctly date as double with all decimals, but doesn't print it
    data = map(path, \(x) read_csv(x, col_types = 'dddddddd'))
  ) %>% 
  unnest(data) %>% 
  select(
    path,
    filename,
    initials,
    position,
    repetition,
    age,
    gender,
    grip_force = `/gdx/grip_force_stream/grip_force`,
    grip_ts = `/gdx/grip_force_stream/header/stamp`,
    emg = `/shimmer0/emg_stream/emg_ch1`,
    emg_ts = `/shimmer0/emg_stream/header/stamp`
  ) %>% 
  mutate(
    # Convert grip_ts and emg_ts in seconds since epoch to POSIXct format
    # The resolution is a bit better then microsecond, enough for application
    grip_ts = as.POSIXct(grip_ts),
    emg_ts = as.POSIXct(emg_ts)
  )

# If rounded to 6 decimals, there is one duplicate entry in the grip_ts column
# Drop one duplicate entry, both are actually 0 N measurements
emg_grip_dat_raw %>% 
  drop_na(grip_ts) %>% 
  get_dupes(grip_ts) %>% 
  select(grip_ts, grip_force) %>% 
  pull(grip_ts) %>% as.numeric() %>% sprintf("%.10f", .)
#     grip_ts                                               grip_force
#     <dttm>                                                   <dbl>
# 1  2024-06-11 10:20:29.346070 (1718094029.3460702896)   -0.0268351060 
# 2  2024-06-11 10:20:29.346070 (1718094029.3460702896)   -0.0091567760

emg_grip_dat_raw <- emg_grip_dat_raw %>%
  filter(
    # Due to something weird, ts entries with NA also get dropped, so include them in or statement
    !(grip_ts == as.POSIXct(1718094029.3460702896) & grip_force == -0.0268351060) | is.na(grip_ts)
    )

# Mark each non NA timestamp in grip_ts with number
emg_grip_dat_raw <- emg_grip_dat_raw %>% 
  group_by(initials, position, repetition) %>% 
  mutate(
    group_grip = as.integer(!is.na(grip_ts)) %>% cumsum()
    ) %>% 
  # Filter out group_grip == 0 since, we need first timestamp to interpolate
  # Data should already be arranged by timestamp
  filter(group_grip != 0)

emg_grip_dat_raw %>% ungroup() %>% str()
  
# Split dataframe to emg and grip data, and join to match ts in grip_ts
emg_dat <- emg_grip_dat_raw %>% 
  ungroup() %>% 
  # Remove grip and grouping data
  select(-c(grip_force, grip_ts, group_grip)) %>%
  drop_na()

emg_dat %>% str()

# Join emg_dat to grip data
emg_grip_dat_aligned <- emg_grip_dat_raw %>%
  ungroup() %>% 
  # Remove grouping column and emg data and drop NA values
  select(-group_grip, -emg, -emg_ts) %>%
  drop_na() %>%
  # Join emg data to grip data using timestamp
  full_join(
    emg_dat %>% rename(grip_ts = emg_ts),
    by = c(
      'path', 'filename', 'initials', 'position', 'repetition', 'age', 'gender', 
      'grip_ts'
      )
    ) %>%
  arrange(grip_ts)
  
emg_grip_dat_aligned %>% View()

# Interpolate grip_force values using grip_ts and grip_force
emg_grip_dat_complete <- emg_grip_dat_aligned %>%
   group_by(initials, position, repetition) %>%
   mutate(
     # Create zoo object with grip_force and grip_ts
     # Interpolate NA values in grip_force using grip_ts (jointd with emg_ts)
     # Extract too core data
     grip_force_interp = coredata(
       na.approx(
         zoo(x = grip_force, order.by = grip_ts),
         na.rm = F
         ))
     ) %>%
   ungroup() %>%
   # Drop observations with NA values either for emg (original dynamometer)
   # or interpolated grip_force_interp (NA values at the end)
   drop_na(emg, grip_force_interp) %>%
   # Drop grip_force column, rename grip_force_interp to grip_force and grip_ts
   # to measure_ts
   select(-grip_force) %>%
   rename(
     grip_force = grip_force_interp,
     measure_ts = grip_ts
     )

# Save complete emg data
emg_grip_dat_complete %>% 
  write_csv('measurements_june_2024/csv/emg_grip_dat_complete.csv')
# Compress the file
zip(
  'measurements_june_2024/csv/emg_grip_dat_complete.zip', 
  'measurements_june_2024/csv/emg_grip_dat_complete.csv'
  )

# Panel plot all the data
emg_grip_dat_complete %>%
  # Filter for one plot
  filter(initials == 'md', position == 3, repetition == 1) %>%
  # Filter only position 3/4
  filter(position == 3) %>%
  # pivot longer using emg and grip_force
  # rename(grip = grip_force) %>%
  rename(value = grip_force) %>%
  # pivot_longer(
  #   cols = c(emg, grip),
  #   names_to = "measure",
  #   values_to = "value"
  #   ) %>%
  # Set each id by pasting initials, position and repetition
  # mutate(id = paste(initials, position, repetition, measure, sep = '-')) %>%
  mutate(id = paste(initials, position, repetition, sep = '-')) %>%
  ggplot(aes(x = measure_ts, y = value, color = initials)) +
  # ggplot(aes(x = measure_ts, y = emg, color = initials)) +
  facet_wrap(~ id, scales = 'free') +
  geom_line(linewidth = 0.4) +
  scale_color_d3(
    palette = "category20b",
  ) +
  labs(
    x = "Timestamp [s]",
    # y = "EMG [mV] / Grip force [N]",
    y = "Grip force [N]",
    # y = "EMG [mV]",
    # title = "Position 4 - EMG and Grip Force Data"
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
    panel.spacing.x = unit(0, 'mm'),
    axis.title = element_text(face="bold", size = 8),
    axis.text = element_text(
      color="black", size = 8, margin = margin(0.0, 0.0, 0.0, 0.0, 'mm')
      ),
    # Remove x axis text and title
    axis.line = element_line(linewidth = 0.1, colour = "black"),
    plot.title = element_text(hjust = 0.5, size = 5),
    panel.grid = element_line(colour = 'grey', size = 0.1),
    panel.border = element_rect(linewidth = 0.1),
    strip.background = element_rect(linewidth = 0.01),
    strip.text = element_text(
      colour = 'black', size = 8.0, margin = margin(b = 0.3, t = 0.3, unit='mm')
    ),
    axis.ticks = element_line(linewidth = 0.1),
    axis.ticks.y = element_blank(),
    axis.ticks.length = unit(0.1, 'lines')
  )
ggsave(
  # 'plots/emg_grip_raw_position_3.png',
  'plots/emg_grip_raw_position_3_mdp3m1_24_m-grip.png',
  # 'plots/emg_grip_raw_position_3_mdp3m1_24_m-emg.png',
  # width = 18, height = 10, units = 'cm', dpi = 320
  width = 3.5, height = 2, units = 'in', dpi = 320
  )

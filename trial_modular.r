library(devtools)
load_all()
library(tidyverse)
library(concaveman)
library(sf)
library(RColorBrewer)
library(colorspace)
library(ggpubr)
library(ggsci)
library(patchwork)
library(rstatix)
library(ggrepel)
library(forcats)
library(gridExtra)
library(scales)
library(ggplot2)

BIDS_DIR = "/Users/edwardclarkson/Downloads/hypercat_slices"

bids = Bids$new(BIDS_DIR)

CUSTOM_THEME <- theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())

bids$print_tree()

fs_labels <- read_fs_labels("FreeSurferColorLUT.txt")

raw_slices <- bids$get_slices(scan_pattern = "reorient", folder_pattern = "aMRI") %>% 
  filter(!grepl("first|brain", scan_name))

padded_slices <- bids$get_slices(scan_pattern = "resizePad", folder_pattern = "aMRI") %>%
  filter(!grepl("mask|brain", scan_name))

motion_slices <- bids$get_slices(scan_pattern = "motion", folder_pattern = "aMRI")

parcel_slices <- bids$get_slices(scan_pattern = "parc", folder_pattern = "aMRI")

# Create scan mapping
session_labels <- unique(raw_slices$session_id)
motion_labels <- str_split(session_labels, "_", simplify = TRUE)[,1:2]  %>% 
  apply(1, function(x) paste(x[1], x[2], sep = " ")) %>%
  str_replace_all(c("HB" = "HK"))

scan_mapping <- data.frame(
    session_id = session_labels,
    scan_label = motion_labels
  ) %>% 
  filter(grepl("12", scan_label)) %>%
  mutate(scan_label = str_replace_all(scan_label, "12", "1.2mm"))

padded_data <- padded_slices %>%
  right_join(scan_mapping, by = "session_id") %>%
  mutate(slice_data = map(data, slice_to_dataframe)) %>%
  unnest(slice_data) %>%
  filter(time %% 3 == 0) %>%
  mutate(time = factor(time, labels = paste("Timepoint", sort(unique(time))))) %>%
  mutate(slice_type = str_to_title(slice_type)) %>%
  mutate(x = ifelse(slice_type == "Axial", -x, x)) %>% 
  mutate(x = ifelse(slice_type == "Axial", x + 256, x))

## 1. Raw Reorient Plot
raw_reorient_plot <- ggplot() +
  scalar_layer(
    raw_slices %>%
  group_by(subject_id, session_id, scan_name) %>%
  mutate(slice_data = map(data, slice_to_dataframe)) %>%
  unnest(slice_data) %>%
  filter(time == 1) %>%
  right_join(scan_mapping, by = "session_id")) +
  facet_grid(slice_type ~ scan_label) +
  scale_fill_gradient(low = "black", high = "white", trans = "sqrt") +
  ggplot2::coord_fixed(ratio = 1) +
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())+
  labs(x = "", y = "", fill = "Signal")
x11()
raw_reorient_plot
ggsave("raw_reorient_plot_hyperkat.png", raw_reorient_plot, width = 10, height = 10)

# motion_data <- motion_slices %>% 
#   filter(session_id == "2024ANON4295Se4") %>%
#   filter(scan_name == "desc-motion_map") %>%
#   mutate(slice_data = map(data, slice_to_dataframe)) %>% 
#   unnest(slice_data) %>%
#   downsample_vector_data(8) %>%
#   filter(time %% 3 == 0) %>%
#   process_vector_data() %>% 
#   mutate(time = factor(time, labels = paste("Timepoint", sort(unique(time))))) %>% 
#   mutate(slice_type = str_to_title(slice_type))

# ## Figure 1: Timeline plot
# plot <- ggplot() + 
#   scalar_layer(padded_data, alpha = 0.3) +
#   scale_fill_gradient(low = "black", high = "white") +
#   vector_layer(motion_data, amplification_factor = 200) +
#   scale_linewidth_continuous(range = c(0, 0.5)) +
#   scale_color_viridis_c(option = "inferno") +
#   ggplot2::coord_fixed(ratio = 1) +
#   CUSTOM_THEME +
#   facet_grid(slice_type ~ time) +
#   labs(x = "", y = "", fill = "Signal", color = "Displacement") +
#   guides(linewidth = "none", fill = "none")

# ggsave("motion_slices_facet_hyperkat.png", plot, width = 10, height = 10)

## Figure 2: Raw vs Motion Slices

raw_slice_data <- padded_slices %>%
  right_join(scan_mapping, by = "session_id") %>%
  mutate(slice_data = map(data, slice_to_dataframe)) %>%
  unnest(slice_data) %>%
  mutate(slice_type = str_to_title(slice_type))

motion_slice_data <- motion_slices %>% 
  filter(grepl("desc-motion_map", scan_name)) %>%
  right_join(scan_mapping, by = "session_id") %>%
  mutate(slice_data = map(data, slice_to_dataframe)) %>%
  unnest(slice_data) %>% 
  downsample_vector_data(8) %>%
  filter(time == 3) %>%
  process_vector_data() %>% 
  mutate(time = factor(time, labels = paste("Timepoint", sort(unique(time))))) %>% 
  mutate(slice_type = str_to_title(slice_type))

plot <- ggplot() +
  scalar_layer(raw_slice_data, alpha = 0.3) +
  vector_layer(motion_slice_data, amplification_factor = 200) +
  scale_linewidth_continuous(range = c(0, 0.5)) +
  scale_color_viridis_c(option = "inferno", begin = 0.3) +
  scale_fill_gradient(low = "black", high = "white") +
  facet_grid(slice_type ~ scan_label) +
  CUSTOM_THEME +
  labs(x = "", y = "", fill = "Signal")+
  ggplot2::coord_fixed(ratio = 1)

ggsave("raw_vs_motion_slices_all_hyperkat.png", plot, width = 15, height = 10)

## Figure 3: Standard deviation in raw slices

raw_sd_data <- padded_slices %>%
  right_join(scan_mapping, by = "session_id") %>%
  mutate(data = map(data, slice_temporal_robust_sd)) %>%
  group_by(subject_id, session_id, scan_name) %>%
  mutate(slice_type = toupper(slice_type)) %>%
  mutate(slice_data = map(data, slice_to_dataframe)) %>%
  unnest(slice_data)

plot <- ggplot() + 
  scalar_layer(raw_sd_data %>% filter(value > 0)) +
  ggplot2::coord_fixed(ratio = 1) +
  facet_grid(slice_type ~ scan_label) +
  scale_fill_viridis_c(option = "inferno", trans = "sqrt") +
  CUSTOM_THEME +
  labs(x = "", y = "", fill = "SD")
ggsave("raw_sd_plot_hyperkat.png", plot, width = 10, height = 10)
 
## Figure 4: How standard deviation within ROIs correlates with motion measurements

# Get parcellation data
parcel_data <- parcel_slices %>%
  left_join(scan_mapping, by = "session_id") %>%
  select(-subject_id, -session_id) %>%
  mutate(slice_data = map(data, slice_to_dataframe)) %>%
  unnest(slice_data) %>%
  mutate(value = as.integer(round(value))) %>%
  left_join(fs_labels %>% filter(parcel_label != "CSF"), by = c("value" = "parcel_id")) %>% 
  filter(value > 0) %>% 
  select(-data) %>% 
  filter(parcel_label != "Csf") %>% 
  mutate(slice_type = str_to_title(slice_type)) %>%
  filter((slice_type == "Sagittal" & !grepl("Cortex|Cerebral|Left|Ventricle|Ventral|Accumbens", parcel_label)) | (slice_type == "Axial" & !grepl("Putamen|Thalamus|White|Caudate|Left", parcel_label))) %>%
  filter(slice_type != "Coronal") %>% 
  mutate(x = ifelse(slice_type == "Axial", -x, x)) %>%
  mutate(x = ifelse(slice_type == "Axial", x + 256, x))


parcel_plot <- ggplot() +
  scalar_layer(padded_data %>% 
    filter(slice_type != "Coronal") %>% 
    filter(scan_label == "1.2mm HK2") %>% 
    filter(x > 20 & x < 230 & y > 20 & y < 240)) +
  scale_fill_gradient(low = "black", high = "white") +
  point_layer(parcel_data,
    variable = "parcel_label",
    size = 0.2
  ) +
  coord_fixed(ratio = 1) +
  facet_grid(slice_type ~ .) +
  scale_color_brewer(palette = "Set1") +
  guides(color = guide_legend(override.aes = list(size = 3)), fill = "none") +
  CUSTOM_THEME +
  theme_void() +
  theme(legend.position = "none", strip.text.y = element_blank())+
  labs(color = "ROI")+
  labs(x = "", y = "")
ggsave("parcel_outline_hyperkat.png", parcel_plot, width = 5, height = 10)
   
## Motion plot

multiple_motion_data_vector <- motion_slices %>%
  filter(grepl("desc-motion_map", scan_name)) %>%
  right_join(scan_mapping, by = "session_id") %>%
  select(-subject_id, -session_id) %>%
  mutate(slice_data = map(data, slice_to_dataframe)) %>%
  select(-data) %>%
  unnest(slice_data) %>%
  process_vector_data()

multiple_motion_data <- multiple_motion_data_vector %>%
  mutate(slice_type = str_to_title(slice_type)) %>%
  right_join(
    parcel_data %>% select(-scan_name, -scan_label),
    by = c("slice_type", "x", "y")
  ) %>%
  filter((slice_type == "Sagittal" & !grepl("Cortex|Cerebral|Left|Ventricle|Ventral|Accumbens|Putamen", parcel_label)) | (slice_type == "Axial" & !grepl("Putamen|Thalamus|White|Caudate|Left", parcel_label)))

motion_plot_data <- multiple_motion_data %>%
  group_by(scan_label, slice_type, parcel_label, time) %>%
  summarise(
    mean_displacement = mean(total_displacement, na.rm = TRUE),
    sd_displacement = sd(total_displacement, na.rm = TRUE),
    n_points = n(),
    .groups = "drop"
  )

parcel_motion_plot <- motion_plot_data %>%
  mutate(parcel_label = str_replace(parcel_label, "White-Matter", "WM")) %>%
  drop_na(parcel_label, scan_label) %>%
  ggplot(aes(x = time, y = mean_displacement, color = parcel_label, fill = parcel_label)) +
  geom_line(linewidth = 0.2, color = 'slategray')+
  geom_ribbon(aes(ymin = mean_displacement - sd_displacement, ymax = mean_displacement + sd_displacement), color = 'slategray', alpha = 0.4, lty = 2,linewidth = 0.2)+
  facet_grid(parcel_label~scan_label, scales = "free_y")+
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Timepoint", y = "Average Displacement   (mm)") +
  theme_bw()+
  theme(legend.position = "none")
ggsave("parcel_motion_plot_hyperkat.png",parcel_motion_plot, width = 10, height = 10)
 

# Create displacement comparison table relative to Reference
displacement_comparison <- multiple_motion_data %>%
  mutate(parcel_label = str_replace(parcel_label, "White-Matter", "WM")) %>%
  mutate(hyperkat = as.numeric(str_extract(scan_label, "(?<=HK)\\d+"))) %>%
  group_by(slice_type, parcel_label, hyperkat) %>%
  summarise(
    mean_displacement = mean(total_displacement, na.rm = TRUE),
    sd_displacement = sd(total_displacement, na.rm = TRUE),
    .groups = "drop"
  )

hyperkat_displacement_plot <- displacement_comparison %>%
  ggplot(aes(x = hyperkat, y = mean_displacement, color = parcel_label)) +
  geom_errorbar(aes(ymin = mean_displacement - sd_displacement, ymax = mean_displacement + sd_displacement), width = 0.2) +
  geom_point(color = "slategray") +
  facet_grid(parcel_label ~ ., scales = "free_y") +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  labs(x = "Hyperkat", y = "Average Displacement (mm)")
ggsave("hyperkat_displacement_plot.png", hyperkat_displacement_plot, width = 10, height = 10)

combined_parcel_plot <- parcel_plot + parcel_motion_plot + hyperkat_displacement_plot + plot_layout(ncol = 3)
ggsave("combined_parcel_plot_hyperkat.png", combined_parcel_plot, width = 20, height = 9)





parcel_plot <- ggplot() +
  scalar_layer(padded_data %>% 
    filter(slice_type != "Coronal") %>% 
    # filter(scan_label == "1.2mm HK2") %>% 
    filter(x > 20 & x < 230 & y > 20 & y < 240)) +
  scale_fill_gradient(low = "black", high = "white") +
  point_layer(padded_data %>%
  left_join(parcel_data, by = c("slice_type", "x", "y")),
    variable = "parcel_label",
    size = 0.2
  ) +
  coord_fixed(ratio = 1) +
  facet_grid(slice_type ~ scan_label) +
  scale_color_brewer(palette = "Set1") +
  guides(color = guide_legend(override.aes = list(size = 3)), fill = "none") +
  CUSTOM_THEME +
  theme_void() +
  theme(legend.position = "none", strip.text.y = element_blank())+
  labs(color = "ROI")+
  labs(x = "", y = "")
ggsave("parcel_outlinegrid_hyperkat.png", parcel_plot, width = 15, height = 10)


# Motion within the cortex across hyperkat factors

# Step 1: Get unique voxel coordinates and sample them
sampled_voxels <- multiple_motion_data %>%
  # filter(grepl("Right-Cerebral-Cortex", parcel_label)) %>%
  distinct(scan_label,  parcel_label, slice_type, x, y) %>%
  group_by(scan_label, parcel_label, slice_type) %>%
  slice_sample(n = 100) %>%
  ungroup()

spatial_color_map <- function(x, y, sat_boost = 1.2, val_range = c(0.6, 1)) {
  x_norm <- (x - min(x)) / (max(x) - min(x))
  y_norm <- (y - min(y)) / (max(y) - min(y))
  dx <- x_norm - 0.5
  dy <- y_norm - 0.5
  angle <- atan2(dy, dx)
  hue <- (angle + pi) / (2 * pi)
  radius <- sqrt(dx^2 + dy^2)
  radius <- pmin(radius * sat_boost, 1)
  value <- scales::rescale(radius, to = val_range)
  hsv(h = hue, s = radius, v = value)
}

## HSV color map regions

colored_points <- padded_data %>% select(-data) %>%
    filter(scan_label == "1.2mm HK2") %>%
    right_join(parcel_data %>% select(-scan_name, -scan_label), by = c("slice_type", "x", "y")) %>%
    group_by(slice_type, scan_label, parcel_label) %>%
    group_modify(~ {
    .x %>%
      mutate(color = spatial_color_map(x, y))
    })

parcel_color_plot <- ggplot() +
  scalar_layer(padded_data %>% select(-data) %>% 
    filter(slice_type != "Coronal") %>% 
    filter(scan_label == "1.2mm HK2") %>%
    filter(x > 20 & x < 230 & y > 20 & y < 240)) +
  scale_fill_gradient(low = "black", high = "white") +
  point_layer(colored_points,
    variable = "color",
    size = 0.2
  ) +
  scale_color_identity() +
  coord_fixed(ratio = 1) +
  facet_grid(slice_type ~ scan_label) +
  guides(color = guide_legend(override.aes = list(size = 3)), fill = "none") +
  CUSTOM_THEME +
  theme_void() +
  theme(legend.position = "none", strip.text.y = element_blank())+
  labs(color = "ROI")+
  labs(x = "", y = "")
ggsave("parcel_color_hyperkat.png", parcel_plot, width = 15, height = 10)

## HSV color map motion 

cortex_motion_data <- multiple_motion_data %>%
  semi_join(sampled_voxels, by = c("scan_label", "slice_type", "x", "y")) %>%
  group_by(scan_label, parcel_label, slice_type) %>%
  drop_na() %>%
  group_modify(~ {
    .x %>%
      mutate(color = spatial_color_map(x, y))
  }) %>%
  ungroup()

x11()
motion_color_plot <- cortex_motion_data %>%
  ggplot(aes(x = time, y = total_displacement, color = color, group = interaction(scan_label, slice_type, x, y))) +
  geom_line(linewidth = 0.2, alpha = 0.5) +
  scale_color_identity() +
  theme_bw() +
  facet_grid(parcel_label~scan_label, scales = "free_y")
ggsave("motion_color_hyperkat.png", motion_color_plot, width = 10, height = 10)

parcel_color_plot + motion_color_plot + plot_layout(ncol = 2)
ggsave("parcel_color_motion_hyperkat.png", width = 15, height = 10)

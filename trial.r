library(devtools)
load_all()
library(tidyverse)
library(concaveman)
library(sf)
library(RColorBrewer)
library(ggpubr)
library(ggsci)
library(patchwork)
library(rstatix)
library(ggrepel)
library(forcats)
library(gridExtra)

BIDS_DIR = "/Users/edwardclarkson/Downloads/motion_slices"

bids = Bids$new(BIDS_DIR)

CUSTOM_THEME <- theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())

bids$print_tree()

# Create a mapping between session IDs and motion labels
motion_labels <- c("Reference", "Brow Motion", "LR 15 Degrees", "SI Motion", "Unconstrained", "Swallowing")
session_labels <- c("2024ANON4295Se4","2024ANON4295Se6", "2024ANON4295Se8", "2024ANON4295Se10_2", "2024ANON4295Se12", "2024ANON4295Se14")
scan_mapping <- data.frame(
  session_id = session_labels,
  scan_label = factor(motion_labels, levels = motion_labels)
)

fs_labels <- read_fs_labels("FreeSurferColorLUT.txt")


raw_slices <- bids$get_slices(scan_pattern = "reorient", folder_pattern = "aMRI") %>% 
  filter(!grepl("first|brain", scan_name))

padded_slices <- bids$get_slices(scan_pattern = "resizePad", folder_pattern = "aMRI") %>%
  filter(!grepl("mask|brain", scan_name))

motion_slices <- bids$get_slices(scan_pattern = "motion", folder_pattern = "aMRI")
 
parcel_slices <- bids$get_slices(scan_pattern = "padded", folder_pattern = "segmentation")

raw_reorient_plot <- raw_slices %>%
  group_by(subject_id, session_id, scan_name) %>%
  mutate(slice_data = map(data, slice_to_dataframe)) %>%
  unnest(slice_data) %>%
  filter(time == 1) %>%
  left_join(scan_mapping, by = "session_id") %>%
  plot_single_slice() +
  facet_grid(slice_type ~ scan_label) +
  scale_fill_gradient(low = "black", high = "white", trans = "sqrt") +
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())+
  labs(x = "", y = "", fill = "Signal")
x11()
raw_reorient_plot
ggsave("raw_reorient_plot.png", raw_reorient_plot, width = 10, height = 10)

raw_sd_plot <- raw_slices %>%
  mutate(data = map(data, slice_temporal_robust_sd)) %>%
  group_by(subject_id, session_id, scan_name) %>%
  mutate(slice_type = toupper(slice_type)) %>%
  mutate(slice_data = map(data, slice_to_dataframe)) %>%
  unnest(slice_data) %>%
  left_join(scan_mapping, by = "session_id") %>%
  plot_single_slice() +
  facet_grid(slice_type ~ scan_label) +
  scale_fill_viridis_c(option = "inferno", trans = "sqrt") +
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())+
  labs(x = "", y = "", fill = "SD")
ggsave("raw_sd_plot.png", raw_sd_plot, width = 10, height = 10)

padded_data <- padded_slices %>%
  filter(session_id == "2024ANON4295Se4") %>%
  filter(scan_name == "desc-aMRI_resizePad") %>%
  mutate(slice_data = map(data, slice_to_dataframe)) %>%
  unnest(slice_data) %>%
  filter(time %% 3 == 0) %>% 
  mutate(time = factor(time, labels = paste("Timepoint", sort(unique(time))))) %>%
  mutate(slice_type = str_to_title(slice_type))

motion_data <- motion_slices %>% 
  filter(session_id == "2024ANON4295Se4") %>%
  filter(scan_name == "desc-motion_map") %>%
  mutate(slice_data = map(data, slice_to_dataframe)) %>% 
  unnest(slice_data) %>%
  downsample_vector_data(8) %>%
  filter(time %% 3 == 0) %>%
  process_vector_data() %>% 
  mutate(time = factor(time, labels = paste("Timepoint", sort(unique(time))))) %>% 
  mutate(slice_type = str_to_title(slice_type))

## Figure 1: Timeline plot
plot <- ggplot() + 
  scalar_layer(padded_data, alpha = 0.3) +
  scale_fill_gradient(low = "black", high = "white") +
  vector_layer(motion_data, amplification_factor = 200) +
  scale_linewidth_continuous(range = c(0, 0.5)) +
  scale_color_viridis_c(option = "inferno") +
  ggplot2::coord_fixed(ratio = 1) +
  CUSTOM_THEME +
  facet_grid(slice_type ~ time) +
  labs(x = "", y = "", fill = "Signal", color = "Displacement") +
  guides(linewidth = "none", fill = "none")

ggsave("motion_slices_facet.png", plot, width = 10, height = 10)

## Figure 2: Raw vs Motion Slices

raw_slice_data <- padded_slices %>% 
  left_join(scan_mapping, by = "session_id") %>%
  mutate(slice_data = map(data, slice_to_dataframe)) %>%
  unnest(slice_data) %>% 
  mutate(slice_type = str_to_title(slice_type))

motion_slice_data <- motion_slices %>% 
  filter(grepl("desc-motion_map", scan_name)) %>%
  mutate(slice_data = map(data, slice_to_dataframe)) %>%
  unnest(slice_data) %>% 
  downsample_vector_data(8) %>%
  filter(time == 3) %>%
  process_vector_data() %>% 
  mutate(time = factor(time, labels = paste("Timepoint", sort(unique(time))))) %>% 
  mutate(slice_type = str_to_title(slice_type)) %>% 
  left_join(scan_mapping, by = "session_id")

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

ggsave("raw_vs_motion_slices_all.png", plot, width = 15, height = 10)

## Figure 3: Standard deviation in raw slices

raw_sd_data <- padded_slices %>%
  left_join(scan_mapping, by = "session_id") %>%
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
ggsave("raw_sd_plot.png", plot, width = 10, height = 10)
 
## Figure 4: How standard deviation within ROIs correlates with motion measurements

parcel_colors <- c("Reference" = "black", setNames(brewer.pal(8, "Set2")[2:8], unique(scan_mapping$scan_label[scan_mapping$scan_label != "Reference"])))

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
  filter(slice_type != "Coronal")


parcel_plot <- ggplot() +
  scalar_layer(padded_data %>% filter(x>20 & x<230 & y>20 & y<240) %>% filter(slice_type != "Coronal")) +
  scale_fill_gradient(low = "black", high = "white") +
  point_layer(parcel_data %>% filter(scan_label == "Reference"),
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
ggsave("parcel_outline.png", parcel_plot, width = 5, height = 10)
  
## Motion plot

multiple_motion_data_vector <- motion_slices %>%
  left_join(scan_mapping, by = "session_id") %>%
  select(-subject_id, -session_id) %>%
  filter(grepl("desc-motion_map", scan_name)) %>%
  mutate(slice_data = map(data, slice_to_dataframe)) %>% 
  unnest(slice_data) %>%
  process_vector_data()

multiple_motion_data <- multiple_motion_data_vector %>%
  mutate(slice_type = str_to_title(slice_type)) %>%
  right_join(
    parcel_data %>% select(-scan_name),
    by = c("scan_label", "slice_type", "x", "y")
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
ggsave("parcel_motion_plot.png",parcel_motion_plot, width = 10, height = 10)
 

# Create displacement comparison table relative to Reference
displacement_comparison <- multiple_motion_data %>%
  mutate(parcel_label = str_replace(parcel_label, "White-Matter", "WM")) %>%
  group_by(slice_type, parcel_label, scan_label) %>%
  summarise(
    mean_displacement = mean(total_displacement, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  drop_na(parcel_label, scan_label) %>%
  # Get the reference displacement for each parcel
  group_by(slice_type, parcel_label) %>%
  mutate(
    reference_displacement = mean_displacement[scan_label == "Reference"],
    displacement_difference = mean_displacement - reference_displacement,
    percent_increase = (displacement_difference / reference_displacement) * 100
  ) %>%
  ungroup() %>%
  filter(scan_label != "Reference") %>%  # Remove reference rows since difference would be 0
  select(parcel_label, scan_label, displacement_difference, percent_increase) %>%
  arrange(scan_label, desc(percent_increase))

# Reshape for faceting
long_displacement <- displacement_comparison %>%
  select(parcel_label, scan_label, displacement_difference, percent_increase) %>%
  pivot_longer(
    cols = c(displacement_difference, percent_increase),
    names_to = "metric_type",
    values_to = "value"
  ) %>%
  mutate(
    metric_type = recode(metric_type,
                         displacement_difference = "Absolute Displacement Difference (mm)",
                         percent_increase = "Percentage Increase from Reference (%)")
  )

# Create unified plot with faceting
combined_displacement_plot <- long_displacement %>%
  filter(metric_type != "Absolute Displacement Difference (mm)") %>%
  ggplot(aes(y = fct_rev(parcel_label), x = value, fill = scan_label)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text_repel(aes(label = scan_label, color = scan_label),
                  position = position_dodge(width = 0.9),
                  direction = "y",
                  hjust = -0.5,
                  size = 3.5,
                  segment.size = 0.2,
                  box.padding = 0.2,
                  data = long_displacement %>% filter(parcel_label == "Brain-Stem" & metric_type == "Percentage Increase from Reference (%)")) +
  facet_grid(parcel_label ~ ., scales = "free_y") +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(
    x = "Percentage Increase in Average Displacement from Reference (%)",
    y = "",
    fill = "Motion Type"
  )

combined_parcel_plot <- parcel_plot + parcel_motion_plot + combined_displacement_plot + plot_layout(ncol = 3)
ggsave("combined_parcel_plot.png", combined_parcel_plot, width = 20, height = 9)



library(tidyverse)
library(devtools)
document()

BIDS_DIR = "/eresearch/qamri-mtbi/ecla535/BIDS_holly_motion/"

CUSTOM_THEME <- theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())

# Create a mapping between session IDs and motion labels
motion_labels <- c("Reference", "Brow Motion", "LR 15 Degrees", "SI Motion", "Unconstrained", "Swallowing")
session_labels <- c("2024ANON4295Se4","2024ANON4295Se6", "2024ANON4295Se8", "2024ANON4295Se10_2", "2024ANON4295Se12", "2024ANON4295Se14")
scan_mapping <- data.frame(
  session_id = session_labels,
  scan_label = factor(motion_labels, levels = motion_labels)
)

fs_labels <- read_fs_labels("FreeSurferColorLUT.txt")

bids = Bids$new(BIDS_DIR)
bids$set_oblique(".*padded.*", ".*segmentation.*", verbose = TRUE)
resized_mask <- bids$get_slices(".*padded.*", ".*segmentation.*", verbose = TRUE)

mask_data <- resized_mask %>%
  group_by(subject_id, session_id, scan_name) %>%
  mutate(data = map(data, crop_2d_slice_with_padding, padding = 10)) %>%
  mutate(slice_data = map(data, slice_to_dataframe)) %>%
  select(-data) %>%
  unnest(slice_data) %>%
  left_join(scan_mapping, by = "session_id")

mask_plot <- ggplot() + 
  point_layer(mask_data %>% mutate(value = as.factor(value)), variable = "value", alpha = 0.5) +
  CUSTOM_THEME +
  facet_grid(scan_label ~ slice_type, scales = "free")
  # coord_fixed(ratio = 1)
ggsave("mask_plot4.png", mask_plot, width = 10, height = 10)

bids$print_tree()

bids$subjects[[1]]$sessions[[1]]$scans[[15]]$path


read_fs_labels <- function(file_path) {
  fs_labels <- read.table(file_path, header = FALSE, sep = "", 
                          comment.char = "#", strip.white = TRUE,
                          colClasses = c("integer", "character", "integer", "integer", "integer", "integer"),
                          col.names = c("parcel_id", "parcel_label", "R", "G", "B", "A"))
  fs_labels %>% clean_fs_labels()
}

clean_fs_labels <- function(fs_labels) {
  if (!all(names(fs_labels) == c("parcel_id", "parcel_label", "R", "G", "B", "A"))) {
    stop("fs_labels must have 6 columns: parcel_id, parcel_label, R, G, B, A")
  }
  fs_labels %>%
    dplyr::filter(!is.na(parcel_id)) %>%
    dplyr::mutate(
      hemisphere = dplyr::case_when(
        grepl("^ctx-lh-", parcel_label) | grepl("^Left-", parcel_label) | grepl("^left", parcel_label, ignore.case = TRUE) ~ "Left",
        grepl("^ctx-rh-", parcel_label) | grepl("^Right-", parcel_label) | grepl("^right", parcel_label, ignore.case = TRUE) ~ "Right",
        grepl("^wm-lh-", parcel_label) ~ "Left",
        grepl("^wm-rh-", parcel_label) ~ "Right",
        TRUE ~ "Unknown"
    ),
    region = dplyr::case_when(
      grepl("^ctx-[lr]h-", parcel_label) ~ gsub("^ctx-[lr]h-", "", parcel_label),
      grepl("^wm-[lr]h-", parcel_label) ~ gsub("^wm-[lr]h-", "", parcel_label),
      grepl("^(Left|Right)-", parcel_label) ~ gsub("^(Left|Right)-", "", parcel_label),
      TRUE ~ parcel_label
    ),
    tissue_type = dplyr::case_when(
      grepl("^ctx-", parcel_label) ~ "Cortical",
      grepl("^wm-", parcel_label) ~ "White Matter",
      grepl("Cortex$", parcel_label) ~ "Cortical",
      grepl("White-Matter", parcel_label) ~ "White Matter",
      grepl("WM", parcel_label) ~ "White Matter",
      TRUE ~ "Subcortical"
    ),
  ) %>%
  dplyr::mutate(
    region = stringr::str_to_title(region),
    parcel_label = stringr::str_to_title(parcel_label)
  ) %>%
  tibble::as_tibble()
}

left_right_parcels <- function(fs_labels, subcortical_included) {
  subcortical_labels <- fs_labels %>%
    dplyr::filter(tissue_type == "Subcortical" & stringr::str_detect(tolower(region), paste(tolower(subcortical_included), collapse = "|")))

  right_side <- subcortical_labels %>%
    dplyr::filter(hemisphere == "Right") %>%
    dplyr::pull(parcel_id)

  left_side <- subcortical_labels %>%
    dplyr::filter(hemisphere == "Left") %>%
    dplyr::pull(parcel_id)

  list(right_side = right_side, left_side = left_side)
}
slice_to_dataframe <- function(slice) {
  dim_names <- switch(as.character(length(dim(slice))),
    "2" = c("x", "y"),
    "3" = c("x", "y", "time"),
    "4" = c("x", "y", "component", "time"),
    stop("Unsupported array dimensionality")
  )
  df <- reshape2::melt(slice, varnames = dim_names, value.name = "value")
  df <- df %>%
    tibble::as_tibble() %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.factor), ~ as.numeric(as.character(.x))))
  df
}

process_vector_data <- function(slice_df) {
  if (!all(c("scan_label", "slice_type", "x", "y", "time", "value") %in% names(slice_df))) {
    stop("slice_df must have 6 columns: scan_label, slice_type, x, y, time, value")
  }
  vector_data <- slice_df %>%
    dplyr::filter(value != 0) %>%
    tidyr::pivot_wider(
      names_from = component,
      values_from = value,
      names_prefix = "displacement_",
      id_cols = c(scan_label, slice_type, x, y, time)
    ) %>%
    dplyr::mutate(
      total_displacement = sqrt(displacement_1^2 + displacement_2^2 + displacement_3^2),
    ) %>%
    dplyr::group_by(slice_type) %>%
    dplyr::group_modify(function(group_data, key) {
      current_slice_type <- key$slice_type
      displacement_info <- get_displacement_components(current_slice_type)
      dplyr::mutate(group_data, !!!displacement_info)
    }) %>%
    dplyr::ungroup() %>%
    dplyr::select(-displacement_1, -displacement_2, -displacement_3)
  vector_data
}

downsample_vector_data <- function(slice_df, factor = 4) {
  slice_df %>%
    dplyr::filter(
      x %% factor == 0,
      y %% factor == 0
    )
}

amplify_vector_data <- function(slice_df, amplification_factor = 200) {
  slice_df %>%
    dplyr::mutate(
      dplyr::across(dplyr::contains("displacement"), ~ . * amplification_factor)
    )
}
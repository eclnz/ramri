DEFAULT_PLOT_PARAMS <- list(
  labs = ggplot2::labs(
    x = "", y = ""
  ),
  thin_factor = 4,
  theme = ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank()
  )
)

scalar_layer <- function(slice_df, alpha = 1) {
  ggplot2::geom_tile(
    data = slice_df,
    mapping = ggplot2::aes(x = x, y = y, fill = value),
    alpha = alpha
  )
}

point_layer <- function(slice_df, variable = "value", alpha = 1, size = 1) {
  ggplot2::geom_point(
    pch=15,
    data = slice_df,
    mapping = ggplot2::aes(x = x, y = y, color = .data[[variable]]),
    alpha = alpha,
    size = size
  )
}

outline_layer <- function(data,
                          label_column = "value",
                          colour = "black",
                          linewidth = 0.5,
                          alpha = 1,
                          ...) {
  label_sym <- rlang::sym(label_column)

  kernel_smooth <- function(x, window_size = 5) {
    n <- length(x)
    if (n < window_size) return(x)

    k <- window_size %/% 2
    weights <- stats::dnorm(-k:k, sd = k/2)
    weights <- weights / sum(weights)

    padded_x <- c(x, x, x)

    smoothed <- stats::filter(padded_x, weights, sides = 2)

    result <- smoothed[(n+1):(2*n)]
    if (any(is.na(result))) result <- x
    result[n] <- result[1]
    result
  }

  existing_cols <- names(data)
  group_cols <- c(label_column, setdiff(existing_cols, c("x", "y", label_column)))

  outlines <- data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::filter(n() >= 4) %>%
    dplyr::group_split() %>%
    purrr::map_dfr(function(group_df) {
      hull_sf <- sf::st_as_sf(group_df, coords = c("x", "y"))
      hull_outline <- concaveman::concaveman(hull_sf)

      coords <- sf::st_coordinates(hull_outline)
      if (nrow(coords) < 4) return(NULL)

      closed_coords <- rbind(coords, coords[1, , drop = FALSE])

      window_size <- min(9, nrow(closed_coords) - 1)
      if (window_size %% 2 == 0) window_size <- window_size - 1

      smoothed_x <- kernel_smooth(closed_coords[, "X"], window_size)
      smoothed_y <- kernel_smooth(closed_coords[, "Y"], window_size)

      result <- tibble::tibble(
        x = smoothed_x,
        y = smoothed_y
      )

      for (col in group_cols) {
        result[[col]] <- unique(group_df[[col]])
      }

      result
    })

  ggplot2::geom_path(
    data = outlines,
    mapping = ggplot2::aes(x = x, y = y, group = !!label_sym, colour = !!label_sym),
    linewidth = linewidth,
    alpha = alpha,
    ...
  )
}

get_displacement_components <- function(slice_type) {
  slice_type <- tolower(slice_type)
  if (slice_type == "axial") {
    list(
      x_displacement = rlang::expr(displacement_1),
      y_displacement = rlang::expr(displacement_2)
    )
  } else if (slice_type == "coronal") {
    list(
      x_displacement = rlang::expr(displacement_1),
      y_displacement = rlang::expr(displacement_3)
    )
  } else if (slice_type == "sagittal") {
    list(
      x_displacement = rlang::expr(displacement_2),
      y_displacement = rlang::expr(displacement_3)
    )
  } else {
    stop("Invalid slice type: '", slice_type, "'. Must be 'axial', 'coronal', or 'sagittal'.")
  }
}

vector_layer <- function(vector_data, amplification_factor = 200) {
  if(!"time" %in% names(vector_data)) {
    warning("Time variable is missing from vector_data")
  }
  if(!"slice_type" %in% names(vector_data)) {
    warning("slice_type variable is missing from vector_data")
  }
  ggarrow::geom_arrow_segment(
    data = vector_data,
    mapping = ggplot2::aes(
      x = x,
      y = y,
      xend = x + x_displacement * amplification_factor,
      yend = y + y_displacement * amplification_factor,
      color = total_displacement,
      linewidth = total_displacement
    )
  )
}

crop_2d_slice_with_padding <- function(img, padding = 0) {
  if (!is.array(img)) stop("Input must be an array")
  nd <- length(dim(img))
  if (nd < 2) stop("Image must be at least 2D")

  spatial_dims <- dim(img)[1:2]

  # Mask of non-zero in 2D space
  nonzero_mask <- apply(img[1:spatial_dims[1], 1:spatial_dims[2], drop = FALSE] != 0, 1:2, any)

  if (!any(nonzero_mask)) return(img)  # Return unchanged if all-zero

  idx <- which(nonzero_mask, arr.ind = TRUE)
  mins <- pmax(apply(idx, 2, min) - padding, 1)
  maxs <- pmin(apply(idx, 2, max) + padding, spatial_dims)

  # Build slicing list: crop first 2 dims, keep rest untouched
  slice_list <- c(
    list(mins[1]:maxs[1], mins[2]:maxs[2]),
    lapply(dim(img)[-(1:2)], seq_len)
  )

  cropped <- do.call(`[`, c(list(img), slice_list, list(drop = FALSE)))
  return(cropped)
}

interpolate_to_isotropic_spacing <- function(img, target_size = c(256, 256)) {
  if (!is.array(img)) stop("Input must be an array")
  nd <- length(dim(img))
  if (nd < 2) stop("Image must be at least 2D")

  spatial_dims <- dim(img)[1:2]
  new_dims <- target_size

  # Generate new coordinates
  new_x <- seq(1, spatial_dims[1], length.out = new_dims[1])
  new_y <- seq(1, spatial_dims[2], length.out = new_dims[2])

  # Setup output array with new dimensions, keeping higher dims unchanged
  output_dims <- c(new_dims, dim(img)[-(1:2)])
  output <- array(0, dim = output_dims)

  # Interpolate each slice across higher dimensions
  higher_dims <- if (nd > 2) expand.grid(lapply(dim(img)[-(1:2)], seq_len)) else data.frame(idx = 1)

  for (i in seq_len(nrow(higher_dims))) {
    idx <- as.integer(higher_dims[i, ])
    if (nd == 2) {
      slice <- img
    } else {
      full_idx <- c(list(quote(expr = )), list(quote(expr = )), as.list(idx))
      slice <- do.call(`[`, c(list(img), full_idx, list(drop = FALSE)))
    }
    slice <- matrix(slice, nrow = spatial_dims[1], ncol = spatial_dims[2])

    # Bilinear interpolation
    interpolated <- outer(
      new_x, new_y,
      Vectorize(function(xi, yi) {
        xi1 <- floor(xi); xi2 <- ceiling(xi)
        yi1 <- floor(yi); yi2 <- ceiling(yi)
        if (xi1 < 1 || xi2 > spatial_dims[1] || yi1 < 1 || yi2 > spatial_dims[2]) return(0)
        Q11 <- slice[xi1, yi1]
        Q12 <- slice[xi1, yi2]
        Q21 <- slice[xi2, yi1]
        Q22 <- slice[xi2, yi2]
        dx <- xi - xi1
        dy <- yi - yi1
        Q11 * (1 - dx) * (1 - dy) + Q21 * dx * (1 - dy) + Q12 * (1 - dx) * dy + Q22 * dx * dy
      })
    )

    # Assign interpolated result
    if (nd == 2) {
      output <- interpolated
    } else {
      assign_idx <- c(list(quote(expr = )), list(quote(expr = )), as.list(idx))
      output <- do.call(`[<-`, c(list(output), assign_idx, list(interpolated)))
    }
  }

  return(output)
}
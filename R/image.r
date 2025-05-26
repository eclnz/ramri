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
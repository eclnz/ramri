arr2slicedf <- function(img) {
  slices <- extract_nonzero_slices(img)
  tibble::tibble(
    slice_type = names(slices),
    data = purrr::map(slices, slice_to_dataframe)
  )
}

extract_nonzero_slices <- function(img, center = NULL) {
  dims <- dim(img)
  if (is.null(center)) {
    center <- floor(dims[1:3] / 2)
  }
  if (length(dims) == 3) {
    sagittal <- img[center[1], , ]
    coronal  <- img[ , center[2], ]
    axial    <- img[ , , center[3]]
  } else if (length(dims) == 4) {
    sagittal <- img[center[1], , , ]
    coronal  <- img[ , center[2], , ]
    axial    <- img[ , , center[3], ]
  } else if (length(dims) == 5) {
    sagittal <- img[center[1], , , , ]
    coronal  <- img[ , center[2], , , ]
    axial    <- img[ , , center[3], , ]
  } else {
    stop("Unsupported number of dimensions")
  }
  list(
    sagittal = sagittal,
    coronal = coronal,
    axial = axial
  )
}


#' Extract slices in canonical anatomical orientation using the affine matrix
extract_oriented_slices <- function(img, rotation = NULL, center = NULL, spacing = 1.0) {
  dims <- dim(img)
  if (length(dims) < 3) stop("Image must be at least 3D")

  if (is.null(center)) {
    center <- floor(dims[1:3] / 2)
  }

  parcellation_extent <- apply(which(img > 0, arr.ind = TRUE), 2, function(x) diff(range(x)))
  slice_size <- ceiling(max(parcellation_extent) * 1.2)

  extract_slice_plane <- function(ax1, ax2, center, img, slice_size, spacing) {
    offset <- (slice_size - 1) / 2
    grid_x <- seq(-offset, offset) * spacing
    grid_y <- seq(-offset, offset) * spacing
    grid <- expand.grid(grid_x, grid_y)

    coords_3d <- t(center + t(grid[,1] %*% t(ax1) + grid[,2] %*% t(ax2)))
    coords_3d <- round(coords_3d)

    dims <- dim(img)
    values <- rep(NA, nrow(coords_3d))
    for (i in seq_len(nrow(coords_3d))) {
      pt <- coords_3d[i, ]
      if (all(pt >= 1 & pt <= dims[1:3])) {
        values[i] <- img[pt[1], pt[2], pt[3]]
      } else {
        values[i] <- 0
      }
    }

    matrix(values, nrow = slice_size, ncol = slice_size, byrow = TRUE) |> t()
  }

  if (is.null(rotation)) {
    slices <- list(
      sagittal = img[center[1], , ],
      coronal  = img[ , center[2], ],
      axial    = img[ , , center[3]]
    )
    return(slices)
  }

  if (is.null(rotation)) stop("Rotation matrix must be provided")

  x_axis <- norm_vector(rotation[, 1])
  y_axis <- norm_vector(rotation[, 2])
  z_axis <- norm_vector(rotation[, 3])

  axial    <- extract_slice_plane(x_axis, y_axis, center, img, slice_size, spacing)
  sagittal <- extract_slice_plane(y_axis, z_axis, center, img, slice_size, spacing)
  coronal  <- extract_slice_plane(x_axis, z_axis, center, img, slice_size, spacing)

  list(
    axial    = axial,
    coronal  = coronal,
    sagittal = sagittal
  )
}
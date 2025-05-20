arr2slicedf <- function(img) {
  slices <- extract_nonzero_slices(img)
  tibble::tibble(
    slice_type = names(slices),
    data = purrr::map(slices, slice_to_dataframe)
  )
}

extract_nonzero_slices <- function(img, threshold = 1) {
  dims <- dim(img)
  indices <- extract_nonzero_slices_cpp(img, dims)
  if (length(dims) == 3) {
    sagittal <- img[indices[1], , ]
    coronal  <- img[ , indices[2], ]
    axial    <- img[ , , indices[3]]
  } else if (length(dims) == 4) {
    sagittal <- img[indices[1], , , ]
    coronal  <- img[ , indices[2], , ]
    axial    <- img[ , , indices[3], ]
  } else if (length(dims) == 5) {
    sagittal <- img[indices[1], , , , ]
    coronal  <- img[ , indices[2], , , ]
    axial    <- img[ , , indices[3], , ]
  } else {
    stop("Unsupported number of dimensions")
  }
  list(
    sagittal = sagittal,
    coronal = coronal,
    axial = axial
  )
}
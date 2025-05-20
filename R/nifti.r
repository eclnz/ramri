load_nifti <- function(filepath) {
  if (!file.exists(filepath)) {
    stop("File does not exist: ", filepath)
  }
  tryCatch({
    volume <- as.array(RNifti::readNifti(filepath))
    return(volume)
  }, error = function(e) {
    stop("Failed to load NIfTI file: ", e$message)
  })
}

nifti_to_slicedf <- function(filepath) {
  volume <- load_nifti(filepath)
  arr2slicedf(volume)
}
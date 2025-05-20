#' @title Compute Principal Direction of a Vector Field
#' @description Computes the principal direction of a vector field using PCA across a 3D array.
#' @param x A 3D array representing vector data with dimensions (m, n, t), where t represents time.
#' @param na.rm A logical indicating whether NA values should be removed. Default is TRUE.
#' @return A 4D array of shape (m, n, 3, 1), where each (m, n) position contains a 3D principal direction vector.
#' @examples
#' set.seed(123)
#' x <- array(stats::rnorm(3 * 4 * 10), dim = c(3, 4, 10)) # Random vector field data
#' result <- vector_principal_direction(x)
#' @importFrom stats prcomp
#' @export
vector_principal_direction <- function(x, na.rm = TRUE) {
  # Store the dimensions of the input array
  dims <- dim(x)

  # Function to compute the principal component for a given position across time
  calc_principal_component <- function(pos_array) {
    tryCatch(
      {
        # Transpose to have time as rows for PCA analysis
        timeseries <- t(pos_array)
        # Perform Principal Component Analysis (PCA)
        pca <- prcomp(timeseries, center = TRUE, scale = FALSE)
        # Extract the first principal component and scale it using the standard deviation
        pca$rotation[, 1] * sqrt(pca$sdev[1]^2)
      },
      error = function(e) c(0, 0, 0)
    ) # Return a zero vector in case of an error
  }

  # Apply the PCA function across spatial dimensions (m, n)
  result <- apply(x, c(1, 2), calc_principal_component)

  # Adjust the array dimensions to maintain the correct shape
  result <- aperm(result, c(2, 3, 1))

  # Return as a 4D array to maintain expected output format
  array(result, dim = c(dims[1], dims[2], 3, 1))
}



#' @title Compute Principal Direction and Standard Deviation
#' @description Computes the principal direction and standard deviation of a vector field using PCA across a 3D array.
#' @param x A 3D array representing vector data with dimensions (m, n, t), where t represents time.
#' @param na.rm A logical indicating whether NA values should be removed. Default is TRUE.
#' @return A list containing two 3D arrays:
#' \item{principal_sd}{A 3D array of shape (m, n, 1), containing the standard deviation of principal directions.}
#' \item{orthogonal_sd}{A 3D array of shape (m, n, 2), containing the standard deviation of orthogonal directions.}
#' @export
vector_principal_direction_with_sd <- function(x, na.rm = TRUE) {
  dims <- dim(x)

  # Inner function to calculate the principal component
  calc_principal_component <- function(pos_array) {
    tryCatch(
      {
        timeseries <- t(pos_array)
        pca <- prcomp(timeseries, center = TRUE, scale = FALSE)
        # Return the first principal component and the rotation matrix
        list(direction = pca$rotation[, 1] * sqrt(pca$sdev[1]^2), rotation = pca$rotation)
      },
      error = function(e) list(direction = c(0, 0, 0), rotation = diag(3))
    ) # Return zero vector and identity matrix on error
  }

  # Calculate the principal direction for each slice
  principal_components <- apply(x, c(1, 2), calc_principal_component)

  # Extract the principal directions and rotation matrices
  principal_directions <- aperm(simplify2array(lapply(principal_components, `[[`, "direction")), c(2, 3, 1))
  rotation_matrices <- aperm(simplify2array(lapply(principal_components, `[[`, "rotation")), c(2, 3, 1, 4))

  # Calculate the standard deviation in the principal direction
  std_dev_principal <- apply(principal_directions, c(1, 2), sd, na.rm = na.rm)

  # Calculate the orthogonal directions using the rotation matrices
  orthogonal_directions <- array(0, dim = c(dims[1], dims[2], 3, 2)) # Prepare array for orthogonal directions
  for (i in 1:(dims[1] * dims[2])) {
    orthogonal_directions[i, , , 1] <- rotation_matrices[i, , , 2] # Second principal component
    orthogonal_directions[i, , , 2] <- rotation_matrices[i, , , 3] # Third principal component
  }

  # Calculate the standard deviation in the orthogonal directions
  std_dev_orthogonal <- apply(orthogonal_directions, c(1, 2, 4), sd, na.rm = na.rm)

  # Return the standard deviations in the same format as principal directions
  std_dev_principal <- array(std_dev_principal, dim = c(dims[1], dims[2], 1))
  std_dev_orthogonal <- array(std_dev_orthogonal, dim = c(dims[1], dims[2], 2))

  return(list(principal_sd = std_dev_principal, orthogonal_sd = std_dev_orthogonal))
}

#' @title Calculate Temporal Standard Deviation for a Slice
#' @description Computes the standard deviation across the time dimension for each voxel in a 3D slice.
#' @param x A 3D array with dimensions (x, y, t), where t represents time.
#' @param na.rm A logical indicating whether NA values should be removed. Default is TRUE.
#' @param mask An optional binary mask of dimensions (x, y) to restrict calculations to specific regions. Default is NULL.
#' @return A 2D array of shape (x, y) containing the standard deviation values for each voxel.
#' @examples
#' set.seed(123)
#' # Create a sample 3D slice with dimensions 10x10x20 (x, y, time)
#' slice <- array(stats::rnorm(10 * 10 * 20), dim = c(10, 10, 20))
#' # Calculate temporal standard deviation
#' sd_map <- slice_temporal_sd(slice)
#' @export
slice_temporal_sd <- function(x, na.rm = TRUE, mask = NULL) {
  # Check input dimensions
  dims <- dim(x)
  if (length(dims) != 3) {
    stop("Input must be a 3D array with dimensions (x, y, t)")
  }
  
  # Apply mask if provided
  if (!is.null(mask)) {
    if (!all(dim(mask) == dims[1:2])) {
      stop("Mask dimensions must match the spatial dimensions of the input array")
    }
    for (t in 1:dims[3]) {
      x[,,t] <- x[,,t] * mask
    }
  }
  
  # Calculate standard deviation along the time dimension (3rd dimension)
  result <- apply(x, c(1, 2), sd, na.rm = na.rm)
  
  # Return 2D array with standard deviations
  return(result)
}

#' @title Calculate Temporal Coefficient of Variation for a Slice
#' @description Computes the coefficient of variation (CV = sd/mean) across the time dimension for each voxel.
#' @param x A 3D array with dimensions (x, y, t), where t represents time.
#' @param na.rm A logical indicating whether NA values should be removed. Default is TRUE.
#' @param mask An optional binary mask of dimensions (x, y) to restrict calculations to specific regions. Default is NULL.
#' @param threshold Minimum mean value to calculate CV; avoids division by small values. Default is 1e-6.
#' @return A 2D array of shape (x, y) containing the coefficient of variation for each voxel.
#' @examples
#' set.seed(123)
#' # Create a sample 3D slice with dimensions 10x10x20 (x, y, time)
#' slice <- array(stats::rnorm(10 * 10 * 20, mean = 10), dim = c(10, 10, 20))
#' # Calculate temporal coefficient of variation
#' cv_map <- slice_temporal_cv(slice)
#' @export
slice_temporal_cv <- function(x, na.rm = TRUE, mask = NULL, threshold = 1e-6) {
  # Check input dimensions
  dims <- dim(x)
  if (length(dims) != 3) {
    stop("Input must be a 3D array with dimensions (x, y, t)")
  }
  
  # Apply mask if provided
  if (!is.null(mask)) {
    if (!all(dim(mask) == dims[1:2])) {
      stop("Mask dimensions must match the spatial dimensions of the input array")
    }
    for (t in 1:dims[3]) {
      x[,,t] <- x[,,t] * mask
    }
  }
  
  # Calculate mean and standard deviation along time dimension
  mean_values <- apply(x, c(1, 2), mean, na.rm = na.rm)
  sd_values <- apply(x, c(1, 2), sd, na.rm = na.rm)
  
  # Calculate coefficient of variation (CV = sd/mean)
  # Apply threshold to avoid division by very small numbers
  cv <- sd_values / pmax(mean_values, threshold)
  
  # If mask provided, set areas outside mask to NA
  if (!is.null(mask)) {
    cv[mask == 0] <- NA
  }
  
  return(cv)
}

#' @title Calculate Robust Temporal Standard Deviation for a Slice
#' @description Computes a robust estimate of standard deviation using median absolute deviation (MAD) across time.
#' @param x A 3D array with dimensions (x, y, t), where t represents time.
#' @param na.rm A logical indicating whether NA values should be removed. Default is TRUE.
#' @param mask An optional binary mask of dimensions (x, y) to restrict calculations to specific regions. Default is NULL.
#' @param scale_factor Scaling factor to make MAD comparable to standard deviation. Default is 1.4826.
#' @return A 2D array of shape (x, y) containing the robust standard deviation values.
#' @examples
#' set.seed(123)
#' # Create a sample 3D slice with dimensions 10x10x20 (x, y, time)
#' slice <- array(rnorm(10 * 10 * 20), dim = c(10, 10, 20))
#' # Add some outliers
#' slice[5, 5, 10] <- 100  # Add an outlier
#' # Calculate robust standard deviation
#' robust_sd <- slice_temporal_robust_sd(slice)
#' @export
slice_temporal_robust_sd <- function(x, na.rm = TRUE, mask = NULL, scale_factor = 1.4826) {
  # Check input dimensions
  dims <- dim(x)
  if (length(dims) != 3) {
    stop("Input must be a 3D array with dimensions (x, y, t)")
  }
  
  # Apply mask if provided
  if (!is.null(mask)) {
    if (!all(dim(mask) == dims[1:2])) {
      stop("Mask dimensions must match the spatial dimensions of the input array")
    }
    for (t in 1:dims[3]) {
      x[,,t] <- x[,,t] * mask
    }
  }
  
  # Function to calculate MAD for a time series
  calc_mad <- function(time_series) {
    if (all(is.na(time_series))) {
      return(NA)
    }
    # Remove NA values if requested
    if (na.rm) {
      time_series <- time_series[!is.na(time_series)]
    }
    # Calculate median absolute deviation and scale it
    mad_value <- stats::mad(time_series, constant = scale_factor, na.rm = na.rm)
    return(mad_value)
  }
  
  # Apply MAD calculation to each voxel time series
  result <- apply(x, c(1, 2), calc_mad)
  
  return(result)
}



# slice.R

# Load necessary libraries
library(RNifti)
library(rgl)
library(misc3d)
library(Rvcg)
library(geometry)
library(dplyr)
library(stringr)
library(pracma)
library(devtools)

# Check if voxels are within specified parcels
is_within_parcels <- function(parc, parcels) {
  result <- array(FALSE, dim = dim(parc))
  for (p in parcels) {
    result <- result | (parc == p)
  }
  result
}

calculate_ellipsoid <- function(voxel_coords) {
  if (nrow(voxel_coords) < 3) {
    stop("Too few points to compute ellipsoid.")
  }

  if (any(!is.finite(voxel_coords))) {
    stop("Non-finite voxel coordinates detected.")
  }

  cov_matrix <- cov(voxel_coords)

  if (any(!is.finite(cov_matrix))) {
    stop("Covariance matrix contains NA, Inf, or NaN.")
  }

  eig <- tryCatch(eigen(cov_matrix), error = function(e) {
    print("Covariance matrix:")
    print(cov_matrix)
    stop("Eigen decomposition failed: ", e$message)
  })

  V <- eig$vectors

  for (i in 1:3) {
    dominant_index <- which.max(abs(V[, i]))
    if (V[dominant_index, i] < 0) {
      V[, i] <- -V[, i]
    }
  }

  if (det(V) < 0) V[,3] <- -V[,3]

  list(
    direction_vectors = V,
    scaling_factors = sqrt(eig$values),
    cov_matrix = cov_matrix
  )
}

normalize <- function(v) {
  norm <- sqrt(sum(v^2))
  if (norm < .Machine$double.eps) {
    warning("Attempted to normalize a near-zero vector.")
    return(rep(NA, length(v)))
  }
  v / norm
}

combine_ellipsoid_orientations <- function(left_center, right_center, left_ellipsoid, right_ellipsoid) {
  x_axis <- normalize(right_center - left_center)

  x1 <- left_ellipsoid$direction_vectors[, 1]
  x2 <- right_ellipsoid$direction_vectors[, 1]

  # Align signs
  if (sum(x1 * x2) < 0) x2 <- -x2
  x_axis <- normalize(x1 + x2)

  y1 <- left_ellipsoid$direction_vectors[, 2]
  y2 <- right_ellipsoid$direction_vectors[, 2]
  if (sum(y1 * y2) < 0) y2 <- -y2
  y_avg <- normalize(y1 + y2)

  y_axis <- y_avg - sum(y_avg * x_axis) * x_axis
  y_axis <- normalize(y_axis)

  z_axis <- normalize(pracma::cross(x_axis, y_axis))
  y_axis <- normalize(pracma::cross(z_axis, x_axis))

  rotation_matrix <- cbind(x_axis, y_axis, z_axis)

  if (abs(det(rotation_matrix) - 1) > 1e-3) {
    warning("Final rotation matrix is not orthonormal.")
  }

  return(rotation_matrix)
}

visualize_rotation_frame <- function(rotation_matrix, origin = c(0, 0, 0), scale = 20, colors = c("red", "green", "blue"), labels = c("X", "Y", "Z")) {
  for (i in 1:3) {
    dir_vec <- rotation_matrix[, i]
    end_point <- origin + dir_vec * scale
    arrow3d(p0 = origin, p1 = end_point, col = colors[i], type = "lines", lwd = 3)
    text3d(end_point, texts = labels[i], col = colors[i], adj = 1.2)
  }
  axes3d()
  box3d()
}

plot_parcel <- function(parcel_voxels, label, color) {
  if (!any(parcel_voxels)) return()
  voxel_coords <- which(parcel_voxels, arr.ind = TRUE)
  if (nrow(voxel_coords) < 3) return()

  hull <- convhulln(voxel_coords)
  mesh <- tmesh3d(vertices = t(voxel_coords), indices = t(hull), homogeneous = FALSE)
  smooth_mesh <- vcgSmooth(mesh, lambda = 0.5)
  shade3d(smooth_mesh, color = color, alpha = 0.2)

  ellipsoid <- calculate_ellipsoid(voxel_coords)
  principal_direction_1 <- ellipsoid$direction_vectors[, 1]
  center <- colMeans(voxel_coords)

  wire3d(ellipse3d(ellipsoid$cov_matrix, centre = center, level = 0.80), col = color, alpha = 0.5)
  arrow3d(center, center + principal_direction_1 * ellipsoid$scaling_factors[1], col = "red", type = "lines", lwd = 2)
}

main <- function(parc_path, lut_path) {
  parc <- readNifti(parc_path)
  fs_labels <- read_fs_labels(lut_path)
  subcortical_labels <- c("thalamus", "caudate", "putamen", "pallidum", "accumbens")
  parcels <- left_right_parcels(fs_labels, subcortical_labels)

  left_mask <- is_within_parcels(parc, parcels$left_side)
  right_mask <- is_within_parcels(parc, parcels$right_side)

  left_coords <- which(left_mask, arr.ind = TRUE)
  right_coords <- which(right_mask, arr.ind = TRUE)

  left_center <- colMeans(left_coords)
  right_center <- colMeans(right_coords)

  left_ellipsoid <- calculate_ellipsoid(left_coords)
  right_ellipsoid <- calculate_ellipsoid(right_coords)

  rotation_matrix <- combine_ellipsoid_orientations(left_center, right_center, left_ellipsoid, right_ellipsoid)

  open3d()
  plot_parcel(left_mask, "Left Subcortical", "lightgreen")
  plot_parcel(right_mask, "Right Subcortical", "lightblue")
  visualize_rotation_frame(rotation_matrix, origin = (left_center + right_center) / 2, scale = 30)

  title3d(main = "Oriented Subcortical Ellipsoids", xlab = "X", ylab = "Y", zlab = "Z")
}

# Example usage:
main(parc_path = files[[6]], lut_path = "FreeSurferColorLUT.txt")
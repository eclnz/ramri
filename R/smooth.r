# slice.R

# Check if voxels are within specified parcels
is_within_parcels <- function(parc, parcels) {
  result <- array(FALSE, dim = dim(parc))
  for (p in parcels) {
    result <- result | (parc == p)
  }
  result
}

calculate_ellipsoid <- function(coords) {
  cov_matrix <- cov(coords)
  eig <- eigen(cov_matrix)
  vectors <- eig$vectors
  for (i in 1:3) {
    dominant_index <- which.max(abs(vectors[, i]))
    if (vectors[dominant_index, i] < 0) {
      vectors[, i] <- -vectors[, i]
    }
  }
  if (det(vectors) < 0) vectors[,3] <- -vectors[,3]
  list(
    direction_vectors = vectors,
    scaling_factors = sqrt(eig$values),
    cov_matrix = cov_matrix
  )
}

norm_vector <- function(v) {
  norm <- sqrt(sum(v^2))
  if (norm < .Machine$double.eps) return(rep(NA, length(v)))
  v / norm
}

combine_ellipsoid_orientations <- function(left_center, right_center, left_ellipsoid, right_ellipsoid) {
  # X axis: average of first principal directions (side-to-side)
  x1 <- left_ellipsoid$direction_vectors[, 1]
  x2 <- right_ellipsoid$direction_vectors[, 1]
  if (sum(x1 * x2) < 0) x2 <- -x2
  x_axis <- norm_vector(x1 + x2)

  # Y axis: average of second principal directions (forward-backward)
  y1 <- left_ellipsoid$direction_vectors[, 2]
  y2 <- right_ellipsoid$direction_vectors[, 2]
  if (sum(y1 * y2) < 0) y2 <- -y2
  y_avg <- norm_vector(y1 + y2)
  y_axis <- norm_vector(y_avg - sum(y_avg * x_axis) * x_axis)

  # Z axis: superior-inferior
  z_axis <- norm_vector(pracma::cross(x_axis, y_axis))
  y_axis <- norm_vector(pracma::cross(z_axis, x_axis))

  rotation_matrix <- cbind(x_axis, y_axis, z_axis)
  if (abs(det(rotation_matrix) - 1) > 1e-3) warning("Rotation matrix is not orthonormal.")

  # Ensure Z (superior) points up and Y (anterior) points forward
  if (rotation_matrix[3, 3] < 0) rotation_matrix[, 3] <- -rotation_matrix[, 3]
  if (rotation_matrix[2, 2] < 0) rotation_matrix[, 2] <- -rotation_matrix[, 2]

  axes <- list(x_axis, y_axis, z_axis)
  axis_scores <- t(sapply(axes, function(a) abs(a)))

  best_x <- which.max(axis_scores[, 1])
  best_y <- which.max(axis_scores[, 2])
  best_z <- which.max(axis_scores[, 3])

  if (length(unique(c(best_x, best_y, best_z))) < 3) {
    warning("Ambiguous axis alignment; unable to assign unique canonical axes.")
  }

  canonical_axes <- list(
    x_axis = axes[[best_x]],
    y_axis = axes[[best_y]],
    z_axis = axes[[best_z]]
  )

  x_axis <- canonical_axes$x_axis
  y_axis <- canonical_axes$y_axis
  z_axis <- canonical_axes$z_axis

  cbind(x_axis, y_axis, z_axis)
}

get_ras_coords <- function(voxel_coords, affine) {
  coords_homogeneous <- cbind(voxel_coords, 1)
  ras_coords <- t(apply(coords_homogeneous, 1, function(v) affine %*% v))[, 1:3]
  ras_coords
}

visualize_ellipsoids_and_frame <- function(left_coords, right_coords, left_ellipsoid, right_ellipsoid, rotation_matrix) {
  open3d()

  left_center <- colMeans(left_coords)
  right_center <- colMeans(right_coords)
  wire3d(ellipse3d(left_ellipsoid$cov_matrix, centre = left_center, level = 0.8), col = "lightgreen", alpha = 0.4)
  wire3d(ellipse3d(right_ellipsoid$cov_matrix, centre = right_center, level = 0.8), col = "lightblue", alpha = 0.4)

  origin <- (left_center + right_center) / 2
  scale <- 30
  colors <- c("red", "green", "blue")
  labels <- c("X", "Y", "Z")
  for (i in 1:3) {
    dir_vec <- rotation_matrix[, i]
    end_point <- origin + dir_vec * scale
    arrow3d(p0 = origin, p1 = end_point, col = colors[i], type = "lines", lwd = 3)
    text3d(end_point, texts = labels[i], col = colors[i], adj = 1.2)
  }
  axes3d()
  box3d()
  title3d(main = "Oriented Ellipsoids and Rotation Matrix", xlab = "X", ylab = "Y", zlab = "Z")
}

compute_combined_rotation <- function(parc, lut_path = "FreeSurferColorLUT.txt") {
  fs_labels <- read_fs_labels(lut_path)
  subcortical_labels <- c("thalamus", "caudate", "putamen", "pallidum", "accumbens")
  parcels <- left_right_parcels(fs_labels, subcortical_labels)

  affine <- RNifti::xform(parc)

  left_mask <- is_within_parcels(parc, parcels$left_side)
  right_mask <- is_within_parcels(parc, parcels$right_side)

  left_coords_vox <- which(left_mask, arr.ind = TRUE)
  right_coords_vox <- which(right_mask, arr.ind = TRUE)

  left_coords <- get_ras_coords(left_coords_vox, affine)
  right_coords <- get_ras_coords(right_coords_vox, affine)
  left_coords_vox <- which(left_mask, arr.ind = TRUE)
  right_coords_vox <- which(right_mask, arr.ind = TRUE)

  left_coords <- get_ras_coords(left_coords_vox, affine)
  right_coords <- get_ras_coords(right_coords_vox, affine)

  left_center <- colMeans(left_coords)
  right_center <- colMeans(right_coords)

  left_ellipsoid <- calculate_ellipsoid(left_coords)
  right_ellipsoid <- calculate_ellipsoid(right_coords)

  rot <- combine_ellipsoid_orientations(left_center, right_center, left_ellipsoid, right_ellipsoid)

  x_axis <- rot[, 1]
  y_axis <- rot[, 2]
  z_axis <- rot[, 3]

  # Enforce RAS orientation convention:
  # X: Right (positive X), Y: Anterior (positive Y), Z: Superior (positive Z)
  if (x_axis[1] < 0) x_axis <- -x_axis
  if (y_axis[2] < 0) y_axis <- -y_axis
  if (z_axis[3] < 0) z_axis <- -z_axis

  # Re-orthonormalize
  y_axis <- norm_vector(pracma::cross(z_axis, x_axis))
  z_axis <- norm_vector(pracma::cross(x_axis, y_axis))

  rotation_matrix <- cbind(x_axis, y_axis, z_axis)
  visualize_ellipsoids_and_frame(left_coords, right_coords, left_ellipsoid, right_ellipsoid, rotation_matrix)

  rotation_matrix
}

main <- function(parc_path, lut_path) {
  parc <- readNifti(parc_path)
  rot <- compute_combined_rotation(parc, lut_path)
  print("Combined rotation matrix (RAS enforced):")
  print(rot)
}

# # Load dependencies
# library(RNifti)
# library(rgl)
# library(misc3d)
# library(Rvcg)
# library(geometry)
# library(dplyr)
# library(stringr)
# library(pracma)

# Usage
# Create a 2x1 subplot layout
# par(mfrow = c(2, 1))

# # # First subject
# path1 = "sub-expANONYMIZED_ses-2024ANON4295Se4_desc-padded_segmentation.nii.gz"
# main(parc_path = path1, lut_path = "FreeSurferColorLUT.txt")

# # # Second subject
# path2 = '/Users/edwardclarkson/git/qaMRI-clone/testData/BIDS4/derivatives/segmentation/sub-control1/ses-2021B/sub-control1_ses-2021B_desc-padded_segmentation.nii.gz'
# main(parc_path = path2, lut_path = "FreeSurferColorLUT.txt")

# # Reset plot layout
# par(mfrow = c(1, 1))

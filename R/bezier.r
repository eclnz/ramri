# Updated R script: nearest‐neighbor ordering + natural spline through all anchors

#––– helper: find extreme point for a given label –––
get_extreme_point <- function(img, label,
                              ext_axes, ext_dirs,
                              weight_axis = NULL) {
  vox <- which(img == label, arr.ind = TRUE)
  if (!nrow(vox)) stop("Label not found: ", label)
  df <- setNames(as.data.frame(vox), c("x","y","z"))
  # apply each “extreme” axis constraint
  for (i in seq_along(ext_axes)) {
    ax  <- ext_axes[i]
    dir <- ext_dirs[i]
    val <- if (dir == "max") max(df[[ax]]) else min(df[[ax]])
    df  <- df[df[[ax]] == val, , drop = FALSE]
  }
  # conditional median weighting
  if (!is.null(weight_axis)) {
    m <- median(df[[weight_axis]])
    df <- df[which.min(abs(df[[weight_axis]] - m)), , drop = FALSE]
  }
  # return single 3D coordinate
  unname(as.numeric(df[1, c("x","y","z")]))
}

#––– nearest‐neighbor ordering of anchor indices –––
reorder_indices <- function(pts) {
  n <- nrow(pts)
  if (n < 2) return(seq_len(n))
  visited <- rep(FALSE, n)
  ord <- integer(n)
  ord[1] <- 1; visited[1] <- TRUE
  for (i in 2:n) {
    prev_idx <- ord[i - 1]
    # compute Euclidean distances from previous
    dists <- sqrt(rowSums((pts - matrix(pts[prev_idx, ], nrow = n, ncol = 3, byrow = TRUE))^2))
    dists[visited] <- Inf
    next_idx <- which.min(dists)
    ord[i] <- next_idx; visited[next_idx] <- TRUE
  }
  ord
}
#––– optimized closed spline curve (minimizing distance between neighbors) –––
make_spline_curve <- function(pts, n = 200) {
  # Ensure closed curve: no start/end discontinuity
  # Find optimal reordering via traveling salesman-like greedy strategy
  n_pts <- nrow(pts)
  dist_matrix <- as.matrix(dist(pts))

  best_order <- NULL
  best_total <- Inf

  for (start in 1:n_pts) {
    visited <- rep(FALSE, n_pts)
    order <- integer(n_pts)
    order[1] <- start
    visited[start] <- TRUE
    total <- 0
    for (i in 2:n_pts) {
      last <- order[i - 1]
      d <- dist_matrix[last, ]
      d[visited] <- Inf
      next_idx <- which.min(d)
      visited[next_idx] <- TRUE
      order[i] <- next_idx
      total <- total + d[next_idx]
    }
    # Add distance to close the loop
    total <- total + dist_matrix[order[n_pts], order[1]]

    if (total < best_total) {
      best_total <- total
      best_order <- order
    }
  }

  pts_ordered <- pts[best_order, , drop = FALSE]
  pts_closed <- rbind(pts_ordered, pts_ordered[1, ])
  t_vals <- seq(0, 1, length.out = nrow(pts_closed))

  fx <- splinefun(t_vals, pts_closed[, 1], method = "fmm")
  fy <- splinefun(t_vals, pts_closed[, 2], method = "fmm")
  fz <- splinefun(t_vals, pts_closed[, 3], method = "fmm")
  
  t_sample <- seq(0, 1, length.out = n)
  cbind(fx(t_sample), fy(t_sample), fz(t_sample))
}

#––– label specifications –––
vertical_specs <- list(
  list(label = 11, ext_axes = "z",              ext_dirs = "max", weight_axis = "y"),  # LR Caudate
  list(label = 50, ext_axes = "z",              ext_dirs = "max", weight_axis = "y"),  # RL Caudate
  list(label = 12, ext_axes = "x",              ext_dirs = "max", weight_axis = "y"),  # LR Putamen
  list(label = 51, ext_axes = "x",              ext_dirs = "min", weight_axis = "y"),  # RL Putamen
  list(label = 14, ext_axes = c("y","z"),     ext_dirs = c("max","min"))                # 3rd-ventricle
)
horizontal_specs <- list(
  list(label = 11, ext_axes = "y",              ext_dirs = "max"),                    # LR Caudate
  list(label = 50, ext_axes = "y",              ext_dirs = "max"),                    # RL Caudate
  list(label = 12, ext_axes = "x",              ext_dirs = "max", weight_axis = "y"),  # LR Putamen
  list(label = 51, ext_axes = "x",              ext_dirs = "min", weight_axis = "y"),  # RL Putamen
  list(label = 10, ext_axes = c("y","x"),     ext_dirs = c("min","max"), weight_axis = "z"),             # LR Thalamus
  list(label = 49, ext_axes = c("y","x"),     ext_dirs = c("min","min"), weight_axis = "z")              # RL Thalamus
)

#––– Compute centroid between two 3D closed curves –––
compute_enclosed_centroid <- function(curve1, curve2) {
  if (nrow(curve1) != nrow(curve2)) stop("Curves must have same number of points.")
  # Combine corresponding point-pairs and interpolate in-between
  midpoints <- (curve1 + curve2) / 2
  # Compute centroid
  colMeans(midpoints)
}

#––– main routine: load, extract, reorder, spline & plot –––
make_and_plot <- function(parc_path, n = 300) {
  # load parcellation
  parc <- readNifti(parc_path)
  parc <- round(parc)
  # extract raw anchor points + labels
  V_raw <- t(sapply(vertical_specs,   function(s)
    get_extreme_point(parc, s$label, s$ext_axes, s$ext_dirs, s$weight_axis)
  ))
  H_raw <- t(sapply(horizontal_specs, function(s)
    get_extreme_point(parc, s$label, s$ext_axes, s$ext_dirs, s$weight_axis)
  ))
  rownames(V_raw) <- rownames(H_raw) <- NULL
  labs_V <- sapply(vertical_specs,   `[[`, "label")
  labs_H <- sapply(horizontal_specs, `[[`, "label")

  # reorder by nearest‐neighbor path
  ordV <- reorder_indices(V_raw)
  ordH <- reorder_indices(H_raw)
  V <- V_raw[ordV, , drop = FALSE]
  H <- H_raw[ordH, , drop = FALSE]
  labs_V <- labs_V[ordV]
  labs_H <- labs_H[ordH]

  # Combine points for ellipsoid fitting
  combined_pts <- rbind(V, H)
  
  # Open a new RGL device to ensure the plot is visible
  open3d()

  # Plot spline curves and points
  curveV <- make_spline_curve(V, n)
  curveH <- make_spline_curve(H, n)
  lines3d(curveV, lwd = 4, col = "red");    points3d(V, size = 8, col = "darkred")
  lines3d(curveH, lwd = 4, col = "blue");   points3d(H, size = 8, col = "darkblue")
  texts3d(V, texts = labs_V, col = "red",     adj = c(1.2,1.2))
  texts3d(H, texts = labs_H, col = "blue",    adj = c(1.2,1.2))

  #––– Compute and plot the centroid between the curves –––
  central_point <- compute_enclosed_centroid(curveV, curveH)
  points3d(t(central_point), col = "black", size = 12)
  texts3d(t(central_point), texts = "Centroid", col = "black", adj = c(0, -1))
  
  # Compute the convex hull of the combined points
  hull_indices <- convhulln(combined_pts)
  
  # Plot the 3D surface of the ROIs and add labels
  for (label in unique(c(labs_V, labs_H))) {
    roi <- parc == label
    if (any(roi)) {
      # Use contour3d to create a surface mesh
      contour3d(roi, level = 0.5, add = TRUE, color = rainbow(length(unique(c(labs_V, labs_H))))[which(unique(c(labs_V, labs_H)) == label)], alpha = 0.5)
      
      # Calculate the centroid of the ROI for labeling
      roi_indices <- which(roi, arr.ind = TRUE)
      centroid <- colMeans(roi_indices)
      
      # Add text label at the centroid
      texts3d(centroid[1], centroid[2], centroid[3], texts = as.character(label), col = "black", adj = c(0.5, 0.5))
    }
  }
  # Add axis labels and a box with ticks
  axes3d(edges = "bbox", labels = TRUE, tick = TRUE)
  box3d()

  title3d(xlab = "i (x)", ylab = "j (y)", zlab = "k (z)")
  legend3d("topright", legend = c("vertical","horizontal"), col = c("red","blue"), lwd = 4)
}

#––– Function to calculate the centroid of a parcellation scheme –––
calculate_parcellation_centroid <- function(parc) {
  # Extract raw anchor points from vertical and horizontal curves
  V_raw <- t(sapply(vertical_specs, function(s)
    get_extreme_point(parc, s$label, s$ext_axes, s$ext_dirs, s$weight_axis)
  ))
  H_raw <- t(sapply(horizontal_specs, function(s)
    get_extreme_point(parc, s$label, s$ext_axes, s$ext_dirs, s$weight_axis)
  ))

  # Reorder points by nearest-neighbor path
  ordV <- reorder_indices(V_raw)
  ordH <- reorder_indices(H_raw)
  V <- V_raw[ordV, , drop = FALSE]
  H <- H_raw[ordH, , drop = FALSE]

  # Create spline curves from the points
  curveV <- make_spline_curve(V)
  curveH <- make_spline_curve(H)

  # Compute the centroid between the two curves
  central_point <- compute_enclosed_centroid(curveV, curveH)

  # Return the x, y, z coordinates of the centroid
  return(round(central_point))
}

# #––– Usage example –––
# file = "/eresearch/qamri-mtbi/ecla535/BIDS_holly_motion/compressed/derivatives/segmentation/sub-expANONYMIZED/ses-2024ANON4295Se10_2/sub-expANONYMIZED_ses-2024ANON4295Se10_2_desc-padded_segmentation.nii.gz"
# make_and_plot(file)
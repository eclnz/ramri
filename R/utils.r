random_slice <- function(dims, noise_sd = 0.1) {
  if (length(dims) < 3 || length(dims) > 5) stop("Input dimensions must be 3, 4, or 5")

  arr <- array(0, dim = dims[1:3])
  mid <- dims[1:3] / 2

  arr[mid[2], , ] <- matrix(abs(stats::rnorm(dims[2] * dims[3], 0, noise_sd)), dims[2], dims[3])
  arr[, mid[1], ] <- matrix(abs(stats::rnorm(dims[1] * dims[3], 0, noise_sd)), dims[1], dims[3])
  arr[, , mid[3]] <- matrix(abs(stats::rnorm(dims[1] * dims[2], 0, noise_sd)), dims[1], dims[2])

  if (length(dims) > 3) {
    arr <- array(rep(arr, prod(dims[4:length(dims)])), dim = dims)
    indices <- as.matrix(expand.grid(lapply(dims[4:length(dims)], seq_len)))

    for (n in seq_len(nrow(indices))) {
      idx <- as.list(indices[n, ])
      arr_subset <- do.call(`[`, c(list(arr), rep(list(TRUE), 3), idx))
      arr_patch <- array(abs(stats::rnorm(prod(dims[1:3]), 0, noise_sd)), dim = dims[1:3])
      do.call(`[<-`, c(list(arr), rep(list(TRUE), 3), idx, list(arr_subset + arr_patch)))
    }
  }
  arr
}
calc_CI <- function(vector, alpha) {

  se <- sd(vector)/sqrt(length(vector))
  t <- qt((1 - alpha)/2 + 0.5, length(vector) - 1)
  CI = t*se
  CI

}

calc_density_over_UMAP <- function(obj, umap_slot_name, data) {

  x <- obj@reductions[[umap_slot_name]]@cell.embeddings[, 1]
  y <- obj@reductions[[umap_slot_name]]@cell.embeddings[, 2]

  z <- calc_density(x, y, data)
  names(z) <- colnames(obj)
  return(z)

}

calc_density <- function(x, y, data) {

  if (NaN %in% data) {

    paste0("NaN in data. Density not calculated")
    exit()

  }
  w <- data / ifelse(sum(data) / length(data) == 0, 1, sum(data) / length(data))
  # h <- c(
  #   ks::hpi(x),
  #   ks::hpi(y)
  # )
  # h <- h * 0.1
  if (all(w == 0)) {

    w <- w + 1

  }
  dens <- ks::kde(x = matrix(data = c(x, y), ncol = 2), w = w, bgridsize = rep(1000, 2))
  ix <- findInterval(x, dens$eval.points[[1]])
  iy <- findInterval(y, dens$eval.points[[2]])
  ii <- cbind(ix, iy)
  z <- dens$estimate[ii]
  return(z)

}

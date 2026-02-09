##--------------------------------------------------------------
## Functions to implement spatial block bootstrap for regression models                                                       
##--------------------------------------------------------------

library(sf)
library(gstat)
library(parallel)

make_block_structure <- function(df, resid_vec, lon_col = "lon", lat_col = "lat",
                                 proj_crs = 5070, cutoff = 500000) {

  df_coords <- data.frame(x = df[[lon_col]], y = df[[lat_col]])
  df_sf <- st_as_sf(df_coords, coords = c("x", "y"), crs = 4326)
  df_sf <- st_transform(df_sf, crs = proj_crs)
  coords_proj <- st_coordinates(df_sf)
  df$x <- coords_proj[,1]
  df$y <- coords_proj[,2]

  resid_df <- data.frame(z = resid_vec, x = df$x, y = df$y)
  v1 <- variogram(z ~ 1, data = resid_df, locations = ~x + y, cutoff = cutoff, width = cutoff/20)
  f1 <- fit.variogram(v1, vgm("Sph"))
  max.dist <- f1$range[2]
  bb <- ceiling(max.dist)

  area <- list(
    xmin = min(df$x, na.rm = TRUE),
    xmax = max(df$x, na.rm = TRUE),
    ymin = min(df$y, na.rm = TRUE),
    ymax = max(df$y, na.rm = TRUE)
  )

  box_centers <- expand.grid(
    x = seq(area$xmin, area$xmax, by = bb) + bb/2,
    y = seq(area$ymin, area$ymax, by = bb) + bb/2
  )

  sample_points <- list(
    xmin = min(box_centers$x),
    xmax = max(box_centers$x),
    ymin = min(box_centers$y),
    ymax = max(box_centers$y),
    num_boxes = nrow(box_centers)
  )

  list(df = df, bb = bb, box_centers = box_centers, sample_points = sample_points)
}

pull_data <- function(xcoord, ycoord, data, bb) {
  data <- data[data$x > (xcoord - (bb/2)) & data$x < (xcoord + (bb/2)), , drop = FALSE]
  data <- data[data$y > (ycoord - (bb/2)) & data$y < (ycoord + (bb/2)), , drop = FALSE]
  return(data)
}

rpull <- function(nbox, sdata, bb, sample_points, mc.cores = 7) {

  xcoords <- runif(nbox, sample_points[["xmin"]], sample_points[["xmax"]])
  ycoords <- runif(nbox, sample_points[["ymin"]], sample_points[["ymax"]])

  bs.data <- mcmapply(FUN = pull_data,
                      xcoords,
                      ycoords,
                      MoreArgs = list(data = sdata, bb = bb),
                      SIMPLIFY = FALSE,
                      mc.cores = mc.cores)
  out <- do.call(rbind, bs.data)
  return(out)
}

run_bootstrap <- function(iter,
                          full_data,
                          bb_info,
                          fields,
                          formula,
                          fitfun = lm,
                          weight_var = "se_pca",
                          nbox_per_draw = 10,
                          cols = NULL,
                          mc.cores = 7) {

  outlist <- list()
  enough <- FALSE
  i <- 1

  while (!enough) {

    outlist[[i]] <- rpull(nbox = nbox_per_draw,
                          sdata = full_data[, fields, drop = FALSE],
                          bb = bb_info$bb,
                          sample_points = bb_info$sample_points,
                          mc.cores = mc.cores)

    combined <- do.call(rbind, outlist)

    cond1 <- nrow(combined) >= nrow(full_data) * 0.95

    # ensure that the bootstrap contains all factor levels of e.g. ecoregion
    cond2 <- TRUE
    if ("ecoregion" %in% fields) {
      cond2 <- length(unique(combined$ecoregion)) == length(unique(full_data$ecoregion))
    }

    if (cond1 & cond2) {
      enough <- TRUE
    } else {
      i <- i + 1
      # continue
    }
  }

  bs.data <- combined

  if (nrow(bs.data) > nrow(full_data) * 1.05) {
    bs.data <- bs.data[sample(1:nrow(bs.data), nrow(full_data)), , drop = FALSE]
  }

  if (!is.null(weight_var) && weight_var %in% names(bs.data)) {
    w <- 1 / (bs.data[[weight_var]]^2)
    call_args <- list(formula = formula,
                      data = bs.data,
                      weights = w)
    bs.model <- do.call(fitfun, call_args)
  } else {
    bs.model <- fitfun(formula = formula, data = bs.data)
  }

  ests <- coef(bs.model)
  
  if (!is.null(cols)) {
    outrow <- as.data.frame(matrix(NA, nrow = 1, ncol = length(cols)))
    colnames(outrow) <- cols
    common <- intersect(names(ests), cols)
    outrow[1, common] <- ests[common]
  } else {
    outrow <- as.data.frame(as.list(ests))
  }

  rm(bs.model)
  gc()

  return(outrow)
}

bootstrap_model_blocks <- function(full_data,
                                   formula,
                                   fitfun = lm,
                                   weight_var = "se_pca",
                                   fields,
                                   filepath = here("data", "bootstraps", "bootstrap_coefs_ecoregion_blocks.csv"),
                                   N = 200,
                                   nbox_per_draw = 10,
                                   mc.cores = 7,
                                   progress = TRUE) {

  if (!is.null(weight_var) && weight_var %in% names(full_data)) {

    weights_vec <- 1 / (full_data[[weight_var]]^2)
    call_args <- list(formula = formula,
                      data = full_data,
                      weights = weights_vec)
    base_fit <- do.call(fitfun, call_args)
    
  } else {
    base_fit <- fitfun(
        formula = formula,
        data = full_data
    )
  }
  
  cols <- names(coef(base_fit))

  message("Estimating scale of spatial autocorrelation")
  bb_info <- make_block_structure(df = full_data, resid_vec = residuals(base_fit),
                                  lon_col = "lon", lat_col = "lat",
                                  proj_crs = 5070, cutoff = 500000)
  
  rm(base_fit); gc()

  bs.coefs <- as.data.frame(matrix(NA, nrow = 0, ncol = length(cols)))
  colnames(bs.coefs) <- cols

  for (j in seq_len(N)) {
    if (progress) message("starting iteration: ", j)
    res_row <- run_bootstrap(iter = j,
                             full_data = bb_info$df,
                             bb_info = bb_info,
                             fields = fields,
                             formula = formula,
                             fitfun = fitfun,
                             weight_var = weight_var,
                             nbox_per_draw = nbox_per_draw,
                             cols = cols,
                             mc.cores = mc.cores)

    bs.coefs <- rbind(bs.coefs, res_row)

    write.csv(bs.coefs, file = filepath,
              row.names = FALSE)

    if (progress) message("completed iteration: ", nrow(bs.coefs))
  }

  return(bs.coefs)
}

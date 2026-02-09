library(ggplot2)

plot_tradeoff_drivers <- function(var, model, data, bootstrap_results,
                                  orig_var = NULL,
                                  n = 100,
                                  ci_mult = 1.97,
                                  xlab = NULL,
                                  ylab = NULL,
                                  ylim = NULL,
                                  xlim = NULL,
                                  delog = FALSE,
                                  fill = "#d0d0d0",
                                  line_color = "#4b4b4b",
                                  span = 0.25,               # LOESS smoothing span for CI curves
                                  quantiles = c(0.025, 0.975)) {

  xseq <- seq(min(data[[var]], na.rm = TRUE),
              max(data[[var]], na.rm = TRUE),
              length.out = n)
  pdata <- data.frame(x = xseq)
  names(pdata)[1] <- var

  all_preds <- all.vars(stats::delete.response(stats::terms(model)))
  others <- setdiff(all_preds, var)
  for (v in others) {
    if (!v %in% names(data)) next
    col <- data[[v]]
    if (is.numeric(col)) {
      pdata[[v]] <- mean(col, na.rm = TRUE)
    } else {
      # factor/character -> modal level
      if (is.character(col)) col <- factor(col)
      modal <- levels(col)[which.max(table(col))]
      pdata[[v]] <- factor(rep(modal, nrow(pdata)), levels = levels(col))
    }
  }


  pf <- tryCatch(
    predict(model, newdata = pdata, type = "response"),
    error = function(e) NULL
  )
  if (!is.null(pf)) {
    if (is.matrix(pf) && ncol(pf) >= 1) {
      pdata$fit <- as.numeric(pf[,1])
    } else {
      pdata$fit <- as.numeric(pf)
    }
  } else {
    warning("predict() failed; using bootstrap mean coefficients to compute fit.")
    pdata$fit <- NA_real_
  }

  br_mat <- as.data.frame(bootstrap_results)
  terms_no_resp <- stats::delete.response(stats::terms(model))
  X <- model.matrix(terms_no_resp, pdata)   # n_grid x p
  needed_coefs <- colnames(X)

  missing_coefs <- setdiff(needed_coefs, colnames(br_mat))
  if (length(missing_coefs) > 0) {
    stop("bootstrap_results is missing coefficient columns: ", paste(missing_coefs, collapse = ", "))
  }

  Bmat <- as.matrix(br_mat[, needed_coefs, drop = FALSE])  # nboot x p
  preds_mat <- Bmat %*% t(X)   # matrix: nboot x n_grid

  qfun <- function(v) stats::quantile(v, probs = quantiles, na.rm = TRUE, names = FALSE)
  qmat <- apply(preds_mat, 2, qfun)  # matrix: length(quantiles) x n_grid
  lower_raw <- as.numeric(qmat[1, ])
  upper_raw <- as.numeric(qmat[2, ])

  if (length(unique(xseq)) > 10) {
    lower_smooth <- stats::predict(stats::loess(lower_raw ~ xseq, span = span), xseq)
    upper_smooth <- stats::predict(stats::loess(upper_raw ~ xseq, span = span), xseq)
    if (any(is.na(lower_smooth))) lower_smooth[is.na(lower_smooth)] <- lower_raw[is.na(lower_smooth)]
    if (any(is.na(upper_smooth))) upper_smooth[is.na(upper_smooth)] <- upper_raw[is.na(upper_smooth)]
  } else {
    lower_smooth <- lower_raw
    upper_smooth <- upper_raw
  }

  pdata$lower <- lower_smooth
  pdata$upper <- upper_smooth

  if (is.na(pdata$fit[1])) {
    pdata$fit <- as.numeric(colMeans(preds_mat, na.rm = TRUE))
  }

  if (delog) {
    pdata$fit <- exp(pdata$fit)
    pdata$lower <- exp(pdata$lower)
    pdata$upper <- exp(pdata$upper)
  }

  if (!is.null(orig_var) && orig_var %in% names(data)) {
    pdata[[var]] <- pdata[[var]] * sd(data[[orig_var]], na.rm = TRUE) +
      mean(data[[orig_var]], na.rm = TRUE)
    if (is.null(xlab)) xlab <- orig_var
  } else {
    if (is.null(xlab)) xlab <- var
  }
  if (is.null(ylab)) ylab <- "pred"

  gg <- ggplot(pdata, aes_string(x = var, y = "fit")) +
    geom_ribbon(aes_string(ymin = "lower", ymax = "upper"), fill = fill, alpha = 0.5) +
    geom_line(size = 1.2, color = line_color) +
    theme_bw() +
    xlab(xlab) + ylab(ylab) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))

  if (!is.null(ylim)) gg <- gg + scale_y_continuous(expand = c(0,0), limits = ylim)
  if (!is.null(xlim)) gg <- gg + scale_x_continuous(expand = c(0,0), limits = xlim)

  return(gg)
}

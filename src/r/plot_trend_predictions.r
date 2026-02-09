##--------------------------------------------------------------
## Functions to plot trend model predictions
##
## author: Jacob Levine; email: jacob.levine@duke.edu
##--------------------------------------------------------------


fit_trend_model <- function(form, data, response = "beta_pca") {

  form_full <- as.formula(paste0(response, form))

  if (response == "beta_pca") {
    mod <- bam(formula = form_full,
               weights = 1 / se_pca^2,
               data = data, method = "fREML",
               discrete = TRUE)
  } else {
    mod <- bam(formula = form_full,
               data = data, method = "fREML",
               discrete = TRUE)
  }

  data$cell_id <- paste(data$lon, data$lat, sep = "_")
  
  ## mark start of each time series
  data <- data[order(data$cell_id, data$year), ]
  data$AR_start <- c(TRUE, diff(as.numeric(as.factor(data$cell_id))) != 0)

  # Get model residuals
  data$resid <- residuals(mod)
  
  cell_ids <- unique(data$cell_id)
  rho_values <- numeric(length(cell_ids))

  for (i in seq_along(cell_ids)) {
    this_resid <- data$resid[data$cell_id == cell_ids[i]]
    if (length(this_resid) > 1) {
      acf_val <- acf(this_resid, lag.max = 1, plot = FALSE)$acf[2]
      rho_values[i] <- acf_val
    } else {
      rho_values[i] <- NA
    }
  }

  # Final estimate
  estimated_rho <- median(rho_values, na.rm = TRUE)

  if (response == "beta_pca") {
    mod <- bam(formula = form_full,
               weights = 1 / se_pca^2,
               data = data,
               rho = estimated_rho,
               AR.start = data$AR_start,
               method = "fREML", discrete = TRUE)
  } else { 
    mod <- bam(formula = form_full,
               data = data,
               rho = estimated_rho,
               AR.start = data$AR_start,
               method = "fREML", discrete = TRUE)
  }

  return(mod)

}



plot_trend_prediction <- function(model, data, variable = "year", pred_length = 100, 
                                  interaction_var = NULL, interaction_levels = 3, full_data = NULL,
                                  ylab = "tradeoff slope (<1 indicates tradeoff)", with_intercept = TRUE) {

  all_vars <- attr(model$terms, "term.labels")
  interaction_terms <- all_vars[grepl(":", all_vars)]
  all_vars <- all_vars[!(all_vars %in% interaction_terms)]
  
  if (variable == "year") {
    var <- "year_centered"
  } else {
    var <- paste0(variable, "_scaled")
  }

  # Detect interactions with the focal variable
  focal_interactions <- interaction_terms[grepl(var, interaction_terms)]
  
  # If no interaction_var specified but interactions exist, use the first one
  if (is.null(interaction_var) && length(focal_interactions) > 0) {
    # Extract the other variable from the interaction term
    interaction_components <- strsplit(focal_interactions[1], ":")[[1]]
    interaction_var <- interaction_components[interaction_components != var][1]
  }

  ## if full data not provided, use supplied data
  if (is.null(full_data)) {
    full_data <- data
  }
  
  # Handle interactions
  if (!is.null(interaction_var) && interaction_var %in% all_vars) {
    # Create levels for the interaction variable
    int_var_values <- quantile(data[[interaction_var]], 
                               probs = seq(0.1, 0.9, length.out = interaction_levels),
                               na.rm = TRUE)

    plot_list <- list()
    
    for (i in seq_along(int_var_values)) {
      pdata <- as.data.frame(matrix(0, pred_length, 0))
      
      # Set all other variables to their means
      for (v in all_vars[!(all_vars %in% c(var, interaction_var))]) {
        pdata[,v] <- rep(mean(data[,v], na.rm = TRUE), times = pred_length)
      }
      
      # Set the focal variable range
      pdata[,var] <- seq(quantile(data[,var], 0.05, na.rm = TRUE),
                         quantile(data[,var], 0.95, na.rm = TRUE),
                         length.out = pred_length)
      
      # Set the interaction variable to this level
      pdata[,interaction_var] <- rep(int_var_values[i], times = pred_length)
      
      # Get predictions
      p <- predict.gam(model, newdata = pdata, se.fit = TRUE, type = "link", unconditional = TRUE,
                       exclude = c("s(lon,lat)", "(Intercept)")) ## exclude variance from spatial term and intercept

      ## add intercept back in:
      if (with_intercept) { 
        p$fit <- p$fit + coef(model)[1]
      }
      
      linkinv <- family(model)$linkinv
      
      pdata$beta_pca <- linkinv(p$fit)
      pdata$beta_pca_lower <- linkinv(p$fit - 1.97 * p$se.fit)
      pdata$beta_pca_upper <- linkinv(p$fit + 1.97 * p$se.fit)

      # Transform back to original scale
      if (variable == "year") {
        pdata$variable <- pdata[,var] + mean(full_data[,variable], na.rm = TRUE)
      } else {
        pdata$variable <- pdata[,var] * sd(full_data[,variable], na.rm = TRUE) + 
          mean(full_data[,variable], na.rm = TRUE)
      }
      
      # Add interaction level info
      if (grepl("_scaled$", interaction_var)) {
        int_var_orig <- gsub("_scaled$", "", interaction_var)
        int_level_orig <- int_var_values[i] * sd(full_data[[int_var_orig]], na.rm = TRUE) + 
          mean(full_data[[int_var_orig]], na.rm = TRUE)
        pdata$interaction_level <- round(int_level_orig, -2)
      } else {
        pdata$interaction_level <- round(int_var_values[i], -2)
      }

      plot_list[[i]] <- pdata
    }
    
    # Combine all data
    plot_data <- do.call(rbind, plot_list)

    plot_data$interaction_level <- as.factor(plot_data$interaction_level)
    
    # create plot with interaction levels
    ggplot(plot_data, aes(x = variable, y = beta_pca)) +
      geom_ribbon(data = plot_data,
                  aes(ymin = beta_pca_lower, ymax = beta_pca_upper, color = interaction_level, fill = interaction_level),
                  color = NA, alpha = 0.15) + 
      geom_line(data = plot_data, aes(color = interaction_level), linewidth = 2.5) +
      theme_bw() +
      xlab(variable) +
      ylab(ylab) +
      scale_x_continuous(expand = c(0,0)) +
      scale_fill_brewer(palette = "Reds") +
      scale_color_brewer(palette = "Reds") + 
      labs(color = interaction_var, fill = interaction_var) +
      theme(legend.position = "bottom")
    
  } else {
    
    # No interactions - original behavior
    pdata <- as.data.frame(matrix(0, pred_length, 0))
    for (v in all_vars[all_vars != var]) {
      pdata[,v] <- rep(mean(data[,v], na.rm = TRUE), times = pred_length)
    }
    pdata[,var] <- seq(quantile(data[,var], 0.05, na.rm = TRUE),
                       quantile(data[,var], 0.95, na.rm = TRUE),
                       length.out = pred_length)

    p <- predict.gam(model, newdata = pdata, se.fit = TRUE, type = "link", unconditional = TRUE,
                     exclude = c("s(lon,lat)", "(Intercept)"))

    if (with_intercept) { 
      p$fit <- p$fit + coef(model)[1]
    }
    
    linkinv <- family(model)$linkinv
    
    pdata$beta_pca <- linkinv(p$fit)
    pdata$beta_pca_lower <- linkinv(p$fit - 1.97 * p$se.fit)
    pdata$beta_pca_upper <- linkinv(p$fit + 1.97 * p$se.fit)

    if (variable == "year") {
      pdata$variable <- pdata[,var] + mean(data[,variable], na.rm = TRUE)
    } else {
      pdata$variable <- pdata[,var] * sd(data[,variable], na.rm = TRUE) + 
        mean(data[,variable], na.rm = TRUE)
    }
    
    ggplot(pdata, aes(x = variable, y = beta_pca)) +
      geom_ribbon(aes(ymin = beta_pca_lower, ymax = beta_pca_upper),
                  fill = "lightgray", alpha = 0.8) + 
      geom_line(linewidth = 2) +
      theme_bw() +
      xlab(variable) +
      ylab(ylab) +
      scale_x_continuous(expand = c(0,0))
  }
}


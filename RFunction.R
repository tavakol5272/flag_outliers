library(move2)
library(sf)
library(terra)
library(dplyr)
library(lubridate)


########## Helper: plot per individual ###################################

plot_per_individual <- function(dat, prob_label, drop_na) {
  trk_col <- mt_track_id_column(dat)
  ids     <- unique(as.character(dat[[trk_col]]))
  n_ind   <- length(ids)
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  
  # 2 panels per individual
  par(mfrow = c(n_ind, 2), mar = c(4, 4, 2, 1))
  pal <- grDevices::hcl.colors(100, "Viridis")
  
  for (id in ids) {
    dat_i <- dat[dat[[trk_col]] == id, ]
    coords <- sf::st_coordinates(dat_i)
    is_na_prob_i <- dat_i$is_na_prob
    is_outlier_i <- dat_i$is_outlier
    n_outliers_i <- sum(is_outlier_i, na.rm = TRUE)
    n_na_i       <- sum(is_na_prob_i)
    
    ## Panel 1: Probability colouring 
    plot(coords, type = "l",
         main = paste0("All locations (", id, ")"),
         col = "black", xlab = "Longitude", ylab = "Latitude", asp = 1)
    
    ok <- !is.na(dat_i$log_prob) & !is_na_prob_i
    if (any(ok)) {
      rng <- range(dat_i$log_prob[ok], na.rm = TRUE)
      idx <- pmax(
        1, pmin(
          100,
          ceiling((dat_i$log_prob[ok] - rng[1]) / max(1e-12, rng[2] - rng[1]) * 99) + 1
        )
      )
      points(coords[ok, 1], coords[ok, 2], col = pal[idx], pch = 19, cex = 0.7)
    }
    if (any(is_na_prob_i)) {
      na_ix <- which(is_na_prob_i)
      points(coords[na_ix, 1], coords[na_ix, 2], col = "gray", pch = 19, cex = 0.7)
    }
    legend("topright", legend = NA,
           col = "gray", pch = 19, bty = "n", title = prob_label)
    
    ##  Panel 2: Kept vs removed
    to_remove_i <- is_outlier_i
    if (isTRUE(drop_na)) to_remove_i <- to_remove_i | is_na_prob_i
    
    kept_ix    <- which(!to_remove_i)
    removed_ix <- which(to_remove_i)
    
    plot(coords, type = "n",
         main = paste0("Kept vs removed (", id, ")"),
         xlab = "Longitude", ylab = "Latitude", asp = 1)
    if (length(kept_ix) > 1) {
      lines(coords[kept_ix, ], col = "black")
    }
    if (length(kept_ix) > 0) {
      points(coords[kept_ix, 1], coords[kept_ix, 2],
             col = "blue", pch = 19, cex = 0.7)
    }
    if (length(removed_ix) > 0) {
      points(coords[removed_ix, 1], coords[removed_ix, 2],
             col = "red", pch = 19, cex = 0.7)
    }
    legend(
      "topright",
      legend = c(
        paste0("Outliers: ", n_outliers_i),
        if (isTRUE(drop_na)) paste0("NAs: ", n_na_i) else NULL
      ),
      col = c("red", if (isTRUE(drop_na)) "gray" else NULL),
      pch = 19, bty = "n"
    )
  }
}

############ Main Function ######################################

rFunction <-  function(data,
                       threshold,
                       prob_type = c("joint", "step_turn", "delta_step", "delta_turn", "custom"),
                       remove = FALSE,
                       plot = TRUE,
                       drop_na = FALSE,
                       .per_track = FALSE) {
  
  # Handle NULL data
  if (is.null(data) || nrow(data) == 0) {
    logger.info("Input is NULL or has 0 rows — returning NULL.")
    return(NULL)
  }
  
  
  # split by individual 
  if (!.per_track) {
    trk_col <- mt_track_id_column(data)
    ids     <- unique(data[[trk_col]])
    
    if (length(ids) > 1L) {
      split_list <- split(data, data[[trk_col]])
      
      # compute per individual 
      res_list <- lapply(split_list, function(tr) {
        rFunction(
          data       = tr,
          threshold  = threshold,
          prob_type  = prob_type,
          remove     = remove,
          plot       = FALSE,
          drop_na    = drop_na,
          .per_track = TRUE
        )
      })
      
      result <- do.call(rbind, res_list)
      
      # label for plots
      prob_label <- switch(
        prob_type,
        joint      = "Joint probability",
        step_turn  = "Step length & turning angle probability",
        delta_step = "Step length change probability",
        delta_turn = "Turning angle change probability",
        custom     = "Custom probability product"
      )
      
      # 2 panels per individual
      if (isTRUE(plot)) {
        plot_per_individual(result, prob_label, drop_na)
      }
      
      return(result)
    }
  }
  
  ####### cleaning per individual #######################
  #remove:empty locations, bad timestamps, extreme speeds (top 1%), very short time lags (< 60 s)
  
  
  data <- data[order(mt_track_id(data)), ]
  
  data <- data %>%
    arrange(timestamp) %>%
    mutate(
      timelag      = as.numeric(difftime(timestamp, lag(timestamp), units = "secs")),
      step_length  = as.numeric(st_distance(geometry, lag(geometry), by_element = TRUE)), 
      ground_speed_ms = step_length / timelag
    ) %>%
    filter(
      !st_is_empty(geometry),
      !is.na(timestamp),
      is.na(timelag) | timelag >= 60        
    ) %>%
    filter(
      ground_speed_ms <= quantile(ground_speed_ms, 0.99, na.rm = TRUE)  
    ) %>%
    ungroup()
  
  
  if (nrow(data) < 3) {
    logger.info("Fewer than 3 locations after cleaning; returning cleaned data without outlier flags.")
    return(data)
  }
  
  ###############################################################
  
  # Function to calculate 2D histogram of turning angle and step length
  TurnStepHist <- function(x, y, stand = TRUE, verbose = FALSE) {
    bwx <- nclass.FD(x[!is.na(x)])
    bwy <- nclass.FD(y[!is.na(y)])
    nx <- max(bwx, 12)
    ny <- max(bwy, 12)
    
    test <- terra::rast(
      ncol = nx, nrow = ny,
      xmin = -pi, xmax = pi,
      ymin = 0, ymax = 1.1 * (max(y, na.rm = TRUE)),
      crs = ""
    )
    
    xyRaster <- terra::rasterize(terra::vect(cbind(x, y)), test, fun = "count")
    xyRaster[is.na(xyRaster)] <- 0
    xyRaster <- xyRaster / sum(terra::values(xyRaster), na.rm = TRUE)
    
    if (verbose) {
      logger.info("x ", nrow(xyRaster), " y ", ncol(xyRaster))
    }
    
    l <- xyRaster
    r <- xyRaster
    terra::ext(l) <- terra::ext(terra::ext(l)[1] - 2*pi, terra::ext(l)[2] - 2*pi,
                                terra::ext(l)[3], terra::ext(l)[4])
    terra::ext(r) <- terra::ext(terra::ext(r)[1] + 2*pi, terra::ext(r)[2] + 2*pi,
                                terra::ext(r)[3], terra::ext(r)[4])
    
    xyRaster <- terra::merge(l, xyRaster, r)
    
    # increase resolution through bilinear interpolation
    rasterXY <- terra::resample(
      xyRaster,
      terra::rast(
        ncol = 150, nrow = 150,
        xmin = -pi, xmax = pi,
        ymin = 0, ymax = max(y, na.rm = TRUE),
        crs = ""
      ),
      method = "bilinear"
    )
    
    rasterXY[rasterXY < 0] <- 0
    rasterXY <- rasterXY / sum(terra::values(rasterXY), na.rm = TRUE)
    
    return(list(rasterXY, xyRaster)[[2 - stand]])
  }
  
  # Function to calculate probability distributions
  get.densities.2d <- function(stepLength, turningAngle, deltaStep, deltaTurn) {
    if (length(turningAngle) != length(stepLength)) {
      stop("Vector lengths of step length and turning angles do not match.")
    }
    rasterTS <- TurnStepHist(x = as.vector(turningAngle), y = as.vector(stepLength))
    
    dS <- stats::density(deltaStep[!is.na(deltaStep)])
    dT <- stats::density(deltaTurn[!is.na(deltaTurn)])
    autoS <- stats::approxfun(dS$x, dS$y, rule = 2)
    autoT <- stats::approxfun(dT$x, dT$y, rule = 2)
    
    return(list(TSRaster = rasterTS, autoT = autoT, autoS = autoS))
  }
  
  # Function to calculate joint movement probabilities
  get_joint_movement_probabilities <- function(stepLength, turningAngle, TwoDHist) {
    n <- length(stepLength)
    
    if (length(turningAngle) != n) {
      stop("stepLength and turningAngle must have the same length")
    }
    
    # Initialize probability vectors
    step_turn_probs <- rep(NA_real_, n)
    delta_step_probs <- rep(NA_real_, n)
    delta_turn_probs <- rep(NA_real_, n)
    joint_probs <- rep(NA_real_, n)
    
    rast <- TwoDHist$TSRaster
    auto_step_func <- TwoDHist$autoS
    auto_turn_func <- TwoDHist$autoT
    
    coords <- cbind(as.numeric(turningAngle), as.numeric(stepLength))
    probs_df <- terra::extract(rast, coords)
    probs_df[is.na(probs_df)] <- 0
    
    step_turn_probs <- probs_df[, 1]
    
    delta_steps <- diff(c(NA, stepLength))
    delta_turns <- diff(c(NA, turningAngle))
    
    for (i in 2:n) {
      if (!is.na(delta_steps[i])) {
        delta_step_probs[i] <- auto_step_func(delta_steps[i])
      }
      
      wrap <- function(x) { (x + pi) %% (2 * pi) - pi }
      if (!is.na(delta_turns[i])) {
        wrapped_delta <- wrap(delta_turns[i])
        delta_turn_probs[i] <- auto_turn_func(wrapped_delta)
      }
      
      if (!is.na(step_turn_probs[i]) && !is.na(delta_step_probs[i]) && !is.na(delta_turn_probs[i])) {
        joint_probs[i] <- step_turn_probs[i] * delta_step_probs[i] * delta_turn_probs[i]
      }
    }
    
    return(list(
      step_turn_probs = step_turn_probs,
      delta_step_probs = delta_step_probs,
      delta_turn_probs = delta_turn_probs,
      joint_probs = joint_probs
    ))
  }
  
  # Calculate initial step lengths
  step_lengths0 <- mt_distance(data)
  if (any(as.vector(step_lengths0) == 0, na.rm = TRUE)) {
    st_geometry(data) <- sf::st_jitter(st_geometry(data), amount = 1e-6)
  }
  
  # Calculate movement metrics
  stepLength   <- mt_distance(data, units = "m")
  turningAngle <- mt_turnangle(data)
  
  if (length(stepLength) < 3) {
    logger.info("Fewer than 3 locations; cannot compute deltas. Returning input.")
    return(data)
  }
  
  deltaStep <- diff(as.numeric(stepLength))
  deltaTurn <- diff(as.numeric(turningAngle))
  
  # logger.info("Calculating probability distributions...")
  TwoDHist <- get.densities.2d(stepLength, turningAngle, deltaStep, deltaTurn)
  
  # logger.info("Calculating joint probabilities...")
  probability_components <- get_joint_movement_probabilities(
    as.numeric(stepLength),
    as.numeric(turningAngle),
    TwoDHist
  )
  
  # Add probability components to the move object
  data$step_turn_prob  <- probability_components$step_turn_probs
  data$delta_step_prob <- probability_components$delta_step_probs
  data$delta_turn_prob <- probability_components$delta_turn_probs
  data$joint_prob      <- probability_components$joint_probs
  
  # choose probability column
  prob_col <- switch(
    prob_type,
    joint      = "joint_prob",
    step_turn  = "step_turn_prob",
    delta_step = "delta_step_prob",
    delta_turn = "delta_turn_prob",
    custom     = { data$custom_prob <- data$step_turn_prob * data$delta_step_prob; "custom_prob" }
  )
  prob_label <- switch(
    prob_type,
    joint      = "Joint probability",
    step_turn  = "Step length & turning angle probability",
    delta_step = "Step length change probability",
    delta_turn = "Turning angle change probability",
    custom     = "Custom probability product"
  )
  
  eps <- 1e-300
  p   <- data[[prob_col]]
  data$log_prob <- log10(pmax(p, eps))
  
  # Identify NA probabilities
  is_na_prob <- is.na(data[[prob_col]])
  n_na <- sum(is_na_prob)
  
  # Identify outliers
  # logger.info("Identifying outliers...")
  probs <- data[[prob_col]]
  ecdf_func <- stats::ecdf(probs[!is.na(probs)])
  data$outlier_percentile <- round((1 - ecdf_func(probs)) * 100, 2)
  is_outlier <- data$outlier_percentile >= (100 - threshold * 100)
  data$is_outlier <- is_outlier
  data$is_na_prob <- is_na_prob
  
  n_outliers <- sum(is_outlier, na.rm = TRUE)
  n_total    <- sum(!is.na(probs))
  
  outlier_message <- paste0(
    "Identified ", n_outliers, " outliers (",
    round(100 * n_outliers / max(1, n_total), 2),
    "% of locations with probabilities) using ", prob_label,
    " @ ", threshold * 100, "th percentile."
  )
  logger.info(outlier_message)
  if (n_na > 0) {
    logger.info(paste0(
      "NA ", prob_label, " in ", n_na, " locations (",
      round(100 * n_na / nrow(data), 2), "% of all). ",
      ifelse(drop_na, "They will be removed.", "They will be kept.")
    ))
  }
  
  # plots for each individual
  if (isTRUE(plot)) {
    plot_per_individual(data, prob_label, drop_na)
  }
  
  ######## step_length_mv, turning_angle_mv, and first/last day ########
  
  # step length and turning angle
  data$step_length_mv   <- as.numeric(stepLength)
  data$turning_angle_mv <- as.numeric(turningAngle)
  
  # time from current location to the next
  data$timelag_tonext <- as.numeric(
    difftime(dplyr::lead(data$timestamp), data$timestamp, units = "secs")
  )
  
  # speed from current to next point
  data$ground_speed_ms_tonext <- data$step_length_mv / data$timelag_tonext
  
  # remove first and last day 
  data <- data %>%
    dplyr::mutate(date = as.Date(timestamp)) %>%
    dplyr::filter(
      date != min(date, na.rm = TRUE),
      date != max(date, na.rm = TRUE)
    )
  
  #################################################################
  
  # print the table of outliers
  trk_col <- mt_track_id_column(data)
  
  outliers_table <- subset(
    cbind(data.frame(track_id = as.character(data[[trk_col]])), as.data.frame(data)),
    is_outlier == TRUE,
    select = c("track_id", "timestamp", "step_turn_prob",
               "joint_prob", "outlier_percentile", "is_outlier")
  )
  
  n_outliers_after <- sum(data$is_outlier, na.rm = TRUE)
  
  logger.info(paste(
    "Preview", n_outliers_after,
    "outliers for individual(s) (not including first/last day):",
    paste(unique(as.character(data[[trk_col]])), collapse = ", ")
  ))
  
  cat("\n")
  print(outliers_table)
  cat("\n\n")
  
  ####################################################
  
  
  # Return filtered or original object
  if (isTRUE(remove)) {
    to_remove <- data$is_outlier | (drop_na & data$is_na_prob)
    filtered_data <- data[!to_remove, ]
    return(filtered_data)
  } else {
    return(data)
  }
}

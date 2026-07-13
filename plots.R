library(plotly)
library(tidyr)
library(dplyr)
library(MASS)
library(expm)
library(ggnewscale)
library(patchwork)


plot_2d_curve <- function(df, par,
                          line_size = 0.2,
                          palette   = "viridis") {
  
  if (!all(c("X1", "X2") %in% names(df))) {
    stop("df must contain columns 'X1' and 'X2'")
  }
  
  df$idx <- seq_len(nrow(df))
  
  cols <- viridis::viridis(256, option = palette)
  
  xs <- seq(-2.3, 2.2, length.out = 400)
  
  # Nullclines
  df_xnull <- data.frame(
    x = xs,
    y = xs - xs^3 + par[2]
  )
  
  df_ynull <- data.frame(
    x = xs,
    y = par[3] * xs + par[4]
  )
  
  ggplot(df, aes(x = X1, y = X2, color = idx)) +
    
    geom_path(
      linewidth = line_size,
      lineend = "round"
    ) +
    
    scale_color_gradientn(
      colors = cols,
      name = "index",
      guide = "colourbar"
    ) +
    
    # x-nullcline
    geom_line(
      data = df_xnull,
      aes(x = x, y = y),
      inherit.aes = FALSE,
      color = "black",
      linetype = "dashed",
      linewidth = 0.5
    ) +
    
    # y-nullcline
    geom_line(
      data = df_ynull,
      aes(x = x, y = y),
      inherit.aes = FALSE,
      color = "black",
      linetype = "dashed",
      linewidth = 0.5
    ) +
    coord_cartesian(xlim = c(-1.4,1.3),ylim=c(-0.3,1.5)) + 
    labs(x = "X", y = "Y") +
    
    theme_bw(base_size = 16, base_family = "Times") +
    
    theme(
      legend.position = "none",
      legend.key.size = unit(1.6, "cm"),
      legend.text = element_text(size = 16),
      
      strip.text.x = element_text(
        size = 16,
        color = "black",
        face = "bold"
      ),
      
      strip.text.y = element_text(
        size = 16,
        color = "black",
        face = "bold"
      ),
      
      strip.text = element_text(face = "bold"),
      
      strip.background = element_rect(fill = "grey")
    )
}

plot_FHN_phase_portrait <- function(par=c(0.1,0.5,1.5,1.4),
                             x_range = c(-1.4, 1.4),
                             y_range = c(-0.3, 1.2),
                             n_points = 40,
                             arrow_scale = 0.12,
                             show_nullclines = TRUE,
                             show_equilibria = TRUE) {
  
  epsilon   <- par[1]; alpha <- par[2]; gamma <- par[3]; beta  <- par[4]
  f <- function(x,y) {  
    (-x^3+x+alpha-y)/epsilon
  }
  g <- function(x,y) {    # y' = (2 - 0.4*x)*x + (0.3 + 0.1*x)
    gamma*x+beta-y
  }
  # build grid
  x <- seq(x_range[1], x_range[2], length.out = n_points)
  y <- seq(y_range[1], y_range[2], length.out = n_points)
  grid <- expand.grid(x = x, y = y, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  
  # evaluate f and g on the grid 
  grid$xdot <- f(grid$x, grid$y)
  grid$ydot <- g(grid$x, grid$y)
  
  # magnitude and normalization
  grid$magnitude <- sqrt(grid$xdot^2 + grid$ydot^2)
  grid$xdot_n <- grid$xdot / grid$magnitude
  grid$ydot_n <- grid$ydot / grid$magnitude
  
  # compute arrow endpoints scaled to grid spacing
  dx <- (x_range[2] - x_range[1]) / max(1, (n_points - 1))
  dy <- (y_range[2] - y_range[1]) / max(1, (n_points - 1))
  base_step <- min(dx, dy)
  len <- arrow_scale * base_step
  grid$xend <- grid$x + len * grid$xdot_n
  grid$yend <- grid$y + len * grid$ydot_n
  
  # plot
  p <- ggplot(grid, aes(x = x, y = y)) +
    geom_raster(aes(fill = magnitude), interpolate = TRUE) +
    geom_segment(aes(xend = xend, yend = yend),
                 arrow = arrow(length = unit(0.08, "inches")), color = "black", alpha = 0.6) +
    labs(x = "X", y = "Y", fill = "||F(X,Y)||") +
    theme_minimal(base_size = 17)
  
  p <- p + scale_fill_viridis_c(option = "plasma",trans='log10')

  return(p)
}

plot_metric_doublewell <- function(est_long,
                                   true_params,
                                   methods,
                                   method_names,
                                   method_palette,
                                   methods_keep = NULL,
                                   title = NULL,
                                   cap_outliers = FALSE) {
  
  metric_df <- compute_metric_df(est_long, "MSE")
  
  # keep only selected methods
  if (!is.null(methods_keep)) {
    
    metric_df <- metric_df %>%
      dplyr::filter(method %in% methods_keep)
    
    methods_used <- methods_keep
    
  } else {
    methods_used <- methods
  }
  if (cap_outliers) {
    est_long <- est_long %>%
      group_by(parameter, method) %>%
      filter(
        estimate >= quantile(estimate, 0.01, na.rm = TRUE),
        estimate <= quantile(estimate, 0.99, na.rm = TRUE)
      ) %>%
      ungroup()
  }
  
  
  ggplot(metric_df,
         aes(sub_h, value,
             colour = method,
             group = method)) +
    
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    
    facet_wrap(~parameter,
               scales = "free_y",
               labeller = label_parsed) +
    
    scale_colour_manual(values = method_palette, drop = FALSE) +
    scale_x_continuous(breaks = c(0.01,0.02,0.04,0.06,0.08))+
    labs(
      title = title,
      x = "Step size",
      y = "Mean squared error",
      colour = NULL
    ) +
    
    theme_bw(base_size = 14, base_family = "Times") +
    
    theme(
      plot.title = element_text(size = 24, face = "bold", hjust = .5),
      axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      strip.text = element_text(size = 22),
      legend.position = "bottom",
      legend.text = element_text(size = 22),
      axis.text.x = element_text(angle = 30, hjust = 1)
    )
}
plot_box_doublewell <- function(est_long,
                                true_params,
                                param_labels,
                                methods,
                                method_names,
                                method_palette,
                                ratio_value = 1000,
                                methods_keep = NULL,
                                cap_outliers = TRUE,
                                title = NULL) {
  
  est_df <- est_long %>%
    dplyr::filter(ratio == ratio_value)
  
  if (!is.null(methods_keep)) {
    
    est_df <- est_df %>%
      dplyr::filter(method %in% methods_keep)
    
    methods_used <- methods_keep
    
  } else {
    methods_used <- methods
  }
  
  truth_df <- tibble(
    parameter = param_labels,
    true_value = true_params
  )
  
  est_df <- est_df %>%
    left_join(truth_df, by = "parameter")
  
  # optional outlier trimming
  if (cap_outliers) {
    
    est_df <- est_df %>%
      group_by(parameter, method) %>%
      filter(
        estimate >= quantile(estimate, 0.05, na.rm = TRUE),
        estimate <= quantile(estimate, 0.95, na.rm = TRUE)
      ) %>%
      ungroup()
  }
  
  est_df <- est_df %>%
    mutate(
      parameter = factor(parameter, levels = param_labels),
      
      method = factor(
        method,
        levels = methods_used,
        labels = method_names[match(methods_used, methods)]
      )
    )
  
  truth_lines <- tibble(
    parameter = factor(param_labels, levels = param_labels),
    true_value = true_params
  )
  
  ggplot(est_df,
         aes(method, estimate, fill = method)) +
    
    geom_boxplot(width = 0.65, outlier.size = 0.25) +
    
    geom_hline(
      data = truth_lines,
      aes(yintercept = true_value),
      linetype = 2,
      linewidth = 0.8,
      inherit.aes = FALSE
    ) +
    
    facet_wrap(~parameter,
               scales = "free_y",
               labeller = label_parsed) +
    
    scale_fill_manual(values = method_palette, drop = FALSE) +
    
    labs(
      title = title,
      x = NULL,
      y = "Estimate",
      fill = NULL
    ) +
    
    theme_bw(base_size = 14, base_family = "Times") +
    
    theme(
      plot.title = element_text(size = 24, face = "bold", hjust = .5),
      strip.text = element_text(size = 24),
      axis.title.y = element_text(size = 22),
      axis.text.y = element_text(size = 18),
      
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      
      legend.position = "bottom",
      legend.text = element_text(size = 18)
    )
}


# Kræver df med d+1 columns X1,X2,...,Xd
plot_variable_vs_time <- function(df,times, to_time = 50) {
  df$times <- times
  h <- times[2]-times[1]
  steps <- to_time / h
  
  df <- df[1:steps,]
  
  df_long <- df %>%
    pivot_longer(cols = starts_with("X"), 
                 names_to = "variable", 
                 values_to = "value")

  ggplot(df_long, aes(x = times, y = value, color = variable)) +
    geom_line() +
    labs(x = "Time", y = "Value", title = "Each Xi vs Time")
}

plot_steps <- function(df_path, EM_steps, df_steps, b,
                       title = "Linearization around fixpoint",
                       add_nullclines = TRUE,
                       xlim = c(-2, 2), n = 400, plane = c('X','Y')) {
  
  # Helper to ensure x,y columns exist
  ensure_xy <- function(df) {
    df <- as.data.frame(df)
    if (all(c("x", "y") %in% names(df))) return(df)
    if (all(c("X1", "X2") %in% names(df))) {
      names(df)[names(df) == "X1"] <- "x"
      names(df)[names(df) == "X2"] <- "y"
      return(df)
    }
    names(df)[1:2] <- c("x", "y")
    return(df)
  }
  
  df_path  <- ensure_xy(df_path)
  EM_steps <- ensure_xy(EM_steps)
  df_steps <- as.data.frame(df_steps)
  
  if (!"label" %in% names(df_steps)) stop("df_steps must contain a 'label' column (e.g. 'S1','S2','S3').")
  if (!all(c("x0_1","x0_2") %in% names(df_steps))) stop("df_steps must contain 'x0_1' and 'x0_2' columns.")
  
  # Prepare plot data
  path_df  <- transform(df_path, which = "Path realization")
  df_S3    <- transform(subset(df_steps, label == "S3"), which = "Distribution of Strang step")
  EM_df    <- transform(EM_steps, which = "True distribution (EM simulated)")
  init_pts <- unique(df_steps[, c("x0_1", "x0_2")])
  names(init_pts) <- c("x", "y"); init_pts$which <- "Initial points"
  
  # ----- handle b (centers) : allow single point or multiple -----
  # Accept: numeric length-2 c(x,y), data.frame/matrix with two columns, or named list/data.frame with x/y
  make_centers_df <- function(b) {
    if (is.null(b)) return(NULL)
    if (is.numeric(b) && length(b) == 2) {
      return(data.frame(x = b[1], y = b[2], which = "Center of linearization"))
    }
    bdf <- as.data.frame(b)
    # If columns named x and y exist, use them; otherwise use first two columns
    if (all(c("x","y") %in% names(bdf))) {
      centers <- bdf[, c("x","y")]
    } else if (all(c("X1","X2") %in% names(bdf))) {
      centers <- bdf[, c("X1","X2")]
      names(centers) <- c("x","y")
    } else {
      centers <- bdf[, 1:2]
      names(centers) <- c("x","y")
    }
    centers$which <- "Center of linearization"
    # drop NA rows (if any)
    centers <- centers[!is.na(centers$x) & !is.na(centers$y), , drop = FALSE]
    rownames(centers) <- NULL
    return(centers)
  }
  
  center_df <- make_centers_df(b)
  
  # color map
  color_map <- c(
    "Path realization" = "grey80",
    "Distribution of Strang step" = "purple",
    "True distribution (EM simulated)" = "gold",
    "Initial points" = "black",
    "Center of linearization" = "red",
    "Nullclines" = "darkgrey"
  )
  
  # base plot
  p <- ggplot() +
    geom_path(data = path_df, aes(x = x, y = y, color = which),
              linewidth = 0.2, alpha = 0.7) +
    geom_point(data = df_S3, aes(x = x, y = y, color = which), size = 0.001) +
    geom_point(data = EM_df, aes(x = x, y = y, color = which), size = 0.001) +
    geom_point(data = init_pts, aes(x = x, y = y, color = which), size = 1.5)
  
  # add centers if provided (can be multiple rows)
  if (!is.null(center_df) && nrow(center_df) > 0) {
    p <- p + geom_point(data = center_df, aes(x = x, y = y, color = which),
                        size = 3, shape = 21, stroke = 0.6)
  }
  
  # Add nullclines (separate dashed lines, same legend)
  if (add_nullclines) {
    # note: 'par' must exist in calling environment (keeps same behavior as original)
    eps   <- par[1]; alpha <- par[2]; gamma <- par[3]; beta <- par[4]
    xs <- seq(xlim[1], xlim[2], length.out = n)
    df_xnull <- data.frame(x = xs, y = xs - xs^3 + alpha, which = "Nullclines")
    df_ynull <- data.frame(x = xs, y = gamma * xs + beta, which = "Nullclines")
    
    p <- p +
      geom_line(data = df_xnull, aes(x = x, y = y, color = which, linetype = which),
                linewidth = 0.4, alpha = 1) +
      geom_line(data = df_ynull, aes(x = x, y = y, color = which, linetype = which),
                linewidth = 0.4, alpha = 1) +
      scale_linetype_manual(values = c("Nullclines" = "dashed"))
  }
  
  p <- p +
    scale_color_manual(
      name = "",
      values = color_map,
      breaks = c(
        "Initial points",
        "Center of linearization",
        "Distribution of Strang step",
        "True distribution (EM simulated)"
      )
    ) +
    labs(x = plane[1], y = plane[2], title = title) +
    theme_bw() +
    guides(
      color = guide_legend(override.aes = list(size = 5)),
      linetype = "none"
    )
  
  return(p)
}

plot_multiple_steps <- function(df_path, EM_steps,
                                df_steps_list, b_list,
                                titles = NULL,
                                ncol = NULL,  # number of columns in the layout; if NULL a sensible default is chosen
                                xlim = NULL, ylim = NULL,
                                legend_pos = "bottom",
                                legend_text_size=20) {
  
  # validate inputs
  if (!is.list(df_steps_list)) stop("df_steps_list must be a list of data.frames (one per linearization).")
  n_plots <- length(df_steps_list)
  
  # Helper to normalize one b-element (can be numeric length2, data.frame/matrix with 2 cols, NULL)
  normalize_b_element <- function(be) {
    if (is.null(be)) return(NULL)
    # numeric vector length 2 -> single center
    if (is.numeric(be) && length(be) == 2) return(as.numeric(be)[1:2])
    # data.frame or matrix -> return 2-col data.frame (rows = centers)
    if (is.data.frame(be) || is.matrix(be)) {
      bdf <- as.data.frame(be)
      if (ncol(bdf) < 2) stop("b element data.frame/matrix must have at least two columns (x,y).")
      bdf <- bdf[, 1:2, drop = FALSE]
      names(bdf) <- c("x", "y")
      # keep as data.frame (possibly multiple rows)
      return(bdf)
    }
    # list: try to convert to data.frame or numeric
    if (is.list(be)) {
      # if it looks like c(x=..., y=...) or list of numeric pairs
      # try to coerce to data.frame
      try_df <- try(as.data.frame(be), silent = TRUE)
      if (!inherits(try_df, "try-error") && ncol(try_df) >= 2) {
        df2 <- try_df[, 1:2, drop = FALSE]
        names(df2) <- c("x", "y")
        return(df2)
      }
      # if it's numeric-ish vector inside list
      if (length(unlist(be)) == 2 && all(sapply(unlist(be), is.numeric))) {
        return(as.numeric(unlist(be))[1:2])
      }
    }
    stop("Unsupported b_list element type. Each element must be NULL, numeric length-2, or a 2-column data.frame/matrix (rows = centers).")
  }
  
  # If b_list provided as matrix/data.frame with nrow == n_plots, treat each row as a single center
  if (is.matrix(b_list) || is.data.frame(b_list)) {
    if (nrow(b_list) != n_plots) {
      stop("When b_list is a matrix/data.frame it must have one row per df_steps_list element (nrow(b_list) == length(df_steps_list)).")
    }
    # convert rows to list of numeric pairs
    b_list <- lapply(seq_len(nrow(as.data.frame(b_list))), function(i) {
      row <- as.numeric(as.data.frame(b_list)[i, 1:2])
      row
    })
  } else {
    # if not a data.frame/matrix, coerce to list if needed
    if (!is.list(b_list)) {
      b_list <- as.list(b_list)
    }
    if (length(b_list) != n_plots) stop("Length of b_list must match length of df_steps_list.")
  }
  
  # Normalize each element (can become numeric length2 or a data.frame of centers or NULL)
  b_list_norm <- vector("list", n_plots)
  for (i in seq_len(n_plots)) {
    b_list_norm[[i]] <- normalize_b_element(b_list[[i]])
  }
  
  if (is.null(titles)) {
    titles <- paste("Linearization", seq_len(n_plots))
  } else {
    if (length(titles) != n_plots) stop("Length of titles must match number of items in df_steps_list.")
  }
  
  # sensible default for ncol
  if (is.null(ncol)) {
    ncol <- if (n_plots == 1) 1 else if (n_plots == 2) 2 else ceiling(sqrt(n_plots))
  }
  
  # build individual plots
  plots <- vector("list", n_plots)
  for (i in seq_len(n_plots)) {
    # pass b_list_norm[[i]] directly to plot_steps(); plot_steps handles single or multiple centers
    plots[[i]] <- plot_steps(df_path, EM_steps, df_steps_list[[i]], b_list_norm[[i]], title = titles[[i]])
    # apply axis limits if provided
    if (!is.null(xlim) || !is.null(ylim)) {
      plots[[i]] <- plots[[i]] + coord_cartesian(xlim = xlim, ylim = ylim)
    }
  }
  
  # combine with patchwork, collect legends and set legend position
  combined <- wrap_plots(plots, ncol = ncol) +
    plot_layout(guides = "collect") &
    theme_bw(base_size=16,base_family = 'Times')+
    theme(legend.position=legend_pos,legend.key.size= unit(1.6, "cm"),legend.text = element_text(size = 16),)+
    
    theme(strip.text.x = element_text(size = 16, color = "black", face = "bold"),
          strip.text.y = element_text(size = 16, color = "black", face = "bold"),
          strip.text = element_text(face="bold"),
          strip.background = element_rect(fill="grey")) # facet strips
  
  return(combined)
}


plot_top_k_b_choices <- function(paths, x0, b_results, par, h, k = 5,
                                 xlim = c(-2, 2), ylim = c(-0.2, 2.1)) {
  library(ggplot2)
  
  # Helper: ensure x,y column names
  ensure_xy <- function(df) {
    df <- as.data.frame(df)
    if (all(c("x", "y") %in% names(df))) return(df)
    if (all(c("X1", "X2") %in% names(df))) {
      names(df)[names(df) == "X1"] <- "x"
      names(df)[names(df) == "X2"] <- "y"
      return(df)
    }
    names(df)[1:2] <- c("x", "y")
    return(df)
  }
  
  # Convert x0 into a data frame with proper column names
  x0_df <- as.data.frame(t(x0))
  names(x0_df) <- c("x", "y")
  x0_df$which <- "Initial point"
  
  # Nullcline data
  xs <- seq(xlim[1], xlim[2], length.out = 400)
  bs <- seq(xlim[1], xlim[2], length.out = 400)
  df_xnull <- data.frame(x = xs, y = xs - xs^3 + par[2], which = "Nullclines")
  df_ynull <- data.frame(x = xs, y = par[3] * xs + par[4], which = "Nullclines")
  
  # Polynomial overlay using x = x0[1]
  x_val <- x0[1]
  y_val <- x0[2]
  eps <- par[1]
  alpha <- par[2]
  gamma <- par[3]
  beta <- par[4]
  DFx <- matrix(c(
      (-3 * x_val^2 + 1) / eps,  -1 / eps,
      gamma,-1
    ),
    nrow = 2,
    byrow = TRUE)
  eig <- eigen(DFx)

  idx_s <- which.min(abs(Re(eig$values)))
  idx_f <- which.max(abs(Re(eig$values)))

  v_s <- Re(eig$vectors[, idx_s])
  v_f <- Re(eig$vectors[, idx_f])

  x_tilde_s <- x_val-v_s[1]*(gamma*x_val+beta-y_val)/(gamma*v_s[1]-v_s[2])
  y_tilde_s <- y_val-v_s[2]*(gamma*x_val+beta-y_val)/(gamma*v_s[1]-v_s[2])
  Fvec <- c(
    (x_val -x_val^3 + alpha - y_val)/eps,
    gamma * x_val+ beta - y_val
  )
  x_tilde_t <- x_val-Fvec[1]*(gamma*x_val+beta-y_val)/(gamma*Fvec[1]-Fvec[2])
  y_tilde_t <- y_val-Fvec[2]*(gamma*x_val+beta-y_val)/(gamma*Fvec[1]-Fvec[2])
  r1_t <- polyroot(c(alpha-beta-(gamma-1)*x_tilde_t, 0, -3*x_tilde_t, 2))
  r2_t <- polyroot(c(alpha*gamma-beta-(gamma-1)*y_tilde_t,0,3*(beta-y_tilde_t),2*gamma))
  real_r1_t <- Re(r1_t[abs(Im(r1_t)) < 1e-8])

  project_to_cubic <- function(x_val, y_val, v_f, a_const) {
    v1 <- v_f[1]
    v2 <- v_f[2]

    coefs <- c(
      y_val - x_val + x_val^3 - a_const,   # t^0
      v1 - v2 - 3 * v1 * x_val^2,          # t^1
      3 * v1^2 * x_val,                    # t^2
      -v1^3                                # t^3
    )

    roots <- polyroot(coefs)
    roots_real <- Re(roots)[abs(Im(roots)) < 1e-8]

    if (length(roots_real) == 0) stop("No real intersection found.")

    t_star <- roots_real[which.min(abs(roots_real))]

    list(
      x_tilde = x_val - t_star * v1,
      y_tilde = y_val - t_star * v2,
      t = t_star
    )
  }

  proj_f <- project_to_cubic(x_val, y_val, v_f, a_const = alpha)
  x_tilde_f <- proj_f$x_tilde
  y_tilde_f <- proj_f$y_tilde
  proj_t <- project_to_cubic(x_val, y_val,Fvec, a_const = alpha)
  x_tilde_t <- proj_t$x_tilde
  y_tilde_t <- proj_t$y_tilde
  # 
  # 
  # 
  # r1_s <- polyroot(c(alpha-beta-(gamma-1)*x_tilde_s, 0, -3*x_tilde_s, 2))
  # r2_s <- polyroot(c(alpha*gamma-beta-(gamma-1)*y_tilde_s,0,3*(beta-y_tilde_s),2*gamma))
  # 
  #r1_f <- polyroot(c(alpha-beta-(gamma-1)*x_tilde_f, 0, -3*x_tilde_f, 2))
  # r2_f <- polyroot(c(alpha*gamma-beta-(gamma-1)*y_tilde_f,0,3*(beta-y_tilde_f),2*gamma))
  # DFx_tilde_f <- matrix(c(
  #     (-3 * x_tilde_f^2 + 1) / eps,  -1 / eps,
  #     gamma,-1
  #   ),
  #   nrow = 2,
  #   byrow = TRUE)
  # 
  # eig_at_f <- eigen(DFx_tilde_f)
  # 
  # idx_s_at_f <- which.min(abs(Re(eig_at_f$values)))
  # idx_f_at_f <- which.max(abs(Re(eig_at_f$values)))
  # 
  # v_s_at_f <- Re(eig_at_f$vectors[, idx_s_at_f])
  # v_f_at_f <- Re(eig_at_f$vectors[, idx_f_at_f])
  # 
  # x_tilde_m <- x_tilde_f-v_s_at_f[1]*(gamma*x_tilde_f+beta-y_tilde_f)/(gamma*v_s_at_f[1]-v_s_at_f[2])
  # y_tilde_m <- y_tilde_f-v_s_at_f[2]*(gamma*x_tilde_f+beta-y_tilde_f)/(gamma*v_s_at_f[1]-v_s_at_f[2])
  # 
  # r1_m <- polyroot(c(alpha-beta-(gamma-1)*x_tilde_m, 0, -3*x_tilde_m, 2))
  # r2_m <- polyroot(c(alpha*gamma-beta-(gamma-1)*y_tilde_m,0,3*(beta-y_tilde_m),2*gamma))
  # 
  # real_r1_s <- Re(r1_s[abs(Im(r1_s)) < 1e-8])
  #real_r1_f <- Re(r1_f[abs(Im(r1_f)) < 1e-8])
  # real_r1_m <- Re(r1_m[abs(Im(r1_m)) < 1e-8])
  # 
  # 
  # # Top k choices
  top_k_df <- b_results
  top_k_df <- top_k_df %>%
    mutate(
      error = ifelse(
        is.finite(error),
        pmin(pmax(error, 1e-8), 1e4),
        1e4
      )
    )
  top_k_df <- top_k_df[order(top_k_df$error, decreasing = TRUE), ]
  best_idx <- which.min(top_k_df$error)
  # bt <- b_tilde(top_k_df$b1)
  # 
  # top_k_bt_df <- data.frame(
  #   b1 = bt[,1],
  #   b2 = bt[,2],
  #   error = top_k_df$error
  # )
  # 
  # pt1_s <- b_tilde(b_inv_roots1_s)
  # pt1_f <- b_tilde(b_inv_roots1_f)
  # 
  # pt2_s <- b_tilde(b_inv_roots2_s)
  # pt2_f <- b_tilde(b_inv_roots2_f)
  # 
  # library(expm)
  # 
  # Jman <- DF(c(x_tilde_f, y_tilde_f))
  # 
  # Fman <- c(
  #   (x_tilde_f -x_tilde_f^3 + alpha - y_tilde_f)/eps,
  #   gamma * x_tilde_f+ beta - y_tilde_f
  # )
  # n_step <- solve(Jman,Fman)
  # 
  # x_tilde_n <- x_tilde_f-n_step[1]*(gamma*x_tilde_f+beta-y_tilde_f)/(gamma*n_step[1]-n_step[2])
  # y_tilde_n <- y_tilde_f-n_step[2]*(gamma*x_tilde_f+beta-y_tilde_f)/(gamma*n_step[1]-n_step[2])
  # 
  # 
  # 
  # 
  # r1_n <- polyroot(c(alpha-beta-(gamma-1)*x_tilde_n, 0, -3*x_tilde_n, 2))
  # r2_n <- polyroot(c(alpha*gamma-beta-(gamma-1)*y_tilde_n,0,3*(beta-y_tilde_n),2*gamma))
  # 
  # real_r1_n <- Re(r1_n[abs(Im(r1_n)) < 1e-8])

  # Build the plot
  p <- ggplot() +
    geom_path(data = ensure_xy(paths[[1]]),
              aes(x = x, y = y),
              color = "grey80", linewidth = 0.4) +
    
    # Nullclines
    geom_line(data = df_xnull, aes(x = x, y = y, linetype = which),
              color = "darkgrey",
              show.legend = FALSE) +
    geom_line(data = df_ynull, aes(x = x, y = y, linetype = which),
              color = "darkgrey",
              show.legend = FALSE) +
    
    # Polynomial overlay
    # geom_line(data = df_poly, aes(x = b1, y = poly_val2),
    #           color = "green", linewidth = 0.6) +
    # geom_line(data = df_poly, aes(x = b1, y = poly_val1),
    #           color = "purple", linewidth = 0.6) +
    # geom_line(data = df_poly, aes(x = b1, y = poly_val3),
    #           color = "red", linewidth = 0.6) +
    # geom_line(data = df_poly, aes(x = b1, y = poly_val4),
    #           color = "red", linewidth = 0.6) +
    # geom_line(data = df_poly, aes(x = b1, y = poly_val5),
    #           color = "red", linewidth = 0.6) +
    # geom_vline(xintercept = b_tilde,color = "red", linewidth = 0.6) +
    # geom_vline(xintercept = b_trust[1],color = "red", linewidth = 0.6) +
    # geom_hline(yintercept = b_trust[2],color = "red", linewidth = 0.6) +
    #geom_vline(xintercept = real_r1_s,color = "red", linewidth = 0.1) +
    #geom_vline(xintercept = x_tilde_f,color = "red", linewidth = 0.3,linetype='dashed') +
    geom_vline(xintercept = x_tilde_t,color = "red", linewidth = 0.3,linetype='dashed') +
    
    #geom_vline(xintercept = real_r1_t,color = "black", linewidth = 0.1) +
    #geom_vline(xintercept = real_r1_n,color = "gold", linewidth = 0.5) +
    # geom_vline(xintercept = real_r1_m,color = "red", linewidth = 0.5,linetype='dashed') +
    #geom_vline(xintercept = real_r1_t,color = "brown", linewidth = 0.5) +
    
    
    

    
    # geom_segment(
    #   aes(
    #     x = x_val,
    #     y = y_val,
    #     xend = x_tilde_f,
    #     yend = y_tilde_f
    #   ),
    #   arrow = arrow(length = unit(0.15, "cm")),
    #   color = "red",
    #   linewidth = 0.5,
    # )+
    geom_segment(
      aes(
        x = x_val,
        y = y_val,
        xend = x_tilde_t,
        yend = y_tilde_t
      ),
      ,
      arrow = arrow(length = unit(0.15, "cm")),
      color = "red",
      linewidth = 0.5,
    )+
    # geom_segment(
    #   aes(
    #     x = x_tilde_f,
    #     y = y_tilde_f,
    #     xend = x_tilde_m,
    #     yend = y_tilde_m
    #   ),
    #   color = "red",
    #   linewidth = 0.5,
    # )+
    # geom_segment(
    #   aes(
    #     x = x_val,
    #     y = y_val,
    #     xend = x_tilde_f,
    #     yend = y_tilde_f
    #   ),
    #   color = "blue",
    #   linewidth = 0.5,
    #   linetype = "dashed"
    # ) +
    # geom_segment(
    #   aes(
    #     x = x_val,
    #     y = y_val,
    #     xend = x_tilde_n,
    #     yend = y_tilde_n
    #   ),
    #   arrow = arrow(length = unit(0.15, "cm")),
    #   color = "gold",
    #   linewidth = 0.6,
    #   linetype='dashed'
    # )+
    # geom_segment(
    #   aes(
    #     x = x_val,
    #     y = y_val,
    #     xend = x_tilde_t,
    #     yend = y_tilde_t
    #   ),
    #   arrow = arrow(length = unit(0.15, "cm")),
    #   color = "brown",
    #   linewidth = 0.3,
    #   linetype = 'dashed'
    # )+
    #geom_hline(yintercept = b_OU[2],color = "red", linewidth = 0.6) +
    # geom_vline(xintercept = N_roots[1],color = "red", linewidth = 0.6) +
    # geom_vline(xintercept = N_roots[2],color = "red", linewidth = 0.6) +
    # geom_vline(xintercept = N_roots[3],color = "red", linewidth = 0.6) +

    #Top k points
    geom_point(data = top_k_df,
               aes(x = b1, y = b2, color = error),
               size = 0.8) +
    #best_point
    geom_point(
      data = top_k_df[best_idx, ],
      aes(x = b1, y = b2),
      shape = 8,
      size = 2,
      color = "yellow",
      stroke = 1.2
    )+
    
    # geom_point(
    #   data = top_k_bt_df,
    #   aes(x = b1, y = b2, color = error),
    #   shape = 1, size = 1
    # )+
    
    # Initial point
    geom_point(data = x0_df,
               aes(x = x, y = y, shape = which),
               color = "black", size = 2) +
    
    scale_color_gradientn(
      colors = c("yellow", "red"),
      trans = "log10",
      name = "Prediction Error"
    )+
    scale_shape_manual(name = "", values = c("Initial point" = 19)) +
    scale_linetype_manual(name = "", values = c("Nullclines" = "dashed")) +
    
    labs(
      #title = paste0("Initial point x = (",round(x0[1], 2), ", ", round(x0[2], 2)),
      x = 'x',
      y = 'y'
    ) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    theme_bw()
  
  return(p)
}


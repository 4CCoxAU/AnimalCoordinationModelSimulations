# ============================================================================
# Functions to create the figures
# ============================================================================

# ============================================================================
# Helper: build vocalization dataframe
# ============================================================================

make_voc_df <- function(sim) {
  bind_rows(
    tibble(Individual = "Individual 1",
           StartTime = sim$voc_times_1, Duration = sim$voc_durations_1),
    tibble(Individual = "Individual 2",
           StartTime = sim$voc_times_2, Duration = sim$voc_durations_2)
  ) %>%
    mutate(EndTime = StartTime + Duration)
}

# ============================================================================
# Row 1: Vocalization timelines (0–30s window)
# ============================================================================

window_low <- 0
window_high <- 30
window <- c(window_low, window_high)

make_timeline <- function(sim, toggle = element_blank(), time = "Time (s)", l_position = c(1.1, 0.5)) {
  voc_df <- make_voc_df(sim)
  
  voc_df %>%
    filter(StartTime >= window[1], StartTime <= window[2]) %>%
    mutate(y_pos = ifelse(Individual == "Individual 1", 0.57, 0.43)) %>%
    mutate(Individual = case_when(Individual == "Individual 1" ~ "Caller A",
                                  Individual == "Individual 2" ~ "Caller B")) %>%
    ggplot() +
    geom_rect(aes(xmin = StartTime, xmax = EndTime,
                  ymin = y_pos - 0.05, ymax = y_pos + 0.05,
                  fill = Individual),
              alpha = 0.9, color = "black", linewidth = 0.3) +
    scale_fill_manual(
      values = c("Caller A" = "#fca636", "Caller B" = "#446455")
    ) +
    scale_y_continuous(breaks = c(0.43, 0.57),
                       labels = c(" ", " "),
                       #limits = c(0.32, 0.68)
    ) +
    scale_x_continuous(breaks = c(seq(0, window_high, 10)),
                       limits = c(0, window_high)) +
    labs(x = time, y = "") +
    theme_combined +
    theme(panel.grid.major.y = element_blank(),
          plot.title = element_blank(),
          axis.text.y = toggle,
          axis.text.x = element_text(size = 12, color = "black"),
          axis.title.x = element_text(size = 12, color = "black"),
          legend.title = element_blank(),
          legend.position = l_position,
          legend.key.size = unit(1, "lines"),
          legend.text = element_text(size = 10),
          legend.background = element_rect(fill = "transparent", color = NA),
          legend.key = element_rect(fill = "transparent", color = NA),
    )
}

# ============================================================================
# Shared theme for all panels
# ============================================================================

theme_combined <- theme_minimal(base_size = 10) +
  theme(
    plot.title       = element_text(size = 9, face = "bold", hjust = 0),
    axis.title       = element_text(size = 8),
    axis.text        = element_text(size = 7),
    legend.position  = "none",
    panel.grid.minor = element_blank(),
    plot.margin      = margin(4, 6, 4, 6)
  )


# ===============================
# Make Inter-Call Interval Plots
# ===============================

make_ici <- function(sim, ref_lines = NULL, ref_labels = NULL,
                     ref_colors = NULL, x_upper = 2.5, x_upper_lim = x_upper, x_lower_lim = 0, breaks = 1, bins = 50, 
                     toggle = element_blank(), toggle_x = "Inter-Call Interval (s)") {
  ici_df <- bind_rows(
    tibble(ici = diff(sim$voc_times_1)),
    tibble(ici = diff(sim$voc_times_2))
  )
  
  p <- ici_df %>%
    filter(ici <= x_upper & ici > 0) %>%
    ggplot(aes(x = ici)) +
    geom_histogram(aes(y = after_stat(density)), bins = bins,
                   fill = "#005CAB", alpha = 0.7,
                   color = "black", linewidth = 0.3) +
    geom_density(linewidth = 0.8) +
    coord_cartesian(xlim = c(x_lower_lim, x_upper_lim)) +
    scale_x_continuous(breaks = seq(x_lower_lim, x_upper_lim, breaks)) +
    labs(x = toggle_x, y = "Density") +
    theme_combined +
    theme(plot.title = element_blank(),
          axis.text.y  = toggle,
          axis.title.y = toggle,
          axis.text.x  = element_text(size = 12, color = "black"),
          axis.title.x = element_text(size = 12, color = "black"),
          panel.grid.major.y = element_blank()
    )
  
  if (!is.null(ref_lines)) {
    for (j in seq_along(ref_lines)) {
      p <- p +
        geom_vline(xintercept = ref_lines[j], linetype = "dashed",
                   color = ref_colors[j], linewidth = 0.8) +
        annotate("text", x = ref_lines[j] + 0.03, y = Inf,
                 label = ref_labels[j], vjust = 2, hjust = 0,
                 color = ref_colors[j], fontface = "bold", size = 3)
    }
  }
  p
}

# ===============================
# Make Phase Dynamics Plots
# ===============================

# --- Model 1: Phase dynamics (oscillatory jitter) ---
make_phase_dynamics_m1 <- function(sim, window) {
  p <- sim$params
  expected_T <- mean(c(p$T1, p$T2))
  
  all_calls <- bind_rows(
    tibble(time = sim$voc_times_1, caller = 1),
    tibble(time = sim$voc_times_2, caller = 2)
  ) %>% arrange(time)
  
  phase_diffs <- numeric()
  phase_times <- numeric()
  
  if (nrow(all_calls) > 4) {
    for (k in 3:nrow(all_calls)) {
      prev_1 <- all_calls %>%
        filter(caller == 1, time < all_calls$time[k]) %>%
        pull(time) %>% tail(1)
      prev_2 <- all_calls %>%
        filter(caller == 2, time < all_calls$time[k]) %>%
        pull(time) %>% tail(1)
      
      if (length(prev_1) > 0 && length(prev_2) > 0) {
        t_now <- all_calls$time[k]
        ph1 <- 2 * pi * (t_now - prev_1) / expected_T
        ph2 <- 2 * pi * (t_now - prev_2) / expected_T
        pd  <- (ph1 - ph2) %% (2 * pi)
        phase_diffs <- c(phase_diffs, pd)
        phase_times <- c(phase_times, t_now)
      }
    }
  }
  
  tibble(Time = phase_times, PhaseDiff = phase_diffs) %>%
    ggplot(aes(x = Time, y = PhaseDiff)) +
    geom_hline(yintercept = pi, linetype = "dashed",
               color = "gray40", linewidth = 0.6, alpha = 0.7) +
    geom_line(linewidth = 0.6, color = "black") +
    scale_x_continuous(breaks = seq(0, window[2], by = 5), limits = window) +
    scale_y_continuous(
      breaks = c(0, pi/2, pi, 3*pi/2, 2*pi),
      #labels = c("0", "\u03c0/2", "\u03c0", "3\u03c0/2", "2\u03c0"),
      labels = expression(0, pi/2, pi, 3*pi/2, 2*pi),
      limits = c(0, 2 * pi)
    ) +
    labs(x = " ", y = "Phase Difference",
         title = "C) Phase Dynamics") +
    theme_combined +
    theme(panel.grid.minor.y = element_blank(),
          plot.title = element_blank(),
          axis.text.y  = element_text(size = 12, color = "black"),
          axis.text.x  = element_text(size = 12, color = "black"),
          axis.title.x = element_text(size = 12, color = "black"),
          axis.title.y = element_blank()
    )
}

# --- Model 2: Phase dynamics (smooth convergence) ---
make_phase_dynamics_m2 <- function(sim, window) {
  p <- sim$params
  
  phase_diff_raw <- sim$phases[, 1] - sim$phases[, 2]
  phase_diff_mod <- phase_diff_raw %% (2 * pi)
  target_phase   <- ifelse(p$K < 0, pi, 0)
  
  tibble(Time = sim$times, PhaseDiff = phase_diff_mod) %>%
    filter(row_number() %% 100 == 1) %>%
    ggplot(aes(x = Time, y = PhaseDiff)) +
    geom_hline(yintercept = target_phase, linetype = "dashed",
               color = "gray40", linewidth = 0.6, alpha = 0.7) +
    geom_line(linewidth = 0.6, color = "black") +
    scale_x_continuous(breaks = seq(0, window[2], by = 5), limits = window) +
    scale_y_continuous(
      breaks = c(0, pi/2, pi, 3*pi/2, 2*pi),
      #labels = c("0", "\u03c0/2", "\u03c0", "3\u03c0/2", "2\u03c0"),
      labels = expression(0, pi/2, pi, 3*pi/2, 2*pi),
      limits = c(0, 2 * pi)
    ) +
    labs(x = "Time (s)", y = "Phase Difference",
         title = "C) Phase Dynamics") +
    theme_combined +
    theme(panel.grid.minor.y = element_blank(),
          plot.title = element_blank(),
          axis.text.y  = element_text(size = 12, color = "black"),
          axis.text.x  = element_text(size = 12, color = "black"),
          axis.title.x = element_text(size = 12, color = "black"),
          axis.title.y = element_blank()
    )
}

# --- Model 3: State occupancy (colored bands) ---
make_phase_dynamics_m3 <- function(sim, window) {
  expected_T <- mean(c(diff(sim$voc_times_1), diff(sim$voc_times_2)), na.rm = TRUE)
  
  all_calls <- bind_rows(
    tibble(time = sim$voc_times_1, caller = 1),
    tibble(time = sim$voc_times_2, caller = 2)
  ) %>% arrange(time)
  
  phase_diffs <- numeric()
  phase_times <- numeric()
  
  if (nrow(all_calls) > 4) {
    for (k in 3:nrow(all_calls)) {
      prev_1 <- all_calls %>%
        filter(caller == 1, time < all_calls$time[k]) %>%
        pull(time) %>% tail(1)
      prev_2 <- all_calls %>%
        filter(caller == 2, time < all_calls$time[k]) %>%
        pull(time) %>% tail(1)
      
      if (length(prev_1) > 0 && length(prev_2) > 0) {
        t_now   <- all_calls$time[k]
        ph1     <- 2 * pi * (t_now - prev_1) / expected_T
        ph2     <- 2 * pi * (t_now - prev_2) / expected_T
        phase_diffs <- c(phase_diffs, (ph1 - ph2) %% (2 * pi))
        phase_times <- c(phase_times, t_now)
      }
    }
  }
  
  tibble(Time = phase_times, PhaseDiff = phase_diffs) %>%
    filter(Time >= window[1], Time <= window[2]) %>%
    ggplot(aes(x = Time, y = PhaseDiff)) +
    geom_hline(yintercept = pi, linetype = "dashed", color = "gray40", linewidth = 0.6, alpha = 0.7) +
    geom_point(size = 2, alpha = 0.9, color = "black") +
    #geom_line(linewidth = 0.6, color = "black") +
    scale_x_continuous(breaks = seq(0, window[2], by = 5), limits = window) +
    scale_y_continuous(
      breaks = c(0, pi/2, pi, 3*pi/2, 2*pi),
      #labels = c("0", "\u03c0/2", "\u03c0", "3\u03c0/2", "2\u03c0"),
      labels = expression(0, pi/2, pi, 3*pi/2, 2*pi),
      limits = c(0, 2 * pi)
    ) +
    labs(x = " ", y = "Phase Difference", title = "C) Phase Dynamics") +
    theme_combined +
    theme(panel.grid.minor.y = element_blank(),
          plot.title = element_blank(),
          axis.text.y  = element_text(size = 12, color = "black"),
          axis.text.x  = element_text(size = 12, color = "black"),
          axis.title.x = element_text(size = 12, color = "black"),
          axis.title.y = element_blank()
    )
}

# =========================
# Add labels for Figure 1:
# =========================
row_labels <- cowplot::plot_grid(
  cowplot::ggdraw() + cowplot::draw_label("Coordination\nTimeline", 
                                          fontface = "bold", size = 20, angle = 90),
  NULL,
  cowplot::ggdraw() + cowplot::draw_label("Inter-Call\nIntervals",      
                                          fontface = "bold", size = 20, angle = 90),
  NULL,
  cowplot::ggdraw() + cowplot::draw_label("Phase\nDynamics",        
                                          fontface = "bold", size = 20, angle = 90),
  NULL,
  ncol = 1,
  rel_heights = c(0.3, 0.05, 0.3, 0.05, 0.3, 0.05)
)

# Column headers with Colour
make_col_header <- function(title, equation) {
  cowplot::ggdraw() +
    cowplot::draw_label(title, fontface = "bold", size = 20, 
                        x = 0.5, y = 0.65, hjust = 0.5)
}

# Column headers with Colour
make_colored_col_header <- function(label, color) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label,
             fontface = "bold", size = 5.5, color = color,
             hjust = 0.5, vjust = 0.5) +
    theme_void()
}

header_row <- cowplot::plot_grid(
  NULL,
  make_col_header("Event-Triggered"),
  NULL,
  make_col_header("Tempo-Tracking"), 
  NULL,
  make_col_header("State-Dependent"),
  NULL,
  ncol = 7,
  rel_widths = c(0.05, 0.3, 0.05, 0.3, 0.05, 0.3, 0.05)
)

# Make panel title
function(label, size = 7) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, hjust = 0.5,
             fontface = "bold", label = label,
             size = 7, color = "black") +
    theme_void()
}


# Figure 2, raw data:
tax_data <- data.frame(
  model = c(
    rep("Event-Triggered", 5),
    rep("Tempo-Tracking",  5),
    rep("State-Dependent", 4)
  ),
  taxon = c(
    "Orthoptera", "Anurans", "Primates", "Herpestids", "Pinnipeds",
    "Orthoptera", "Anurans", "Aves", "Primates", "Pinnipeds",
    "Anurans", "Aves", "Primates", "Cetaceans"
  ),
  n = c(
    9, 4, 2, 1, 1,
    3, 6, 3, 5, 2,
    1, 2, 1, 1
  )
)

# Factor ordering
tax_data$model <- factor(tax_data$model,
                         levels = c("State-Dependent", "Tempo-Tracking", "Event-Triggered"))
tax_data$taxon <- factor(tax_data$taxon,
                         levels = c("Orthoptera", "Anurans", "Aves", "Primates",
                                    "Cetaceans", "Pinnipeds", "Herpestids"))

# Full grid with zeros
full_grid <- expand.grid(
  model = levels(tax_data$model),
  taxon = levels(tax_data$taxon)
)

tax_full <- full_grid %>%
  left_join(tax_data, by = c("model", "taxon")) %>%
  mutate(n = replace_na(n, 0))


compute_8_stats <- function(voc_times_1, voc_times_2) {
  
  # ICI statistics
  ici <- c(diff(voc_times_1), diff(voc_times_2))
  ici <- ici[ici > 0.05 & ici < 8]
  if (length(ici) < 20) return(rep(NA, 8))
  
  n   <- length(ici)
  m2  <- mean((ici - mean(ici))^2)
  m3  <- mean((ici - mean(ici))^3)
  m4  <- mean((ici - mean(ici))^4)
  skw <- m3 / (m2^1.5)
  krt <- m4 / (m2^2) - 3
  bc  <- (skw^2 + 1) / (krt + 3 * ((n-1)^2 / ((n-2)*(n-3))))
  
  s1 <- mean(ici)
  s2 <- sd(ici)
  s3 <- bc
  s4 <- mean(ici < mean(ici) * 0.5)
  s5 <- dip.test(ici)$statistic
  
  # Phase statistics
  expected_T <- mean(ici)
  all_calls  <- data.frame(
    time   = c(voc_times_1, voc_times_2),
    caller = c(rep(1, length(voc_times_1)), rep(2, length(voc_times_2)))
  )
  all_calls   <- all_calls[order(all_calls$time), ]
  phase_diffs <- numeric()
  
  for (k in 3:nrow(all_calls)) {
    prev_1 <- tail(all_calls$time[all_calls$caller == 1 &
                                    all_calls$time < all_calls$time[k]], 1)
    prev_2 <- tail(all_calls$time[all_calls$caller == 2 &
                                    all_calls$time < all_calls$time[k]], 1)
    if (length(prev_1) > 0 && length(prev_2) > 0) {
      t_now <- all_calls$time[k]
      ph1   <- 2 * pi * (t_now - prev_1) / expected_T
      ph2   <- 2 * pi * (t_now - prev_2) / expected_T
      phase_diffs <- c(phase_diffs, (ph1 - ph2) %% (2 * pi))
    }
  }
  
  if (length(phase_diffs) < 10) return(rep(NA, 8))
  
  s6 <- mean(abs(phase_diffs - pi))
  s7 <- var(phase_diffs)
  s8 <- cor(phase_diffs[-length(phase_diffs)], phase_diffs[-1])
  
  c(s1, s2, s3, s4, s5, s6, s7, s8)
}

compute_dist <- function(ref_stats, obs_stats) {
  denom <- pmax(abs(obs_stats), 1e-6)
  apply(ref_stats, 1, function(row)
    sqrt(sum(((row - obs_stats) / denom)^2)))
}

simulate_ici_event <- function(best, n=n_ppc, x_upper=5) {
  out <- list()
  for (i in 1:n) {
    sim <- simulate_event_triggered(
      duration=1000, T1=best$T[i], T2=best$T[i],
      rho=best$rho[i], d_call=best$d_call[i],
      noise_sd=best$noise_sd[i], duration_noise_sd=best$noise_sd[i]
    )
    ici <- c(diff(sim$voc_times_1), diff(sim$voc_times_2))
    out[[i]] <- ici[ici > 0.05 & ici < x_upper]
  }
  out
}

simulate_ici_tempo <- function(best, n=n_ppc, x_upper=3) {
  out <- list()
  for (i in 1:n) {
    sim <- simulate_coupled_oscillator(
      duration=50, freq1_init=best$freq1[i], freq2_init=best$freq2[i],
      K=-best$K[i], epsilon=best$epsilon[i], dt=0.02,
      phase_noise_sd=best$phase_noise_sd[i], freq_noise_sd=0.001,
      d_call=0.2, duration_noise_sd=0.1, ici_noise_sd=0.05,
      initial_phase_diff=runif(1,0,pi)
    )
    ici <- c(diff(sim$voc_times_1), diff(sim$voc_times_2))
    out[[i]] <- ici[ici > 0.05 & ici < x_upper]
  }
  out
}

simulate_ici_state <- function(best, n=n_ppc, x_upper=8) {
  out <- list()
  for (i in 1:n) {
    sim <- simulate_state_dependent(
      duration=300, lambda_active=best$lambda_active[i],
      lambda_listening=0.05, lambda_refractory=0.001,
      p_active_to_listening_base=best$p_off[i],
      p_listening_to_active_base=best$p_on[i],
      p_active_to_listening_partner=0.3,
      p_listening_to_active_partner=0.005,
      refrac_dur=0.3, partner_window=0.5,
      d_call=0.2, duration_noise_sd=0.1, dt=0.01
    )
    ici <- c(diff(sim$voc_times_1), diff(sim$voc_times_2))
    out[[i]] <- ici[ici > 0.05 & ici < x_upper]
  }
  out
}

make_ppc_panel <- function(ppc_list, obs_voc1, obs_voc2,
                           color, title, x_upper, show_y = FALSE,
                           is_correct = FALSE) {
  
  obs_ici <- c(diff(obs_voc1), diff(obs_voc2))
  obs_ici <- obs_ici[obs_ici > 0.05 & obs_ici < x_upper]
  
  x_grid <- seq(0.05, x_upper, length.out=200)
  
  dens_list <- lapply(ppc_list, function(ici) {
    ici <- ici[ici > 0.05 & ici < x_upper]
    if (length(ici) < 5) return(rep(NA, 200))
    d <- density(ici, from=0.05, to=x_upper, n=200)
    approx(d$x, d$y, xout=x_grid)$y
  })
  
  dens_mat <- do.call(rbind, dens_list)
  df_band  <- data.frame(
    x   = x_grid,
    med = apply(dens_mat, 2, median,   na.rm=TRUE),
    lo  = apply(dens_mat, 2, quantile, 0.05, na.rm=TRUE),
    hi  = apply(dens_mat, 2, quantile, 0.95, na.rm=TRUE)
  )
  
  alpha_ribbon <- if (is_correct) 0.5 else 0.1
  alpha_line   <- if (is_correct) 1.0  else 0.1
  width_line   <- if (is_correct) 1.0  else 0.5
  
  df_band$lo_smooth  <- predict(loess(lo  ~ x, data=df_band, span=0.1))
  df_band$hi_smooth  <- predict(loess(hi  ~ x, data=df_band, span=0.1))
  df_band$med_smooth <- predict(loess(med ~ x, data=df_band, span=0.1))
  
  ggplot() +
    geom_ribbon(data=df_band, aes(x=x, ymin=lo_smooth, ymax=hi_smooth),
                fill=color, alpha=alpha_ribbon) +
    #geom_line(data=df_band, aes(x=x, y=med), color=color, linewidth=1.0, alpha=alpha_line) +
    #geom_point(data=df_band, aes(x=x, y=med), color=color, linewidth=1.0, alpha=alpha_line) +
    geom_smooth(data=df_band, aes(x=x, y=med), method="loess", span=0.1, color=color, linewidth=width_line, se=FALSE, alpha = alpha_line) +
    geom_density(data=data.frame(ici=obs_ici), aes(x=ici), color=col_dark, linewidth=width_line, linetype="dashed", alpha = alpha_line) +
    #annotate("text", x=x_upper*0.98, y=Inf, label=if(is_correct) "\u2713 correct" else "\u2717 wrong", hjust=1, vjust=1.5, size=5, fontface="bold", color=if(is_correct) color else "grey50") +
    labs(x="ICI (s)", y="Density", title=title) +
    scale_x_continuous(limits=c(0.0, x_upper), expand=c(0,0)) +
    scale_y_continuous(expand=c(0.01,0)) +
    theme_classic() +
    theme(panel.grid.minor=element_blank(),
          plot.title=element_text(size=20, face="bold", hjust = 0.5),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 15, color = "black"),
          axis.text.y = element_blank(),
          axis.title.y = element_text(size = 15, color = "black"),
          axis.title.x = if(show_y) element_text(size=12, color = "black") else element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.border = element_blank(),
    )
}

make_panel_title <- function(label, size = 7) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, hjust = 0.5,
             fontface = "bold", label = label,
             size = 7, color = "black") +
    theme_void()
}

make_model_card <- function(color,
                            title,
                            mechanism,
                            equations,   # list of TeX strings
                            eq_y,        # y positions for equations
                            predictions, # named list of label = value
                            subtitle,
                            size_t = 5) {
  
  p <- ggplot() +
    
    # Background
    annotate("rect", xmin=0, xmax=1, ymin=0, ymax=1,
             fill=col_light, color=color, linewidth=1.5) +
    
    # Title bar
    annotate("rect", xmin=0, xmax=1, ymin=0.88, ymax=1,
             fill=color, color=color, linewidth=0) +
    annotate("text", x=0.5, y=0.94,
             label=title, hjust=0.5, vjust=0.5,
             fontface="bold", size=size_t, color="white") +
    
    # Subtitle — mechanism description
    annotate("text", x=0.5, y=0.84,
             label=subtitle, hjust=0.5, vjust=1,
             fontface="italic", size=size_t, color=col_dark,
             lineheight=1.3)
  
  # Equations
  for (i in seq_along(equations)) {
    p <- p + annotate("text", x=0.5, y=eq_y[i],
                      label=equations[[i]],
                      hjust=0.5, vjust=0.5,
                      size=size_t, color=col_dark,
                      parse=TRUE)
  }
  
  # Divider before predictions
  p <- p +
    annotate("segment", x=0.05, xend=0.95, y=0.38, yend=0.38,
             color="#cccccc", linewidth=0.5) +
    annotate("text", x=0.5, y=0.37,
             label="Diagnostic Predictions",
             hjust=0.5, vjust=1,
             fontface="bold", size=size_t, color=col_dark)
  
  # Prediction rows
  pred_names <- names(predictions)
  pred_vals  <- unlist(predictions)
  n_pred     <- length(pred_names)
  y_starts   <- seq(0.3, 0.04, length.out=n_pred)
  
  for (i in seq_len(n_pred)) {
    p <- p +
      annotate("text", x=0.08, y=y_starts[i],
               label=pred_names[i],
               hjust=0, vjust=0.5,
               size=size_t, color="black", fontface="italic") +
      annotate("text", x=0.92, y=y_starts[i],
               label=pred_vals[i],
               hjust=1, vjust=0.5,
               size=size_t, color=color, fontface="bold")
  }
  
  p +
    scale_x_continuous(limits=c(0,1), expand=c(0,0)) +
    scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
    theme_void() +
    theme(plot.margin=margin(5,8,5,8))
}

# Scatter helper
make_scatter <- function(x_var, y_var, x_label, y_label,
                         subtitle = NULL, show_legend = FALSE) {
  
  p <- ggplot() +
    geom_point(data = ref_sample,
               aes(x = .data[[x_var]], y = .data[[y_var]],
                   color = model),
               alpha = 0.3, size = 0.9) +
    geom_point(data = obs_df,
               aes(x = .data[[x_var]], y = .data[[y_var]]), color = "black",
               size = 5, shape = 18) +
    geom_point(data = obs_df,
               aes(x = .data[[x_var]], y = .data[[y_var]],
                   color = model),
               size = 4.5, shape = 18) +
    scale_color_manual(
      values = c("Event-Triggered" = col_event,
                 "Tempo-Tracking"  = col_tempo,
                 "State-Dependent" = col_state),
      name   = NULL
    ) +
    labs(x = x_label, y = y_label, subtitle = subtitle) +
    theme_classic(base_size = 10) +
    theme(
      plot.title   = element_text(face = "bold", size = 20, hjust = 0.5, color = "black"),
      axis.text.x  = element_text(size = 15, color = "black"),
      axis.text.y  = element_text(size = 15, color = "black"),
      axis.title.y = element_text(face = "bold", size = 15, color = "black"),
      axis.title.x = element_text(face = "bold", size = 15, color = "black", margin = margin(t = 15)),
      panel.grid   = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 16, face = "bold"),  # Bigger title
      legend.text = element_text(size = 15),  # Bigger numbers
    )
  p
}

col_event <- "#D55E00"
col_tempo <- "#0072B2"
col_state <- "#009E73"
col_dark  <- "#2C2C2C"
col_light <- "white"

# Checkpoint helper: create file with headers on first write, append thereafter
write_or_append <- function(row, path) {
  if (!file.exists(path)) {
    write.table(row, path, sep=",", row.names=FALSE, col.names=TRUE)
  } else {
    write.table(row, path, append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
  }
}


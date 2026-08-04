# from code_2.R file 
# Linear Counting Simulation 
# Baseline table comparison of tablet variability 

library(tibble)
library(ggplot2)
set.seed(42)

n_targets        <- 24978  # fossil grains (targets) per slide
n_markers        <- 12489   # marker grains per slide
strip_height     <- 0.05 # height of the counting strip
stop_threshold <- 100 # k: stop once this many markers are in the counted window
V  <- 1 # V in the concentration formula (unit volume)
num_slides       <- 10000 # number of simulations

# Marker dose parameters
marker_label <- n_markers                       # what the formula assumes
tablet_cv    <- 491 / 12489                     # mahers measured CV  0.0393
marker_sd    <- marker_label * tablet_cv        # SD at our scale  

# Storage
estimates_on  <- numeric(num_slides) #concentration estimate per slide with tablet variability
x_counts_on   <- numeric(num_slides) # fossils counted per slide with tablet variability
n_counts_on   <- numeric(num_slides)
actual_doses  <- numeric(num_slides)


# TABLET ERROR ON 
# marker dose varies slide to slide 

cat("TABLET ERROR ON (", num_slides, "slides) \n")
start_t <- Sys.time()

for (s in 1:num_slides) {
  
  # marker dose drawn from the tablet distribution
  n_markers_actual <- round(rnorm(1, mean = marker_label, sd = marker_sd))
  if (n_markers_actual < 1) n_markers_actual <- 1
  actual_doses[s] <- n_markers_actual
  
  targets <- tibble(x = runif(n_targets, 0, 1),
                    y = runif(n_targets, 0, 1),
                    type = "target")
  markers <- tibble(x = runif(n_markers_actual, 0, 1),
                    y = runif(n_markers_actual, 0, 1),
                    type = "marker")
  slide <- rbind(targets, markers)
  # Keep only markers inside the strip height 
  markers_in_strip <- markers[markers$y <= strip_height, ]   #Filter the particles to the strip height
  markers_in_strip <- markers_in_strip[order(markers_in_strip$x), ] # sort the filtered particles 
  # if not enough markers then na 
  if (nrow(markers_in_strip) < stop_threshold) {
    estimates_on[s] <- NA; x_counts_on[s] <- NA; n_counts_on[s] <- NA
    next
  }
  
  x_end <- markers_in_strip$x[stop_threshold]
  in_window  <- slide$y <= strip_height & slide$x <= x_end
  x_count <- sum(slide$type == "target" & in_window)
  n_count <- sum(slide$type == "marker" & in_window)
  
  
  estimates_on[s] <- (x_count * marker_label) / (n_count * V)
  x_counts_on[s]  <- x_count
  n_counts_on[s]  <- n_count
  
  if (s %% 10000 == 0) cat("  slide", s, "of", num_slides, "\n")
}

cat("Run A finished in:", format(Sys.time() - start_t), "\n\n")


# Analysis 
true_concentration <- n_targets / V
# Empirical measured from simulation
mean_on <- mean(estimates_on, na.rm = TRUE) 
sd_on   <- sd(estimates_on, na.rm = TRUE)
cv_on   <- 100 * sd_on / mean_on

mean_x <- mean(x_counts_on, na.rm = TRUE)
mean_n <- stop_threshold
s1P    <- tablet_cv
# Predicted from mays formula 

T_term <- s1P / sqrt(1)                  # tablet variability
F_term <- sqrt(mean_x) / mean_x           # fossil counting noise
M_term <- sqrt(mean_n) / mean_n           # marker counting noise
predicted_cv_on <- 100 * sqrt(T_term^2 + F_term^2 + M_term^2)

rel_error_on <- 100 * (mean_on - true_concentration) / true_concentration
#acc_conc_on  <- 100 * abs(mean_on - true_concentration) / true_concentration
acc_cv_on    <- 100 * abs(predicted_cv_on - cv_on) / cv_on


###################
# TABLET ERROR not included
# marker dose fixed 

estimates_off <- numeric(num_slides)
x_counts_off  <- numeric(num_slides)
n_counts_off  <- numeric(num_slides)

cat(" TABLET ERROR OFF (", num_slides, "slides) \n")
start_t <- Sys.time()

for (s in 1:num_slides) {
  
  n_markers_fixed <- marker_label   
  
  targets <- tibble(x = runif(n_targets, 0, 1),
                    y = runif(n_targets, 0, 1),
                    type = "target")
  markers <- tibble(x = runif(n_markers_fixed, 0, 1),
                    y = runif(n_markers_fixed, 0, 1),
                    type = "marker")
  slide <- rbind(targets, markers)
  
  markers_in_strip <- markers[markers$y <= strip_height, ]
  markers_in_strip <- markers_in_strip[order(markers_in_strip$x), ]
  
  if (nrow(markers_in_strip) < stop_threshold) {
    estimates_off[s] <- NA; x_counts_off[s] <- NA; n_counts_off[s] <- NA
    next
  }
  
  x_end <- markers_in_strip$x[stop_threshold]
  in_window  <- slide$y <= strip_height & slide$x <= x_end
  x_count <- sum(slide$type == "target" & in_window)
  n_count <- sum(slide$type == "marker" & in_window)
  
  estimates_off[s] <- (x_count * marker_label) / (n_count * V)
  x_counts_off[s]  <- x_count
  n_counts_off[s]  <- n_count
  
  if (s %% 10000 == 0) cat("  slide", s, "of", num_slides, "\n")
}

cat("Run B finished in:", format(Sys.time() - start_t), "\n\n")


# Analysis
mean_off <- mean(estimates_off, na.rm = TRUE)
sd_off   <- sd(estimates_off, na.rm = TRUE)
cv_off   <- 100 * sd_off / mean_off

mean_x_off <- mean(x_counts_off, na.rm = TRUE)
predicted_cv_off <- 100 * sqrt(
  (0)^2 +
    (sqrt(mean_x_off) / mean_x_off)^2 +
    (sqrt(stop_threshold) / stop_threshold)^2
)

rel_error_off <- 100 * (mean_off - true_concentration) / true_concentration
#acc_conc_off  <- 100 * abs(mean_off - true_concentration) / true_concentration
acc_cv_off    <- 100 * abs(predicted_cv_off - cv_off) / cv_off



# COMPARISON TABLE


comparison1 <- tibble(
  quantity = c("True concentration",
               "Mean estimate",
               "Signed relative error of mean (%)",
               "Empirical CV (%)",
               "Predicted CV (%)",
               "CV Difference ",
               "Absolute accuracy of CV (%)",
               "Empirical SD"),
  without_tablet = c(true_concentration,
                     round(mean_off, 1),
                     round(rel_error_off, 4),
                     #round(acc_conc_off, 4),
                     round(cv_off, 2),
                     round(predicted_cv_off, 2),
                     round(abs(cv_off - predicted_cv_off), 2),
                     round(acc_cv_off, 2),
                     round(sd_off, 1)),
  with_tablet    = c(true_concentration,
                     round(mean_on, 1),
                     round(rel_error_on, 4),
                     #round(acc_conc_on, 4),
                     round(cv_on, 2),
                     round(predicted_cv_on, 2),
                     round(abs(cv_on - predicted_cv_on), 2),
                     round(acc_cv_on, 2),
                     round(sd_on, 1))
)

print(comparison1)
write.csv(comparison1, "comparison_baseline_1.csv",  row.names=FALSE)


####################################
# THRESHOLD SWEEP WITH TABLET ERROR

# Thresholds to test
thresholds <- c(10,25,50,100,200,400,500,600,700,800)

sweep_off <- tibble()
sweep_on  <- tibble()

# TABLET ERROR not included 
cat("THRESHOLD SWEEP — TABLET ERROR OFF\n")
start_t <- Sys.time()

for (thresh in thresholds) {
  
  cat("  threshold =", thresh, "\n")
  
  estimates <- numeric(num_slides)
  x_counts  <- numeric(num_slides)
  
  for (s in 1:num_slides) {
    targets <- tibble(x = runif(n_targets), y = runif(n_targets), type = "target")
    markers <- tibble(x = runif(marker_label), y = runif(marker_label), type = "marker")
    slide <- rbind(targets, markers)
    
    markers_in_strip <- markers[markers$y <= strip_height, ]
    markers_in_strip <- markers_in_strip[order(markers_in_strip$x), ]
    
    if (nrow(markers_in_strip) < thresh) { estimates[s] <- NA; x_counts[s] <- NA; next }
    
    x_end <- markers_in_strip$x[thresh]
    in_window  <- slide$y <= strip_height & slide$x <= x_end
    x_count <- sum(slide$type == "target" & in_window)
    n_count <- sum(slide$type == "marker" & in_window)
    
    estimates[s] <- (x_count * marker_label) / (n_count * V)
    x_counts[s]  <- x_count
  }
  
  true_c   <- n_targets / V
  emp_mean <- mean(estimates, na.rm = TRUE)
  emp_sd   <- sd(estimates, na.rm = TRUE)
  emp_cv   <- 100 * emp_sd / emp_mean
  
  mean_x <- mean(x_counts, na.rm = TRUE)
  pred_cv <- 100 * sqrt(
    (0)^2 +
      (sqrt(mean_x) / mean_x)^2 +
      (sqrt(thresh) / thresh)^2
  )
  
  rel_err  <- 100 * (emp_mean - true_c) / true_c
  acc_conc <- 100 * abs(emp_mean - true_c) / true_c
  acc_cv   <- 100 * abs(pred_cv - emp_cv) / emp_cv
  
  sweep_off <- rbind(sweep_off, tibble(
    threshold    = thresh,
    mean_x       = round(mean_x, 1),
    empirical_cv = round(emp_cv, 2),
    predicted_cv = round(pred_cv, 2),
    cv_agreement = round(abs(emp_cv - pred_cv), 2),
    acc_cv       = round(acc_cv, 2),
    rel_error    = round(rel_err, 4),
    acc_conc     = round(acc_conc, 4),
    scenario     = "Without tablet error"
  ))
}

cat("Run 1 finished in:", format(Sys.time() - start_t), "\n\n")
write.csv(sweep_off, "sweep_off_10to800.csv", row.names = FALSE)



# TABLET ERROR included 
cat("THRESHOLD SWEEP — TABLET ERROR ON \n")
start_t <- Sys.time()

for (thresh in thresholds) {
  
  cat("  threshold =", thresh, "\n")
  
  estimates <- numeric(num_slides)
  x_counts  <- numeric(num_slides)
  
  for (s in 1:num_slides) {
    
    n_markers_actual <- round(rnorm(1, mean = marker_label, sd = marker_sd))
    if (n_markers_actual < 1) n_markers_actual <- 1
    
    targets <- tibble(x = runif(n_targets), y = runif(n_targets), type = "target")
    markers <- tibble(x = runif(n_markers_actual), y = runif(n_markers_actual), type = "marker")
    slide <- rbind(targets, markers)
    
    markers_in_strip <- markers[markers$y <= strip_height, ]
    markers_in_strip <- markers_in_strip[order(markers_in_strip$x), ]
    
    if (nrow(markers_in_strip) < thresh) { estimates[s] <- NA; x_counts[s] <- NA; next }
    
    x_end <- markers_in_strip$x[thresh]
    in_window  <- slide$y <= strip_height & slide$x <= x_end
    x_count <- sum(slide$type == "target" & in_window)
    n_count <- sum(slide$type == "marker" & in_window)
    
    estimates[s] <- (x_count * marker_label) / (n_count * V)
    x_counts[s]  <- x_count
  }
  
  true_c   <- n_targets / V
  emp_mean <- mean(estimates, na.rm = TRUE)
  emp_sd   <- sd(estimates, na.rm = TRUE)
  emp_cv   <- 100 * emp_sd / emp_mean
  
  mean_x <- mean(x_counts, na.rm = TRUE)
  pred_cv <- 100 * sqrt(
    (tablet_cv / sqrt(1))^2 +
      (sqrt(mean_x) / mean_x)^2 +
      (sqrt(thresh) / thresh)^2
  )
  
  rel_err  <- 100 * (emp_mean - true_c) / true_c
  acc_conc <- 100 * abs(emp_mean - true_c) / true_c
  acc_cv   <- 100 * abs(pred_cv - emp_cv) / emp_cv
  
  sweep_on <- rbind(sweep_on, tibble(
    threshold    = thresh,
    mean_x       = round(mean_x, 1),
    empirical_cv = round(emp_cv, 2),
    predicted_cv = round(pred_cv, 2),
    cv_agreement = round(abs(emp_cv - pred_cv), 2),
    acc_cv       = round(acc_cv, 2),
    rel_error    = round(rel_err, 4),
    acc_conc     = round(acc_conc, 4),
    scenario     = "With tablet error"
  ))
}

write.csv(sweep_on, "sweep_on_10to800.csv", row.names = FALSE)
cat("Run 2 finished in:", format(Sys.time() - start_t), "\n\n")



# COMPARISON TABLE
comparison2 <- tibble(
  threshold     = sweep_off$threshold,
  cv_off        = sweep_off$empirical_cv,
  cv_on         = sweep_on$empirical_cv,
  cv_difference = round(sweep_on$empirical_cv - sweep_off$empirical_cv, 2),
  pred_off      = sweep_off$predicted_cv,
  pred_on       = sweep_on$predicted_cv,
  acc_cv_off    = sweep_off$acc_cv,
  acc_cv_on     = sweep_on$acc_cv
)

print(comparison2)
write.csv(comparison2, "comparison_threshold_10to800.csv",  row.names=FALSE)


# PRECISION-VS-stopping threshold PLOT
both <- rbind(sweep_off, sweep_on)

ggplot(both, aes(x = threshold, colour = scenario)) +
  geom_line(aes(y = empirical_cv, linetype = "Empirical"), linewidth = 1.2) +
  geom_point(aes(y = empirical_cv), size = 3) +
  geom_line(aes(y = predicted_cv, linetype = "Predicted"), linewidth = 0.8) +
  geom_point(aes(y = predicted_cv), size = 2, shape = 1) +
  scale_colour_manual(values = c("steelblue", "orange")) +
  scale_linetype_manual(values = c(Empirical = "solid", Predicted = "dashed")) +
  scale_x_log10(breaks = thresholds) +
  labs(title    = "Precision vs counting effort, with and without tablet variability",
       subtitle = paste0("Solid = empirical, Dashed = Mays formula prediction. ",
                         num_slides, " slides per point."),
       x = "Markers counted (stopping threshold, log scale)",
       y = "Total error (%)",
       colour = NULL, linetype = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")



# VARIANCE BREAKDOWN by error source Bar Chart 
# Shows what fraction of the total variance each formula term contributes
# at each threshold: tablet (T), fossil (F), marker (M).

breakdown <- tibble()

for (i in 1:nrow(sweep_on)) {
  thresh <- sweep_on$threshold[i]
  mean_x_here <- sweep_on$mean_x[i]
  
  T_sq <- (tablet_cv / sqrt(1))^2
  F_sq <- 1 / mean_x_here
  M_sq <- 1 / thresh
  total <- T_sq + F_sq + M_sq
  
  breakdown <- rbind(breakdown, tibble(
    threshold = rep(thresh, 3),
    source    = c("T (tablet)", "F (fossil)", "M (marker)"),
    pct       = 100 * c(T_sq, F_sq, M_sq) / total
  ))
}

breakdown$source <- factor(breakdown$source,
                           levels = c("M (marker)", "F (fossil)", "T (tablet)"))

ggplot(breakdown, aes(x = factor(threshold), y = pct, fill = source)) +
  geom_col(colour = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("orange", "steelblue", "red")) +
  labs(title    = "Variance contribution from each error source",
       subtitle = "At low thresholds, marker noise dominates. Tablet error grows with threshold.",
       x = "Stopping threshold",
       y = "Contribution to total variance (%)",
       fill = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

###################################
# TARGET-TO-MARKER RATIO SWEEP
# Ratios to test
ratios <- c(1, 2, 3, 6, 10, 30, 60)


ratio_off <- tibble()
ratio_on  <- tibble()

# TABLET ERROR Not included 
cat("RATIO SWEEP — TABLET ERROR OFF\n")
start_t <- Sys.time()

for (r in ratios) {
  
  n_targets_r <- r * n_markers
  true_c    <- n_targets_r / V
  
  cat("  ratio =", r, " (n_targets =", n_targets_r, ")\n")
  
  estimates <- numeric(num_slides)
  x_counts  <- numeric(num_slides)
  
  for (s in 1:num_slides) {
    targets <- tibble(x = runif(n_targets_r), y = runif(n_targets_r), type = "target")
    markers <- tibble(x = runif(n_markers), y = runif(n_markers), type = "marker")
    slide <- rbind(targets, markers)
    
    markers_in_strip <- markers[markers$y <= strip_height, ]
    markers_in_strip <- markers_in_strip[order(markers_in_strip$x), ]
    
    if (nrow(markers_in_strip) < stop_threshold) { 
      estimates[s] <- NA; x_counts[s] <- NA; next 
    }
    
    x_end <- markers_in_strip$x[stop_threshold]
    in_window  <- slide$y <= strip_height & slide$x <= x_end
    x_count <- sum(slide$type == "target" & in_window)
    n_count <- sum(slide$type == "marker" & in_window)
    
    estimates[s] <- (x_count * marker_label) / (n_count * V)
    x_counts[s]  <- x_count
  }
  
  emp_mean <- mean(estimates, na.rm = TRUE)
  emp_sd   <- sd(estimates, na.rm = TRUE)
  emp_cv   <- 100 * emp_sd / emp_mean
  
  mean_x <- mean(x_counts, na.rm = TRUE)
  pred_cv <- 100 * sqrt(
    (0)^2 +
      (sqrt(mean_x) / mean_x)^2 +
      (sqrt(stop_threshold) / stop_threshold)^2
  )
  
  rel_error       <- 100 * (emp_mean - true_c) / true_c
  cv_signed_error <- 100 * (pred_cv - emp_cv) / emp_cv
  
  ratio_off <- rbind(ratio_off, tibble(
    ratio           = r,
    n_targets       = n_targets_r,
    mean_x          = round(mean_x, 1),
    empirical_cv    = round(emp_cv, 2),
    predicted_cv    = round(pred_cv, 2),
    cv_signed_error = round(cv_signed_error, 2),
    rel_error       = round(rel_error, 4),
    failed_slides   = sum(is.na(estimates)),
    scenario        = "Without tablet error"
  ))
}

cat("Run 1 finished in:", format(Sys.time() - start_t), "\n\n")
write.csv(ratio_off, "ratio_off_2.csv",  row.names=FALSE)


# TABLET ERROR Included 

cat(" RATIO SWEEP — TABLET ERROR ON\n")
start_t <- Sys.time()

for (r in ratios) {
  
  n_targets_r <- r * n_markers
  true_c    <- n_targets / V
  
  cat("  ratio =", r, " (n_targets =", n_targets_r, ")\n")
  
  estimates <- numeric(num_slides)
  x_counts  <- numeric(num_slides)
  
  for (s in 1:num_slides) {
    
    n_markers_actual <- round(rnorm(1, mean = marker_label, sd = marker_sd))
    if (n_markers_actual < 1) n_markers_actual <- 1
    
    targets <- tibble(x = runif(n_targets_r), y = runif(n_targets_r), type = "target")
    markers <- tibble(x = runif(n_markers_actual), y = runif(n_markers_actual), type = "marker")
    slide <- rbind(targets, markers)
    
    markers_in_strip <- markers[markers$y <= strip_height, ]
    markers_in_strip <- markers_in_strip[order(markers_in_strip$x), ]
    
    if (nrow(markers_in_strip) < stop_threshold) { 
      estimates[s] <- NA; x_counts[s] <- NA; next 
    }
    
    x_end <- markers_in_strip$x[stop_threshold]
    in_window  <- slide$y <= strip_height & slide$x <= x_end
    x_count <- sum(slide$type == "target" & in_window)
    n_count <- sum(slide$type == "marker" & in_window)
    
    estimates[s] <- (x_count * marker_label) / (n_count * V)
    x_counts[s]  <- x_count
  }
  
  emp_mean <- mean(estimates, na.rm = TRUE)
  emp_sd   <- sd(estimates, na.rm = TRUE)
  emp_cv   <- 100 * emp_sd / emp_mean
  
  mean_x <- mean(x_counts, na.rm = TRUE)
  pred_cv <- 100 * sqrt(
    (tablet_cv / sqrt(1))^2 +
      (sqrt(mean_x) / mean_x)^2 +
      (sqrt(stop_threshold) / stop_threshold)^2
  )
  
  rel_error       <- 100 * (emp_mean - true_c) / true_c
  cv_signed_error <- 100 * (pred_cv - emp_cv) / emp_cv
  
  ratio_on <- rbind(ratio_on, tibble(
    ratio           = r,
    n_targets       = n_targets_r,
    mean_x          = round(mean_x, 1),
    empirical_cv    = round(emp_cv, 2),
    predicted_cv    = round(pred_cv, 2),
    cv_signed_error = round(cv_signed_error, 2),
    rel_error       = round(rel_error, 4),
    failed_slides   = sum(is.na(estimates)),
    scenario        = "With tablet error"
  ))
}

cat("Run 2 finished in:", format(Sys.time() - start_t), "\n\n")

write.csv(ratio_on, "ratio_on_2.csv",  row.names=FALSE)

# COMPARISON TABLE
ratio_comparison <- tibble(
  ratio            = ratio_off$ratio,
  cv_off           = ratio_off$empirical_cv,
  cv_on            = ratio_on$empirical_cv,
  cv_tablet_effect = round(ratio_on$empirical_cv - ratio_off$empirical_cv, 2),
  pred_off         = ratio_off$predicted_cv,
  pred_on          = ratio_on$predicted_cv,
  cv_error_off     = ratio_off$cv_signed_error,
  cv_error_on      = ratio_on$cv_signed_error,
  mean_x_off       = ratio_off$mean_x
)

print(ratio_comparison)

write.csv(ratio_comparison, "comparison_ratio_2.csv",  row.names=FALSE)

# PRECISION-VS-TARGET-TO-MARKER RATIO PLOT
both_ratio <- rbind(ratio_off, ratio_on)

ggplot(both_ratio, aes(x = ratio, colour = scenario)) +
  geom_line(aes(y = empirical_cv, linetype = "Empirical"), linewidth = 1.2) +
  geom_point(aes(y = empirical_cv), size = 3) +
  geom_line(aes(y = predicted_cv, linetype = "Predicted"), linewidth = 0.8) +
  geom_point(aes(y = predicted_cv), size = 2, shape = 1) +
  scale_colour_manual(values = c("steelblue", "orange")) +
  scale_linetype_manual(values = c(Empirical = "solid", Predicted = "dashed")) +
  scale_x_log10(breaks = ratios) +
  labs(title    = "Precision vs target-to-marker ratio, with and without tablet variability",
       subtitle = paste0("Solid = empirical, Dashed = Mays formula prediction. ",
                         num_slides, " slides per point."),
       x = "Target-to-marker ratio (log scale)",
       y = "Total error (%)",
       colour = NULL, linetype = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")
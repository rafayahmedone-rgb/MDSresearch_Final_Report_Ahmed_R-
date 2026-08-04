
# threshold × ratio sweep
# Both with and without tablet error


library(tibble)
library(ggplot2)


n_markers        <- 12489          
strip_height     <- 0.05
V                <- 1
num_slides       <- 10000

marker_label <- n_markers
tablet_cv    <- 491 / 12489
marker_sd    <- marker_label * tablet_cv

# Grid axes 
ratios     <- c(1, 2, 3, 6, 10, 30, 60) # target-to-marker ratios
thresholds <- c(10, 25, 50, 100, 200, 400) # stopping thresholds k

grid_off <- tibble()
grid_on  <- tibble()

# TABLET ERROR not included 
cat("SCENARIO 1 of 2: TABLET ERROR OFF\n")
cat("===", length(ratios) * length(thresholds), 
    "cells x", num_slides, "slides each \n")


scenario_start <- Sys.time()

for (r in ratios) {
  
  n_targets <- r * n_markers
  true_c    <- n_targets / V
  ratio_start <- Sys.time()
  
  cat(paste0("\n RATIO ", r, " (n_targets=", n_targets, ") \n"))
  
  for (thresh in thresholds) {
    
    cell_start <- Sys.time()
    cat(paste0("  k=", thresh, " "))
    
    set.seed(42)
    estimates <- numeric(num_slides)
    x_counts  <- numeric(num_slides)
    n_failed  <- 0
    
    for (s in 1:num_slides) {
      
      targets <- tibble(x = runif(n_targets), y = runif(n_targets), type = "target")
      markers <- tibble(x = runif(n_markers), y = runif(n_markers), type = "marker")
      slide   <- rbind(targets, markers)
      
      ms <- markers[markers$y <= strip_height, ]
      ms <- ms[order(ms$x), ]
      
      if (nrow(ms) < thresh) {
        estimates[s] <- NA; x_counts[s] <- NA
        n_failed <- n_failed + 1
        next
      }
      
      x_end   <- ms$x[thresh]
      in_w    <- slide$y <= strip_height & slide$x <= x_end
      x_count <- sum(slide$type == "target" & in_w)
      n_count <- sum(slide$type == "marker" & in_w)
      
      estimates[s] <- (x_count * marker_label) / (n_count * V)
      x_counts[s]  <- x_count
    }
    
    # Analysis
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
    
    grid_off <- rbind(grid_off, tibble(
      ratio          = r,
      threshold      = thresh,
      n_targets      = n_targets,
      mean_x         = round(mean_x, 1),
      empirical_cv   = round(emp_cv, 2),
      predicted_cv   = round(pred_cv, 2),
      cv_agreement   = round(abs(emp_cv - pred_cv), 2),
      acc_cv         = round(100 * abs(pred_cv - emp_cv) / emp_cv, 2),
      rel_error      = round(rel_err, 4),
      failed_slides  = n_failed,
      cell_time_sec  = round(as.numeric(Sys.time() - cell_start), 1)
    ))
    
    fail_warn <- if (n_failed > num_slides * 0.01) "  MANY FAILED" else ""
    cat(paste0(": CV=", round(emp_cv, 2), "%, ",
               n_failed, " failed, ",
               round(as.numeric(Sys.time() - cell_start), 1), "s",
               fail_warn, "\n"))
  }
  
  # Save intermediate after each ratio
  saveRDS(grid_off, "grid_off_intermediate.rds")
  cat(paste0("  >> ratio ", r, " done in ",
             format(Sys.time() - ratio_start), ". Saved intermediate.\n"))
}

cat(paste0("\n  SCENARIO 1 (OFF) finished in ",
           format(Sys.time() - scenario_start), " \n"))


##################
# TABLET ERROR included 



cat(" SCENARIO 2 of 2: TABLET ERROR ON  \n")
scenario_start <- Sys.time()

for (r in ratios) {
  
  n_targets <- r * n_markers
  true_c    <- n_targets / V
  ratio_start <- Sys.time()
  
  cat(paste0("\n--- RATIO ", r, " (n_targets=", n_targets, ") ---\n"))
  
  for (thresh in thresholds) {
    
    cell_start <- Sys.time()
    cat(paste0("  k=", thresh, " "))
    
    set.seed(42)
    estimates <- numeric(num_slides)
    x_counts  <- numeric(num_slides)
    n_failed  <- 0
    
    for (s in 1:num_slides) {
      
      # KEY DIFFERENCE: random marker dose per slide
      n_markers_actual <- round(rnorm(1, mean = marker_label, sd = marker_sd))
      if (n_markers_actual < 1) n_markers_actual <- 1
      
      targets <- tibble(x = runif(n_targets), y = runif(n_targets), type = "target")
      markers <- tibble(x = runif(n_markers_actual), y = runif(n_markers_actual), type = "marker")
      slide   <- rbind(targets, markers)
      
      ms <- markers[markers$y <= strip_height, ]
      ms <- ms[order(ms$x), ]
      
      if (nrow(ms) < thresh) {
        estimates[s] <- NA; x_counts[s] <- NA
        n_failed <- n_failed + 1
        next
      }
      
      x_end   <- ms$x[thresh]
      in_w    <- slide$y <= strip_height & slide$x <= x_end
      x_count <- sum(slide$type == "target" & in_w)
      n_count <- sum(slide$type == "marker" & in_w)
      
      # Formula uses LABEL value, not actual dose
      estimates[s] <- (x_count * marker_label) / (n_count * V)
      x_counts[s]  <- x_count
    }
    
    # Analysis
    emp_mean <- mean(estimates, na.rm = TRUE)
    emp_sd   <- sd(estimates, na.rm = TRUE)
    emp_cv   <- 100 * emp_sd / emp_mean
    
    mean_x <- mean(x_counts, na.rm = TRUE)
    pred_cv <- 100 * sqrt(
      (tablet_cv / sqrt(1))^2 +     # T term now non-zero
        (sqrt(mean_x) / mean_x)^2 +
        (sqrt(thresh) / thresh)^2
    )
    
    rel_err  <- 100 * (emp_mean - true_c) / true_c
    
    grid_on <- rbind(grid_on, tibble(
      ratio          = r,
      threshold      = thresh,
      n_targets      = n_targets,
      mean_x         = round(mean_x, 1),
      empirical_cv   = round(emp_cv, 2),
      predicted_cv   = round(pred_cv, 2),
      cv_agreement   = round(abs(emp_cv - pred_cv), 2),
      acc_cv         = round(100 * abs(pred_cv - emp_cv) / emp_cv, 2),
      rel_error      = round(rel_err, 4),
      failed_slides  = n_failed,
      cell_time_sec  = round(as.numeric(Sys.time() - cell_start), 1)
    ))
    
    fail_warn <- if (n_failed > num_slides * 0.01) "  MANY FAILED" else ""
    cat(paste0(": CV=", round(emp_cv, 2), "%, ",
               n_failed, " failed, ",
               round(as.numeric(Sys.time() - cell_start), 1), "s",
               fail_warn, "\n"))
  }
  
  # Save intermediate after each ratio
  saveRDS(grid_on, "grid_on_intermediate.rds")
  cat(paste0("  >> ratio ", r, " done in ",
             format(Sys.time() - ratio_start), ". Saved intermediate.\n"))
}

cat(paste0("\n SCENARIO 2 (ON) finished in ",
           format(Sys.time() - scenario_start), " \n"))



# FINAL OUTPUT


cat(" FINAL RESULTS  \n")
# Save final results
saveRDS(grid_off, "grid_off_final_2.rds")
saveRDS(grid_on,  "grid_on_final_2.rds")
write.csv(grid_off, "grid_off_final_2.csv", row.names = FALSE)
write.csv(grid_on,  "grid_on_final_2.csv",  row.names = FALSE)

# Summary table — empirical CV across the grid (without tablet)
cat("EMPIRICAL CV (%) — without tablet error:\n")
cv_matrix_off <- matrix(grid_off$empirical_cv,
                        nrow = length(ratios),
                        ncol = length(thresholds),
                        byrow = TRUE)
rownames(cv_matrix_off) <- paste0("r=", ratios)
colnames(cv_matrix_off) <- paste0("k=", thresholds)
print(cv_matrix_off)

cat("\nEMPIRICAL CV (%) — with tablet error:\n")
cv_matrix_on <- matrix(grid_on$empirical_cv,
                       nrow = length(ratios),
                       ncol = length(thresholds),
                       byrow = TRUE)
rownames(cv_matrix_on) <- paste0("r=", ratios)
colnames(cv_matrix_on) <- paste0("k=", thresholds)
print(cv_matrix_on)

cat("\nDIFFERENCE (cv_on - cv_off):\n")
print(round(cv_matrix_on - cv_matrix_off, 2))

cat("\nFiles saved:\n")
cat("  grid_off_final_2.rds / .csv\n")
cat("  grid_on_final_2.rds  / .csv\n")
cat("  Plus intermediate saves: grid_off_intermediate.rds, grid_on_intermediate.rds\n")


##############
#Practical Guidance finding best results  

library(tibble)
library(dplyr)

# Load 
on <- read.csv("grid_on_final.csv")

# Add the total particles column
on$total <- on$threshold + on$mean_x

# For each target precision, find the minimum-effort cell
targets <- c(50, 30, 20, 15, 10, 8)

for (tgt in targets) {
  qualifying <- on[on$empirical_cv <= tgt, ]
  
  if (nrow(qualifying) == 0) {
    cat(paste0(tgt, "% target: NO cell achieves this\n"))
    
    next
  }
  
  best <- qualifying[which.min(qualifying$total), ]
  cat(paste0(tgt, "% target: ratio=", best$ratio,
             ", k=", best$threshold,
             ", x.bar=", round(best$mean_x),
             ", total=", round(best$total),
             ", CV=", best$empirical_cv, "%\n"))
}




##############################
# Heatmap of empirical CV
# with tablet variability

on <- read.csv("grid_on_final.csv")

# Add a "rounded label" column for nice annotations
on$cv_label <- as.character(round(on$empirical_cv, 1))

# Plot
ggplot(on, aes(x = factor(threshold), y = factor(ratio), fill = empirical_cv)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = cv_label,
                colour = ifelse(empirical_cv > 35, "white", "black")),
            size = 4) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1,
                       name = "Total error\n(% CV)") +
  scale_colour_identity() +     
  labs(
    title    = "Empirical total error across the 6×7 grid",
    subtitle = "10,000 simulated slides per cell, with tablet variability)",
    x = "Stopping threshold k",
    y = "Target-to-marker ratio"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid       = element_blank(),
    legend.position  = "right",
    plot.title       = element_text(face = "bold")
  )


ggsave("fig8_grid_heatmap.png", width = 9, height = 5.5, dpi = 150,
       bg = "white")





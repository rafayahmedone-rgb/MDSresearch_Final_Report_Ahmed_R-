library(tibble)
library(ggplot2)
n_markers    <- 12489
ratio        <- 1
n_targets    <- ratio * n_markers     
threshold    <- 5                     
strip_height <- 0.05
marker_label <- n_markers             
V            <- 1
true_c       <- n_targets / V         
num_slides   <- 20000                 

set.seed(42)
estimates <- numeric(num_slides)      
x_counts  <- numeric(num_slides)

for (s in 1:num_slides) {
  
  targets <- tibble(x = runif(n_targets), y = runif(n_targets), type = "target")
  markers <- tibble(x = runif(n_markers), y = runif(n_markers), type = "marker")
  slide   <- rbind(targets, markers)
  
  markers_in_strip <- markers[markers$y <= strip_height, ]
  markers_in_strip <- markers_in_strip[order(markers_in_strip$x), ]
  
  if (nrow(markers_in_strip) < threshold) {
    estimates[s] <- NA
    x_counts[s]  <- NA
    next
  }
  
  x_end     <- markers_in_strip$x[threshold]
  in_window <- slide$y <= strip_height & slide$x <= x_end
  x_count   <- sum(slide$type == "target" & in_window)
  n_count   <- sum(slide$type == "marker" & in_window)
  
  estimates[s] <- (x_count * marker_label) / (n_count * V)
  x_counts[s]  <- x_count
}


summary(estimates)      
mean(estimates, na.rm = TRUE)   
sd(estimates, na.rm = TRUE)


sim_results <- tibble(iteration = 1:num_slides, estimate = estimates)
write.csv(sim_results, "sim_estimates_2.csv", row.names = FALSE)
cat("Saved sim_estimates_2.csv\n")

# Sim data 
sim_results <- read.csv("./sim_estimates_2.csv")
est <- sim_results$estimate

est
##
# RMSE 
running_rmse <- numeric(num_slides)

for (n in 1:num_slides) {
  errors          <- estimates[1:n] - true_c     
  running_rmse[n] <- sqrt(mean(errors^2))        # RMSE
}

running_pct <- 100 * running_rmse / true_c       # as % of true c

# Plot
conv <- tibble(iteration = 1:num_slides, rmse_pct = running_pct)

ggplot(conv, aes(x = iteration, y = rmse_pct)) +
  geom_line()
##






sq_err       <- (est - true_c)^2
running_rmse <- sqrt(cumsum(sq_err) / seq_along(sq_err))
running_pct  <- 100 * running_rmse / true_c




sq_err      <- (estimates - true_c)^2
running_pct <- 100 * sqrt(cumsum(sq_err) / seq_along(sq_err)) / true_c
conv        <- tibble(iteration = 1:num_slides, rmse_pct = running_pct)

ggplot(conv, aes(x = iteration  , y =rmse_pct  )) +
  geom_line()




sim_est <- read.csv("sim_estimates_2.csv")
true_c_conv <- 12489          # ratio 1 x 5000 markers
est <- sim_est$estimate
sq_err      <- (est - true_c_conv)^2
running_pct <- 100 * sqrt(cumsum(sq_err) / seq_along(sq_err)) / true_c_conv
conv <- tibble(iteration = seq_along(est), rmse_pct = running_pct)

ggplot(conv, aes(x = iteration, y = rmse_pct)) +
  geom_line(colour = col_off, linewidth = 0.4) +
  geom_vline(xintercept = 10000, linetype = "dotted", colour = "grey40") +
  labs(x = "Number of simulated slides", y = "Running RMSE (% of true value)") +
  theme_minimal(base_size = 11)



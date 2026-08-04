# From code.R file 

library(tibble)
library(ggplot2)
set.seed(42)
# manual inter evrnt distance 
interevent_dist <- function(points) {
  n <- nrow(points)  
  
  # number of unique pairs is n*(n-1)/2
  num_pairs <- n * (n - 1) / 2              #for eg if 4 points then 6 distances  12 13 14 23 24 34 cant go to same point and cant repeat 
  
  # preallocate the output vector (faster than growing with c())
  d <- numeric(num_pairs) # distances vector of length n_rows
  k <- 1  # counter for which pair we're on
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      dx <- points$x[i] - points$x[j]
      dy <- points$y[i] - points$y[j]
      d[k] <- sqrt(dx^2 + dy^2)
      k <- k + 1
    }
  }
  
  d
}

interevent_dist(points)
dist(points)


# Diggle

# Diggle's theoretical CDF (Eqn 2.2 for unit square)
diggle_cdf <- function(t) {
  H <- numeric(length(t))
  
  # First piece: 0 <= t <= 1
  i1 <- t >= 0 & t <= 1  # gives like a logical vector of true false 
  H[i1] <- pi * t[i1]^2 - 8 * t[i1]^3 / 3 + t[i1]^4 / 2 # takes only the values which is true 
  
  # Second piece: 1 < t <= sqrt(2)
  i2 <- t > 1 & t <= sqrt(2)
  tt <- t[i2]  # added for convinence 
  arg <- pmax(pmin(2 / tt^2 - 1, 1), -1)   # round off  for asin so values remain between -1 1 CAN TRY TO REMOVE ....
  H[i2] <- 1/3 - 2 * tt^2 - tt^4 / 2 +
    4 * sqrt(tt^2 - 1) * (2 * tt^2 + 1) / 3 +
    2 * tt^2 * asin(arg)
  
  # Beyond sqrt(2), CDF is 1
  H[t > sqrt(2)] <- 1
  
  H
}



# USING geom line 
# build a tibble of t-values and the corresponding H(t) values
diggle_curve <- tibble(
  t = seq(0, sqrt(2), length.out = 200),
  H = diggle_cdf(seq(0, sqrt(2), length.out = 200))
)

# then plot
ggplot(diggle_curve, aes(x = t, y = H)) +
  geom_line(colour = "red", linewidth = 1) +
  labs(title = "Diggle's theoretical CDF for inter-event distances",
       x = "Distance t",
       y = "H(t)") +
  theme_minimal()




# Pool all the distances into one big vector
pooled_distances <- unlist(all_distances)
length(pooled_distances)   # should be 50 * 100*99/2 = 247,500


unlist(all_distances)
all_distances

#####################################################
# Generate many simulated Poisson slides


set.seed(42)
n_points        <- 100   # points per slide
num_simulations <- 50    # number of slides

all_distances <- vector("list", num_simulations)   # Empty list with num_simulation slots 1-50 

for (i in 1:num_simulations) {
  pts <- tibble(x = runif(n_points, 0, 1),     # points tiblle containg x and y 
                y = runif(n_points, 0, 1))
  all_distances[[i]] <- as.vector(dist(pts))    # Note [[i]] are how we access the individual element otherwise [] will give a sublist 
}

# Pool all the distances into one big vector
pooled_distances <- unlist(all_distances)
length(pooled_distances)   # should be 50 * 100*99/2 = 247,500



# Build empirical CDF and pair it with theoretical

# tibble contaiing x distances y empericala nd theoretical values 
comparison <- tibble( 
  distance    = sort(pooled_distances),
  empirical   = (1:length(pooled_distances)) / length(pooled_distances)
)
comparison$theoretical <- diggle_cdf(comparison$distance)



# OR PLOTTING DIFFERENT IN SAME WINDOW 

library(patchwork)

p1 <- ggplot(comparison, aes(x = distance, y = empirical)) +
  geom_step(colour = "blue") +
  labs(title = "Empirical inter-event CDF",
       x = "Distance t", y = "H(t)") +
  theme_minimal()

p2 <- ggplot(comparison, aes(x = distance, y = theoretical)) +
  geom_line(colour = "red", linewidth = 1) +
  labs(title = "Theoretical inter-event CDF ",
       x = "Distance t", y = "H(t)") +
  theme_minimal()

# side by side
p1 + p2

# or stacked vertically
p1 / p2






# NEAREST-NEIGHBOUR DISTANCE VERIFICATION
# compare the theoretical diggle with the simulation nearest neighbour


library(tibble)
library(ggplot2)
library(spatstat)


#  Manual nearest-neighbour distance
#min distance for each point 
nn_dist_by_hand <- function(points) {
  n <- nrow(points)
  nn <- numeric(n)   
  
  for (i in 1:n) {
    min_d <- Inf   # start with infinity so any real distance beats it
    
    for (j in 1:n) {
      if (i != j) {   # skip comparing a point to itself
        dx <- points$x[i] - points$x[j]
        dy <- points$y[i] - points$y[j]
        d  <- sqrt(dx^2 + dy^2)
        
        if (d < min_d) min_d <- d   # keep smallest seen so far
      }
    }
    
    nn[i] <- min_d
  }
  
  nn
}


# Verify 
set.seed(1)
check_pts <- tibble(x = runif(20), y = runif(20))

manual_result   <- nn_dist_by_hand(check_pts)
spatstat_result <- nndist(check_pts$x, check_pts$y)

all(abs(manual_result - spatstat_result) < 1e-10)   # should be TRUE



# Run many simulations and pool NN distances

set.seed(42)
n_points        <- 100
num_simulations <- 50

all_nn <- vector("list", num_simulations)

for (i in 1:num_simulations) {
  pts <- tibble(x = runif(n_points, 0, 1),
                y = runif(n_points, 0, 1))
  all_nn[[i]] <- nndist(pts$x, pts$y)
}

pooled_nn <- unlist(all_nn)
length(pooled_nn)   



# Diggle's Nearest neighbour 

diggle_nn_cdf <- function(y, n) {
  1 - exp(-n * pi * y^2)
}


diggle_nn_cdf(0, n_points)       # should be 0
diggle_nn_cdf(0.2, n_points)     # should be close to 1 (for n = 100)



# Build empirical CDF and compare


nn_comparison <- tibble(
  distance    = sort(pooled_nn),
  empirical   = (1:length(pooled_nn)) / length(pooled_nn)
)
nn_comparison$theoretical <- diggle_nn_cdf(nn_comparison$distance, n_points)



# Plot empirical vs theoretical CDF


ggplot(nn_comparison, aes(x = distance)) +
  geom_step(aes(y = empirical),   colour = "blue") +
  geom_line(aes(y = theoretical), colour = "red", linewidth = 1) +
  labs(title    = "Empirical vs theoretical CDF of nearest-neighbour distances",
       subtitle = "Blue = simulation, Red = Diggle Eqn 2.8",
       x        = "Nearest-neighbour distance y",
       y        = "G(y)") +
  theme_minimal()



# Empirical histogram of NN distances 
# will use this 


ggplot(tibble(nn = pooled_nn), aes(x = nn)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 40, fill = "lightblue", colour = "black") +
  labs(title = "Empirical nearest-neighbour distance distribution",
       x     = "Nearest-neighbour distance y",
       y     = "Density") +
  theme_minimal()  

########################################

#add the markers




# Generate one slide with both targets and markers

n_targets <- 200   #  targets 
n_markers <- 100    #  markers

targets <- tibble(x = runif(n_targets, 0, 1),
                  y = runif(n_targets, 0, 1),
                  type = "target")

markers <- tibble(x = runif(n_markers, 0, 1),
                  y = runif(n_markers, 0, 1),
                  type = "marker")

slide <- rbind(targets, markers)

# PLOT SLIDE  
ggplot(slide, aes(x, y, colour = type, shape = type, size = type)) +
  geom_point() +
  scale_colour_manual(values = c(target = "grey50", marker = "orange")) +
  scale_shape_manual(values = c(target = 16, marker = 18)) +   # shape 
  scale_size_manual(values = c(target = 1.5, marker = 3)) +
  coord_equal() +
  theme_minimal() +
  labs(title = "A simulated slide with targets and markers",
       x = "x", y = "y")



################################
# LINEAR COUNTING PLOT 

# Setup (same slide as before)

n_targets <- 2000
n_markers <- 1000

targets <- tibble(x = runif(n_targets, 0, 1),
                  y = runif(n_targets, 0, 1),
                  type = "target")
markers <- tibble(x = runif(n_markers, 0, 1),
                  y = runif(n_markers, 0, 1),
                  type = "marker")
slide <- rbind(targets, markers)

strip_height <- 0.05

# Find the x-position at each milestone 
markers_in_strip <- markers[markers$y <= strip_height, ]
markers_in_strip <- markers_in_strip[order(markers_in_strip$x), ]

# Helper function: plot the slide with a transect reaching up to a given marker
plot_at_marker <- function(k) {
  x_end <- markers_in_strip$x[k]
  
  # Count what's inside the transect so far
  in_transect <- slide$y <= strip_height & slide$x <= x_end
  x_count <- sum(slide$type == "target" & in_transect)
  n_count <- sum(slide$type == "marker" & in_transect)
  
  
  ggplot(slide, aes(x, y, colour = type, shape = type, size = type)) +
    geom_point() +
    scale_colour_manual(values = c(target = "grey60", marker = "orange"),
                        guide = "none") +
    scale_shape_manual(values = c(target = 16, marker = 18),
                       guide = "none") +
    scale_size_manual(values = c(target = 1.2, marker = 2.5),
                      guide = "none") +
    annotate("rect",                                       # add a rectangle immitating a transect window 
             xmin = 0, xmax = x_end,
             ymin = 0, ymax = strip_height,
             colour = "red", fill = "red", alpha = 0.15, linewidth = 0.8) +
    coord_cartesian(xlim = c(0, 0.3), ylim = c(0, 0.15)) +
    theme_minimal(base_size = 9) +
    labs(title    = paste0("After ", k, " markers"),
         subtitle = paste0("x = ", round(x_end, 3),
                           ", targets = ", x_count))
}

# Build four snapshots and combine them
p1 <- plot_at_marker(3)
p2 <- plot_at_marker(6)
p3 <- plot_at_marker(9)
p4 <- plot_at_marker(10)

(p1 + p2) / (p3 + p4) +
  plot_annotation(title = "Transect growing until the 10th marker",
                  subtitle = "Red rectangle = counted region; grey dots = targets; orange diamonds = markers")




###############################

#Linear Count 




set.seed(42)
n_targets <- 2000
n_markers <- 1000

targets <- tibble(x = runif(n_targets, 0, 1),
                  y = runif(n_targets, 0, 1),
                  type = "target")
markers <- tibble(x = runif(n_markers, 0, 1),
                  y = runif(n_markers, 0, 1),
                  type = "marker")
slide <- rbind(targets, markers)

# fix a dimension of the window 
strip_height <- 0.05

# stopping threshold 
target_n_markers <- 10

# the window grows from x = 0 rightward 
# We find the x-coordinate where the 10th marker sits within the strip.
# That's the rightmost edge of the final window.

# Keep only markers inside the strip height
markers_in_strip <- markers[markers$y <= strip_height, ]    #Filter the particles to the strip height

# sort the filtered particles 
markers_in_strip <- markers_in_strip[order(markers_in_strip$x), ] 

# if not enough markers 
if (nrow(markers_in_strip) < target_n_markers) {                                 # if for eg count is 7 and we needed to stop at 10 not enough so change parameters 
  stop("Strip did not contain enough markers to reach the stopping threshold. ",
       "Increase n_markers or strip_height.")
}

# The x-coordinate of the 10th marker is where the window stops growing
x_end <- markers_in_strip$x[target_n_markers]      

# count everything inside the final window 
in_window <- slide$y <= strip_height & slide$x <= x_end    # particles where the slide is less than strip and the slide is less than the strip end 

x_count <- sum(slide$type == "target" & in_window)
n_count <- sum(slide$type == "marker" & in_window)

# should exactly equal the stopping threshold
stopifnot(n_count == target_n_markers)

# plug into the concentration formula
N1 <- 1000   # known marker dose
Y1 <- 1      # no tablet scaling for unit square
V  <- 1      # unit volume

c_est <- (x_count * N1 * Y1) / (n_count * V)



results <- tibble(
  quantity = c("Window height",
               "Window length",
               "Targets counted",
               "Markers counted",
               "Concentration estimate",
               "True concentration"),
  value    = c(strip_height,
               round(x_end, 4),
               x_count,
               n_count,
               c_est,
               n_targets / V)
)

print(results)




###############################


# THE LOOP (generates 1000 estimates)


n_targets        <- 2000
n_markers        <- 1000
strip_height     <- 0.05
target_n_markers <- 10
N1 <- 1000
Y1 <- 1
V  <- 1
num_slides <- 1000

estimates   <- numeric(num_slides) # concentration estimates
x_counts    <- numeric(num_slides) # fossils count
n_counts    <- numeric(num_slides) # markers count



for (s in 1:num_slides) {
  
  targets <- tibble(x = runif(n_targets, 0, 1),
                    y = runif(n_targets, 0, 1),
                    type = "target")
  markers <- tibble(x = runif(n_markers, 0, 1),
                    y = runif(n_markers, 0, 1),
                    type = "marker")
  slide <- rbind(targets, markers)
  
  
  markers_in_strip <- markers[markers$y <= strip_height, ]
  markers_in_strip <- markers_in_strip[order(markers_in_strip$x), ]
  
  # if not enough markers, store NA and skip
  if (nrow(markers_in_strip) < target_n_markers) {
    estimates[s] <- NA
    x_counts[s]  <- NA
    n_counts[s]  <- NA
    next
  }
  
  x_end <- markers_in_strip$x[target_n_markers]
  
  in_window <- slide$y <= strip_height & slide$x <= x_end
  x_count   <- sum(slide$type == "target" & in_window)
  n_count   <- sum(slide$type == "marker" & in_window)
  
  estimates[s] <- (x_count * N1 * Y1) / (n_count * V)
  x_counts[s]  <- x_count
  n_counts[s]  <- n_count
}


# Emperical vs predicted 

# true value
true_concentration <- n_targets / V   # 2000

# Empirical (measured from simulation) 
empirical_mean <- mean(estimates, na.rm = TRUE)
empirical_sd   <- sd(estimates, na.rm = TRUE)
empirical_cv   <- 100 * empirical_sd / empirical_mean

# Predicted from Mays formula
mean_x <- mean(x_counts, na.rm = TRUE)   # average fossils counted
mean_n <- target_n_markers                # always 10
s1P    <- 0                               # tablet error OFF

T_term <- s1P / sqrt(N1)                  # tablet variability
F_term <- sqrt(mean_x) / mean_x           # fossil counting noise
M_term <- sqrt(mean_n) / mean_n           # marker counting noise

predicted_cv <- 100 * sqrt(T_term^2 + F_term^2 + M_term^2)



# Mean estimate vs true concentration

rel_error_concentration <- 100 * (empirical_mean - true_concentration) / true_concentration
acc_concentration       <- 100 * abs(empirical_mean - true_concentration) / true_concentration



# Predicted CV vs empirical CV


rel_error_cv <- 100 * (predicted_cv - empirical_cv) / empirical_cv
acc_cv       <- 100 * abs(predicted_cv - empirical_cv) / empirical_cv



# RESULTS TABLE


checks <- tibble(
  quantity = c(
    "True concentration",
    "Mean of 1000 estimates",
    "Signed relative error of mean (%)",
    "Absolute accuracy of mean (%)",
    "",
    "Empirical CV (%)",
    "Predicted CV from Mays formula (%)",
    "Difference between CVs (pct points)",
    "Signed relative error of CV (%)",
    "Absolute accuracy of CV (%)"
  ),
  value = c(
    round(true_concentration, 1),
    round(empirical_mean, 1),
    round(rel_error_concentration, 3),
    round(acc_concentration, 3),
    NA,
    round(empirical_cv, 2),
    round(predicted_cv, 2),
    round(abs(empirical_cv - predicted_cv), 2),
    round(rel_error_cv, 2),
    round(acc_cv, 2)
  )
)

print(checks)


# If x is Poisson, mean ≈ variance
mean(x_counts, na.rm = TRUE)   # ~19.85
var(x_counts, na.rm = TRUE)    # should be ~19.85 if Poisson

# Equivalently, the CV of x should be ~ 1/sqrt(mean)
sd(x_counts, na.rm = TRUE) / mean(x_counts, na.rm = TRUE)   # observed CV of x
1 / sqrt(mean(x_counts, na.rm = TRUE))                      # Poisson prediction
##############################


# Jensens 

library(tibble)


n_targets        <- 2000
n_markers        <- 1000
strip_height     <- 0.05
N1 <- 1000
Y1 <- 1
V  <- 1
num_slides <- 10000

# Variable stopping threshold

threshold_min <- 5
threshold_max <- 15


estimates <- numeric(num_slides)
x_counts  <- numeric(num_slides)
n_counts  <- numeric(num_slides)
thresholds_used <- integer(num_slides)



for (s in 1:num_slides) {
  
  # each slide uses a random stopping threshold
  this_threshold <- sample(2:20, 1)
  thresholds_used[s] <- this_threshold
  
  targets <- tibble(x = runif(n_targets, 0, 1),
                    y = runif(n_targets, 0, 1),
                    type = "target")
  markers <- tibble(x = runif(n_markers, 0, 1),
                    y = runif(n_markers, 0, 1),
                    type = "marker")
  slide <- rbind(targets, markers)
  
  markers_in_strip <- markers[markers$y <= strip_height, ]
  markers_in_strip <- markers_in_strip[order(markers_in_strip$x), ]
  
  if (nrow(markers_in_strip) < this_threshold) {
    estimates[s] <- NA
    x_counts[s]  <- NA
    n_counts[s]  <- NA
    next
  }
  
  x_end <- markers_in_strip$x[this_threshold]
  
  in_window <- slide$y <= strip_height & slide$x <= x_end
  x_count   <- sum(slide$type == "target" & in_window)
  n_count   <- sum(slide$type == "marker" & in_window)
  
  estimates[s] <- (x_count * N1 * Y1) / (n_count * V)
  x_counts[s]  <- x_count
  n_counts[s]  <- n_count
}



true_concentration <- n_targets / V

# Two ways to compute the mean
naive_mean  <- mean(estimates, na.rm = TRUE)                   # mean of ratios
jensen_mean <- (mean(x_counts, na.rm = TRUE) * N1 * Y1) /      # ratio of means
  (mean(n_counts, na.rm = TRUE) * V)

# Difference between the two methods 
jensen_difference <- naive_mean - jensen_mean
jensen_pct        <- 100 * jensen_difference / jensen_mean

# Accuracy
acc_naive  <- 100 * abs(naive_mean - true_concentration) / true_concentration
acc_jensen <- 100 * abs(jensen_mean - true_concentration) / true_concentration

# errors 
rel_err_naive  <- 100 * (naive_mean - true_concentration) / true_concentration
rel_err_jensen <- 100 * (jensen_mean - true_concentration) / true_concentration

cat("Stopping thresholds used:\n")
cat("  Range: ", range(n_counts, na.rm = TRUE), "\n")
cat("  Mean:  ", round(mean(n_counts, na.rm = TRUE), 2), "\n")
cat("  SD:    ", round(sd(n_counts, na.rm = TRUE), 2), "\n\n")


# Results table

results <- tibble(
  quantity = c(
    "True concentration",
    "",
    "Naive mean (mean of ratios)",
    "Jensen-corrected mean (ratio of means)",
    "Difference (naive - corrected)",
    "Difference as % of corrected mean",
    "",
    "Signed error - naive (%)",
    "Signed error - corrected (%)",
    "Absolute accuracy - naive (%)",
    "Absolute accuracy - corrected (%)"
  ),
  value = c(
    round(true_concentration, 1),
    NA,
    round(naive_mean, 2),
    round(jensen_mean, 2),
    round(jensen_difference, 2),
    round(jensen_pct, 3),
    NA,
    round(rel_err_naive, 3),
    round(rel_err_jensen, 3),
    round(acc_naive, 3),
    round(acc_jensen, 3)
  )
)

print(results)

# Verify the claim directly
mean(x_counts, na.rm = TRUE)              
mean(n_counts, na.rm = TRUE)              
mean(1/n_counts, na.rm = TRUE)            
1/mean(n_counts, na.rm = TRUE)            


mean(1/n_counts, na.rm = TRUE) * mean(n_counts, na.rm = TRUE)


# Tried varying hte tthresholds still no considerable cjhange observed .....



###########################################
#Varying target marker ratio 
library(tibble)
library(ggplot2)


n_markers        <- 5000
strip_height     <- 0.05
target_n_markers <- 10     # fixed stopping threshold
N1 <- 5000                 # must match n_markers
Y1 <- 1
V  <- 1
num_slides <- 1000

# Ratios to test 
ratios <- c(1, 3, 6, 10, 30, 60)


ratio_results <- tibble(
  ratio        = numeric(),
  n_targets    = integer(),
  mean_est     = numeric(),
  true_conc    = numeric(),
  rel_error    = numeric(),
  empirical_cv = numeric(),
  predicted_cv = numeric(),
  cv_agreement = numeric(),
  mean_x       = numeric(),
  failed_slides = integer()
)




for (r in ratios) {
  
  n_targets <- r * n_markers    # this gives the right ratio
  true_c    <- n_targets / V
  
  cat("Running ratio =", r,
      "(n_targets =", n_targets, ", n_markers =", n_markers, ")...\n")
  
  estimates <- numeric(num_slides)
  x_counts  <- numeric(num_slides)
  n_counts  <- numeric(num_slides)
  n_failed  <- 0
  
  for (s in 1:num_slides) {
    
    targets <- tibble(x = runif(n_targets, 0, 1),
                      y = runif(n_targets, 0, 1),
                      type = "target")
    markers <- tibble(x = runif(n_markers, 0, 1),
                      y = runif(n_markers, 0, 1),
                      type = "marker")
    slide <- rbind(targets, markers)
    
    markers_in_strip <- markers[markers$y <= strip_height, ]
    markers_in_strip <- markers_in_strip[order(markers_in_strip$x), ]
    
    if (nrow(markers_in_strip) < target_n_markers) {
      estimates[s] <- NA
      x_counts[s]  <- NA
      n_counts[s]  <- NA
      n_failed     <- n_failed + 1
      next
    }
    
    x_end <- markers_in_strip$x[target_n_markers]
    
    in_window <- slide$y <= strip_height & slide$x <= x_end
    x_count   <- sum(slide$type == "target" & in_window)
    n_count   <- sum(slide$type == "marker" & in_window)
    
    estimates[s] <- (x_count * N1 * Y1) / (n_count * V)
    x_counts[s]  <- x_count
    n_counts[s]  <- n_count
  }
  
  # Analysis
  emp_mean <- mean(estimates, na.rm = TRUE)
  emp_sd   <- sd(estimates, na.rm = TRUE)
  emp_cv   <- 100 * emp_sd / emp_mean
  
  mean_x <- mean(x_counts, na.rm = TRUE)
  mean_n <- target_n_markers
  s1P    <- 0
  
  pred_cv <- 100 * sqrt(
    (s1P / sqrt(N1))^2 +
      (sqrt(mean_x) / mean_x)^2 +
      (sqrt(mean_n) / mean_n)^2
  )
  
  rel_err <- 100 * (emp_mean - true_c) / true_c
  
  # --- Store ---
  ratio_results <- rbind(ratio_results, tibble(
    ratio        = r,
    n_targets    = n_targets,
    mean_est     = round(emp_mean, 1),
    true_conc    = true_c,
    rel_error    = round(rel_err, 3),
    empirical_cv = round(emp_cv, 1),
    predicted_cv = round(pred_cv, 1),
    cv_agreement = round(abs(emp_cv - pred_cv), 1),
    mean_x       = round(mean_x, 1),
    failed_slides = n_failed
  ))
}

print(ratio_results)



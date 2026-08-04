library(tibble)
library(ggplot2)


# SLIDE 
n <- 10   

x <- runif(n, min = 0, max = 1)
y <- runif(n, min = 0, max = 1)

points <- tibble(x = x, y = y)

print(points)

ggplot(points, aes(x = x, y = y)) +
  geom_point() 


# Distances 
# SHOWS a matrix form ??

# INTER EVENT DISTANCES 

distances <- dist(points)
distances
#convert to vector 
dist_vector <- as.vector(distances)
dist_vector

hist(dist_vector,
     breaks = 20,
     main = "Distribution of Pairwise Distances",
     xlab = "Distance")

# Number of simulations
# Slide 

num_simulations <- 5

# Create an empty list to store results
all_simulations_loop <- list()

# Run the loop
for (i in 1:num_simulations) {
  # Generate data and store in the list
  all_simulations_loop[[i]] <- runif(n, min = 0, max = 1)
}


# Number of simulations
num_simulations <- 5

# Empty list
all_simulations_loop <- list()

for (i in 1:num_simulations) {
  
  x <- runif(n, min = 0, max = 1)
  y <- runif(n, min = 0, max = 1)
  
  sim_points <- tibble(x = x, y = y)
  
  all_simulations_loop[[i]] <- sim_points
  
}


all_simulations_loop[[1]]

ggplot(all_simulations_loop[[1]], aes(x = x, y = y)) +
  geom_point(size = 3)



# adding list for inter event distance
all_distances <- list()

num_simulations <- 5

for (i in 1:num_simulations) {
  
  x <- runif(n, min = 0, max = 1)
  y <- runif(n, min = 0, max = 1)
  
  sim_points <- tibble(x = x, y = y)
  

  all_simulations_loop[[i]] <- sim_points
  

  distances <- dist(sim_points)
  dist_vector <- as.vector(distances)
  all_distances[[i]] <- dist_vector
  
}
all_simulations_loop[[1]]  
all_distances[[1]]         
all_simulations_loop[[2]]  
all_distances[[2]]



############

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

# --- Diggle's theoretical CDF (Eqn 2.2 for unit square) ---
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


# i DOnt think will use this ....
ggplot() +
  geom_function(fun = diggle_cdf, xlim = c(0, sqrt(2))) +
  labs(title = "Diggle's theoretical CDF for inter-event distances",
       subtitle = "Unit square, Eqn 2.2",
       x = "Distance t",
       y = "H(t)") +
  theme_minimal()



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
# STEP 1: Generate many simulated Poisson slides


set.seed(42)
n_points        <- 100   # points per slide
num_simulations <- 50    # number of slides

all_distances <- vector("list", num_simulations)   # Empty list with num_simulation 

for (i in 1:num_simulations) {
  pts <- tibble(x = runif(n_points, 0, 1),     # points tiblle containg x and y 
                y = runif(n_points, 0, 1))
  all_distances[[i]] <- as.vector(dist(pts))    # Note [[i]] are how we access the individual element otherwise [] will give a sublist 
}

# Pool all the distances into one big vector
pooled_distances <- unlist(all_distances)
length(pooled_distances)   # should be 50 * 100*99/2 = 247,500


# ============================================================
# STEP 3: Build empirical CDF and pair it with theoretical

# tibble contaiing x distances y empericala nd theoretical values 
comparison <- tibble( 
  distance    = sort(pooled_distances),
  empirical   = (1:length(pooled_distances)) / length(pooled_distances)
)
comparison$theoretical <- diggle_cdf(comparison$distance)


ggplot(comparison, aes(x = distance)) +
  geom_step(aes(y = empirical),   colour = "blue") +
  geom_line(aes(y = theoretical), colour = "red", linewidth = 1) +
  labs(title    = "Empirical vs theoretical CDF of inter-event distances",
       subtitle = "Blue = simulation, Red = Diggle Eqn 2.2",
       x        = "Distance t",
       y        = "H(t)") +
  theme_minimal()

# OR PLOTTING DIFFERENT IN SAME WINDOW 

library(patchwork)

p1 <- ggplot(comparison, aes(x = distance, y = empirical)) +
  geom_step(colour = "blue") +
  labs(title = "our simulation",
       x = "Distance t", y = "H(t)") +
  theme_minimal()

p2 <- ggplot(comparison, aes(x = distance, y = theoretical)) +
  geom_line(colour = "red", linewidth = 1) +
  labs(title = "Diggle eq ",
       x = "Distance t", y = "H(t)") +
  theme_minimal()

# side by side
p1 + p2

# or stacked vertically
p1 / p2


# NEAREST NEIGHBOUR 
# ============================================================
# STEP 1: Compute nearest-neighbour distances for each slide
# ============================================================

# Manual version 
nn_dist_by_hand <- function(points) {
  n <- nrow(points)  # total points 
  nn <- numeric(n)   # vector of size total points
  for (i in 1:n) {
    min_d <- Inf
    for (j in 1:n) {
      if (i != j) {   # skip the point comapring to itself 
        dx <- points$x[i] - points$x[j]
        dy <- points$y[i] - points$y[j]
        d  <- sqrt(dx^2 + dy^2)
        if (d < min_d) min_d <- d  #If the distance we just computed is smaller than the smallest we've seen so far, update min_d. Otherwise leave it alone.
      }
    }
    nn[i] <- min_d
  }
  nn
}
#############
# Built in function spatstat 
#library(spatstat)

# Your points as a tibble
pts <- tibble(x = runif(10), y = runif(10))

# Get nearest-neighbour distances
nn_spatstat <- nndist(pts$x, pts$y)

#######################



# NOW WONT NEED THIS AS WE HAVE NNDIST 
# Faster version using dist(), for production runs
nn_dist_fast <- function(points) {
  dmat <- as.matrix(dist(points))
  diag(dmat) <- Inf
  apply(dmat, 1, min)
}

# Verify they agree
set.seed(1)
check_pts <- tibble(x = runif(20), y = runif(20))
all(abs(
  sort(nn_dist_by_hand(check_pts)) - sort(nn_dist_fast(check_pts))
) < 1e-10)   # should be TRUE


# ============================================================
# STEP 2: Run many simulations and pool NN distances
# ============================================================

set.seed(42)
n_points        <- 100
num_simulations <- 50

all_nn <- vector("list", num_simulations)

for (i in 1:num_simulations) {
  pts <- tibble(x = runif(n_points, 0, 1),
                y = runif(n_points, 0, 1))
  all_nn[[i]] <- nn_dist_fast(pts)
}

pooled_nn <- unlist(all_nn)
length(pooled_nn)   # should be num_simulations * n_points = 5000


# ============================================================
# STEP 3: Diggle's theoretical CDF for NN distances (Eqn 2.8)
# ============================================================

# G(y) = 1 - exp(-lambda * pi * y^2)
# For a unit square with n points, lambda = n / 1 = n
diggle_nn_cdf <- function(y, n) {
  1 - exp(-n * pi * y^2)
}

# Quick sanity checks
diggle_nn_cdf(0, n_points)       # should be 0
diggle_nn_cdf(0.2, n_points)     # should be close to 1 (for n = 100)
# --- Theoretical density (derivative of G) ---


# ============================================================
# STEP 4: Build empirical CDF and compare
# ============================================================

nn_comparison <- tibble(
  distance    = sort(pooled_nn),
  empirical   = (1:length(pooled_nn)) / length(pooled_nn)
)
nn_comparison$theoretical <- diggle_nn_cdf(nn_comparison$distance, n_points)


# ============================================================
# STEP 5: Plot empirical vs theoretical
# ============================================================

ggplot(nn_comparison, aes(x = distance)) +
  geom_step(aes(y = empirical),   colour = "blue") +
  geom_line(aes(y = theoretical), colour = "red", linewidth = 1) +
  labs(title    = "Empirical vs theoretical CDF of nearest-neighbour distances",
       subtitle = "Blue = simulation, Red = Diggle Eqn 2.8",
       x        = "Nearest-neighbour distance y",
       y        = "G(y)") +
  theme_minimal()



# MANUAL NEAREST NEIGHBOUR 

library(tibble)
library(ggplot2)

# --- Manual nearest-neighbour distance calculation ---
nn_dist_by_hand <- function(points) {
  n <- nrow(points)
  nn <- numeric(n)   # preallocate: one NN distance per point
  
  for (i in 1:n) {
    min_d <- Inf   # start with "infinity" so any real distance beats it
    
    for (j in 1:n) {
      if (i != j) {   # skip comparing a point to itself
        dx <- points$x[i] - points$x[j]
        dy <- points$y[i] - points$y[j]
        d  <- sqrt(dx^2 + dy^2)
        
        if (d < min_d) min_d <- d   # keep track of smallest so far
      }
    }
    
    nn[i] <- min_d   # store the nearest-neighbour distance for point i
  }
  
  nn
}


demo_points <- tibble(x = runif(5), y = runif(5))
print(demo_points)

nn_demo <- nn_dist_by_hand(demo_points)
print(nn_demo)


# Fast version using dist() for comparison
nn_dist_fast <- function(points) {
  dmat <- as.matrix(dist(points))
  diag(dmat) <- Inf
  apply(dmat, 1, min)
}

# Verify the two agree
set.seed(1)
check_pts <- tibble(x = runif(20), y = runif(20))

manual_result <- nn_dist_by_hand(check_pts)
fast_result   <- nn_dist_fast(check_pts)

all(abs(manual_result - fast_result) < 1e-10)   # should be TRUE

set.seed(42)
n_points        <- 100
num_simulations <- 50

all_nn <- vector("list", num_simulations)

for (i in 1:num_simulations) {
  pts <- tibble(x = runif(n_points, 0, 1),
                y = runif(n_points, 0, 1))
  all_nn[[i]] <- nn_dist_fast(pts)   # using the fast version for speed
}

pooled_nn <- unlist(all_nn)
length(pooled_nn)   # should be 50 * 100 = 5000

ggplot(tibble(nn = pooled_nn), aes(x = nn)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 40, fill = "lightblue", colour = "black") +
  labs(title = "Empirical nearest-neighbour distances (simulation only)",
       x = "Nearest-neighbour distance y",
       y = "Density") +
  theme_minimal()


########################################

#add the markers

# --- Generate one slide with both targets and markers ---

n_targets <- 200   #  targets 
n_markers <- 100    #  markers

targets <- tibble(x = runif(n_targets, 0, 1),
                  y = runif(n_targets, 0, 1),
                  type = "target")

markers <- tibble(x = runif(n_markers, 0, 1),
                  y = runif(n_markers, 0, 1),
                  type = "marker")

slide <- rbind(targets, markers)

# --- Plot both together ---
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
#ALTERNATIVE METHOD IMPLementation
library(tibble)
library(ggplot2)

# --- REQUIREMENT: a slide with targets and markers ---
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

# --- REQUIREMENT 2: fix a dimension of the window ---
strip_height <- 0.05

# --- REQUIREMENT 4: the stopping condition ---
target_n_markers <- 10

# --- REQUIREMENT 3: the window grows from x = 0 rightward ---
#    We find the x-coordinate where the 10th marker sits within the strip.
#    That's the rightmost edge of the final window.

# Keep only markers inside the strip height
markers_in_strip <- markers[markers$y <= strip_height, ]    #Filter the particles to the strip height

# Sort them left to right
markers_in_strip <- markers_in_strip[order(markers_in_strip$x), ] # sort the filtered particles 

# --- REQUIREMENT 6: handle "not enough markers" case ---
if (nrow(markers_in_strip) < target_n_markers) {                                 # if for eg count is 7 and we needed to stop at 10 not enough so change parameters 
  stop("Strip did not contain enough markers to reach the stopping threshold. ",
       "Increase n_markers or strip_height.")
}

# The x-coordinate of the 10th marker is where the window stops growing
x_end <- markers_in_strip$x[target_n_markers]      

# --- REQUIREMENT 5: count everything inside the final window ---
in_window <- slide$y <= strip_height & slide$x <= x_end    # particles where the slide is less than strip and the slide is less than the strip end 

x_count <- sum(slide$type == "target" & in_window)
n_count <- sum(slide$type == "marker" & in_window)

# Sanity check — n_count should exactly equal the stopping threshold
stopifnot(n_count == target_n_markers)

# --- REQUIREMENT 7: plug into the concentration formula ---
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


# ============================================================
# PART 1: THE LOOP (generates 1000 estimates)
# ============================================================

library(tibble)

n_targets        <- 2000
n_markers        <- 1000
strip_height     <- 0.05
target_n_markers <- 10
N1 <- 1000
Y1 <- 1
V  <- 1
num_slides <- 1000

estimates   <- numeric(num_slides) # concentrstion estimates 
x_counts    <- numeric(num_slides) # fossils count 
n_counts    <- numeric(num_slides) # markers count 

set.seed(42)

for (s in 1:num_slides) {
  # new slide 
  targets <- tibble(x = runif(n_targets, 0, 1),
                    y = runif(n_targets, 0, 1),
                    type = "target")
  markers <- tibble(x = runif(n_markers, 0, 1),
                    y = runif(n_markers, 0, 1),
                    type = "marker")
  slide <- rbind(targets, markers)
  # finding the stopping marker 
  markers_in_strip <- markers[markers$y <= strip_height, ]
  markers_in_strip <- markers_in_strip[order(markers_in_strip$x), ]
  # if the slide doesnt have enough markers in the strip then we store na values and skip the slide
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
  # save results to the slide s 
  estimates[s] <- (x_count * N1 * Y1) / (n_count * V)
  x_counts[s]  <- x_count
  n_counts[s]  <- n_count
}

# ============================================================
# PART 2: AFTER THE LOOP — compute both errors and compare
# ============================================================
# true concentration 
true_concentration <- n_targets / V   # 2000

# --- The ACTUAL error (measured from the simulation) ---
empirical_mean <- mean(estimates, na.rm = TRUE)
empirical_sd   <- sd(estimates, na.rm = TRUE)
empirical_cv   <- 100 * empirical_sd / empirical_mean

# --- The PREDICTED error (from the Mays formula) ---
mean_x <- mean(x_counts, na.rm = TRUE)   # average fossils counted
mean_n <- target_n_markers                # always 10
s1P    <- 0                               # tablet error OFF for now

T_term <- s1P / sqrt(N1)                  # tablet variability
F_term <- sqrt(mean_x) / mean_x           # fossil counting noise
M_term <- sqrt(mean_n) / mean_n           # marker counting noise

predicted_cv <- 100 * sqrt(T_term^2 + F_term^2 + M_term^2)

# --- Accuracy check: is the mean close to the truth? ---
relative_error <- 100 * (empirical_mean - true_concentration) / true_concentration

# --- Put everything in a tibble ---
checks <- tibble(
  quantity = c("True concentration",
               "Mean of 1000 estimates",
               "Relative error of mean (%)",
               "Empirical CV (%)",
               "Predicted CV from Mays formula (%)",
               "Agreement between CVs (%)"),
  value    = c(true_concentration,
               round(empirical_mean, 1),
               round(relative_error, 3),
               round(empirical_cv, 1),
               round(predicted_cv, 1),
               round(abs(empirical_cv - predicted_cv), 1))
)

print(checks)




##############################
# ============================================================
# FULL CODE: loop + accuracy checks + Jensen's comparison
# ============================================================

library(tibble)

n_targets        <- 2000
n_markers        <- 1000
strip_height     <- 0.05
target_n_markers <- 10
N1 <- 1000
Y1 <- 1
V  <- 1
num_slides <- 1000

estimates <- numeric(num_slides)
x_counts  <- numeric(num_slides)
n_counts  <- numeric(num_slides)

set.seed(42)

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

# ============================================================
# ANALYSIS
# ============================================================

true_concentration <- n_targets / V

# --- Two ways to compute the mean ---
naive_mean  <- mean(estimates, na.rm = TRUE)
jensen_mean <- (mean(x_counts, na.rm = TRUE) * N1 * Y1) /       # jensens 
  (mean(n_counts, na.rm = TRUE) * V)

# --- Accuracy ---
relative_error_naive  <- 100 * (naive_mean - true_concentration) / true_concentration
relative_error_jensen <- 100 * (jensen_mean - true_concentration) / true_concentration

# --- Precision: empirical vs predicted ---
empirical_sd <- sd(estimates, na.rm = TRUE)
empirical_cv <- 100 * empirical_sd / naive_mean

mean_x <- mean(x_counts, na.rm = TRUE)
mean_n <- target_n_markers
s1P    <- 0

predicted_cv <- 100 * sqrt(
  (s1P / sqrt(N1))^2 +
    (sqrt(mean_x) / mean_x)^2 +
    (sqrt(mean_n) / mean_n)^2
)

# --- True SE (method: from deviations around truth) ---
true_se <- sqrt(sum((estimates - true_concentration)^2, na.rm = TRUE) / 
                  (sum(!is.na(estimates)) - 1))
true_cv <- 100 * true_se / true_concentration

# --- Results ---
results <- tibble(
  quantity = c("True concentration",
               "Naive mean (mean of ratios)",
               "Jensen-corrected mean (ratio of means)",
               "Relative error - naive (%)",
               "Relative error - corrected (%)",
               "Empirical CV (%)",
               "Predicted CV - Mays formula (%)",
               "True CV - deviations from truth (%)",
               "CV agreement (empirical vs predicted)"),
  value    = c(true_concentration,
               round(naive_mean, 1),
               round(jensen_mean, 1),
               round(relative_error_naive, 3),
               round(relative_error_jensen, 3),
               round(empirical_cv, 1),
               round(predicted_cv, 1),
               round(true_cv, 1),
               paste0(round(abs(empirical_cv - predicted_cv), 1), " pct points"))
)

print(results)



############################################
# Varying  the stopping thresholds 

library(tibble)
library(ggplot2)

# --- Parameters ---
n_targets    <- 2000
strip_height <- 0.05
N1 <- 1000
Y1 <- 1
V  <- 1
num_slides   <- 1000

# Thresholds to test
thresholds <- c(5, 10, 20, 50, 100, 200)
n_markers <- 5000

# --- Storage for summary results ---
sweep_results <- tibble(
  threshold    = integer(),
  mean_est     = numeric(),
  rel_error    = numeric(),
  empirical_cv = numeric(),
  predicted_cv = numeric(),
  cv_agreement = numeric(),
  failed_slides = integer()
)

# --- Run the sweep ---
set.seed(42)

for (thresh in thresholds) {
  
  cat("Running threshold =", thresh, "...\n")
  
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
    
    if (nrow(markers_in_strip) < thresh) {
      estimates[s] <- NA
      x_counts[s]  <- NA
      n_counts[s]  <- NA
      n_failed     <- n_failed + 1
      next
    }
    
    x_end <- markers_in_strip$x[thresh]
    
    in_window <- slide$y <= strip_height & slide$x <= x_end
    x_count   <- sum(slide$type == "target" & in_window)
    n_count   <- sum(slide$type == "marker" & in_window)
    
    estimates[s] <- (x_count * N1 * Y1) / (n_count * V)
    x_counts[s]  <- x_count
    n_counts[s]  <- n_count
  }
  
  # --- Analysis for this threshold ---
  true_c <- n_targets / V
  
  emp_mean <- mean(estimates, na.rm = TRUE)
  emp_sd   <- sd(estimates, na.rm = TRUE)
  emp_cv   <- 100 * emp_sd / emp_mean
  
  mean_x <- mean(x_counts, na.rm = TRUE)
  mean_n <- thresh
  s1P    <- 0
  
  pred_cv <- 100 * sqrt(
    (s1P / sqrt(N1))^2 +
      (sqrt(mean_x) / mean_x)^2 +
      (sqrt(mean_n) / mean_n)^2
  )
  
  rel_err <- 100 * (emp_mean - true_c) / true_c
  
  # --- Store ---
  sweep_results <- rbind(sweep_results, tibble(
    threshold    = thresh,
    mean_est     = round(emp_mean, 1),
    rel_error    = round(rel_err, 3),
    empirical_cv = round(emp_cv, 1),
    predicted_cv = round(pred_cv, 1),
    cv_agreement = round(abs(emp_cv - pred_cv), 1),
    failed_slides = n_failed
  ))
}

print(sweep_results)




# precision effort curve 
ggplot(sweep_results, aes(x = threshold, y = empirical_cv)) +
  geom_line(colour = "steelblue", linewidth = 1.2) +
  geom_point(colour = "steelblue", size = 3) +
  geom_line(aes(y = predicted_cv), colour = "red", 
            linewidth = 1, linetype = "dashed") +
  geom_point(aes(y = predicted_cv), colour = "red", size = 2) +
  labs(title    = "Precision vs effort: how many markers to count?",
       subtitle = "Blue = empirical CV, Red dashed = Mays formula prediction",
       x = "Markers counted (stopping threshold)",
       y = "Coefficient of variation (%)") +
  theme_minimal(base_size = 13)




###########################################
#Varying target marker ratio 
library(tibble)
library(ggplot2)

# --- Fixed parameters ---
n_markers        <- 5000
strip_height     <- 0.05
target_n_markers <- 10     # fixed stopping threshold
N1 <- 5000                 # must match n_markers
Y1 <- 1
V  <- 1
num_slides <- 1000

# --- Ratios to test ---
ratios <- c(1, 3, 6, 10, 30, 60)

# --- Storage ---
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

# --- Run the sweep ---
set.seed(42)

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
  
  # --- Analysis for this ratio ---
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



# precision vs ratio 
ggplot(ratio_results, aes(x = ratio)) +
  geom_line(aes(y = empirical_cv), colour = "steelblue", linewidth = 1.2) +
  geom_point(aes(y = empirical_cv), colour = "steelblue", size = 3) +
  geom_line(aes(y = predicted_cv), colour = "red",
            linewidth = 1, linetype = "dashed") +
  geom_point(aes(y = predicted_cv), colour = "red", size = 2) +
  scale_x_log10(breaks = ratios) +
  labs(title    = "How target-to-marker ratio affects precision",
       subtitle = "Blue = empirical CV, Red dashed = Mays formula prediction",
       x = "Target-to-marker ratio (log scale)",
       y = "Coefficient of variation (%)") +
  theme_minimal(base_size = 13)


########################################################

# ============================================================
# TABLET ERROR EXTENSION
# ============================================================

library(tibble)
library(ggplot2)

# --- Fixed parameters ---
n_targets        <- 2000
strip_height     <- 0.05
target_n_markers <- 10
V  <- 1
num_slides <- 1000

# --- Marker dose parameters ---
# mean = 12,489 spores, SD = 491

marker_label     <- 1000                         
tablet_cv        <- 491 / 12489                   
marker_sd        <- marker_label * tablet_cv      

# --- Storage ---
estimates_on  <- numeric(num_slides)   # tablet error ON
x_counts_on   <- numeric(num_slides)
n_counts_on   <- numeric(num_slides)
actual_doses  <- numeric(num_slides)   # track what each slide actually got

# --- The loop ---
set.seed(42)

for (s in 1:num_slides) {
  
  # THE KEY CHANGE: actual marker count varies each slide
  n_markers_actual <- round(rnorm(1, mean = marker_label, sd = marker_sd))
  
  
  if (n_markers_actual < 1) n_markers_actual <- 1
  
  actual_doses[s] <- n_markers_actual
  
  # Generate the slide with the ACTUAL (random) number of markers
  targets <- tibble(x = runif(n_targets, 0, 1),
                    y = runif(n_targets, 0, 1),
                    type = "target")
  markers <- tibble(x = runif(n_markers_actual, 0, 1),
                    y = runif(n_markers_actual, 0, 1),
                    type = "marker")
  slide <- rbind(targets, markers)
  
  # Transect counting (same as before)
  markers_in_strip <- markers[markers$y <= strip_height, ]
  markers_in_strip <- markers_in_strip[order(markers_in_strip$x), ]
  
  if (nrow(markers_in_strip) < target_n_markers) {
    estimates_on[s] <- NA
    x_counts_on[s]  <- NA
    n_counts_on[s]  <- NA
    next
  }
  
  x_end <- markers_in_strip$x[target_n_markers]
  
  in_window <- slide$y <= strip_height & slide$x <= x_end
  x_count   <- sum(slide$type == "target" & in_window)
  n_count   <- sum(slide$type == "marker" & in_window)
  
  # THE FORMULA USES THE LABEL VALUE, NOT THE ACTUAL DOSE
  estimates_on[s] <- (x_count * marker_label) / (n_count * V)
  x_counts_on[s]  <- x_count
  n_counts_on[s]  <- n_count
}


# ============================================================
# ANALYSIS: tablet error ON
# ============================================================

true_concentration <- n_targets / V   # 2000

# --- Empirical results ---
mean_on  <- mean(estimates_on, na.rm = TRUE)
sd_on    <- sd(estimates_on, na.rm = TRUE)
cv_on    <- 100 * sd_on / mean_on

# --- Predicted CV WITH tablet error ---
mean_x <- mean(x_counts_on, na.rm = TRUE)
mean_n <- target_n_markers   # still 10

# The T term is now NON-ZERO
s1P    <- tablet_cv          #  0.0393
N1_tablets <- 1              # one tablet used

T_term <- s1P / sqrt(N1_tablets)      #  0.0393
F_term <- sqrt(mean_x) / mean_x      # fossil noise
M_term <- sqrt(mean_n) / mean_n      # marker noise

predicted_cv_on <- 100 * sqrt(T_term^2 + F_term^2 + M_term^2)

# --- Relative error ---
rel_error_on <- 100 * (mean_on - true_concentration) / true_concentration

# ============================================================
# BASELINE: tablet error OFF 
# ============================================================

estimates_off <- numeric(num_slides)
x_counts_off  <- numeric(num_slides)
n_counts_off  <- numeric(num_slides)

set.seed(42)

for (s in 1:num_slides) {
  
  # Fixed number of markers (no tablet error)
  n_markers_fixed <- marker_label   # always exactly 1000
  
  targets <- tibble(x = runif(n_targets, 0, 1),
                    y = runif(n_targets, 0, 1),
                    type = "target")
  markers <- tibble(x = runif(n_markers_fixed, 0, 1),
                    y = runif(n_markers_fixed, 0, 1),
                    type = "marker")
  slide <- rbind(targets, markers)
  
  markers_in_strip <- markers[markers$y <= strip_height, ]
  markers_in_strip <- markers_in_strip[order(markers_in_strip$x), ]
  
  if (nrow(markers_in_strip) < target_n_markers) {
    estimates_off[s] <- NA
    x_counts_off[s]  <- NA
    n_counts_off[s]  <- NA
    next
  }
  
  x_end <- markers_in_strip$x[target_n_markers]
  
  in_window <- slide$y <= strip_height & slide$x <= x_end
  x_count   <- sum(slide$type == "target" & in_window)
  n_count   <- sum(slide$type == "marker" & in_window)
  
  estimates_off[s] <- (x_count * marker_label) / (n_count * V)
  x_counts_off[s]  <- x_count
  n_counts_off[s]  <- n_count
}

# --- Baseline analysis ---
mean_off <- mean(estimates_off, na.rm = TRUE)
sd_off   <- sd(estimates_off, na.rm = TRUE)
cv_off   <- 100 * sd_off / mean_off

mean_x_off <- mean(x_counts_off, na.rm = TRUE)
predicted_cv_off <- 100 * sqrt(
  (0 / sqrt(1))^2 +                                 # T = 0
    (sqrt(mean_x_off) / mean_x_off)^2 +
    (sqrt(target_n_markers) / target_n_markers)^2
)

rel_error_off <- 100 * (mean_off - true_concentration) / true_concentration


# ============================================================
# COMPARISON: with vs without tablet error
# ============================================================

comparison <- tibble(
  quantity = c("True concentration",
               "Mean estimate",
               "Relative error (%)",
               "Empirical CV (%)",
               "Predicted CV (%)",
               "CV agreement (pct points)",
               "Empirical SD"),
  without_tablet = c(true_concentration,
                     round(mean_off, 1),
                     round(rel_error_off, 3),
                     round(cv_off, 1),
                     round(predicted_cv_off, 1),
                     round(abs(cv_off - predicted_cv_off), 1),
                     round(sd_off, 1)),
  with_tablet    = c(true_concentration,
                     round(mean_on, 1),
                     round(rel_error_on, 3),
                     round(cv_on, 1),
                     round(predicted_cv_on, 1),
                     round(abs(cv_on - predicted_cv_on), 1),
                     round(sd_on, 1))
)

print(comparison)


# ============================================================
# VISUAL: overlaid histograms
# ============================================================

plot_data <- tibble(
  estimate = c(estimates_off, estimates_on),
  scenario = factor(
    c(rep("Without tablet error", num_slides),
      rep("With tablet error", num_slides)),
    levels = c("Without tablet error", "With tablet error")
  )
)

ggplot(plot_data[!is.na(plot_data$estimate), ], 
       aes(x = estimate, fill = scenario)) +
  geom_histogram(alpha = 0.5, position = "identity",
                 bins = 40, colour = "black", linewidth = 0.2) +
  geom_vline(xintercept = 2000, colour = "red", linewidth = 1) +
  scale_fill_manual(values = c("steelblue", "orange")) +
  labs(title    = "Effect of tablet variability on concentration estimates",
       subtitle = "Red line = true concentration (2000)",
       x = "Estimated concentration",
       y = "Number of slides",
       fill = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")



###############################################################
# ============================================================
# TASK 9: Threshold sweep WITH tablet error
# ============================================================

library(tibble)
library(ggplot2)

# --- Fixed parameters ---
n_targets    <- 2000
strip_height <- 0.05
V  <- 1
num_slides   <- 1000

# --- Marker dose parameters  ---
marker_label <- 5000                          
tablet_cv    <- 491 / 12489                    
marker_sd    <- marker_label * tablet_cv      

# --- Thresholds to test ---
thresholds <- c(5, 10, 20, 50, 100, 200)

# --- Storage for BOTH scenarios ---
sweep_off <- tibble(
  threshold = integer(), empirical_cv = numeric(),
  predicted_cv = numeric(), rel_error = numeric(),
  scenario = character()
)

sweep_on <- tibble(
  threshold = integer(), empirical_cv = numeric(),
  predicted_cv = numeric(), rel_error = numeric(),
  scenario = character()
)

# ============================================================
# RUN 1: tablet error OFF
# ============================================================

cat("=== TABLET ERROR OFF ===\n")
set.seed(42)

for (thresh in thresholds) {
  
  cat("  Threshold =", thresh, "...\n")
  
  estimates <- numeric(num_slides)
  x_counts  <- numeric(num_slides)
  
  for (s in 1:num_slides) {
    targets <- tibble(x = runif(n_targets), y = runif(n_targets), type = "target")
    markers <- tibble(x = runif(marker_label), y = runif(marker_label), type = "marker")
    slide <- rbind(targets, markers)
    
    ms <- markers[markers$y <= strip_height, ]
    ms <- ms[order(ms$x), ]
    
    if (nrow(ms) < thresh) { estimates[s] <- NA; x_counts[s] <- NA; next }
    
    x_end <- ms$x[thresh]
    in_w <- slide$y <= strip_height & slide$x <= x_end
    x_count <- sum(slide$type == "target" & in_w)
    n_count <- sum(slide$type == "marker" & in_w)
    
    estimates[s] <- (x_count * marker_label) / (n_count * V)
    x_counts[s]  <- x_count
  }
  
  true_c   <- n_targets / V
  emp_mean <- mean(estimates, na.rm = TRUE)
  emp_sd   <- sd(estimates, na.rm = TRUE)
  emp_cv   <- 100 * emp_sd / emp_mean
  
  mean_x <- mean(x_counts, na.rm = TRUE)
  s1P    <- 0
  pred_cv <- 100 * sqrt(
    (s1P / sqrt(1))^2 +
      (sqrt(mean_x) / mean_x)^2 +
      (sqrt(thresh) / thresh)^2
  )
  
  rel_err <- 100 * (emp_mean - true_c) / true_c
  
  sweep_off <- rbind(sweep_off, tibble(
    threshold = thresh, empirical_cv = round(emp_cv, 1),
    predicted_cv = round(pred_cv, 1), rel_error = round(rel_err, 3),
    scenario = "Without tablet error"
  ))
}

# ============================================================
# RUN 2: tablet error ON
# ============================================================

cat("\n=== TABLET ERROR ON ===\n")
set.seed(42)

for (thresh in thresholds) {
  
  cat("  Threshold =", thresh, "...\n")
  
  estimates <- numeric(num_slides)
  x_counts  <- numeric(num_slides)
  
  for (s in 1:num_slides) {
    
    # THE KEY CHANGE: actual markers varies each slide
    n_markers_actual <- round(rnorm(1, mean = marker_label, sd = marker_sd))
    if (n_markers_actual < 1) n_markers_actual <- 1
    
    targets <- tibble(x = runif(n_targets), y = runif(n_targets), type = "target")
    markers <- tibble(x = runif(n_markers_actual), y = runif(n_markers_actual), type = "marker")
    slide <- rbind(targets, markers)
    
    ms <- markers[markers$y <= strip_height, ]
    ms <- ms[order(ms$x), ]
    
    if (nrow(ms) < thresh) { estimates[s] <- NA; x_counts[s] <- NA; next }
    
    x_end <- ms$x[thresh]
    in_w <- slide$y <= strip_height & slide$x <= x_end
    x_count <- sum(slide$type == "target" & in_w)
    n_count <- sum(slide$type == "marker" & in_w)
    
    # Formula still uses LABEL value, not actual dose
    estimates[s] <- (x_count * marker_label) / (n_count * V)
    x_counts[s]  <- x_count
  }
  
  true_c   <- n_targets / V
  emp_mean <- mean(estimates, na.rm = TRUE)
  emp_sd   <- sd(estimates, na.rm = TRUE)
  emp_cv   <- 100 * emp_sd / emp_mean
  
  mean_x <- mean(x_counts, na.rm = TRUE)
  s1P    <- tablet_cv    # NOW NON-ZERO
  pred_cv <- 100 * sqrt(
    (s1P / sqrt(1))^2 +
      (sqrt(mean_x) / mean_x)^2 +
      (sqrt(thresh) / thresh)^2
  )
  
  rel_err <- 100 * (emp_mean - true_c) / true_c
  
  sweep_on <- rbind(sweep_on, tibble(
    threshold = thresh, empirical_cv = round(emp_cv, 1),
    predicted_cv = round(pred_cv, 1), rel_error = round(rel_err, 3),
    scenario = "With tablet error"
  ))
}

# ============================================================
# COMPARISON TABLE
# ============================================================

comparison <- tibble(
  threshold     = sweep_off$threshold,
  cv_off        = sweep_off$empirical_cv,
  cv_on         = sweep_on$empirical_cv,
  cv_difference = round(sweep_on$empirical_cv - sweep_off$empirical_cv, 1),
  pred_off      = sweep_off$predicted_cv,
  pred_on       = sweep_on$predicted_cv
)

print(comparison)



# ============================================================
# PLOT: precision vs effort, with and without tablet error
# ============================================================

both <- rbind(sweep_off, sweep_on)

ggplot(both, aes(x = threshold, colour = scenario)) +
  # Empirical (solid lines)
  geom_line(aes(y = empirical_cv, linetype = "Empirical"), linewidth = 1.2) +
  geom_point(aes(y = empirical_cv), size = 3) +
  # Predicted (dashed lines)
  geom_line(aes(y = predicted_cv, linetype = "Predicted"), linewidth = 0.8) +
  geom_point(aes(y = predicted_cv), size = 2, shape = 1) +
  scale_colour_manual(values = c("steelblue", "orange")) +
  scale_linetype_manual(values = c(Empirical = "solid", Predicted = "dashed")) +
  scale_x_log10(breaks = thresholds) +
  labs(title    = "Precision vs effort: effect of tablet variability",
       subtitle = "Solid = empirical, Dashed = Mays formula prediction",
       x = "Markers counted (stopping threshold, log scale)",
       y = "Total error (%)",
       colour = NULL, linetype = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")
# ============================================================
# PLOT: error breakdown showing T, F, M contributions
# ============================================================

breakdown <- tibble()

for (i in 1:nrow(sweep_on)) {
  thresh <- sweep_on$threshold[i]
  
  # Approximate mean_x from the ratio: at our parameters,
  # mean_x ≈ (n_targets / marker_label) * thresh
  approx_mean_x <- (n_targets / marker_label) * thresh
  
  T_sq <- (tablet_cv / sqrt(1))^2
  F_sq <- 1 / approx_mean_x
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
  labs(title    = "Error breakdown: which source dominates at each threshold?",
       subtitle = "At low thresholds, marker noise dominates. At high thresholds, tablet error grows.",
       x = "Stopping threshold",
       y = "Contribution to total error (%)",
       fill = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

# ============================================================
# TASK 10: Ratio sweep WITH tablet error
# ============================================================

library(tibble)
library(ggplot2)

# --- Fixed parameters ---
n_markers        <- 5000
strip_height     <- 0.05
target_n_markers <- 10
V  <- 1
num_slides <- 1000

# --- Marker dose parameters ---
marker_label <- n_markers
tablet_cv    <- 491 / 12489
marker_sd    <- marker_label * tablet_cv

# --- Ratios to test ---
ratios <- c(1, 3, 6, 10, 30, 60)

# --- Storage ---
ratio_off <- tibble(
  ratio = numeric(), empirical_cv = numeric(),
  predicted_cv = numeric(), rel_error = numeric(),
  mean_x = numeric(), scenario = character()
)

ratio_on <- tibble(
  ratio = numeric(), empirical_cv = numeric(),
  predicted_cv = numeric(), rel_error = numeric(),
  mean_x = numeric(), scenario = character()
)

# ============================================================
# RUN 1: tablet error OFF
# ============================================================

cat("=== TABLET ERROR OFF ===\n")
set.seed(42)

for (r in ratios) {
  
  n_targets <- r * n_markers
  true_c    <- n_targets / V
  
  cat("  Ratio =", r, "(n_targets =", n_targets, ")...\n")
  
  estimates <- numeric(num_slides)
  x_counts  <- numeric(num_slides)
  
  for (s in 1:num_slides) {
    targets <- tibble(x = runif(n_targets), y = runif(n_targets), type = "target")
    markers <- tibble(x = runif(n_markers), y = runif(n_markers), type = "marker")
    slide <- rbind(targets, markers)
    
    ms <- markers[markers$y <= strip_height, ]
    ms <- ms[order(ms$x), ]
    
    if (nrow(ms) < target_n_markers) {
      estimates[s] <- NA; x_counts[s] <- NA; next
    }
    
    x_end <- ms$x[target_n_markers]
    in_w <- slide$y <= strip_height & slide$x <= x_end
    x_count <- sum(slide$type == "target" & in_w)
    n_count <- sum(slide$type == "marker" & in_w)
    
    estimates[s] <- (x_count * marker_label) / (n_count * V)
    x_counts[s]  <- x_count
  }
  
  emp_mean <- mean(estimates, na.rm = TRUE)
  emp_sd   <- sd(estimates, na.rm = TRUE)
  emp_cv   <- 100 * emp_sd / emp_mean
  
  mean_x <- mean(x_counts, na.rm = TRUE)
  s1P    <- 0
  pred_cv <- 100 * sqrt(
    (s1P / sqrt(1))^2 +
      (sqrt(mean_x) / mean_x)^2 +
      (sqrt(target_n_markers) / target_n_markers)^2
  )
  
  rel_err <- 100 * (emp_mean - true_c) / true_c
  
  ratio_off <- rbind(ratio_off, tibble(
    ratio = r, empirical_cv = round(emp_cv, 1),
    predicted_cv = round(pred_cv, 1), rel_error = round(rel_err, 3),
    mean_x = round(mean_x, 1), scenario = "Without tablet error"
  ))
}

# ============================================================
# RUN 2: tablet error ON
# ============================================================

cat("\n=== TABLET ERROR ON ===\n")
set.seed(42)

for (r in ratios) {
  
  n_targets <- r * n_markers
  true_c    <- n_targets / V
  
  cat("  Ratio =", r, "(n_targets =", n_targets, ")...\n")
  
  estimates <- numeric(num_slides)
  x_counts  <- numeric(num_slides)
  
  for (s in 1:num_slides) {
    
    # Tablet error: actual markers varies
    n_markers_actual <- round(rnorm(1, mean = marker_label, sd = marker_sd))
    if (n_markers_actual < 1) n_markers_actual <- 1
    
    targets <- tibble(x = runif(n_targets), y = runif(n_targets), type = "target")
    markers <- tibble(x = runif(n_markers_actual), y = runif(n_markers_actual), type = "marker")
    slide <- rbind(targets, markers)
    
    ms <- markers[markers$y <= strip_height, ]
    ms <- ms[order(ms$x), ]
    
    if (nrow(ms) < target_n_markers) {
      estimates[s] <- NA; x_counts[s] <- NA; next
    }
    
    x_end <- ms$x[target_n_markers]
    in_w <- slide$y <= strip_height & slide$x <= x_end
    x_count <- sum(slide$type == "target" & in_w)
    n_count <- sum(slide$type == "marker" & in_w)
    
    # Formula uses LABEL value
    estimates[s] <- (x_count * marker_label) / (n_count * V)
    x_counts[s]  <- x_count
  }
  
  emp_mean <- mean(estimates, na.rm = TRUE)
  emp_sd   <- sd(estimates, na.rm = TRUE)
  emp_cv   <- 100 * emp_sd / emp_mean
  
  mean_x <- mean(x_counts, na.rm = TRUE)
  s1P    <- tablet_cv
  pred_cv <- 100 * sqrt(
    (s1P / sqrt(1))^2 +
      (sqrt(mean_x) / mean_x)^2 +
      (sqrt(target_n_markers) / target_n_markers)^2
  )
  
  rel_err <- 100 * (emp_mean - true_c) / true_c
  
  ratio_on <- rbind(ratio_on, tibble(
    ratio = r, empirical_cv = round(emp_cv, 1),
    predicted_cv = round(pred_cv, 1), rel_error = round(rel_err, 3),
    mean_x = round(mean_x, 1), scenario = "With tablet error"
  ))
}

# ============================================================
# COMPARISON TABLE
# ============================================================

ratio_comparison <- tibble(
  ratio         = ratio_off$ratio,
  cv_off        = ratio_off$empirical_cv,
  cv_on         = ratio_on$empirical_cv,
  cv_difference = round(ratio_on$empirical_cv - ratio_off$empirical_cv, 1),
  pred_off      = ratio_off$predicted_cv,
  pred_on       = ratio_on$predicted_cv,
  mean_x_off    = ratio_off$mean_x
)

print(ratio_comparison)

# ============================================================
# PLOT: precision vs ratio, with and without tablet error
# ============================================================

both_ratio <- rbind(ratio_off, ratio_on)

ggplot(both_ratio, aes(x = ratio, colour = scenario)) +
  geom_line(aes(y = empirical_cv, linetype = "Empirical"), linewidth = 1.2) +
  geom_point(aes(y = empirical_cv), size = 3) +
  geom_line(aes(y = predicted_cv, linetype = "Predicted"), linewidth = 0.8) +
  geom_point(aes(y = predicted_cv), size = 2, shape = 1) +
  scale_colour_manual(values = c("steelblue", "orange")) +
  scale_linetype_manual(values = c(Empirical = "solid", Predicted = "dashed")) +
  scale_x_log10(breaks = ratios) +
  labs(title    = "Precision vs target-to-marker ratio: effect of tablet variability",
       subtitle = "Solid = empirical, Dashed = Mays formula prediction",
       x = "Target-to-marker ratio (log scale)",
       y = "Total error (%)",
       colour = NULL, linetype = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")


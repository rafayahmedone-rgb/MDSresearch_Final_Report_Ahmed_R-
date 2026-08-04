

library(tibble)
library(ggplot2)

n <- 10   

x <- runif(n, min = 0, max = 1)
y <- runif(n, min = 0, max = 1)

points <- tibble(x = x, y = y)

print(points)

ggplot(points, aes(x = x, y = y)) +
  geom_point() 


# distances 
# SHOWS a matrix form ??

distances <- dist(points)

distances
#convert to vector 
dist_vector <- as.vector(distances)
dist_vector


num_simulations <- 5

all_simulations_loop <- list()

for (i in 1:num_simulations) {
  all_simulations_loop[[i]] <- runif(n, min = 0, max = 1)}



num_simulations <- 5


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



# adding list for distance
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
# Compute nearest-neighbour distances for each slide


# Manual version (the by-hand approach, for demonstration)
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
#pts <- tibble(x = runif(10), y = runif(10))

# Get nearest-neighbour distances
#nn_spatstat <- nndist(pts$x, pts$y)

#######################

#  Run many simulations and pool NN distances


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



# Diggle's theoretical CDF for NN distances (Eqn 2.8)


# G(y) = 1 - exp(-lambda * pi * y^2)
# For a unit square with n points, lambda = n / 1 = n
diggle_nn_cdf <- function(y, n) {
  1 - exp(-n * pi * y^2)
}

# Quick sanity checks
diggle_nn_cdf(0, n_points)       # should be 0
diggle_nn_cdf(0.2, n_points)     # should be close to 1 (for n = 100)
# --- Theoretical density (derivative of G) ---



# Build empirical CDF and compare


nn_comparison <- tibble(
  distance    = sort(pooled_nn),
  empirical   = (1:length(pooled_nn)) / length(pooled_nn)
)
nn_comparison$theoretical <- diggle_nn_cdf(nn_comparison$distance, n_points)



# Plot empirical vs theoretical


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

#  Manual nearest-neighbour distance calculation 
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
#  Generate one slide with both targets and markers

n_targets <- 200   #  targets 
n_markers <- 100    #  markers

targets <- tibble(x = runif(n_targets, 0, 1),
                  y = runif(n_targets, 0, 1),
                  type = "target")

markers <- tibble(x = runif(n_markers, 0, 1),
                  y = runif(n_markers, 0, 1),
                  type = "marker")

slide <- rbind(targets, markers)

#  Plot both together 
ggplot(slide, aes(x, y, colour = type, shape = type, size = type)) +
  geom_point() +
  scale_colour_manual(values = c(target = "grey50", marker = "orange")) +
  scale_shape_manual(values = c(target = 16, marker = 18)) +   # shape 
  scale_size_manual(values = c(target = 1.5, marker = 3)) +
  coord_equal() +
  theme_minimal() +
  labs(title = "A simulated slide with targets and markers",
       x = "x", y = "y")

###########################

#ALTERNATIVE METHOD IMPLementation
library(tibble)
library(ggplot2)

#  a slide with targets and markers ---
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

#  fix a dimension of the window 
strip_height <- 0.05

#  the stopping condition 
target_n_markers <- 10

#the window grows from x = 0 rightward 
#    We find the x-coordinate where the 10th marker sits within the strip.
#    That's the rightmost edge of the final window.

# Keep only markers inside the strip height
markers_in_strip <- markers[markers$y <= strip_height, ]    #Filter the particles to the strip height

# Sort them left to right
markers_in_strip <- markers_in_strip[order(markers_in_strip$x), ] # sort the filtered particles 

#  handle "not enough markers" case 
if (nrow(markers_in_strip) < target_n_markers) {                                 # if for eg count is 7 and we needed to stop at 10 not enough so change parameters 
  stop("Strip did not contain enough markers to reach the stopping threshold. ",
       "Increase n_markers or strip_height.")
}

# The x-coordinate of the 10th marker is where the window stops growing
x_end <- markers_in_strip$x[target_n_markers]      

#  count everything inside the final window 
in_window <- slide$y <= strip_height & slide$x <= x_end    # particles where the slide is less than strip and the slide is less than the strip end 

x_count <- sum(slide$type == "target" & in_window)
n_count <- sum(slide$type == "marker" & in_window)

# Sanity check — n_count should exactly equal the stopping threshold
stopifnot(n_count == target_n_markers)

#  plug into the concentration formula 
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


#####################

# THE LOOP 


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


# AFTER THE LOOP — compute both errors and compare

# true concentration 
true_concentration <- n_targets / V   # 2000

# The ACTUAL error (measured from the simulation)
empirical_mean <- mean(estimates, na.rm = TRUE)
empirical_sd   <- sd(estimates, na.rm = TRUE)
empirical_cv   <- 100 * empirical_sd / empirical_mean

# The PREDICTED error (from the Mays formula)
mean_x <- mean(x_counts, na.rm = TRUE)   # average fossils counted
mean_n <- target_n_markers                # always 10
s1P    <- 0                               # tablet error OFF for now

T_term <- s1P / sqrt(N1)                  # tablet variability
F_term <- sqrt(mean_x) / mean_x           # fossil counting noise
M_term <- sqrt(mean_n) / mean_n           # marker counting noise

predicted_cv <- 100 * sqrt(T_term^2 + F_term^2 + M_term^2)

# Accuracy check: is the mean close to the truth? 
relative_error <- 100 * (empirical_mean - true_concentration) / true_concentration

# Put everything in a tibble 
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

#  + Jensen's comparison


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


# ANALYSIS


true_concentration <- n_targets / V

# Two ways to compute the mean  
naive_mean  <- mean(estimates, na.rm = TRUE)
jensen_mean <- (mean(x_counts, na.rm = TRUE) * N1 * Y1) /       # jensens 
  (mean(n_counts, na.rm = TRUE) * V)

#  Accuracy 
relative_error_naive  <- 100 * (naive_mean - true_concentration) / true_concentration
relative_error_jensen <- 100 * (jensen_mean - true_concentration) / true_concentration

# Precision: empirical vs predicted 
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

#  True SE (Anthony's method: from deviations around truth) 
true_se <- sqrt(sum((estimates - true_concentration)^2, na.rm = TRUE) / 
                  (sum(!is.na(estimates)) - 1))
true_cv <- 100 * true_se / true_concentration

# Results 
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









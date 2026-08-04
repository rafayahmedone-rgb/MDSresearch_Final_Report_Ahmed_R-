# All plots 
grid_on       <- read.csv("./data/grid_on_final_2.csv")
comp_baseline <- read.csv("./data/comparison_baseline_1.csv")
comp_threshold <- read.csv("./data/comparison_threshold_10to800.csv")
comp_ratio    <- read.csv("./data/comparison_ratio_2.csv")
sim_est        <- read.csv("./data/sim_estimates_2.csv")
sweep_on      <- read.csv("./data/sweep_on_10to800.csv")
ied_df    <- read.csv("./data/ied_df.csv")
pooled_nn <- read.csv("./data/pooled_nn.csv")
theory_df <- read.csv("./data/theory_df.csv")

#Linear count visualisation
n_t2 <- 2000; n_m2 <- 1000; sh <- 0.05
tg2 <- tibble(x = runif(n_t2), y = runif(n_t2), type = "Target")
mk2 <- tibble(x = runif(n_m2), y = runif(n_m2), type = "Marker")
sl2 <- rbind(tg2, mk2)
ms2 <- mk2[mk2$y <= sh, ]; ms2 <- ms2[order(ms2$x), ]

snap <- function(k) {
  xe <- ms2$x[k]
  xc <- sum(tg2$y <= sh & tg2$x <= xe)
  ggplot(sl2, aes(x, y, colour = type, shape = type, size = type)) +
    geom_point(alpha = 0.5) +
    scale_colour_manual(values = c(Target = "grey60", Marker = col_on), guide = "none") +
    scale_shape_manual(values  = c(Target = 16, Marker = 18), guide = "none") +
    scale_size_manual(values   = c(Target = 0.8, Marker = 2), guide = "none") +
    annotate("rect", xmin = 0, xmax = xe, ymin = 0, ymax = sh,
             colour = "red", fill = "red", alpha = 0.15, linewidth = 0.6) +
    coord_cartesian(xlim = c(0, 0.25), ylim = c(0, 0.12)) +
    theme_minimal(base_size = 9) +
    labs(title    = sprintf("After %d markers", k),
         subtitle = sprintf("fossils = %d", xc), x = "x", y = "y")
}
(snap(3) + snap(6)) / (snap(9) + snap(10)) +
  plot_annotation(subtitle = "Red = counted region; grey = targets; orange = markers")




###############
#The virtual slide

sl <- rbind(
  tibble(x = runif(200), y = runif(200), type = "Target"),
  tibble(x = runif(100), y = runif(100), type = "Marker"))
ggplot(sl, aes(x, y, colour = type, shape = type, size = type)) +
  geom_point() +
  scale_colour_manual(values = c(Target = "grey50", Marker = col_on)) +
  scale_shape_manual(values  = c(Target = 16, Marker = 18)) +
  scale_size_manual(values   = c(Target = 1.5, Marker = 2.5)) +
  coord_equal() +
  theme_minimal(base_size = 11) +
  labs(x = "x", y = "y", colour = NULL, shape = NULL, size = NULL) +
  theme(legend.position = "bottom")



###########################
# Poisson Verification 

p1 <- ggplot(ied_df, aes(x = d)) +
  geom_step(aes(y = emp), colour = col_off, linewidth = 0.9) +
  geom_line(aes(y = th),  colour = "red",   linewidth = 1) +
  labs(title    = "Inter-event distances",
       subtitle = "Blue = empirical CDF, Red = Theoretical",
       x = "Distance t", y = "H(t)") +
  theme_minimal(base_size = 10)

p3 <- ggplot(pooled_nn, aes(x = nn)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 40, fill = "lightblue",
                 colour = "black", linewidth = 0.2) +
  geom_line(data = theory_df, aes(x = y, y = density),
            colour = "red", linewidth = 1) +
  labs(title    = "Nearest-neighbour distance distribution",
       subtitle = "Histogram = empirical, Red = Theoretical Poisson density ",
       x = "Nearest-neighbour distance y", y = "Density") +
  theme_minimal(base_size = 10)

p1 / p3


#########################

# threshold vary 

both_thresh <- bind_rows(
  comp_threshold |> mutate(scenario = "Without tablet error", emp = cv_off, pred = pred_off),
  comp_threshold |> mutate(scenario = "With tablet error",    emp = cv_on,  pred = pred_on)
) |> select(threshold, scenario, emp, pred)

ggplot(both_thresh, aes(x = threshold, colour = scenario)) +
  geom_vline(xintercept = 500, linetype = "dotted", colour = "grey40") +
  geom_line(aes(y = pred), linewidth = 0.9) +
  geom_point(aes(y = emp), size = 2.6) +
  scale_colour_manual(values = c("Without tablet error" = col_off,
                                 "With tablet error"    = col_on)) +
  scale_x_log10(breaks = c(10, 25, 50, 100, 200, 400, 800)) +
  labs(x = "Stopping threshold k (log scale)",
       y = "Error (%)",
       colour = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

#############################

#Ratio vary 

both_ratio <- bind_rows(
  comp_ratio |> mutate(scenario = "Without tablet error", emp = cv_off, pred = pred_off),
  comp_ratio |> mutate(scenario = "With tablet error",    emp = cv_on,  pred = pred_on)
) |> select(ratio, scenario, emp, pred)

ggplot(both_ratio, aes(x = ratio, colour = scenario)) +
  geom_hline(yintercept = 100 / sqrt(100), linetype = "dotted",
             colour = "grey40", linewidth = 0.8) +
  geom_line(aes(y = pred), linewidth = 0.9) +
  geom_point(aes(y = emp), size = 2.6) +
  scale_colour_manual(values = c("Without tablet error" = col_off,
                                 "With tablet error"    = col_on)) +
  scale_x_log10(breaks = c(1, 2, 3, 6, 10, 30, 60)) +
  labs(x = "Target-to-marker ratio (log scale)",
       y = "Error (%)",
       colour = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")



###############################
# grid 

grid_on$cv_label <- sprintf("%.1f", grid_on$empirical_cv)

ggplot(grid_on, aes(x = factor(threshold), y = factor(ratio), fill = empirical_cv)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label  = cv_label,
                colour = ifelse(empirical_cv > 35, "white", "black")),
            size = 3.5) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1,
                       name = "Simulated error\n(%)") +
  scale_colour_identity() +
  labs(x = "Stopping threshold k", y = "Target-to-marker ratio") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), legend.position = "right")
#########################################

#simulation stability 

sim_est <- read.csv("./data/sim_estimates_2.csv")
true_c_conv <- 12489     
est <- sim_est$estimate
sq_err      <- (est - true_c_conv)^2
running_pct <- 100 * sqrt(cumsum(sq_err) / seq_along(sq_err)) / true_c_conv
conv <- tibble(iteration = seq_along(est), rmse_pct = running_pct)

ggplot(conv, aes(x = iteration, y = rmse_pct)) +
  geom_line(colour = col_off, linewidth = 0.8) +
  geom_vline(xintercept = 10000, linetype = "dotted", colour = "grey40") +
  labs(x = "Number of simulated slides", y = "RMSE (% of true value)") +
  theme_minimal(base_size = 11)

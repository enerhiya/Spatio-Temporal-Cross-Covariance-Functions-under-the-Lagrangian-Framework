directory <- '/home/salvanmo/Desktop/'

root <- paste(directory, 'Spatio-Temporal-Cross-Covariance-Functions-under-the-Lagrangian-Framework/', sep = '')

source(file = paste(root, "Functions/load_packages.R", sep = ''))
source(file = paste(root, "Functions/cov_func.R", sep = ''))

N <- 30
t <- 3
grid_x <- seq(from = 0, to = 1, length.out = N)
grid_y <- seq(from = 0, to = 1, length.out = N)
sim_grid_locations <- expand.grid(grid_x, grid_y) %>% as.matrix()

wind <- c(0.1001, 0.1001)
theta <- c(1, 1, 0.23, 0.5, 1, 0.5)

cov1 <- matern_cov(theta, wind, max_time_lag = t - 1, LOCS = sim_grid_locations)

# call function to plot the covariance

jpeg(file = paste(root, 'Figures/heatmap_covariance.jpg', sep = ''), width = 800, height = 800)

image.plot(cov1)

dev.off()

set.seed(1234)
r1 <- rmvn(1, rep(0, ncol(cov1)), cov1, ncores = 25)


directory <- '/home/salvanmo/Desktop/'

root <- paste(directory, 'Spatio-Temporal-Cross-Covariance-Functions-under-the-Lagrangian-Framework/', sep = '')

source(file = paste(root, "Functions/load_packages.R", sep = ''))
source(file = paste(root, "Functions/cov_func.R", sep = ''))

n <- 121
N <- sqrt(n)
t <- 3
grid_x <- seq(from = 0, to = 1, length.out = N)
grid_y <- seq(from = 0, to = 1, length.out = N)
sim_grid_locations <- expand.grid(grid_x, grid_y) %>% as.matrix()

wind <- c(0.1001, 0.1001)
theta <- c(1, 1, 0.23, 0.5, 1)

cov1 <- matern_cov(theta, wind, max_time_lag = t - 1, LOCS = sim_grid_locations)

set.seed(1234)
r1 <- mvrnorm(1, mu = rep(0, ncol(cov1)), Sigma = cov1)

root <- '/home/salvanmo/Desktop/Spatio-Temporal-Cross-Covariance-Functions-under-the-Lagrangian-Framework/'

n <- 121
N <- sqrt(n)
t <- 3
grid_x <- seq(from = 0, to = 1, length.out = N)
grid_y <- seq(from = 0, to = 1, length.out = N)
sim_grid_locations <- expand.grid(grid_x, grid_y) %>% as.matrix()

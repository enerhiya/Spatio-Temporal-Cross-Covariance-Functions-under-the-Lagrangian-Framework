directory <- '/home/salvanmo/Desktop/'

root <- paste(directory, 'Spatio-Temporal-Cross-Covariance-Functions-under-the-Lagrangian-Framework/', sep = '')

source(file = paste(root, "Functions/load_packages.R", sep = ''))
source(file = paste(root, "Functions/auxiliary_functions.R", sep = ''))
source(file = paste(root, "Functions/cov_func.R", sep = ''))

load_default_data = F

if(load_default_data){

	r1 <- read.table(paste(root, "Data/frozen_matern_realizations", sep = ''), header = FALSE, sep = ",")

}else{

	N <- 50
	n <- N^2
	t <- 5
	grid_x <- seq(from = 0, to = 1, length.out = N)
	grid_y <- seq(from = 0, to = 1, length.out = N)
	sim_grid_locations <- expand.grid(grid_x, grid_y) %>% as.matrix()

	wind <- c(0.1001, 0.1001)
	theta <- c(1, 1, 0.23, 0.5, 1, 0.8)

	cov1 <- frozen_matern_cov(theta, wind, max_time_lag = t - 1, LOCS = sim_grid_locations)

	write_to_txt(cov1, paste(root, "Data/frozen_matern_cov", sep = ''))

	# call function to plot the covariance

	jpeg(file = paste(root, 'Figures/heatmap_covariance.jpg', sep = ''), width = 800, height = 800)

	image.plot(cov1)

	dev.off()

	set.seed(1234)
	r1 <- rmvn(1, rep(0, ncol(cov1)), cov1, ncores = 25)

	write_to_txt(cbind(r1[1, 1:(n * t)], r1[1, n * t + 1:(n * t)]), paste(root, "Data/frozen_matern_realizations", sep = ''))
}

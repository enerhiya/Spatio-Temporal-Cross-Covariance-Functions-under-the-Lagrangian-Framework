directory <- '/home/salvanmo/Desktop/'

root <- paste(directory, 'Spatio-Temporal-Cross-Covariance-Functions-under-the-Lagrangian-Framework/', sep = '')

source(file = paste(root, "Functions/load_packages.R", sep = ''))
source(file = paste(root, "Functions/auxiliary_functions.R", sep = ''))
source(file = paste(root, "Functions/cov_func.R", sep = ''))
source(file = paste(root, "Functions/plot_func.R", sep = ''))

load_default_data = F

if(load_default_data){

	cov1 <- read.table(paste(root, "Data/frozen_matern_cov_for_heatmap", sep = ''), header = FALSE, sep = " ") %>% as.matrix() 
	r1 <- read.table(paste(root, "Data/frozen_matern_realizations", sep = ''), header = FALSE, sep = ",") %>% as.matrix() 
	locs <- read.table(paste(root, "Data/simulation_locs", sep = ''), header = FALSE, sep = " ") %>% as.matrix() 
	TT <- max(unique(locs[, 3]))
	n <- nrow(locs)/TT	

	covariances <- list()

	covariances[[1]] <- cov1
	covariances[[2]] <- cov1
	covariances[[3]] <- cov1
	covariances[[4]] <- cov1

	realizations <- list()

	realizations[[1]] <- r1
	realizations[[2]] <- r1
	realizations[[3]] <- r1
	realizations[[4]] <- r1

	N <- 51
        n <- N^2
        TT <- 3
        grid_x <- seq(from = -0.5, to = 0.5, length.out = N)
        sim_grid_locations <- expand.grid(grid_x, grid_x) %>% as.matrix()
	
	fig1(COVARIANCES = covariances, file_name = paste(root, 'Figures/simulation_covariances.pdf', sep = ''), LOCS = sim_grid_locations)
	fig2(REALIZATIONS = realizations, file_name = paste(root, 'Figures/simulation_realizations.pdf', sep = ''), LOCS = locs[1:n, 1:2])

}else{

	N <- 50
	n <- N^2
	t <- 5
	grid_x <- seq(from = 0, to = 1, length.out = N)
	grid_y <- seq(from = 0, to = 1, length.out = N)
	sim_grid_locations <- expand.grid(grid_x, grid_y) %>% as.matrix()

	write_to_txt(rbind(cbind(sim_grid_locations, 1), cbind(sim_grid_locations, 2), cbind(sim_grid_locations, 3), cbind(sim_grid_locations, 4), cbind(sim_grid_locations, 5)), paste(root, "Data/simulation_locs", sep = ''))

	wind <- c(0.1001, 0.1001)
	theta <- c(1, 1, 0.23, 0.5, 1, 0.8)

	cov1 <- frozen_matern_cov(theta, wind, max_time_lag = t - 1, LOCS = sim_grid_locations)
	cov1_for_heatmap <- frozen_matern_cov_for_heatmap(theta, wind)

	write_to_txt(cov1, paste(root, "Data/frozen_matern_cov", sep = ''))
	write_to_txt(cov1_for_heatmap, paste(root, "Data/frozen_matern_cov_for_heatmap", sep = ''))

	# call function to plot the covariance

	jpeg(file = paste(root, 'Figures/heatmap_covariance.jpg', sep = ''), width = 800, height = 800)

	image.plot(cov1)

	dev.off()

	set.seed(1234)
	r1 <- rmvn(1, rep(0, ncol(cov1)), cov1, ncores = 25)

	write_to_txt(cbind(r1[1, 1:(n * t)], r1[1, n * t + 1:(n * t)]), paste(root, "Data/frozen_matern_realizations", sep = ''))
}

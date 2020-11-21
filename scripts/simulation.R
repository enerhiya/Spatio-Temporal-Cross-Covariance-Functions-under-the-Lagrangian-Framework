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

	N <- 30
	n <- N^2
	TT <- 3
	grid_x <- seq(from = 0, to = 1, length.out = N)
	grid_y <- seq(from = 0, to = 1, length.out = N)
	sim_grid_locations <- expand.grid(grid_x, grid_y) %>% as.matrix()

	locs <- cbind(sim_grid_locations, 1)

	for(tt in 2:TT){

		locs <- rbind(locs, cbind(sim_grid_locations, tt))

	}

	#write_to_txt(locs, paste(root, "Data/simulation_locs", sep = ''))


	#-------------------------------     SIMULATE FROM FROZEN MODEL     -------------------------------# 

	wind <- c(0.1001, 0.1001)
	theta <- c(1, 1, 2, 1, 1, 0.8)

	cov1 <- frozen_matern_cov(theta, wind, max_time_lag = TT - 1, LOCS = sim_grid_locations)
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

	#-------------------------------     SIMULATE FROM NONFROZEN MODEL VIA MONTECARLO SIMULATION    -------------------------------# 



	#-------------------------------     SIMULATE UNIVARIATE FROM NONFROZEN VELOCITY FIELD MODEL     -------------------------------# 

	w_cov <- frozen_matern_cov(theta = c(0.1, 0.1, 0.03, 0.5, 1, 0.8), max_time_lag = 0, LOCS = sim_grid_locations)

	#~~~~~~~~~~~~~~~~~~~~~~~~~~     SIMULATE SAMPLE REALIZATIONS FOR THE VELOCITY FIELD     ~~~~~~~~~~~~~~~~~~~~~~~~~~#

	set.seed(1235)
	w <- rmvn(1, rep(0, ncol(w_cov)), w_cov, ncores = 20)

	jpeg(file = paste(root, 'Figures/velocity_field.jpg', sep = ''), width = 1000, height = 1000)

	split.screen( rbind(c(0.1,0.91,0.1,0.95), c(0.92,0.99,0.1,0.95)))

	screen(1)
	par(pty = 's')
	par(mai=c(0.3, 0.3, 0.3, 0.3))

	plot(sim_grid_locations, ylab = '', xlab = '', cex.lab = 4, cex.axis = 2, yaxt = 'n')
	arrows(sim_grid_locations[, 1], sim_grid_locations[, 2], sim_grid_locations[, 1] + 0.1 * w[1:n], sim_grid_locations[, 2] + 0.1 * w[n + 1:n], length = 0.1)
	mtext('Longitude', side = 1, line = 3.5, cex = 3, font = 2)
	mtext('Latitude', side = 2, line = 3.5, cex = 3, font = 2)

	mtext('Advection Velocity Vectors', side = 3, line = 2, adj = 0.5, cex = 4, font = 2)

	close.screen( all=TRUE)
	dev.off()

	#~~~~~~~~~~~~~~~~~~~~~~~~~~     SIMULATE RANDOM FIELD GIVEN THE VELOCITY FIELD     ~~~~~~~~~~~~~~~~~~~~~~~~~~#

	Sigma <- matrix(0, ncol = n * TT, nrow = n * TT)

	for(sim in 1:100){
		cat(sim, '\n')
		w <- rmvn(1, rep(0, ncol(w_cov)), w_cov, ncores = 20)
		locs <- sim_grid_locations 
		for(tt in 1:(TT - 1)){
			locs <- rbind(locs, cbind(sim_grid_locations[, 1] - 0.1 * w[1:n] * tt, sim_grid_locations[, 2] - 0.1 * w[n + 1:n] * tt))
		}	
		dist0 <- dist(locs, diag = TRUE, upper = TRUE) %>% as.matrix()

		Sigma <- Sigma + Matern(dist0, range = 0.23, nu = 1)

	}

	Sigma <- Sigma / 100

	r1 <- rmvn(1, rep(0, n * TT), Sigma, ncores = 20)

	jpeg(file = paste(root, 'Figures/sim_uni.jpg', sep = ''), width = 1800, height = 700)
	par(mfrow = c(1,3))
	for(tt in 1:3){
		par(pty = 's')
		quilt.plot(sim_grid_locations[, 1], sim_grid_locations[, 2], r1[1, (tt - 1) * n + 1:n], nx = 25, ny = 25, ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', cex.lab = 4, add.legend = F, cex.axis = 2)
	}
	dev.off()

}

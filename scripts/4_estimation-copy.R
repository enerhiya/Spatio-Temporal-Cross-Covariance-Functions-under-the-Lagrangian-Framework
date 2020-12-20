
workstation = T

if(workstation) directory <- '/home/salvanmo/Desktop/' 		else 		directory <- '/ibex/scratch/salvanmo/'

root <- paste(directory, 'Spatio-Temporal-Cross-Covariance-Functions-under-the-Lagrangian-Framework/', sep = '')

source(file = paste(root, "Functions/load_packages.R", sep = ''))
source(file = paste(root, "Functions/auxiliary_functions.R", sep = ''))

yr <- 2015

Z1 <- read.table(paste(root, 'Data/ncdf/layer1_residuals_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
Z2 <- read.table(paste(root, 'Data/ncdf/layer2_residuals_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
locs <- LOCS <- read.table(paste(root, 'Data/ncdf/LOCS-3D-dataset', sep = ''), header = FALSE, sep = " ") %>% as.matrix()

locs_cartesian <- cbind(6371 * cos(locs[, 2] / 180.0 * pi) * cos(locs[, 1] / 180.0 * pi), 6371 * cos(locs[, 2] / 180.0 * pi) * sin(locs[, 1] / 180.0 * pi))
locs_demean <- cbind(locs_cartesian[, 1] - mean(locs_cartesian[, 1]), locs_cartesian[, 2] - mean(locs_cartesian[, 2]))
std_locs <- cbind((locs_cartesian[, 1] - mean(locs_cartesian[, 1])) / sd(locs_cartesian[, 1]), (locs_cartesian[, 2] - mean(locs_cartesian[, 2])) / sd(locs_cartesian[, 2]))

set.seed(1234)
subs <- sample(1:550, 200)
sim_grid_locations <- locs_demean[subs, ]

n <- nrow(sim_grid_locations)

TT <- 3

Z1 <- Z1[, subs]
Z2 <- Z2[, subs]

############# CHANGE DIRECTORY #############

workstation = T

if(workstation) directory <- '/home/salvanmo/Desktop/'          else            directory <- '/ibex/scratch/salvanmo/'

root <- paste(directory, 'thesis/', sep = '')

source(file = paste(root, "R_codes/Functions/load_packages.R", sep = ''))
source(file = paste(root, "R_codes/Functions/auxiliary_functions.R", sep = ''))
source(file = paste(root, "R_codes/Functions/cov_func.R", sep = ''))

sourceCpp(file=paste(root, "R_codes/Functions/spatially_varying_parameters2.cpp",sep=''))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  FITTING M1  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

NEGLOGLIK_ST <- function(p){

	if(p[6] <= -1 | p[6] >= 1){
		return(Inf)
	}else{

		theta <- c(exp(p[1:5]), p[6])
		mu <- p[7:8]
		wind_var_chol <- matrix(c(p[9], p[10], 0, p[11]), ncol = 2, byrow = T)
		wind_var <- t(wind_var_chol) %*% wind_var_chol

		Sigma <- nonfrozen_matern_cov_multi_small_scale(theta, wind_mu = mu, wind_var = wind_var, max_time_lag = TT - 1, LOCS = sim_grid_locations)
		cholmat <- t(cholesky(Sigma, parallel = TRUE))
		if( length(cholmat) == 0){
			return(Inf)
		}else{
			z <- forwardsolve(cholmat, t(Z_rand_sample))
			logsig  <- 2 * sum(log(diag(cholmat))) * nrow(Z_rand_sample)
			out  <- 1/2 * logsig + 1/2 * sum(z^2)
			
			return(out)
		}
	}

}

init <- c(0, 0, log(100), 0, 0, 0, 0.001, 0.001, 1, 0, 1)

r1 <- c(Z1[1, ], Z1[2, ], Z1[3, ], Z2[1, ], Z2[2, ], Z2[3, ])

Z_rand_sample <- matrix( r1 - mean(r1), nrow = 1)

fit1 <- optim(par = init, fn = NEGLOGLIK_ST, control = list(trace = 5, maxit = 500)) #

for(tt in 1:10){
	write.table(fit1$par, file = paste(root, 'Results/multivariate_stationary_real_data_parameter_estimates_M1-small-scale', sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
	fit1 <- optim(par = fit1$par, fn = NEGLOGLIK_ST, control = list(trace = 5, maxit = 500)) #
}

p <- fit1$par
theta <- c(exp(p[1:5]), p[6])
mu <- p[7:8]
wind_var_chol <- matrix(c(p[9], p[10], 0, p[11]), ncol = 2, byrow = T)
wind_var <- t(wind_var_chol) %*% wind_var_chol

#raw params: -1.2511348 -1.8163147  4.9434041  0.3562039  0.3421842  0.1720987 -0.9255224  5.0388504  2.0005501  2.8191134  2.4391513
#transformed params: 0.2861799   0.1626240 140.2468528   1.4278986   1.4080196   0.1720987 -0.9255224  5.0388504 4.002201  5.639777 13.896859

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  FITTING M2  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

NEGLOGLIK_ST <- function(p){
	cat(p, '\n')

	if(p[6] <= -1 | p[6] >= 1){
		return(Inf)
	}else{

		theta <- c(exp(p[1:5]), p[6])
		mu1 <- p[7:8]
		mu2 <- p[9:10]
		wind_var_chol <- matrix(c(p[11], p[12], p[13], p[14], 0, p[15], p[16], p[17], 0, 0, p[18], p[19], 0, 0, 0, p[20]), ncol = 4, byrow = T)
		wind_var <- t(wind_var_chol) %*% wind_var_chol

		Sigma <- nonfrozen_matern_cov_multi_advec_small_scale(theta, wind_mu1 = mu1, wind_mu2 = mu2, wind_var1 = wind_var[1:2, 1:2], wind_var2 = wind_var[3:4, 3:4], wind_var12 = wind_var[1:2, 3:4], max_time_lag = TT - 1, LOCS = sim_grid_locations)	

		cholmat <- t(cholesky(Sigma, parallel = TRUE))
		if( length(cholmat) == 0){
			return(Inf)
		}else{
			z <- forwardsolve(cholmat, t(Z_rand_sample))
			logsig  <- 2 * sum(log(diag(cholmat))) * nrow(Z_rand_sample)
			out  <- 1/2 * logsig + 1/2 * sum(z^2)
			
			return(out)
		}
	}

}

init <- c(0, 0, log(100), 0, 0, 0, 0.001, 0.001, 0.001, 0.001, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1)

r1 <- c(Z1[1, ], Z1[2, ], Z1[3, ], Z2[1, ], Z2[2, ], Z2[3, ])

Z_rand_sample <- matrix( r1 - mean(r1), nrow = 1)

fit1 <- optim(par = init, fn = NEGLOGLIK_ST, control = list(trace = 5, maxit = 500)) #

for(tt in 1:5){
	write.table(fit1$par, file = paste(root, 'Results/multivariate_stationary_real_data_parameter_estimates_M2', sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
	fit1 <- optim(par = fit1$par, fn = NEGLOGLIK_ST, control = list(trace = 5, maxit = 500)) #
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  FITTING M4  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

NEGLOGLIK_ST <- function(p){

	theta <- exp(p[1:4])
	mu1 <- p[5:6]
	mu2 <- p[7:8]
	wind_var_chol <- matrix(c(p[9], p[10], 0, p[11]), ncol = 2, byrow = T)
	wind_var1 <- t(wind_var_chol) %*% wind_var_chol
	wind_var_chol <- matrix(c(p[12], p[13], 0, p[14]), ncol = 2, byrow = T)
	wind_var2 <- t(wind_var_chol) %*% wind_var_chol

	coef_chol <- matrix(c(p[15], p[16], 0, p[17]), ncol = 2, byrow = T)
	A1 <- t(coef_chol) %*% coef_chol
	coef_chol <- matrix(c(p[18], p[19], 0, p[20]), ncol = 2, byrow = T)
	A2 <- t(coef_chol) %*% coef_chol

	Sigma <- nonfrozen_lmc_cov_small_scale(theta = c(1, 1, theta), wind_mu1 = mu1, wind_mu2 = mu2, wind_var1 = wind_var1, wind_var2 = wind_var2, coef_mat1 = A1, coef_mat2 = A2, max_time_lag = TT - 1, LOCS = sim_grid_locations)

	cholmat <- t(cholesky(Sigma, parallel = TRUE))
	if( length(cholmat) == 0){
		return(Inf)
	}else{
		z <- forwardsolve(cholmat, t(Z_rand_sample))
		logsig  <- 2 * sum(log(diag(cholmat))) * nrow(Z_rand_sample)
		out  <- 1/2 * logsig + 1/2 * sum(z^2)
		
		return(out)
	}
}

init <- c(0, 0, 0, 0, 0.001, 0.001, 0.001, 0.001, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1)

r1 <- c(Z1[1, ], Z1[2, ], Z2[1, ], Z2[2, ])

Z_rand_sample <- matrix( r1 - mean(r1), nrow = 1)

fit1 <- optim(par = init, fn = NEGLOGLIK_ST, control = list(trace = 5, maxit = 500)) #

for(tt in 1:5){
	write.table(fit1$par, file = paste(root, 'Results/multivariate_stationary_real_data_parameter_estimates_M4', sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
	fit1 <- optim(par = fit1$par, fn = NEGLOGLIK_ST, control = list(trace = 5, maxit = 500)) #
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  FITTING M5  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

locs_s <- locs_t <- NULL

for(tt in 1:TT){

	locs_s <- rbind(locs_s, locs_demean)
	locs_t <- c(locs_t, rep(tt, n))

}

dist_s <- as.matrix(dist(locs_s, diag = TRUE, upper = TRUE))	
dist_t <- as.matrix(dist(locs_t, diag = TRUE, upper = TRUE))	

NEGLOGLIK_ST <- function(p){

	if(p[6] <= -1 | p[6] >= 1 | p[9] < 0 | p[9] > 1){
		return(Inf)
	}else{

		theta <- c(exp[1:5], p[6])

		dist_temp1 <- (exp(p[7]) * abs(dist_t)^(2 * exp(p[8])) + 1)
		DIST <- dist_s / (dist_temp1)^(p[9] / 2)

		#Sigma <- exp(p[1]) * exp(- DIST / exp(p[2])) / dist_temp1      + 0.00001 * diag(TT * n)

		Sigma <- purely_spatial_parsimonious_matern(theta, DIST)  / dist_temp1 

		cholmat <- tryCatch(t(cholesky(Sigma, parallel = TRUE)) , error = function(a) numeric(0) )

		if( length(cholmat) == 0){
			return(Inf)
		}else{
			z <- forwardsolve(cholmat, t(Z_rand_sample))
			logsig  <- 2 * sum(log(diag(cholmat))) * nrow(Z_rand_sample)
			out  <- 1/2 * logsig + 1/2 * sum(z^2)
			return(out)
		}
	}
}


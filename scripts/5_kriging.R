
workstation = T

if(workstation) directory <- '/home/salvanmo/Desktop/' 		else 		directory <- '/ibex/scratch/salvanmo/'

root <- paste(directory, 'Spatio-Temporal-Cross-Covariance-Functions-under-the-Lagrangian-Framework/', sep = '')

source(file = paste(root, "Functions/load_packages.R", sep = ''))
source(file = paste(root, "Functions/auxiliary_functions.R", sep = ''))

Z1 <- read.table(paste(root, 'Data/ncdf/layer1_residuals_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
Z2 <- read.table(paste(root, 'Data/ncdf/layer2_residuals_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
locs <- LOCS <- read.table(paste(root, 'Data/ncdf/LOCS-3D-dataset', sep = ''), header = FALSE, sep = " ") %>% as.matrix()

locs_cartesian <- cbind(6371 * cos(locs[, 2] / 180.0 * pi) * cos(locs[, 1] / 180.0 * pi), 6371 * cos(locs[, 2] / 180.0 * pi) * sin(locs[, 1] / 180.0 * pi))
locs_demean <- cbind(locs_cartesian[, 1] - mean(locs_cartesian[, 1]), locs_cartesian[, 2] - mean(locs_cartesian[, 2]))
std_locs <- cbind((locs_cartesian[, 1] - mean(locs_cartesian[, 1])) / sd(locs_cartesian[, 1]), (locs_cartesian[, 2] - mean(locs_cartesian[, 2])) / sd(locs_cartesian[, 2]))

n <- nrow(locs)

sim_grid_locations <- locs_demean
TT <- 5

############# CHANGE DIRECTORY #############

workstation = T

if(workstation) directory <- '/home/salvanmo/Desktop/'          else            directory <- '/ibex/scratch/salvanmo/'

root <- paste(directory, 'thesis/', sep = '')

source(file = paste(root, "R_codes/Functions/load_packages.R", sep = ''))
source(file = paste(root, "R_codes/Functions/auxiliary_functions.R", sep = ''))
source(file = paste(root, "R_codes/Functions/cov_func.R", sep = ''))

sourceCpp(file=paste(root, "R_codes/Functions/spatially_varying_parameters2.cpp",sep=''))

r1 <- c(Z1[1, ], Z1[2, ], Z1[3, ], Z2[1, ], Z2[2, ], Z2[3, ])
Z_rand_sample_in <- matrix( r1 - mean(r1), nrow = 1)

r1 <- c(Z1[4, ], Z1[5, ], Z2[4, ], Z2[5, ])
Z_rand_sample_out <- matrix( r1 - mean(r1), nrow = 1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  KRIGING M1  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

p <- read.table(paste(root, 'Results/multivariate_stationary_real_data_parameter_estimates_M1', sep = ''), header = FALSE, sep = " ") %>% as.matrix()

theta <- c(exp(p[1:5]), p[6])
mu <- p[7:8]
wind_var_chol <- matrix(c(p[9], p[10], 0, p[11]), ncol = 2, byrow = T)
wind_var <- t(wind_var_chol) %*% wind_var_chol

Sigma <- nonfrozen_matern_cov_multi_small_scale(theta, wind_mu = mu, wind_var = wind_var, max_time_lag = TT - 1, LOCS = sim_grid_locations)

index_in <- c(1:(n * (TT - 2)), n * TT + 1:(n * (TT - 2)))
index_out <- c(n * (TT - 2) + 1:(n * 2), n * TT + n * (TT - 2) + 1:(n * 2))

C11 <- Sigma[index_in, index_in]
C22 <- Sigma[index_out, index_out]
C12 <- Sigma[index_in, index_out]

Z_pred <- t(C12) %*% solve(C11) %*% t(Z_rand_sample_in)

mean((c(Z_rand_sample_out) - c(Z_pred))^2)
#0.1836779

#write.table(c(Z_pred), file = paste(root, "Results/real-data-predictions-M1", sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  FITTING M2  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

NEGLOGLIK_ST <- function(p){

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

r1 <- c(Z1[1, ], Z1[2, ], Z2[1, ], Z2[2, ])


fit1 <- optim(par = init, fn = NEGLOGLIK_ST, control = list(trace = 5, maxit = 500)) #

for(tt in 1:5){
	write.table(fit1$par, file = paste(root, 'Results/multivariate_stationary_real_data_parameter_estimates_M2', sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
	fit1 <- optim(par = fit1$par, fn = NEGLOGLIK_ST, control = list(trace = 5, maxit = 500)) #
}


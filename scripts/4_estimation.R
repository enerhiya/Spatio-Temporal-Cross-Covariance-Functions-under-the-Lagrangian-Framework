
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
TT <- 3

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

for(tt in 1:5){
	fit1 <- optim(par = fit1$par, fn = NEGLOGLIK_ST, control = list(trace = 5, maxit = 500)) #
}

p <- fit1$par
theta <- c(exp(p[1:5]), p[6])
mu <- p[7:8]
wind_var_chol <- matrix(c(p[9], p[10], 0, p[11]), ncol = 2, byrow = T)
wind_var <- t(wind_var_chol) %*% wind_var_chol

#raw params: -1.2511348 -1.8163147  4.9434041  0.3562039  0.3421842  0.1720987 -0.9255224  5.0388504  2.0005501  2.8191134  2.4391513
#transformed params: 0.2861799   0.1626240 140.2468528   1.4278986   1.4080196   0.1720987 -0.9255224  5.0388504 4.002201  5.639777 13.896859

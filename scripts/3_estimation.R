#fit models

parallel = T

workstation = F

if(workstation) directory <- '/home/salvanmo/Desktop/' 		else 		directory <- '/ibex/scratch/salvanmo/'

root <- paste(directory, 'Spatio-Temporal-Cross-Covariance-Functions-under-the-Lagrangian-Framework/', sep = '')

source(file = paste(root, "Functions/load_packages.R", sep = ''))
source(file = paste(root, "Functions/auxiliary_functions.R", sep = ''))
source(file = paste(root, "Functions/cov_func.R", sep = ''))
source(file = paste(root, "Functions/plot_func.R", sep = ''))

Z1 <- read.table(paste(root, 'Data/ncdf/DUSMASS25_residuals_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
Z2 <- read.table(paste(root, 'Data/ncdf/BCSMASS_residuals_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
locs <- LOCS <- read.table(paste(root, 'Data/ncdf/LOCS', sep = ''), header = FALSE, sep = " ") %>% as.matrix()

locs_cartesian <- cbind(6371 * cos(locs[, 2] / 180.0 * pi) * cos(locs[, 1] / 180.0 * pi), 6371 * cos(locs[, 2] / 180.0 * pi) * sin(locs[, 1] / 180.0 * pi))
locs_demean <- cbind(locs_cartesian[, 1] - mean(locs_cartesian[, 1]), locs_cartesian[, 2] - mean(locs_cartesian[, 2]))
std_locs <- cbind((locs_cartesian[, 1] - mean(locs_cartesian[, 1])) / sd(locs_cartesian[, 1]), (locs_cartesian[, 2] - mean(locs_cartesian[, 2])) / sd(locs_cartesian[, 2]))

n <- nrow(locs)

dist0 <- as.matrix(dist(locs_demean, diag = TRUE, upper = TRUE))	

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  FITTING ST STATIONARY SCHLATHER MODEL  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# --------------------  STEP 1: Fit purely spatial stationary Schlather (squared exponential)  -------------------- #

NEGLOGLIK1 <- function(p){

	theta <- c(exp(p[1:5]), p[6])

	if(theta[6] > 1 | theta[6] < -1){

		return(Inf)

	}else{

		Sigma <- frozen_matern_cov(theta = theta, max_time_lag = 0, LOCS = locs_demean)

		cholmat <- t(cholesky(Sigma, parallel = TRUE))
		out <- 0
		for(ss in 1:744){
			Z_rand_sample <- matrix(c(Z1[ss, ] - mean(Z1[ss, ]), Z2[ss, ] - mean(Z2[ss, ])), nrow = 1)
			z <- forwardsolve(cholmat, t(Z_rand_sample))
			logsig  <- 2 * sum(log(diag(cholmat))) * nrow(Z_rand_sample)
			out  <- out + 1/2 * logsig + 1/2 * sum(z^2)
		}
		return(out)
	
	}
}

init <- c(0, 0, 0, 0, 0, 0.5)
fit1 <- optim(par = init, fn = NEGLOGLIK1, control = list(trace = 5, maxit = 500)) #
fit1 <- optim(par = fit1$par, fn = NEGLOGLIK1, control = list(trace = 5, maxit = 500)) #
theta <- c(exp(fit1$par[1:5]), fit1$par[6])
cat('theta: ', theta, '\n')
#est: 0.6423467  0.1691413 79.4845434  1.3730585  1.3424806  0.1621744

# --------------------  STEP 2: Fit space-time stationary Schlather  -------------------- #

if(!parallel){

	NEGLOGLIK2 <- function(p, space_params){
				
		wind_var_chol <- matrix(c(p[3], p[4], 0, p[5]), ncol = 2, byrow = T)
		wind_var <- t(wind_var_chol) %*% wind_var_chol		
		cat(p[1:2], wind_var[1,1], wind_var[2,2], wind_var[1, 2], '\n')

		Sigma <- matrix(0, ncol = n * TT * 2, nrow = n * TT * 2)

		for(sim in 1:100){
			w <- mvrnorm(1, p[1:2], wind_var)
			Sigma_temp <- frozen_matern_cov(theta = space_params, wind = w, max_time_lag = 1, LOCS = locs_demean)
			Sigma <- Sigma + Sigma_temp
		}

		Sigma <- Sigma / 100

		cholmat <- t(cholesky(Sigma, parallel = TRUE))
		z <- forwardsolve(cholmat, t(Z_rand_sample))
		logsig  <- 2 * sum(log(diag(cholmat))) * nrow(Z_rand_sample)
		out  <- 1/2 * logsig + 1/2 * sum(z^2)

		return(out)
	}

	TT <- 2
	space_params <- theta
	Z1_rand_sample <- matrix(c(Z1[1, ], Z1[2, ]), nrow = 1)
	Z2_rand_sample <- matrix(c(Z2[1, ], Z2[2, ]), nrow = 1)
	Z_rand_sample <- cbind(Z1_rand_sample - mean(Z1_rand_sample), Z2_rand_sample - mean(Z2_rand_sample))

	init <- c(0.1, 0.1, 0.1, 0, 0.1)
	fit2 <- optim(par = init, fn = NEGLOGLIK2, space_params = theta, control = list(trace = 5, maxit = 500)) #
	fit2 <- optim(par = fit2$par, fn = NEGLOGLIK2, space_params = theta, control = list(trace = 5, maxit = 500)) #
	cat('w_mean: ', fit2$par[1:2], '; w_var: ', fit2$par[3:5], '\n')

}else{
	NEGLOGLIK2 <- function(p, space_params){

		clusterExport(cl, c("p", "space_params"), envir = environment())

		output <- foreach(i=1:number_of_cores_to_use, .combine = '+', .packages = "MASS") %dopar% {

			wind_var_chol <- matrix(c(p[3], p[4], 0, p[5]), ncol = 2, byrow = T)
			wind_var <- t(wind_var_chol) %*% wind_var_chol		

			S <- matrix(0, ncol = n * TT * 2, nrow = n * TT * 2)

			for(sim in 1:100){
				w <- mvrnorm(1, p[1:2], wind_var)
				Sigma_temp <- frozen_matern_cov(theta = space_params, wind = w, max_time_lag = 1, LOCS = locs_demean)
				S <- S + Sigma_temp
			}
			return(S)
		}

		Sigma <- output / (number_of_cores_to_use * 100)

		cholmat <- t(cholesky(Sigma, parallel = TRUE))
		z <- forwardsolve(cholmat, t(Z_rand_sample))
		logsig  <- 2 * sum(log(diag(cholmat))) * nrow(Z_rand_sample)
		out  <- 1/2 * logsig + 1/2 * sum(z^2) #+  90 * sum(abs(c(beta1, beta2)))

		return(out)
				
	}

	TT <- 2
	space_params <- theta
	Z1_rand_sample <- matrix(c(Z1[1, ], Z1[2, ]), nrow = 1)
	Z2_rand_sample <- matrix(c(Z2[1, ], Z2[2, ]), nrow = 1)
	Z_rand_sample <- cbind(Z1_rand_sample - mean(Z1_rand_sample), Z2_rand_sample - mean(Z2_rand_sample))

	init <- c(0.1, 0.1, 0.1, 0, 0.1)

	cores=detectCores()
	number_of_cores_to_use = cores[1]-1 # not to overload the computer
	cl <- makeCluster(number_of_cores_to_use) 
	registerDoParallel(cl)

	clusterExport(cl, c("root", "TT", "locs_demean", "n", "theta"), envir = environment())
	clusterEvalQ(cl, source(paste(root, "../Functions/load_packages.R", sep = '')))
	clusterEvalQ(cl, source(paste(root, "../Functions/cov_func.R", sep = '')))

	fit2 <- optim(par = init, fn = NEGLOGLIK2, space_params = theta, control = list(trace = 5, maxit = 500)) #
	fit2 <- optim(par = fit2$par, fn = NEGLOGLIK2, space_params = theta, control = list(trace = 5, maxit = 500)) #
	cat('w_mean: ', fit2$par[1:2], '; w_var: ', fit2$par[3:5], '\n')
	stopCluster(cl)

}




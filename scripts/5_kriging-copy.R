
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

Z1 <- Z1[, subs]
Z2 <- Z2[, subs]

n <- nrow(sim_grid_locations)

TT <- 5

############# CHANGE DIRECTORY #############

workstation = T

if(workstation) directory <- '/home/salvanmo/Desktop/'          else            directory <- '/ibex/scratch/salvanmo/'

root <- paste(directory, 'thesis/', sep = '')

source(file = paste(root, "R_codes/Functions/load_packages.R", sep = ''))
source(file = paste(root, "R_codes/Functions/auxiliary_functions.R", sep = ''))
source(file = paste(root, "R_codes/Functions/cov_func.R", sep = ''))

sourceCpp(file=paste(root, "R_codes/Functions/spatially_varying_parameters2.cpp",sep=''))

r1 <- c(Z1[6, ], Z1[7, ], Z1[8, ], Z1[9, ], Z2[6, ], Z2[7, ], Z2[8, ], Z2[9, ])
Z_rand_sample_in <- matrix( r1 - mean(r1), nrow = 1)

r1 <- c(Z1[10, ], Z2[10, ])
Z_rand_sample_out <- matrix( r1 - mean(r1), nrow = 1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  FITTING M1  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
index_in <- c(1:(n * (TT - 1)), n * TT + 1:(n * (TT - 1)))
index_out <- c(n * (TT - 1) + 1:n, n * TT + n * (TT - 1) + 1:n)
#index_in <- c(1:(n * (TT - 2)), n * TT + 1:(n * (TT - 2)))
#index_out <- c(n * (TT - 2) + 1:(n * 2), n * TT + n * (TT - 2) + 1:(n * 2))

NEGLOGLIK_SPACE <- function(p){

	theta <- exp(p[1:5])

	out <- 0

	for(ii in 1:2){

		Sigma <- matern_multi_step_one(theta[c(ii, 3, ii + 3)], DIST = dist0)
		cholmat <- t(cholesky(Sigma, parallel = TRUE))
		if( length(cholmat) == 0){
			return(Inf)
		}else{
			z <- forwardsolve(cholmat, t(Z_rand_sample0[[ii]]))
			logsig  <- 2 * sum(log(diag(cholmat))) * nrow(Z_rand_sample0[[ii]])
			out  <- out + 1/2 * logsig + 1/2 * sum(z^2)
			
		}
	}
	return(out)
}

dist0 <- as.matrix(dist(sim_grid_locations, diag = TRUE, upper = TRUE))	

Z_rand_sample0 <- list()

for(ii in 1:2){
	if(ii == 1) r1 <- Z1 	else	r1 <- Z2
	Z_rand_sample0_temp <- matrix(, ncol = n, nrow = TT)

	for(tt in 1:TT){
		Z_rand_sample0_temp[tt, ] <- r1[5 + tt, ] - mean(r1[5 + tt, ])
	}

	Z_rand_sample0[[ii]] <- Z_rand_sample0_temp
}

init <- c(0, 0, 0, 0, 0)

fit1 <- optim(par = init, fn = NEGLOGLIK_SPACE, control = list(trace = 5, maxit = 500)) #

for(tt in 1:5){
	fit1 <- optim(par = fit1$par, fn = NEGLOGLIK_SPACE, control = list(trace = 5, maxit = 500)) #
}

p_space <- fit1$par

#

NEGLOGLIK_ST <- function(p){

	if(p[1] <= -1 | p[1] >= 1){
		return(Inf)
	}else{

		theta <- c(exp(p_space), p[1])
		mu <- p[2:3]
		wind_var_chol <- matrix(c(p[4], p[5], 0, p[6]), ncol = 2, byrow = T)
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

init <- c(0, 0.001, 0.001, 1, 0, 1)

r1 <- c(Z1[6, ], Z1[7, ], Z1[8, ], Z1[9, ], Z2[6, ], Z2[7, ], Z2[8, ], Z1[9, ])

Z_rand_sample <- matrix( r1 - mean(r1), nrow = 1)

fit1 <- optim(par = init, fn = NEGLOGLIK_ST, control = list(trace = 5, maxit = 500)) #

for(tt in 1:10){
	write.table(c(p_space, fit1$par), file = paste(root, 'Results/multivariate_stationary_real_data_parameter_estimates_M1-two-step', sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
	fit1 <- optim(par = fit1$par, fn = NEGLOGLIK_ST, control = list(trace = 5, maxit = 500)) #
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  FITTING M2  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

NEGLOGLIK_ST <- function(p){
	cat(p, '\n')

	if(p[1] <= -1 | p[1] >= 1){
		return(Inf)
	}else{

		theta <- c(exp(p_space[1:5]), p[1])
		mu1 <- p[2:3]
		mu2 <- p[4:5]
		wind_var_chol <- matrix(c(p[6], p[7], p[8], p[9], 0, p[10], p[11], p[12], 0, 0, p[13], p[14], 0, 0, 0, p[15]), ncol = 4, byrow = T)
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

init <- c(0, 0.001, 0.001, 0.001, 0.001, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1)
init <- c(0.3044761, -1.7, 5.9, -1.7, 5.9, 0.67, 0.67, 0, 0, 0.35, 0, 0, 0.67, 0.67, 0.35)
#0.17644981 4.33303062 5.25706830 4.56330587 1.58028973 0.44942522 0.40738188 0.18691291 0.19320501 0.44648839 0.27776018 0.50376127 0.76536680 0.73309070 0.08729368

fit1 <- optim(par = init, fn = NEGLOGLIK_ST, control = list(trace = 5, maxit = 500)) #

for(tt in 1:10){
	write.table(c(p_space, fit1$par), file = paste(root, 'Results/multivariate_stationary_real_data_parameter_estimates_M2-two-step', sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
	fit1 <- optim(par = fit1$par, fn = NEGLOGLIK_ST, control = list(trace = 5, maxit = 500)) #
}

######################################  FAST M2  #######################################

NEGLOGLIK_ST_FAST <- function(p){

	if(p[1] <= -1 | p[1] >= 1){
		return(Inf)
	}else{
		cat(p, '\n')

		theta <- c(exp(p_space[1:5]), p[1])
		wind_mu1 <- mu1 <- p[2:3]
		wind_mu2 <- mu2 <- p[4:5]
		wind_var_chol <- matrix(c(p[6], p[7], p[8], p[9], 0, p[10], p[11], p[12], 0, 0, p[13], p[14], 0, 0, 0, p[15]), ncol = 4, byrow = T)
		wind_var <- t(wind_var_chol) %*% wind_var_chol
		wind_var1 = wind_var[1:2, 1:2]
		wind_var2 = wind_var[3:4, 3:4]
		wind_var12 = wind_var[1:2, 3:4]

		rho_temp <- theta[6]
		nu <- theta[4:5]
		beta <- theta[3]
		sigma2 <- theta[1:2]

		nu1 <- nu[1]
		nu2 <- nu[2]
		nu3 <- (nu1 + nu2)/2

		rho <- rho_temp * (gamma(nu1 + 3/2) / gamma(nu1))^(1/2) * (gamma(nu2 + 3/2) / gamma(nu2))^(1/2) * gamma(nu3) / (gamma(nu3 + 3/2))
		clusterExport(cl, c("theta", "sigma2", "beta", "nu", "nu3", "rho", "wind_mu1", "wind_mu2", "wind_var1", "wind_var2", "wind_var12"), envir = environment())

		output11 <- foreach(aa = 1:length(locs_cbind), .packages = "MASS", .noexport = "nonfrozen_rcpp") %dopar% {
			locs1 <- locs_cbind[[aa]][, 1:3]
			locs2 <- locs_cbind[[aa]][, 3 + 1:3]

			Sigma_temp <- tryCatch(nonfrozen_rcpp_multi_cross(Loc1 = locs1, Loc2 = locs2, param = c(sigma2[1], beta, nu[1]), v_mean = c(wind_mu1, wind_mu1), v_var = rbind(cbind(wind_var1, wind_var12), cbind(t(wind_var12), wind_var1))), error = function(a) numeric(0) )
			if(length(Sigma_temp) == 0){
				return(matrix(, ncol = 550, nrow = 550))
			}else{
				return(Sigma_temp)
			}

		}
		output22 <- foreach(aa = 1:length(locs_cbind), .packages = "MASS", .noexport = "nonfrozen_rcpp") %dopar% {
			locs1 <- locs_cbind[[aa]][, 1:3]
			locs2 <- locs_cbind[[aa]][, 3 + 1:3]

			Sigma_temp <- tryCatch(nonfrozen_rcpp_multi_cross(Loc1 = locs1, Loc2 = locs2, param = c(sigma2[2], beta, nu[2]), v_mean = c(wind_mu2, wind_mu2), v_var = rbind(cbind(wind_var2, wind_var12), cbind(t(wind_var12), wind_var2))), error = function(a) numeric(0) )

			if(length(Sigma_temp) == 0){
				return(matrix(, ncol = 550, nrow = 550))
			}else{
				return(Sigma_temp)
			}

		}
		output12 <- foreach(aa = 1:length(locs_cbind), .packages = "MASS", .noexport = "nonfrozen_rcpp") %dopar% {
			locs1 <- locs_cbind[[aa]][, 1:3]
			locs2 <- locs_cbind[[aa]][, 3 + 1:3]

			Sigma_temp <- tryCatch(nonfrozen_rcpp_multi_cross(Loc1 = locs1, Loc2 = locs2, param = c( sqrt(sigma2[1] * sigma2[2]) * rho, beta, nu3), v_mean = c(wind_mu1, wind_mu2), v_var = rbind(cbind(wind_var1, wind_var12), cbind(t(wind_var12), wind_var2))), error = function(a) numeric(0) )
			if(length(Sigma_temp) == 0){
				return(matrix(, ncol = 550, nrow = 550))
			}else{
				return(Sigma_temp)
			}

		}

		C11 <- C22 <- C12 <- matrix(, ncol = n * TT, nrow = n * TT)
		for(l1 in 1:TT){
			for(l2 in 1:TT){
				C11[(l1 - 1) * n + 1:n, (l2 - 1) * n + 1:n] <- output11[[(l1 - 1) * TT + l2]]
				C22[(l1 - 1) * n + 1:n, (l2 - 1) * n + 1:n] <- output22[[(l1 - 1) * TT + l2]]
				C12[(l1 - 1) * n + 1:n, (l2 - 1) * n + 1:n] <- output12[[(l1 - 1) * TT + l2]]
			}
		}

		Sigma <- rbind(cbind(C11, C12), cbind(t(C12), C22))

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

LOCS = sim_grid_locations
locs1 <- cbind(LOCS, rep(0, n))
for(tt in 1:(TT - 1)){
	locs1 <- rbind(locs1, cbind(LOCS, rep(tt, n)))
}

locs_list <- list()
for(ll in 1:TT){
	locs_list[[ll]] <- locs1[(ll - 1) * n + 1:n, ]
}

locs_cbind <- list()
for(l1 in 1:TT){
	for(l2 in 1:TT){
		locs_cbind[[(l1 - 1) * TT + l2]] <- cbind(locs_list[[l1]], locs_list[[l2]])
	}
}

p_space <- c(-2.2068347 -2.5709620  4.2753554  0.2801407  0.4386500)

cores=detectCores()
number_of_cores_to_use = cores[1]-1 # not to overload the computer
cl <- makeCluster(number_of_cores_to_use) 
registerDoParallel(cl)

clusterExport(cl, c("root", "TT", "n", "p_space", "locs_cbind"), envir = environment())
clusterEvalQ(cl, source(file = paste(root, "R_codes/Functions/load_packages.R", sep = '')))
clusterEvalQ(cl, sourceCpp(file = paste(root, "R_codes/Functions/spatially_varying_parameters2.cpp", sep = '')))
clusterEvalQ(cl, source(paste(root, "R_codes/Functions/cov_func.R", sep = '')))

fit1 <- optim(par = fit1$par, fn = NEGLOGLIK_ST_FAST, control = list(trace = 5, maxit = 500)) #

for(tt in 1:10){
	write.table(c(p_space, fit1$par), file = paste(root, 'Results/multivariate_stationary_real_data_parameter_estimates_M2-two-step', sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
	fit1 <- optim(par = fit1$par, fn = NEGLOGLIK_ST_FAST, control = list(trace = 5, maxit = 1500)) #
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  KRIGING M2  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

#write.table(c(Z_pred), file = paste(root, "Results/real-data-predictions-M1", sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)

p <- read.table(paste(root, 'Results/multivariate_stationary_real_data_parameter_estimates_M1-two-step', sep = ''), header = FALSE, sep = " ") %>% as.matrix()

TT <- 5

theta <- c(exp(p[1:5]), p[6])
mu1 <- p[7:8]
mu2 <- p[9:10]
wind_var_chol <- matrix(c(p[11], p[12], p[13], p[14], 0, p[15], p[16], p[17], 0, 0, p[18], p[19], 0, 0, 0, p[20]), ncol = 4, byrow = T)
wind_var <- t(wind_var_chol) %*% wind_var_chol

Sigma <- nonfrozen_matern_cov_multi_advec_small_scale(theta, wind_mu1 = mu1, wind_mu2 = mu2, wind_var1 = wind_var[1:2, 1:2], wind_var2 = wind_var[3:4, 3:4], wind_var12 = wind_var[1:2, 3:4], max_time_lag = TT - 1, LOCS = sim_grid_locations)	

C11 <- Sigma[index_in, index_in]
C22 <- Sigma[index_out, index_out]
C12 <- Sigma[index_in, index_out]

Z_pred <- t(C12) %*% solve(C11) %*% t(Z_rand_sample_in)

mean((c(Z_rand_sample_out) - c(Z_pred))^2)
#0.1917722


pdf(file = paste(root, 'thesis-defense/Figures/kriging-plots.pdf', sep = ''), width = 20, height = 5.5)

split.screen( rbind(c(0.04, 0.265, 0.1, 0.9), c(0.275, 0.485, 0.1, 0.9), c(0.515, 0.725, 0.1, 0.9), c(0.735, 0.945, 0.1, 0.9),
		      c(0.96,0.99,0.05,0.95)))
  
split.screen( figs = c( 2, 2 ), screen = 1 )
split.screen( figs = c( 2, 2 ), screen = 2 )
split.screen( figs = c( 2, 2 ), screen = 3 )
split.screen( figs = c( 2, 2 ), screen = 4 )

for(aa in 1:2){
	
	screen(5 + aa)
	par(pty = 's')
	par(mai=c(0.2,0.2,0.2,0.2))
	#quilt.plot(sim_grid_locations[, 1], sim_grid_locations[, 2], Z_pred[(aa - 1) * n + 1:n, ], nx = 25, ny = 25, ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', cex.lab = 4, add.legend = F, cex.axis = 2)
	quilt.plot(locs[, 1], locs[, 2], Z_pred[(aa - 1) * 550 + 1:550, ], nx = 25, ny = 25, ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', cex.lab = 4, add.legend = F, cex.axis = 2)
	map("worldHires", xlim = c(26.719, 85.078), ylim = c(5.625, 42.188), lwd = 0.75, add = T)

}

screen(5)
x1 <- c(0.025,0.09,0.09,0.025) + 0.12
y1 <- c(0.3,0.3,0.7,0.7)
  
legend.gradient2(cbind(x1,y1), title = "", limits = round(seq(0, 1, length.out = 3), 1), CEX = 1)
	  
close.screen( all=TRUE)
dev.off()

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

fit1 <- optim(par = init, fn = NEGLOGLIK_ST, control = list(trace = 5, maxit = 500)) #

for(tt in 1:10){
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

	if(p[1] <= -1 | p[1] >= 1 | p[4] < 0 | p[4] > 1){
		return(Inf)
	}else{

		theta <- c(exp(p_space[1:5]), p[1])

		dist_temp1 <- (exp(p[2]) * abs(dist_t)^(2 * exp(p[3])) + 1)
		DIST0 <- dist_s / (dist_temp1)^(p[4] / 2)

		Sigma <- purely_spatial_parsimonious_matern(theta, max_time_lag = TT - 1, DIST0)  / dist_temp1 

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

init <- c(0, 0, 0, 0.1)

fit1 <- optim(par = init, fn = NEGLOGLIK_ST, control = list(trace = 5, maxit = 500)) #

for(tt in 1:10){
	write.table(fit1$par, file = paste(root, 'Results/multivariate_stationary_real_data_parameter_estimates_M5', sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
	fit1 <- optim(par = fit1$par, fn = NEGLOGLIK_ST, control = list(trace = 5, maxit = 500)) #
}

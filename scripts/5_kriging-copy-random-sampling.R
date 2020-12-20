
workstation = T

if(workstation) directory <- '/home/salvanmo/Desktop/' 		else 		directory <- '/ibex/scratch/salvanmo/'

root <- paste(directory, 'Spatio-Temporal-Cross-Covariance-Functions-under-the-Lagrangian-Framework/', sep = '')

source(file = paste(root, "Functions/load_packages.R", sep = ''))
source(file = paste(root, "Functions/auxiliary_functions.R", sep = ''))

yr <- 2015

Z1 <- dat <- read.table(paste(root, 'Data/ncdf/layer1_residuals_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
Z2 <- dat2 <- read.table(paste(root, 'Data/ncdf/layer2_residuals_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
locs <- LOCS <- dat3 <- read.table(paste(root, 'Data/ncdf/LOCS-3D-dataset', sep = ''), header = FALSE, sep = " ") %>% as.matrix()

locs_cartesian <- cbind(6371 * cos(locs[, 2] / 180.0 * pi) * cos(locs[, 1] / 180.0 * pi), 6371 * cos(locs[, 2] / 180.0 * pi) * sin(locs[, 1] / 180.0 * pi))
locs_demean <- cbind(locs_cartesian[, 1] - mean(locs_cartesian[, 1]), locs_cartesian[, 2] - mean(locs_cartesian[, 2]))
std_locs <- cbind((locs_cartesian[, 1] - mean(locs_cartesian[, 1])) / sd(locs_cartesian[, 1]), (locs_cartesian[, 2] - mean(locs_cartesian[, 2])) / sd(locs_cartesian[, 2]))

set.seed(1234)
subs <- sample(1:550, 385)
sim_grid_locations <- locs_demean[subs, ]

Z1 <- Z1[, subs]
Z2 <- Z2[, subs]

n <- nrow(sim_grid_locations)

TT <- 4

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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  KRIGING M1  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
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

#-2.5128371 -2.3847379  4.2463860  0.4564371  0.5135349
# -2.2068347 -2.5709620  4.2753554  0.2801407  0.4386500

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

subs_out <- seq(1, 550, 1)[-subs]
index_in <- c(subs, 550 + subs, 550 * 2 + subs, 500 * 3 + subs, 550 * TT + subs, 550 * TT + 550 + subs, 550 * TT + 550 * 2 + subs, 550 * TT + 500 * 3 + subs)
index_out <- c(subs_out, 550 + subs_out, 550 * 2 + subs_out, 500 * 3 + subs_out, 550 * TT + subs_out, 550 * TT + 550 + subs_out, 550 * TT + 550 * 2 + subs_out, 550 * TT + 500 * 3 + subs_out)

theta <- c(exp(p[1:5]), p[6])
mu <- p[7:8]
wind_var_chol <- matrix(c(p[9], p[10], 0, p[11]), ncol = 2, byrow = T)
wind_var <- t(wind_var_chol) %*% wind_var_chol

Sigma <- nonfrozen_matern_cov_multi_small_scale(theta, wind_mu = mu, wind_var = wind_var, max_time_lag = TT - 1, LOCS = sim_grid_locations)
C11 <- Sigma[index_in, index_in]
C22 <- Sigma[index_out, index_out]
C12 <- Sigma[index_in, index_out]

Z_pred <- t(C12) %*% solve(C11) %*% t(Z_rand_sample_in)

DAT <- c(dat[6, ], dat[7, ], dat[8, ], dat[9, ], dat2[6, ], dat2[7, ], dat2[8, ], dat2[9, ])

DAT[-index_in] <- c(Z_pred)


pdf(file = paste(root, 'thesis-defense/Figures/kriging-plots.pdf', sep = ''), width = 20, height = 5.5)

zlim_range1 <- zlim_range2 <- range(DAT)

split.screen( rbind(c(0.08,0.95,0.1,0.95), c(0.95,0.99,0.1,0.95)))
split.screen( figs = c( 2, 5 ), screen = 1 )

hr_count <- 0
hr_label <- c('0:00', '3:00', '6:00', '9:00', '12:00', '15:00', '18:00', '21:00')
for(hr in 1:4){
	
	for(variable in 1:2){
		
		screen((variable - 1) * 5 + 2 + hr_count)
		par(pty = 's')
		par(mai=c(0.2,0.2,0.2,0.2))
		
		if(hr_count == 1 & variable == 2){
		quilt.plot(dat3[, 1], dat3[, 2], DAT[550 * TT + (hr - 1) * 550 + 1:550], nx = 25, ny = 25, zlim = zlim_range2, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 2)
		mtext('985 hPa', side = 2, line = 7, adj = 0.5, cex = 3, font = 2, col = 'blue')
		}else if(hr_count == 1 & variable == 1){
		quilt.plot(dat3[, 1], dat3[, 2], DAT[(hr - 1) * 550 + 1:550], nx = 25, ny = 25, zlim = zlim_range1, ylab = '', xlab = '', xaxt = 'n', cex.lab = 4, add.legend = F, cex.axis = 2)
		mtext('850 hPa', side = 2, line = 7, adj = 0.5, cex = 3, font = 2, col = 'blue')
		}else if(variable == 2){
		quilt.plot(dat3[, 1], dat3[, 2], DAT[550 * TT + (hr - 1) * 550 + 1:550], nx = 25, ny = 25, zlim = zlim_range2, ylab = '', xlab = '', yaxt = 'n', cex.lab = 4, add.legend = F, cex.axis = 2)
		}else{
		quilt.plot(dat3[, 1], dat3[, 2], DAT[(hr - 1) * 550 + 1:550], nx = 25, ny = 25, zlim = zlim_range1, ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', cex.lab = 4, add.legend = F, cex.axis = 2)
		}
		map("worldHires", xlim = c(26.719, 85.078), ylim = c(5.625, 42.188), lwd = 0.75, add = T)
	}				
}

screen(2)
x1 <- c(0.025,0.12,0.12,0.025) + 0.1
y1 <- c(0.15,0.15,0.8,0.8)
legend.gradient2(cbind(x1,y1), title = "", limits = round(seq(zlim_range2[1], zlim_range2[2], length.out = 5), 1), CEX = 2)

close.screen( all=TRUE)

dev.off()




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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  KRIGING M2  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

p <- read.table(paste(root, 'Results/multivariate_stationary_real_data_parameter_estimates_M2', sep = ''), header = FALSE, sep = " ") %>% as.matrix()

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

for(aa in 1:4){
	
	screen(5 + aa)
	par(pty = 's')
	par(mai=c(0.2,0.2,0.2,0.2))
	quilt.plot(locs[, 1], locs[, 2], Z_pred[(aa - 1) * 550 + 1:550, ], nx = 25, ny = 25, ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', cex.lab = 4, add.legend = F, cex.axis = 2)
	map("worldHires", xlim = c(26.719, 85.078), ylim = c(5.625, 42.188), lwd = 0.75, add = T)

}

screen(5)
x1 <- c(0.025,0.09,0.09,0.025) + 0.12
y1 <- c(0.3,0.3,0.7,0.7)
  
legend.gradient2(cbind(x1,y1), title = "", limits = round(seq(0, 1, length.out = 3), 1), CEX = 1)
	  
close.screen( all=TRUE)
dev.off()



#remove mean

directory <- '/home/salvanmo/Desktop/'

root <- paste(directory, 'Spatio-Temporal-Cross-Covariance-Functions-under-the-Lagrangian-Framework/', sep = '')

source(file = paste(root, "Functions/load_packages.R",sep=''))
sourceCpp(file=paste("/home/salvanmo/Desktop/thesis/R_codes/Functions/distR.cpp",sep=''))

yr = 2019 #2008	2015

#ln <- NULL

#for(yr in 1980:2019){

cat(yr, '\n')

dat_temp <- read.table(paste(root, 'Data/ncdf/layer1_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
dat2_temp <- read.table(paste(root, 'Data/ncdf/layer2_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
dat3 <- locs <- read.table(paste(root, 'Data/ncdf/LOCS-3D-dataset', sep = ''), header = FALSE, sep = " ") %>% as.matrix()

humid <- humid_temp <- read.table(paste(root, 'Data/ncdf/humid_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
temp <- temp_temp <- read.table(paste(root, 'Data/ncdf/temp_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()

load(paste(root, 'Data/ncdf/covariates_' , yr, '.RData', sep = ''))

locs_cartesian <- cbind(6371 * cos(locs[, 2] / 180.0 * pi) * cos(locs[, 1] / 180.0 * pi), 6371 * cos(locs[, 2] / 180.0 * pi) * sin(locs[, 1] / 180.0 * pi))
locs_demean <- cbind(locs_cartesian[, 1] - mean(locs_cartesian[, 1]), locs_cartesian[, 2] - mean(locs_cartesian[, 2]))
std_locs <- cbind((locs_cartesian[, 1] - mean(locs_cartesian[, 1])) / sd(locs_cartesian[, 1]), (locs_cartesian[, 2] - mean(locs_cartesian[, 2])) / sd(locs_cartesian[, 2]))

std_locs_with_time <- cbind(rep((locs[, 1] - mean(locs[, 1])) / sd(locs[, 1]), 6), rep((locs[, 2] - mean(locs[, 2])) / sd(locs[, 2]), 6), c(rep(0, 550), rep(1, 550), rep(2, 550), rep(3, 550), rep(4, 550), rep(5, 550)))
std_locs_with_time[, 3] <- std_locs_with_time[, 3] / 6

#ln <- c(ln, length(which(dat_temp < -25)))
#}

dat_temp[which(dat_temp < -27)] <- -27
dat2_temp[which(dat2_temp < -27)] <- -27

DAT <- dat_temp
DAT2 <- dat2_temp

args <- commandArgs(trailingOnly = TRUE)
#start_hr <- 125 #130
start_hr <- as.numeric(args[1])

Z_rand_sample <- matrix(c(DAT[start_hr + 1, ], DAT[start_hr + 2, ], DAT[start_hr + 3, ], DAT[start_hr + 4, ], DAT[start_hr + 5, ], DAT[start_hr + 6, ]), ncol = 1)

NN <- 10
grid_x <- seq(from = min(std_locs_with_time[, 1]), to = max(std_locs_with_time[, 1]), length.out = NN)
grid_y <- seq(from = min(std_locs_with_time[, 2]), to = max(std_locs_with_time[, 2]), length.out = NN)
X <- expand.grid(grid_x, grid_y) %>% as.matrix()

htarg <- distR_C(std_locs_with_time, cbind(X, rep(0, nrow(X))))
sigma <- htarg^2 * log(htarg)
sigma[htarg == 0] <- 0

Y <- cbind(rep(1, nrow(std_locs_with_time)), std_locs_with_time, sigma)
beta_hat <- solve(t(Y) %*% Y + 100 * diag(ncol(Y))) %*% t(Y) %*% Z_rand_sample
range(beta_hat)

mu_pred <- Y %*% beta_hat
sum((mu_pred - Z_rand_sample)^2)

res <- mu_pred - Z_rand_sample

jpeg(file = paste('/home/salvanmo/Desktop/thesis/thesis-defense/Figures/residuals.jpg', sep = ''), width = 1500, height = 1000)
par(mfrow = c(2, 3))
for(ll in 1:6){
	quilt.plot(dat3[, 1], dat3[, 2], res[(ll - 1) * 550 + 1:550, ], zlim = range(res), nx = 25, ny = 25, ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', cex.lab = 4, add.legend = F, cex.axis = 2)
}
dev.off()

###########################################################################
start_hr <- 140 #167, 182

Z_rand_sample <- matrix(c(DAT[start_hr + 1, ], DAT[start_hr + 2, ], DAT[start_hr + 3, ], DAT[start_hr + 4, ], DAT[start_hr + 5, ], DAT[start_hr + 6, ]), ncol = 1)

X1 <- X2 <- X3 <- X4 <- X5 <- NULL
for(ll in 1:6){
	X1 <- c(X1, data_array[start_hr + ll, , 1])
	X2 <- c(X2, data_array[start_hr + ll, , 2])
	X3 <- c(X3, data_array[start_hr + ll, , 3])
	X4 <- c(X4, data_array[start_hr + ll, , 4])
	X5 <- c(X5, rep(ll, 550))
}

X <- cbind(rep(1, length(X1)), (X1 - mean(X1)) / sd(X1), (X2 - mean(X2)) / sd(X2), (X3 - mean(X3)) / sd(X3), (X4 - mean(X4)) / sd(X4), rep(std_locs[, 1], 6), rep(std_locs[, 2], 6), X5)

data_full <- data.frame(Z_rand_sample, X)
colnames(data_full) <- c('Y', 'intercept', 'humidity', 'temperature', 'U', 'V', 'locsx', 'locsy', 'time')

mod_lm2 <- gam(Y ~ humidity + temperature + U + V + locsx + locsy + time, data=data_full)
summary(mod_lm2)

beta_hat <- matrix(mod_lm2[[1]], ncol = 1)
mu_pred <- X %*% beta_hat
sum((mu_pred - Z_rand_sample)^2)
res <- mu_pred - Z_rand_sample


jpeg(file = paste('/home/salvanmo/Desktop/thesis/thesis-defense/Figures/residuals.jpg', sep = ''), width = 1000, height = 1000)

quilt.plot(dat3[, 1], dat3[, 2], res, nx = 25, ny = 25, ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', cex.lab = 4, add.legend = F, cex.axis = 2)

dev.off()








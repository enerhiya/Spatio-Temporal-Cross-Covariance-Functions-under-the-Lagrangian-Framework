#remove mean

directory <- '/home/salvanmo/Desktop/'

root <- paste(directory, 'Spatio-Temporal-Cross-Covariance-Functions-under-the-Lagrangian-Framework/', sep = '')

source(file = paste(root, "Functions/load_packages.R",sep=''))

yr = 2019

cat(yr, '\n')

dat_temp <- read.table(paste(root, 'Data/ncdf/layer1_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
dat2_temp <- read.table(paste(root, 'Data/ncdf/layer2_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
dat3 <- locs <- read.table(paste(root, 'Data/ncdf/LOCS-3D-dataset', sep = ''), header = FALSE, sep = " ") %>% as.matrix()

load(paste(root, 'Data/ncdf/covariates_' , yr, '.RData', sep = ''))

locs_cartesian <- cbind(6371 * cos(locs[, 2] / 180.0 * pi) * cos(locs[, 1] / 180.0 * pi), 6371 * cos(locs[, 2] / 180.0 * pi) * sin(locs[, 1] / 180.0 * pi))
locs_demean <- cbind(locs_cartesian[, 1] - mean(locs_cartesian[, 1]), locs_cartesian[, 2] - mean(locs_cartesian[, 2]))
std_locs <- cbind((locs_cartesian[, 1] - mean(locs_cartesian[, 1])) / sd(locs_cartesian[, 1]), (locs_cartesian[, 2] - mean(locs_cartesian[, 2])) / sd(locs_cartesian[, 2]))

dat_temp[which(dat_temp < -27)] <- -27
dat2_temp[which(dat2_temp < -27)] <- -27

DAT <- dat_temp
DAT2 <- dat2_temp

start_hr <- 140

VAR <- 2

if(VAR == 1){
	Z_rand_sample <- matrix(c(DAT[start_hr + 1, ], DAT[start_hr + 2, ], DAT[start_hr + 3, ], DAT[start_hr + 4, ], DAT[start_hr + 5, ], DAT[start_hr + 6, ]), ncol = 1)
}else{
	Z_rand_sample <- matrix(c(DAT2[start_hr + 1, ], DAT2[start_hr + 2, ], DAT2[start_hr + 3, ], DAT2[start_hr + 4, ], DAT2[start_hr + 5, ], DAT2[start_hr + 6, ]), ncol = 1)
}

X1 <- X2 <- X3 <- X4 <- X5 <- NULL
for(ll in 1:6){
	X1 <- c(X1, data_array[start_hr + ll, , 1, VAR])
	X2 <- c(X2, data_array[start_hr + ll, , 2, VAR])
	X3 <- c(X3, data_array[start_hr + ll, , 3, VAR])
	X4 <- c(X4, data_array[start_hr + ll, , 4, VAR])
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

jpeg(file = paste('/home/salvanmo/Desktop/thesis/thesis-defense/Figures/residuals.jpg', sep = ''), width = 1500, height = 1000)
par(mfrow = c(2, 3))
for(ll in 1:6){
	quilt.plot(dat3[, 1], dat3[, 2], res[(ll - 1) * 550 + 1:550, ], zlim = range(res), nx = 25, ny = 25, ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', cex.lab = 4, add.legend = F, cex.axis = 2)
}
dev.off()








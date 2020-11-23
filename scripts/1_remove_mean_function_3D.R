#remove mean

directory <- '/home/salvanmo/Desktop/'

root <- paste(directory, 'Spatio-Temporal-Cross-Covariance-Functions-under-the-Lagrangian-Framework/', sep = '')

source(file = paste(root, "Functions/load_packages.R",sep=''))

for(yr in 1980:2019){

	cat(yr, '\n')

	dat <- read.table(paste(root, 'Data/ncdf/layer1_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()

	dat2 <- read.table(paste(root, 'Data/ncdf/layer2_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()

	DAT <- dat
	DAT2 <- dat2

	Yhat1 <- res_mat1 <- DAT <- DAT - matrix(colMeans(DAT), ncol = ncol(DAT), nrow = nrow(DAT), byrow = T)
	Yhat2 <- res_mat2 <- DAT2 <- DAT2 - matrix(colMeans(DAT2), ncol = ncol(DAT2), nrow = nrow(DAT2), byrow = T)

	obs_mat_SVD1 <- svd(DAT)
	obs_mat_SVD2 <- svd(DAT2)

	variance.explained1 = prop.table(obs_mat_SVD1$d^2)
	variance.explained2 = prop.table(obs_mat_SVD2$d^2)

	percent_sum_squared_variation1 <- cumsum(variance.explained1)[cumsum(variance.explained1) >= 0.85]
	min_percent_sum_squared_variation1 <- min(percent_sum_squared_variation1)
	num_singular_vec1 <- which(cumsum(variance.explained1) == min_percent_sum_squared_variation1) 
	#num_singular_vec1 <- 7 

	percent_sum_squared_variation2 <- cumsum(variance.explained2)[cumsum(variance.explained2) >= 0.85]
	min_percent_sum_squared_variation2 <- min(percent_sum_squared_variation2)
	num_singular_vec2 <- which(cumsum(variance.explained2) == min_percent_sum_squared_variation2) 
	#num_singular_vec2 <- 7

	X1 <- cbind(rep(1, nrow(obs_mat_SVD1$u)), obs_mat_SVD1$u)
	X2 <- cbind(rep(1, nrow(obs_mat_SVD2$u)), obs_mat_SVD2$u)

	for(nn in 1:ncol(DAT)){

		beta_hat <- solve(t(X1[, 1:num_singular_vec1]) %*% X1[, 1:num_singular_vec1]) %*% t(X1[, 1:num_singular_vec1]) %*% DAT[, nn]

		Yhat1[, nn] <- Yhat_temp <- X1[, 1:num_singular_vec1] %*% beta_hat

		err <- Yhat_temp - DAT[, nn]

		res_mat1[, nn] <- err

		beta_hat <- solve(t(X2[, 1:num_singular_vec2]) %*% X2[, 1:num_singular_vec2]) %*% t(X2[, 1:num_singular_vec2]) %*% DAT2[, nn]

		Yhat2[, nn] <- Yhat_temp <- X2[, 1:num_singular_vec2] %*% beta_hat

		err <- Yhat_temp - DAT2[, nn]

		res_mat2[, nn] <- err

	}
	
	#res_mat1 <- DAT
	#res_mat2 <- DAT2

	write.table(res_mat1, file = paste(root, "Data/ncdf/layer1_residuals_", yr, sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
	write.table(res_mat2, file = paste(root, "Data/ncdf/layer2_residuals_", yr, sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)

	write.table(Yhat1, file = paste(root, "Results/estimated_mean/layer1_trend_", yr, sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
	write.table(Yhat2, file = paste(root, "Results/estimated_mean/layer2_trend_", yr, sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)

}

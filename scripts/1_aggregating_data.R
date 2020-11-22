#remove mean

directory <- '/home/salvanmo/Desktop/'

root <- paste(directory, 'Spatio-Temporal-Cross-Covariance-Functions-under-the-Lagrangian-Framework/', sep = '')

source(file = paste(root, "Functions/load_packages.R",sep=''))
source(file = paste(root, "Functions/auxiliary_functions.R",sep=''))

DAT_array <- DAT2_array <- array(, dim = c(248, 550, 40))

for(yr in 1980:2019){

	cat(yr, '\n')

	dat_temp <- read.table(paste(root, 'Data/ncdf/layer1_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()

	dat2_temp <- read.table(paste(root, 'Data/ncdf/layer2_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()

	DAT_array[,, yr - 1979] <- dat_temp
	DAT2_array[,, yr - 1979] <- dat2_temp

}

DAT <- apply(DAT_array, c(1, 2), mean)
DAT2 <- apply(DAT2_array, c(1, 2), mean)

Yhat1 <- res_mat1 <- DAT <- DAT - matrix(colMeans(DAT), ncol = ncol(DAT), nrow = nrow(DAT), byrow = T)
Yhat2 <- res_mat2 <- DAT2 <- DAT2 - matrix(colMeans(DAT2), ncol = ncol(DAT2), nrow = nrow(DAT2), byrow = T)

obs_mat_SVD1 <- svd(DAT)
obs_mat_SVD2 <- svd(DAT2)

variance.explained1 = prop.table(obs_mat_SVD1$d^2)
variance.explained2 = prop.table(obs_mat_SVD2$d^2)

percent_sum_squared_variation1 <- cumsum(variance.explained1)[cumsum(variance.explained1) >= 0.85]
min_percent_sum_squared_variation1 <- min(percent_sum_squared_variation1)
#num_singular_vec1 <- which(cumsum(variance.explained1) == min_percent_sum_squared_variation1) 
num_singular_vec1 <- 3 

percent_sum_squared_variation2 <- cumsum(variance.explained2)[cumsum(variance.explained2) >= 0.85]
min_percent_sum_squared_variation2 <- min(percent_sum_squared_variation2)
#num_singular_vec2 <- which(cumsum(variance.explained2) == min_percent_sum_squared_variation2) 
num_singular_vec2 <- 3

#lm_model <- lm(DAT[, 1] ~ obs_mat_SVD1$u[, 1:11])
#summary(lm_model)

#anova(lm(DAT[, 1] ~ obs_mat_SVD1$u[, 1:4]), lm(DAT[, 1] ~ obs_mat_SVD1$u[, 1:5]))

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

dat <- res_mat1	
dat2 <- res_mat2	


pdf(file = paste(root, 'Figures/spacetime-maps-residuals-supplementary.pdf', sep = ''), width = 25, height = 10)

day_count <- 0

for(start_hr in seq(1, 200, by = 5)){

	cat(start_hr, '\n')

	hr_index <- seq(start_hr, start_hr + 4, by = 1)

	dat3 <- read.table(paste(root, 'Data/ncdf/LOCS-3D-dataset', sep = ''), header = FALSE, sep = " ") %>% as.matrix()

	#zlim_range1 <- range(dat)
	#zlim_range2 <- range(dat2)

	zlim_range1 <- c(-0.4, 0.4)
	zlim_range2 <- c(-0.4, 0.4)

	split.screen( rbind(c(0.08,0.95,0.1,0.95), c(0.95,0.99,0.1,0.95)))
	split.screen( figs = c( 2, 5 ), screen = 1 )

	hr_count <- 0
	for(hr in hr_index){
		
		if(mod(hr - 1, 3) == 0){
			hr_label <- 0
		}else{
			hr_label <- hr_label + 1
		}

		if(mod(hr - 1, 3) == 0){
			day_count <- day_count + 1
		}

		if(day_count < 10)	day_label <- paste('0', day_count, sep = '')	else 	day_label <- day_count

		hr_count <- hr_count + 1
		
		for(variable in 1:2){
			
			screen((variable - 1) * 5 + 2 + hr_count)
			par(pty = 's')
			par(mai=c(0.2,0.2,0.2,0.2))
			
			if(hr_count == 1 & variable == 2){
			quilt.plot(dat3[, 1], dat3[, 2], dat2[hr, ], zlim = zlim_range2, nx = 25, ny = 25, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 2)
			mtext('1000 hPa', side = 2, line = 7, adj = 0.5, cex = 3, font = 2, col = 'blue')
			}else if(hr_count == 1 & variable == 1){
			quilt.plot(dat3[, 1], dat3[, 2], dat[hr, ], zlim = zlim_range1, nx = 25, ny = 25, ylab = '', xlab = '', xaxt = 'n', cex.lab = 4, add.legend = F, cex.axis = 2)
			mtext('850 hPa', side = 2, line = 7, adj = 0.5, cex = 3, font = 2, col = 'blue')
			}else if(variable == 2){
			quilt.plot(dat3[, 1], dat3[, 2], dat2[hr, ], zlim = zlim_range2, nx = 25, ny = 25, ylab = '', xlab = '', yaxt = 'n', cex.lab = 4, add.legend = F, cex.axis = 2)
			}else{
			quilt.plot(dat3[, 1], dat3[, 2], dat[hr, ], zlim = zlim_range1, nx = 25, ny = 25, ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', cex.lab = 4, add.legend = F, cex.axis = 2)
			}
			map("worldHires", xlim = c(26.719, 85.078), ylim = c(5.625, 42.188), lwd = 0.75, add = T)
			
			if(hr_count == 1){
				mtext('Latitude', side = 2, line = 4, adj = 0.5, cex = 2.5, font = 2)
			}

			if(variable == 1){
				mtext(paste('January ', day_label, ', ', hr_label * 8, ':00', sep = ''), side = 3, line = 1, adj = 0.5, cex = 3, font = 2)
			}else{
				mtext('Longitude', side = 1, line = 4, adj = 0.5,  cex = 2.5, font = 2)
			}
		}				
	}
	screen(2)
	x1 <- c(0.025,0.12,0.12,0.025) + 0.1
	y1 <- c(0.05,0.05,0.35,0.35)
	legend.gradient2(cbind(x1,y1), title = "", limits = round(seq(zlim_range2[1], zlim_range2[2], length.out = 5), 1), cex = 2)

	x1 <- c(0.025,0.12,0.12,0.025) + 0.1
	y1 <- c(0.6,0.6,0.9,0.9)
	legend.gradient2(cbind(x1,y1), title = "", limits = round(seq(zlim_range1[1], zlim_range1[2], length.out = 5), 1), cex = 2)

	close.screen( all=TRUE)

}

dev.off()

pdf(file = paste(root, 'Figures/spacetime-maps-residuals-manuscript.pdf', sep = ''), width = 25, height = 10)

day_count <- 0

for(start_hr in 1:1){

	cat(start_hr, '\n')

	hr_index <- seq(start_hr, start_hr + 4, by = 1)

	dat3 <- read.table(paste(root, 'Data/ncdf/LOCS-3D-dataset', sep = ''), header = FALSE, sep = " ") %>% as.matrix()

	zlim_range1 <- zlim_range2 <- c(-0.3, 0.3)

	split.screen( rbind(c(0.08,0.95,0.1,0.95), c(0.95,0.99,0.1,0.95)))
	split.screen( figs = c( 2, 5 ), screen = 1 )

	hr_count <- 0
	for(hr in hr_index){
		
		if(mod(hr - 1, 3) == 0){
			hr_label <- 0
		}else{
			hr_label <- hr_label + 1
		}

		if(mod(hr - 1, 3) == 0){
			day_count <- day_count + 1
		}

		if(day_count < 10)	day_label <- paste('0', day_count, sep = '')	else 	day_label <- day_count

		hr_count <- hr_count + 1
		
		for(variable in 1:2){
			
			screen((variable - 1) * 5 + 2 + hr_count)
			par(pty = 's')
			par(mai=c(0.2,0.2,0.2,0.2))
			
			if(hr_count == 1 & variable == 2){
			quilt.plot(dat3[, 1], dat3[, 2], dat2[hr, ], zlim = zlim_range2, nx = 25, ny = 25, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 2)
			mtext('1000 hPa', side = 2, line = 7, adj = 0.5, cex = 3, font = 2, col = 'blue')
			}else if(hr_count == 1 & variable == 1){
			quilt.plot(dat3[, 1], dat3[, 2], dat[hr, ], zlim = zlim_range1, nx = 25, ny = 25, ylab = '', xlab = '', xaxt = 'n', cex.lab = 4, add.legend = F, cex.axis = 2)
			mtext('850 hPa', side = 2, line = 7, adj = 0.5, cex = 3, font = 2, col = 'blue')
			}else if(variable == 2){
			quilt.plot(dat3[, 1], dat3[, 2], dat2[hr, ], zlim = zlim_range2, nx = 25, ny = 25, ylab = '', xlab = '', yaxt = 'n', cex.lab = 4, add.legend = F, cex.axis = 2)
			}else{
			quilt.plot(dat3[, 1], dat3[, 2], dat[hr, ], zlim = zlim_range1, nx = 25, ny = 25, ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', cex.lab = 4, add.legend = F, cex.axis = 2)
			}
			map("worldHires", xlim = c(26.719, 85.078), ylim = c(5.625, 42.188), lwd = 0.75, add = T)
			
			if(hr_count == 1){
				mtext('Latitude', side = 2, line = 4, adj = 0.5, cex = 2.5, font = 2)
			}

			if(variable == 1){
				mtext(paste('January ', day_label, ', ', hr_label * 8, ':00', sep = ''), side = 3, line = 1, adj = 0.5, cex = 3, font = 2)
			}else{
				mtext('Longitude', side = 1, line = 4, adj = 0.5,  cex = 2.5, font = 2)
			}
		}				
	}
	screen(2)
	x1 <- c(0.025,0.12,0.12,0.025) + 0.1
	y1 <- c(0.15,0.15,0.8,0.8)
	legend.gradient2(cbind(x1,y1), title = "", limits = round(seq(zlim_range2[1], zlim_range2[2], length.out = 5), 1), CEX = 2)

	close.screen( all=TRUE)

}

dev.off()


pdf(file = paste(root, 'Figures/spacetime-maps-mean-manuscript.pdf', sep = ''), width = 25, height = 10)

day_count <- 0

for(start_hr in 1:1){

	cat(start_hr, '\n')

	hr_index <- seq(start_hr, start_hr + 4, by = 1)

	dat3 <- read.table(paste(root, 'Data/ncdf/LOCS-3D-dataset', sep = ''), header = FALSE, sep = " ") %>% as.matrix()

	zlim_range1 <- zlim_range2 <- c(-0.4, 0.4)
	#zlim_range2 <- range(Yhat2[start_hr:(start_hr + 4), ])

	split.screen( rbind(c(0.08,0.95,0.1,0.95), c(0.95,0.99,0.1,0.95)))
	split.screen( figs = c( 2, 5 ), screen = 1 )

	hr_count <- 0
	for(hr in hr_index){
		
		if(mod(hr - 1, 3) == 0){
			hr_label <- 0
		}else{
			hr_label <- hr_label + 1
		}

		if(mod(hr - 1, 3) == 0){
			day_count <- day_count + 1
		}

		if(day_count < 10)	day_label <- paste('0', day_count, sep = '')	else 	day_label <- day_count

		hr_count <- hr_count + 1
		
		for(variable in 1:2){
			
			screen((variable - 1) * 5 + 2 + hr_count)
			par(pty = 's')
			par(mai=c(0.2,0.2,0.2,0.2))
			
			if(hr_count == 1 & variable == 2){
			quilt.plot(dat3[, 1], dat3[, 2], Yhat2[hr, ], zlim = zlim_range2, nx = 25, ny = 25, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 2)
			mtext('1000 hPa', side = 2, line = 7, adj = 0.5, cex = 3, font = 2, col = 'blue')
			}else if(hr_count == 1 & variable == 1){
			quilt.plot(dat3[, 1], dat3[, 2], Yhat1[hr, ], zlim = zlim_range1, nx = 25, ny = 25, ylab = '', xlab = '', xaxt = 'n', cex.lab = 4, add.legend = F, cex.axis = 2)
			mtext('850 hPa', side = 2, line = 7, adj = 0.5, cex = 3, font = 2, col = 'blue')
			}else if(variable == 2){
			quilt.plot(dat3[, 1], dat3[, 2], Yhat2[hr, ], zlim = zlim_range2, nx = 25, ny = 25, ylab = '', xlab = '', yaxt = 'n', cex.lab = 4, add.legend = F, cex.axis = 2)
			}else{
			quilt.plot(dat3[, 1], dat3[, 2], Yhat1[hr, ], zlim = zlim_range1, nx = 25, ny = 25, ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', cex.lab = 4, add.legend = F, cex.axis = 2)
			}
			map("worldHires", xlim = c(26.719, 85.078), ylim = c(5.625, 42.188), lwd = 0.75, add = T)
			
			if(hr_count == 1){
				mtext('Latitude', side = 2, line = 4, adj = 0.5, cex = 2.5, font = 2)
			}

			if(variable == 1){
				mtext(paste('January ', day_label, ', ', hr_label * 8, ':00', sep = ''), side = 3, line = 1, adj = 0.5, cex = 3, font = 2)
			}else{
				mtext('Longitude', side = 1, line = 4, adj = 0.5,  cex = 2.5, font = 2)
			}
		}				
	}
	screen(2)
	x1 <- c(0.025,0.12,0.12,0.025) + 0.1
	y1 <- c(0.15,0.15,0.8,0.8)
	legend.gradient2(cbind(x1,y1), title = "", limits = round(seq(zlim_range2[1], zlim_range2[2], length.out = 5), 1), CEX = 2)

	close.screen( all=TRUE)

}

dev.off()

write.table(res_mat1, file = paste(root, "Data/ncdf/layer1_residuals_aggregate", sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
write.table(res_mat2, file = paste(root, "Data/ncdf/layer2_residuals_aggregate", sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)

write.table(Yhat1, file = paste(root, "Results/estimated_mean/layer1_trend_aggregate", sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
write.table(Yhat2, file = paste(root, "Results/estimated_mean/layer2_trend_aggregate", sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)


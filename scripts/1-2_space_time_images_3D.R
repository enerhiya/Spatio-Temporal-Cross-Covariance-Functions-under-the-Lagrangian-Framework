
directory <- '/home/salvanmo/Desktop/'

root <- paste(directory, 'Spatio-Temporal-Cross-Covariance-Functions-under-the-Lagrangian-Framework/', sep = '')

source(file = paste(root, "Functions/load_packages.R",sep=''))
source(file = paste(root, "Functions/auxiliary_functions.R",sep=''))

yr <- 2015 #2010, 2016 is good

#dat <- read.table(paste(root, 'Data/ncdf/layer1_residuals_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
#dat2 <- read.table(paste(root, 'Data/ncdf/layer2_residuals_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
#dat <- read.table(paste(root, 'Data/ncdf/layer1_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
#dat2 <- read.table(paste(root, 'Data/ncdf/layer2_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
#dat <- read.table(paste(root, 'Data/ncdf/TEST-layer1_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
#dat2 <- read.table(paste(root, 'Data/ncdf/TEST-layer2_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
dat <- read.table(paste(root, 'Data/ncdf/TEST-layer1_residuals_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
dat2 <- read.table(paste(root, 'Data/ncdf/TEST-layer2_residuals_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
dat3 <- read.table(paste(root, 'Data/ncdf/LOCS-3D-dataset', sep = ''), header = FALSE, sep = " ") %>% as.matrix()

Yhat1 <- read.table(paste(root, "Results/estimated_mean/layer1_trend_", yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
Yhat2 <- read.table(paste(root, "Results/estimated_mean/layer2_trend_", yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()

pdf(file = paste(root, 'Figures/spacetime-maps-residuals-supplementary.pdf', sep = ''), width = 25, height = 10)

day_count <- 0

for(start_hr in 1:1){
#for(start_hr in seq(1, 200, by = 5)){

	cat(start_hr, '\n')

	hr_index <- seq(start_hr, start_hr + 4, by = 1)

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

pdf(file = paste(root, 'Figures/spacetime-maps-residuals-manuscript-NEW.pdf', sep = ''), width = 25, height = 10)
#pdf(file = paste(root, 'Figures/spacetime-maps-residuals-manuscript.pdf', sep = ''), width = 25, height = 10)

day_count <- 0

for(start_hr in seq(1, 240, by = 5)){

	cat(start_hr, '\n')

	dat2 <- dat2 - matrix(colMeans(dat2), nrow = nrow(dat2), ncol = ncol(dat2), byrow = T) 
	dat <- dat - matrix(colMeans(dat), nrow = nrow(dat), ncol = ncol(dat), byrow = T) 

	hr_index <- seq(start_hr, start_hr + 4, by = 1)

	zlim_range1 <- range(dat[start_hr:(start_hr + 4), ])
	zlim_range2 <- range(dat2[start_hr:(start_hr + 4), ])

	split.screen( rbind(c(0.08,0.95,0.1,0.95), c(0.95,0.99,0.1,0.95)))
	split.screen( figs = c( 2, 5 ), screen = 1 )

	hr_count <- 0
	hr_label <- c('0:00', '3:00', '6:00', '9:00', '12:00', '15:00', '18:00', '21:00')
	for(hr in hr_index){
		
		hr_count <- hr_count + 1
		
		for(variable in 1:2){
			
			screen((variable - 1) * 5 + 2 + hr_count)
			par(pty = 's')
			par(mai=c(0.2,0.2,0.2,0.2))
			
			if(hr_count == 1 & variable == 2){
			quilt.plot(dat3[, 1], dat3[, 2], dat2[hr, ], zlim = zlim_range2, nx = 25, ny = 25, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 2)
			points(cbind(dat3[320, 1], dat3[320, 2]), col = 'black', pch = 4, cex = 4, lwd = 4)
			mtext('985 hPa', side = 2, line = 7, adj = 0.5, cex = 3, font = 2, col = 'blue')
			}else if(hr_count == 1 & variable == 1){
			quilt.plot(dat3[, 1], dat3[, 2], dat[hr, ], zlim = zlim_range1, nx = 25, ny = 25, ylab = '', xlab = '', xaxt = 'n', cex.lab = 4, add.legend = F, cex.axis = 2)
			points(cbind(dat3[412, 1], dat3[412, 2]), col = 'black', pch = 4, cex = 4, lwd = 4)
			mtext('850 hPa', side = 2, line = 7, adj = 0.5, cex = 3, font = 2, col = 'blue')
			}else if(variable == 2){
			quilt.plot(dat3[, 1], dat3[, 2], dat2[hr, ], zlim = zlim_range2, nx = 25, ny = 25, ylab = '', xlab = '', yaxt = 'n', cex.lab = 4, add.legend = F, cex.axis = 2)
			points(cbind(dat3[320, 1], dat3[320, 2]), col = 'black', pch = 4, cex = 4, lwd = 4)
			}else{
			quilt.plot(dat3[, 1], dat3[, 2], dat[hr, ], zlim = zlim_range1, nx = 25, ny = 25, ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', cex.lab = 4, add.legend = F, cex.axis = 2)
			points(cbind(dat3[412, 1], dat3[412, 2]), col = 'black', pch = 4, cex = 4, lwd = 4)
			}
			map("worldHires", xlim = c(26.719, 85.078), ylim = c(5.625, 42.188), lwd = 0.75, add = T)
			
			if(hr_count == 1){
				mtext('Latitude', side = 2, line = 4, adj = 0.5, cex = 2.5, font = 2)
			}

			if(variable == 1){
				if(mod(hr, 8) == 0){
					mtext(hr_label[1], side = 3, line = 1, adj = 0.5, cex = 3, font = 2)
				}else{
					mtext(hr_label[mod(hr, 8)], side = 3, line = 1, adj = 0.5, cex = 3, font = 2)
				}
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

for(start_hr in 130:130){

	cat(start_hr, '\n')

	hr_index <- seq(start_hr, start_hr + 4, by = 1)

	zlim_range1 <- zlim_range2 <- c(-1.9, 1.9)
	#zlim_range1 <- range(Yhat1[start_hr:(start_hr + 4), ])
	#zlim_range2 <- range(Yhat2[start_hr:(start_hr + 4), ])

	split.screen( rbind(c(0.08,0.95,0.1,0.95), c(0.95,0.99,0.1,0.95)))
	split.screen( figs = c( 2, 5 ), screen = 1 )

	hr_count <- 0
	hr_label <- 0
	for(hr in hr_index){
		
		if(mod(hr - 1, 8) == 0){
			hr_label <- 0
		}else{
			hr_label <- hr_label + 1
		}

		if(mod(hr - 1, 8) == 0){
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
			mtext('985 hPa', side = 2, line = 7, adj = 0.5, cex = 3, font = 2, col = 'blue')
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
				mtext(paste(hr_label * 3, ':00', sep = ''), side = 3, line = 1, adj = 0.5, cex = 3, font = 2)
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


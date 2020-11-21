
directory <- '/home/salvanmo/Desktop/'

root <- paste(directory, 'Spatio-Temporal-Cross-Covariance-Functions-under-the-Lagrangian-Framework/', sep = '')

source(file = paste(root, "Functions/load_packages.R",sep=''))
source(file = paste(root, "Functions/auxiliary_functions.R",sep=''))

#jpeg(file = paste(root, 'Figures/spacetime-maps.jpg', sep = ''), width = 1600, height = 700)
pdf(file = paste(root, 'Figures/spacetime-maps.pdf', sep = ''), width = 25, height = 10)

yr <- 2018

for(start_hr in seq(1, 200, by = 5)){

	cat(start_hr, '\n')
	#cat(yr, '\n')

	hr_index <- seq(start_hr, start_hr + 4, by = 1)

	dat <- read.table(paste(root, 'Data/ncdf/layer1_residuals_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
	dat2 <- read.table(paste(root, 'Data/ncdf/layer2_residuals_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
	#dat <- read.table(paste(root, 'Data/ncdf/layer1_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
	#dat2 <- read.table(paste(root, 'Data/ncdf/layer2_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
	dat3 <- read.table(paste(root, 'Data/ncdf/LOCS-3D-dataset', sep = ''), header = FALSE, sep = " ") %>% as.matrix()

	zlim_range1 <- range(dat[start_hr:(start_hr + 4),])
	zlim_range2 <- range(dat2[start_hr:(start_hr + 4),])

	split.screen( rbind(c(0.08,0.95,0.1,0.95), c(0.95,0.99,0.1,0.95)))
	split.screen( figs = c( 2, 5 ), screen = 1 )

	hr_count <- 0
	for(hr in hr_index){
		
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
				mtext(paste((hr - 1) * 8, ':00', sep = ''), side = 3, line = 1, adj = 0.5, cex = 3, font = 2)
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


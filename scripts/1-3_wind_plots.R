
directory <- '/home/salvanmo/Desktop/'

root <- paste(directory, 'Spatio-Temporal-Cross-Covariance-Functions-under-the-Lagrangian-Framework/', sep = '')

source(file = paste(root, "Functions/load_packages.R",sep=''))
source(file = paste(root, "Functions/auxiliary_functions.R",sep=''))

saudi<- map("world", "Saudi", fill = TRUE)
IDs <- sapply(strsplit(saudi$names, ":"), function(x) x[1])
saudi <- map2SpatialPolygons(saudi, IDs=IDs, proj4string=CRS("+proj=longlat +datum=WGS84"))

WIND_MAT <- matrix(, ncol = 4, nrow = 40)

for (yr in 1980:2019){

	wind_U_component <- wind_V_component <- list()

	for(VAR in 1:2){
		
		data_array1 <- data_array2 <- data_array3 <- NULL

		if(yr < 1992){
			merra_ind <- 100
		}else if(yr >= 1992 & yr < 2001){
			merra_ind <- 200
		}else if (yr >= 2001 & yr < 2011){
			merra_ind <- 300
		}else{
			merra_ind <- 400
		}

		mnth = 1
		if(mnth == 2){
			mnth_end <- 28
		}else if(mnth %in% c(1, 3, 5, 7, 8, 10, 12)){
			mnth_end <- 31
		}else{
			mnth_end <- 30
		}

		if(mnth < 10){
			mo <- paste("0", mnth, sep='')
		}else{
			mo <- mnth
		}

		ncname <- paste("/home/salvanmo/Downloads/MERRA2_", merra_ind, ".inst3_3d_asm_Nv.", yr, "0101.SUB.nc", sep='')   
		#ncname <- paste("/home/salvanmo/Downloads/MERRA/wind/MERRA2_", merra_ind, ".inst3_3d_asm_Nv.", yr, "0101.SUB.nc", sep='')   

		ncin <- nc_open(ncname)
		
		dname1 <- "U"
		dname2 <- "V"

		u_array <- ncvar_get(ncin,dname1)
		v_array <- ncvar_get(ncin,dname2)

		# get longitude and latitude
		lon <- ncvar_get(ncin,"lon")
		lat <- ncvar_get(ncin,"lat")

		nc_close(ncin)

		if(VAR == 1)	lev <- 65	else	lev <- 68

		U <- u_array[,, lev, ]
		V <- v_array[,, lev, ]

		lon.lat <- expand.grid(lon,lat)
		lon_new <- matrix(lon.lat[, 1], ncol = length(lat))
		lat_new <- matrix(lon.lat[, 2], ncol = length(lat))

		test1 <- data.frame(rep(lon.lat[,1], 8), rep(lon.lat[,2], 8),  c(U),  c(V))

		for(day in 2:31){

			cat('READING NETCDF DATA ===> year: ', yr, 'month: ', mnth, 'day: ', day, '\n')
			if(day > 9){
				ncname <- paste("/home/salvanmo/Downloads/MERRA2_", merra_ind, ".inst3_3d_asm_Nv.", yr, mo, day,".SUB.nc", sep='')
				#ncname <- paste("/home/salvanmo/Downloads/MERRA/wind/MERRA2_", merra_ind, ".inst3_3d_asm_Nv.", yr, mo, day,".SUB.nc", sep='')
			}else{
				ncname <- paste("/home/salvanmo/Downloads/MERRA2_", merra_ind, ".inst3_3d_asm_Nv.", yr, mo, "0",day,".SUB.nc", sep='')
				#ncname <- paste("/home/salvanmo/Downloads/MERRA/wind/MERRA2_", merra_ind, ".inst3_3d_asm_Nv.", yr, mo, "0",day,".SUB.nc", sep='')
			}
			ncin <- nc_open(ncname)

			u_array <- ncvar_get(ncin,dname1)
			v_array <- ncvar_get(ncin,dname2)

			nc_close(ncin)

			U <- u_array[,, lev, ]
			V <- v_array[,, lev, ]

			test1 <- rbind(test1, data.frame(rep(lon.lat[,1], 8), rep(lon.lat[,2], 8),  c(U),  c(V)))
		}

		colnames(test1) <- c('lon', 'lat', 'Y1', 'Y2')
		spdf <- SpatialPointsDataFrame(coords = test1[, c("lon", "lat")], data = test1, proj4string = CRS("+proj=longlat +datum=WGS84"))
		saudi_data_orig <- data.frame(spdf[!is.na(over(spdf, as(saudi, "SpatialPolygons"))), ])

		N <- nrow(saudi_data_orig)/(8 * mnth_end)

		data_temp <- matrix(saudi_data_orig[, 3], ncol = N, byrow = T)
		data_array1 <- rbind(data_array1, data_temp)

		data_temp <- matrix(saudi_data_orig[, 4], ncol = N, byrow = T)
		data_array2 <- rbind(data_array2, data_temp)

		wind_U_component[[VAR]] <- data_array1
		wind_V_component[[VAR]] <- data_array2

	}
	WIND_MAT[yr - 1979, 1] <- mean(wind_U_component[[1]][140:144, ])
	WIND_MAT[yr - 1979, 2] <- mean(wind_V_component[[1]][140:144, ])
	WIND_MAT[yr - 1979, 3] <- mean(wind_U_component[[2]][140:144, ])
	WIND_MAT[yr - 1979, 4] <- mean(wind_V_component[[2]][140:144, ])
}



pdf(file = paste('/home/salvanmo/Desktop/thesis/thesis-defense/Figures/wind_data.pdf', sep = ''), width = 15, height = 15)

split.screen( rbind(c(0.06, 0.98, 0.05, 0.98), c(0.99,0.99,0.05,0.95)))
split.screen( figs = c( 3, 3 ), screen = 1 )

screen(3)
par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))

plot(WIND_MAT[, 1:2], xlab = "", ylab = "", cex.axis = 1.5, pch = 20, cex = 1, col = "#808080", ylim = range(WIND_MAT[, c(2, 4)]), xlim = range(WIND_MAT[, c(1, 3)]))
abline(h = 0, v = 0, col = 1, lwd = 2, lty = 2)
mtext(bquote(paste(V[y], " (880 hPa) ")), side = 3, line = 0, adj = 0.5, cex = 2)
mtext(bquote(paste(V[x], " (880 hPa) ")) , side = 2, line = 3, adj = 0.5, cex = 2)

screen(11)
par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))

plot(WIND_MAT[, 3:4], xlab = "", ylab = "", cex.axis = 1.5, pch = 20, cex = 1, col = "#808080", ylim = range(WIND_MAT[, c(2, 4)]), xlim = range(WIND_MAT[, c(1, 3)]))
abline(h = 0, v = 0, col = 1, lwd = 2, lty = 2)
mtext(bquote(paste(V[x], " (925 hPa) ")) , side = 2, line = 3, adj = 0.5, cex = 2)

screen(4)
par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))

plot(WIND_MAT[, c(1, 3)], xlab = "", ylab = "", cex.axis = 1.5, pch = 20, cex = 1, col = "#808080", ylim = range(WIND_MAT[, c(2, 4)]), xlim = range(WIND_MAT[, c(1, 3)]), yaxt = 'n', xaxt = 'n')
abline(h = 0, v = 0, col = 1, lwd = 2, lty = 2)
mtext(bquote(paste(V[x], " (925 hPa) ")), side = 3, line = 0, adj = 0.5, cex = 2)
  
screen(8)
par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))

plot(WIND_MAT[, c(2, 4)], xlab = "", ylab = "", cex.axis = 1.5, pch = 20, cex = 1, col = "#808080", ylim = range(WIND_MAT[, c(2, 4)]), xlim = range(WIND_MAT[, c(1, 3)]), yaxt = 'n', xaxt = 'n')
abline(h = 0, v = 0, col = 1, lwd = 2, lty = 2)

screen(5)
par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))

plot(WIND_MAT[, c(1, 4)], xlab = "", ylab = "", cex.axis = 1.5, pch = 20, cex = 1, col = "#808080", ylim = range(WIND_MAT[, c(2, 4)]), xlim = range(WIND_MAT[, c(1, 3)]), yaxt = 'n', xaxt = 'n')
abline(h = 0, v = 0, col = 1, lwd = 2, lty = 2)
mtext(bquote(paste(V[y], " (925 hPa) ")), side = 3, line = 0, adj = 0.5, cex = 2)

screen(7)
par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))

plot(WIND_MAT[, c(2, 3)], xlab = "", ylab = "", cex.axis = 1.5, pch = 20, cex = 1, col = "#808080", ylim = range(WIND_MAT[, c(2, 4)]), xlim = range(WIND_MAT[, c(1, 3)]))
abline(h = 0, v = 0, col = 1, lwd = 2, lty = 2)
mtext(bquote(paste(V[y], " (880 hPa) ")) , side = 2, line = 3, adj = 0.5, cex = 2)

close.screen( all=TRUE)
dev.off()



pdf(file = paste('/home/salvanmo/Desktop/thesis/thesis-defense/Figures/wind_data_OLD.pdf', sep = ''), width = 11, height = 6)

split.screen( rbind(c(0.08, 0.52, 0.03, 0.98), c(0.55, 0.99, 0.03, 0.98), c(0.99,0.99,0.05,0.95)))

screen(1)
par(pty = 's')
par(mai=c(0.3,0.3,0.3,0.3))

plot(cbind(c(wind_U_component[[1]][141:145, ]), c(wind_V_component[[1]][141:145, ])), xlab = "", ylab = "", cex.axis = 1.5, pch = 20, cex = 0.5, col = "#808080", ylim = c(-10, 10), xlim = c(-10, 10))
abline(h = 0, v = 0, col = 1, lwd = 2, lty = 2)
mtext(expression(V[y]), side = 2, line = 3, adj = 0.5, cex = 2)
mtext(expression(V[x]), side = 1, line = 3, adj = 0.5, cex = 2)
mtext('850 hPa', side = 3, line = 0, adj = 0.5, cex = 2, font = 2, col = 'blue')
points(matrix(c(0.2485011, -0.6025650), ncol = 2), pch = 8, lwd = 4, cex = 2, col = "#1159b9")
points(matrix(c(0.1221316, -1.6079946), ncol = 2), pch = 8, lwd = 4, cex = 2, col = "#e54457")
points(matrix(c(0.81010468, -0.09895062), ncol = 2), pch = 8, lwd = 4, cex = 2, col = "#32CD32")
#points(matrix(c(0.008669427, 0.049909815), ncol = 2), pch = 8, lwd = 4, cex = 2, col = "#9400D3")
#points(matrix(c(0.248501062, -0.602564994), ncol = 2), pch = 8, lwd = 4, cex = 2, col = "#9400D3")

screen(2)
par(pty = 's')
par(mai=c(0.3,0.3,0.3,0.3))

plot(cbind(c(wind_U_component[[2]][141:145, ]), c(wind_V_component[[2]][141:145, ])), xlab = "", ylab = "", cex.axis = 1.5, pch = 20, cex = 0.5, col = "#808080", ylim = c(-10, 10), xlim = c(-10, 10), yaxt = 'n')
abline(h = 0, v = 0, col = 1, lwd = 2, lty = 2)
mtext('985 hPa', side = 3, line = 0, adj = 0.5, cex = 2, font = 2, col = 'blue')
mtext(expression(V[x]), side = 1, line = 3, adj = 0.5, cex = 2)
points(matrix(c(-0.3167604, -2.3540665), ncol = 2), pch = 8, lwd = 4, cex = 2, col = "#1159b9")
points(matrix(c(0.1221316, -1.6079946), ncol = 2), pch = 8, lwd = 4, cex = 2, col = "#e54457")
points(matrix(c(-0.72022097, -0.84885145), ncol = 2), pch = 8, lwd = 4, cex = 2, col = "#32CD32")
#points(matrix(c(0.008669427, 0.049909815), ncol = 2), pch = 8, lwd = 4, cex = 2, col = "#9400D3")
#points(matrix(c(0.248501062, -0.602564994), ncol = 2), pch = 8, lwd = 4, cex = 2, col = "#9400D3")
#ellipse(mu = c(-2, -2), sigma = 2 * diag(2), alpha = .05, npoints = 250, col = "#1159b9", lwd = 2) 
  
close.screen( all=TRUE)
dev.off()


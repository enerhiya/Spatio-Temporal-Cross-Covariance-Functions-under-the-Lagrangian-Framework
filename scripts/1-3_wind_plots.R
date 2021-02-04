
directory <- '/home/salvanmo/Desktop/'

root <- paste(directory, 'Spatio-Temporal-Cross-Covariance-Functions-under-the-Lagrangian-Framework/', sep = '')

source(file = paste(root, "Functions/load_packages.R",sep=''))
source(file = paste(root, "Functions/auxiliary_functions.R",sep=''))

saudi<- map("world", "Saudi", fill = TRUE)
IDs <- sapply(strsplit(saudi$names, ":"), function(x) x[1])
saudi <- map2SpatialPolygons(saudi, IDs=IDs, proj4string=CRS("+proj=longlat +datum=WGS84"))

WIND_ARRAY <- array(, dim = c(550 * 5, 4, 20))

for (yr in 2000:2019){

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

	WIND_ARRAY[, 1, yr - 1999] <- c(wind_U_component[[1]][140:144, ])
	WIND_ARRAY[, 2, yr - 1999] <- c(wind_V_component[[1]][140:144, ])
	WIND_ARRAY[, 3, yr - 1999] <- c(wind_U_component[[2]][140:144, ])
	WIND_ARRAY[, 4, yr - 1999] <- c(wind_V_component[[2]][140:144, ])
}

test <- round(cbind(c(WIND_ARRAY[, 2, ]), c(WIND_ARRAY[, 1, ])), 2)

#pdf(file = paste('/home/salvanmo/Desktop/thesis/thesis-defense/Figures/wind_data.pdf', sep = ''), width = 15, height = 15)
jpeg(file = paste('/home/salvanmo/Desktop/thesis/thesis-defense/Figures/wind_data.jpg', sep = ''), width = 1500, height = 1500)

split.screen( rbind(c(0.06, 0.98, 0.05, 0.98), c(0.99,0.99,0.05,0.95)))
split.screen( figs = c( 3, 3 ), screen = 1 )

screen(3)
par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))

plot(unique(round(cbind(c(WIND_ARRAY[, 2, ]), c(WIND_ARRAY[, 1, ])), 0)), xlab = "", ylab = "", cex.axis = 1.5, pch = 20, cex = 0.5, col = "#808080", ylim = c(-max(WIND_ARRAY), max(WIND_ARRAY)), xlim = c(-max(WIND_ARRAY), max(WIND_ARRAY)))
abline(h = 0, v = 0, col = 1, lwd = 2, lty = 2)
mtext(bquote(paste(v[y], " (880 hPa) ")), side = 3, line = 0, adj = 0.5, cex = 2)
mtext(bquote(paste(v[x], " (880 hPa) ")) , side = 2, line = 3, adj = 0.5, cex = 2)
#ellipse(mu = c(0.13178367, 0.00721022), sigma = matrix(c(27.160213, 4.917272, 4.917272, 5.554407), 2, 2), alpha = .05, npoints = 250, col = "#3fb1e2", lwd = 4) 
ellipse(mu = c(mu1_ms[2], mu1_ms[1]), sigma = wind_var1[c(2,1), c(2,1)], alpha = .05, npoints = 250, col = "#3fb1e2", lwd = 2) 
points(matrix(c(mu1_ms[2], mu1_ms[1]), ncol = 2), pch = 8, lwd = 4, cex = 2, col = "#3fb1e2")

screen(11)
par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))

plot(unique(round(cbind(c(WIND_ARRAY[, 4, ]), c(WIND_ARRAY[, 3, ])), 0)), xlab = "", ylab = "", cex.axis = 1.5, pch = 20, cex = 0.5, col = "#808080", ylim = c(-max(WIND_ARRAY), max(WIND_ARRAY)), xlim = c(-max(WIND_ARRAY), max(WIND_ARRAY)))
abline(h = 0, v = 0, col = 1, lwd = 2, lty = 2)
mtext(bquote(paste(v[x], " (925 hPa) ")) , side = 2, line = 3, adj = 0.5, cex = 2)
#ellipse(mu = c(0.139535500, -0.006487013), sigma = matrix(c(12.8859287, -0.2185663, -0.2185663, 4.4855505), 2, 2), alpha = .05, npoints = 250, col = "#3fb1e2", lwd = 4) 
#ellipse(mu = colMeans(cbind(c(WIND_ARRAY[, 4, ]), c(WIND_ARRAY[, 3, ]))), sigma = var(cbind(c(WIND_ARRAY[, 4, ]), c(WIND_ARRAY[, 3, ]))), alpha = .05, npoints = 250, col = "#3fb1e2", lwd = 4) 
ellipse(mu = c(mu2_ms[2], mu2_ms[1]), sigma = wind_var2[c(2,1), c(2,1)], alpha = .05, npoints = 250, col = "#3fb1e2", lwd = 2) 
points(matrix(c(mu2_ms[2], mu2_ms[1]), ncol = 2), pch = 8, lwd = 4, cex = 2, col = "#3fb1e2")

screen(4)
par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))

plot(unique(round(cbind(c(WIND_ARRAY[, 3, ]), c(WIND_ARRAY[, 1, ])), 0)), xlab = "", ylab = "", cex.axis = 1.5, pch = 20, cex = 0.5, col = "#808080", ylim = c(-max(WIND_ARRAY), max(WIND_ARRAY)), xlim = c(-max(WIND_ARRAY), max(WIND_ARRAY)), yaxt = 'n', xaxt = 'n')
abline(h = 0, v = 0, col = 1, lwd = 2, lty = 2)
mtext(bquote(paste(v[x], " (925 hPa) ")), side = 3, line = 0, adj = 0.5, cex = 2)
#ellipse(mu = colMeans(cbind(c(WIND_ARRAY[, 3, ]), c(WIND_ARRAY[, 1, ]))), sigma = var(cbind(c(WIND_ARRAY[, 3, ]), c(WIND_ARRAY[, 1, ]))), alpha = .05, npoints = 250, col = "#3fb1e2", lwd = 4) 
#ellipse(mu = c(mu2_ms[1], mu1_ms[1]), sigma = wind_var_ms2[c(3, 1), c(3, 1)], alpha = .05, npoints = 250, col = "#3fb1e2", lwd = 2) 
points(matrix(c(mu2_ms[1], mu1_ms[1]), ncol = 2), pch = 8, lwd = 4, cex = 2, col = "#3fb1e2")
  
screen(8)
par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))

plot(unique(round(cbind(c(WIND_ARRAY[, 4, ]), c(WIND_ARRAY[, 2, ])), 0)), xlab = "", ylab = "", cex.axis = 1.5, pch = 20, cex = 0.5, col = "#808080", ylim = c(-max(WIND_ARRAY), max(WIND_ARRAY)), xlim = c(-max(WIND_ARRAY), max(WIND_ARRAY)), yaxt = 'n', xaxt = 'n')
abline(h = 0, v = 0, col = 1, lwd = 2, lty = 2)
#ellipse(mu = colMeans(cbind(c(WIND_ARRAY[, 4, ]), c(WIND_ARRAY[, 2, ]))), sigma = var(cbind(c(WIND_ARRAY[, 4, ]), c(WIND_ARRAY[, 2, ]))), alpha = .05, npoints = 250, col = "#3fb1e2", lwd = 4) 
#ellipse(mu = c(mu2_ms[2], mu1_ms[2]), sigma = wind_var_ms2[c(4, 2), c(4, 2)], alpha = .05, npoints = 250, col = "#3fb1e2", lwd = 2) 
points(matrix(c(mu2_ms[2], mu1_ms[2]), ncol = 2), pch = 8, lwd = 4, cex = 2, col = "#3fb1e2")

screen(5)
par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))

plot(unique(round(cbind(c(WIND_ARRAY[, 4, ]), c(WIND_ARRAY[, 1, ])), 0)), xlab = "", ylab = "", cex.axis = 1.5, pch = 20, cex = 0.5, col = "#808080", ylim = c(-max(WIND_ARRAY), max(WIND_ARRAY)), xlim = c(-max(WIND_ARRAY), max(WIND_ARRAY)), yaxt = 'n', xaxt = 'n')
abline(h = 0, v = 0, col = 1, lwd = 2, lty = 2)
mtext(bquote(paste(v[y], " (925 hPa) ")), side = 3, line = 0, adj = 0.5, cex = 2)
#ellipse(mu = colMeans(cbind(c(WIND_ARRAY[, 4, ]), c(WIND_ARRAY[, 1, ]))), sigma = var(cbind(c(WIND_ARRAY[, 4, ]), c(WIND_ARRAY[, 1, ]))), alpha = .05, npoints = 250, col = "#3fb1e2", lwd = 4) 
#ellipse(mu = c(mu2_ms[2], mu1_ms[1]), sigma = wind_var_ms2[c(4, 1), c(4, 1)], alpha = .05, npoints = 250, col = "#3fb1e2", lwd = 2) 
points(matrix(c(mu2_ms[2], mu1_ms[1]), ncol = 2), pch = 8, lwd = 4, cex = 2, col = "#3fb1e2")

screen(7)
par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))

plot(unique(round(cbind(c(WIND_ARRAY[, 3, ]), c(WIND_ARRAY[, 2, ])), 0)), xlab = "", ylab = "", cex.axis = 1.5, pch = 20, cex = 0.5, col = "#808080", ylim = c(-max(WIND_ARRAY), max(WIND_ARRAY)), xlim = c(-max(WIND_ARRAY), max(WIND_ARRAY)))
abline(h = 0, v = 0, col = 1, lwd = 2, lty = 2)
mtext(bquote(paste(v[y], " (880 hPa) ")) , side = 2, line = 3, adj = 0.5, cex = 2)
#ellipse(mu = colMeans(cbind(c(WIND_ARRAY[, 3, ]), c(WIND_ARRAY[, 2, ]))), sigma = var(cbind(c(WIND_ARRAY[, 3, ]), c(WIND_ARRAY[, 2, ]))), alpha = .05, npoints = 250, col = "#3fb1e2", lwd = 4) 
#ellipse(mu = c(mu2_ms[1], mu1_ms[2]), sigma = wind_var_ms2[c(3, 2), c(3, 2)], alpha = .05, npoints = 250, col = "#3fb1e2", lwd = 2) 
points(matrix(c(mu2_ms[1], mu1_ms[2]), ncol = 2), pch = 8, lwd = 4, cex = 2, col = "#3fb1e2")

close.screen( all=TRUE)
dev.off()



pdf(file = paste('/home/salvanmo/Desktop/thesis/thesis-defense/Figures/wind_data_OLDOLD.pdf', sep = ''), width = 15, height = 15)

split.screen( rbind(c(0.06, 0.98, 0.05, 0.98), c(0.99,0.99,0.05,0.95)))
split.screen( figs = c( 3, 3 ), screen = 1 )

screen(3)
par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))

plot(WIND_MAT[, c(2, 1)], xlab = "", ylab = "", cex.axis = 1.5, pch = 20, cex = 1, col = "#808080", ylim = c(-max(WIND_MAT), max(WIND_MAT)), xlim = c(-max(WIND_MAT), max(WIND_MAT)))
abline(h = 0, v = 0, col = 1, lwd = 2, lty = 2)
mtext(bquote(paste(v[y], " (880 hPa) ")), side = 3, line = 0, adj = 0.5, cex = 2)
mtext(bquote(paste(v[x], " (880 hPa) ")) , side = 2, line = 3, adj = 0.5, cex = 2)
ellipse(mu = c(mu1_ms[2], mu1_ms[1]), sigma = wind_var_ms2[c(2,1), c(2,1)], alpha = .05, npoints = 250, col = "#3fb1e2", lwd = 2) 
points(matrix(c(mu1_ms[2], mu1_ms[1]), ncol = 2), pch = 8, lwd = 4, cex = 2, col = "#3fb1e2")

screen(11)
par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))

plot(WIND_MAT[, c(4, 3)], xlab = "", ylab = "", cex.axis = 1.5, pch = 20, cex = 1, col = "#808080", ylim = c(-max(WIND_MAT), max(WIND_MAT)), xlim = c(-max(WIND_MAT), max(WIND_MAT)))
abline(h = 0, v = 0, col = 1, lwd = 2, lty = 2)
mtext(bquote(paste(v[x], " (925 hPa) ")) , side = 2, line = 3, adj = 0.5, cex = 2)
ellipse(mu = c(mu2_ms[2], mu2_ms[1]), sigma = wind_var_ms2[c(4, 3), c(4, 3)], alpha = .05, npoints = 250, col = "#3fb1e2", lwd = 2) 
points(matrix(c(mu2_ms[2], mu2_ms[1]), ncol = 2), pch = 8, lwd = 4, cex = 2, col = "#3fb1e2")

screen(4)
par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))

plot(WIND_MAT[, c(3, 1)], xlab = "", ylab = "", cex.axis = 1.5, pch = 20, cex = 1, col = "#808080", ylim = c(-max(WIND_MAT), max(WIND_MAT)), xlim = c(-max(WIND_MAT), max(WIND_MAT)), yaxt = 'n', xaxt = 'n')
abline(h = 0, v = 0, col = 1, lwd = 2, lty = 2)
mtext(bquote(paste(v[x], " (925 hPa) ")), side = 3, line = 0, adj = 0.5, cex = 2)
ellipse(mu = c(mu2_ms[1], mu1_ms[1]), sigma = wind_var_ms2[c(3, 1), c(3, 1)], alpha = .05, npoints = 250, col = "#3fb1e2", lwd = 2) 
points(matrix(c(mu2_ms[1], mu1_ms[1]), ncol = 2), pch = 8, lwd = 4, cex = 2, col = "#3fb1e2")
  
screen(8)
par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))

plot(WIND_MAT[, c(4, 2)], xlab = "", ylab = "", cex.axis = 1.5, pch = 20, cex = 1, col = "#808080", ylim = c(-max(WIND_MAT), max(WIND_MAT)), xlim = c(-max(WIND_MAT), max(WIND_MAT)), yaxt = 'n', xaxt = 'n')
abline(h = 0, v = 0, col = 1, lwd = 2, lty = 2)
ellipse(mu = c(mu2_ms[2], mu1_ms[2]), sigma = wind_var_ms2[c(4, 2), c(4, 2)], alpha = .05, npoints = 250, col = "#3fb1e2", lwd = 2) 
points(matrix(c(mu2_ms[2], mu1_ms[2]), ncol = 2), pch = 8, lwd = 4, cex = 2, col = "#3fb1e2")

screen(5)
par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))

plot(WIND_MAT[, c(4, 1)], xlab = "", ylab = "", cex.axis = 1.5, pch = 20, cex = 1, col = "#808080", ylim = c(-max(WIND_MAT), max(WIND_MAT)), xlim = c(-max(WIND_MAT), max(WIND_MAT)), yaxt = 'n', xaxt = 'n')
abline(h = 0, v = 0, col = 1, lwd = 2, lty = 2)
mtext(bquote(paste(v[y], " (925 hPa) ")), side = 3, line = 0, adj = 0.5, cex = 2)
ellipse(mu = c(mu2_ms[2], mu1_ms[1]), sigma = wind_var_ms2[c(4, 1), c(4, 1)], alpha = .05, npoints = 250, col = "#3fb1e2", lwd = 2) 
points(matrix(c(mu2_ms[2], mu1_ms[1]), ncol = 2), pch = 8, lwd = 4, cex = 2, col = "#3fb1e2")

screen(7)
par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))

plot(WIND_MAT[, c(3, 2)], xlab = "", ylab = "", cex.axis = 1.5, pch = 20, cex = 1, col = "#808080", ylim = c(-max(WIND_MAT), max(WIND_MAT)), xlim = c(-max(WIND_MAT), max(WIND_MAT)))
abline(h = 0, v = 0, col = 1, lwd = 2, lty = 2)
mtext(bquote(paste(v[y], " (880 hPa) ")) , side = 2, line = 3, adj = 0.5, cex = 2)
ellipse(mu = c(mu2_ms[1], mu1_ms[2]), sigma = wind_var_ms2[c(3, 2), c(3, 2)], alpha = .05, npoints = 250, col = "#3fb1e2", lwd = 2) 
points(matrix(c(mu2_ms[1], mu1_ms[2]), ncol = 2), pch = 8, lwd = 4, cex = 2, col = "#3fb1e2")

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





####### CONVERTING ADVECTION VELOCITY FROM DEG/3HRS TO KM/3HRS

mu1 <- c( -0.01173841, -0.02089627)
mu2 <- c(0.02873192, -0.03053685)

wind_var <- matrix(c(0.044225813, 0.003904809, 0.05286808, 0.01294508, 0.003904809, 0.011216205, 0.01922954, 0.02109589, 0.052868084, 0.019229541, 0.09889986, 0.04511405, 0.012945081, 0.021095889, 0.04511405, 0.06660899), 4, 4)

####### CONVERTING ADVECTION VELOCITY FROM KM/3HRS TO M/S

est_vals <- c(median(p_vals[, 7]), median(p_vals[, 8]), median(p_vals[, 10]), median(p_vals[, 11]))
est_vals <- c(median(p_vals[, 7]), median(p_vals[, 8]), median(p_vals[, 9]), median(p_vals[, 10]))
est_vals <- c(median(p_vals[, 5]), median(p_vals[, 6]), median(p_vals[, 7]), median(p_vals[, 8]))
(est_vals / 3) * (0.2778)

#reverse the sign of the following numbers
#M2: -0.2485011  0.6025650  0.3167604  2.3540665
#M3: -0.81010468  0.09895062  0.72022097  0.84885145
#M4: 0.008669427  0.049909815 -0.248501062  0.602564994

####### CONVERTING ADVECTION VARIANCES FROM (KM/3HRS)^2 TO (M/S)^2

set=85
p <- read.table(paste(root, 'Results/multivariate_stationary_real_data_parameter_estimates_M2_set_', set, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
wind_var_chol <- matrix(c(p[11], p[12], p[13], p[14], 0, p[15], p[16], p[17], 0, 0, p[18], p[19], 0, 0, 0, p[20]), ncol = 4, byrow = T)
wind_var <- t(wind_var_chol) %*% wind_var_chol

Sig11 <- wind_var[1:2, 1:2]
Sig12 <- wind_var[1:2, 3:4]
Sig22 <- wind_var[3:4, 3:4]
cond_cov <- Sig11 - Sig12 %*% solve(Sig22) %*% t(Sig12)

est_vals <- c(cond_cov[1, 1], cond_cov[1, 2], cond_cov[2, 2])
(est_vals / 3^2) * 0.07716049

#########################################################

mu1 <- c(0.01222045, -0.07345084)
mu2 <- c(0.06025701, 0.07327084)

mu1 <- c(0.04082589, 0.06223707)
mu2 <- c(-0.02917227, 0.04565329)

mu1 <- c(-0.073211462, -0.005373963)
mu2 <- c(0.02514516, 0.01563370)

mu1_ms <- (mu1 * 100 / 3) * (0.2778)
mu2_ms <- (mu2 * 100 / 3) * (0.2778) 


wind_var_ms2 <- (wind_var * 10000 / 3^2) * 0.07716049


wind_var_ms2 <- matrix(c(5.477262,  2.517276, 5.274911, 2.486674, 2.517276, 10.220525, 2.186165, 10.11981, 5.274911, 2.186165, 5.200563,  2.164414, 2.486674, 10.11981, 2.164414, 10.196487), 4, 4)

p <- read.table(paste(root, 'Results/multivariate_stationary_real_M1_A_FULL', sep = ''), header = FALSE, sep = " ") %>% as.matrix()

mu1 <- p[9:10]
mu2 <- p[11:12]

mu1_ms <- (mu1 * 100 / 3) * (0.2778)
mu2_ms <- (mu2 * 100 / 3) * (0.2778) 

wind_var_chol <- matrix(c(p[13], p[14], 0, p[15]), ncol = 2, byrow = T)
wind_var1 <- t(wind_var_chol) %*% wind_var_chol

wind_var_chol <- matrix(c(p[16], p[17], 0, p[18]), ncol = 2, byrow = T)
wind_var2 <- t(wind_var_chol) %*% wind_var_chol




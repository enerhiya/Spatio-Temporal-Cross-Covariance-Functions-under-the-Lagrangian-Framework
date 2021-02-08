
directory <- '/home/salvanmo/Desktop/'

root <- paste(directory, 'Spatio-Temporal-Cross-Covariance-Functions-under-the-Lagrangian-Framework/', sep = '')

source(file = paste(root, "Functions/load_packages.R",sep=''))
source(file = paste(root, "Functions/auxiliary_functions.R",sep=''))

saudi<- map("world", "Saudi", fill = TRUE)
IDs <- sapply(strsplit(saudi$names, ":"), function(x) x[1])
saudi <- map2SpatialPolygons(saudi, IDs=IDs, proj4string=CRS("+proj=longlat +datum=WGS84"))

WIND_ARRAY <- array(, dim = c(550 * 5, 4, 40))

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

	WIND_ARRAY[, 1, yr - 1999] <- c(wind_U_component[[1]][140:144, ])
	WIND_ARRAY[, 2, yr - 1999] <- c(wind_V_component[[1]][140:144, ])
	WIND_ARRAY[, 3, yr - 1999] <- c(wind_U_component[[2]][140:144, ])
	WIND_ARRAY[, 4, yr - 1999] <- c(wind_V_component[[2]][140:144, ])
}

pdf(file = paste('/home/salvanmo/Desktop/thesis/thesis-defense/Figures/wind_data_NEW.pdf', sep = ''), width = 20, height = 20)
#jpeg(file = paste('/home/salvanmo/Desktop/thesis/thesis-defense/Figures/wind_data_NEW.jpg', sep = ''), width = 1500, height = 1500)

split.screen( rbind(c(0.06, 0.98, 0.05, 0.98), c(0.99,0.99,0.05,0.95)))
split.screen( figs = c( 4, 4 ), screen = 1 )

xvals <- c(-max(WIND_ARRAY), max(WIND_ARRAY))

screen(3)
par(pty = 's')
par(mai=c(0.4, 0.4, 0.4, 0.8))

hist(unique(c(WIND_ARRAY[, 1, ])), freq = FALSE, xlab = "", ylab = "", main = NULL, cex.axis = 1.5, pch = 20, cex = 0.5, xlim = xvals)
mtext(bquote(paste(v[x], " (880 hPa) ")), side = 3, line = 2, adj = 0.5, cex = 2)
mtext(bquote(paste(v[x], " (880 hPa) ")) , side = 2, line = 3, adj = 0.5, cex = 2)

yvals <- dnorm(seq(min(xvals), max(xvals), length.out = 100), mu1_ms_M1[1], sqrt(wind_var1_M1[1,1]))
lines(seq(min(xvals), max(xvals), length.out = 100), yvals, col = "#d9534f", lwd = 2)

yvals <- dnorm(seq(min(xvals), max(xvals), length.out = 100), mu_ms_M2[1], sqrt(wind_var_M2[1,1]))
lines(seq(min(xvals), max(xvals), length.out = 100), yvals, col = "#93c44b", lwd = 2)

yvals <- dnorm(seq(min(xvals), max(xvals), length.out = 100), mu1_ms_M3[1], sqrt(wind_var_M3[1,1]))
lines(seq(min(xvals), max(xvals), length.out = 100), yvals, col = "#3fb1e2", lwd = 2)

yvals <- dnorm(seq(min(xvals), max(xvals), length.out = 100), mu1_ms_M4[1], sqrt(wind_var_M4[1,1]))
lines(seq(min(xvals), max(xvals), length.out = 100), yvals, col = "#916fcd", lwd = 2)

screen(8)
par(pty = 's')
par(mai=c(0.4, 0.4, 0.4, 0.8))

hist(unique(c(WIND_ARRAY[, 2, ])), freq = FALSE, xlab = "", ylab = "", main = NULL, cex.axis = 1.5, pch = 20, cex = 0.5, xlim = xvals)
mtext(bquote(paste(v[y], " (880 hPa) ")) , side = 2, line = 3, adj = 0.5, cex = 2)

yvals <- dnorm(seq(min(xvals), max(xvals), length.out = 100), mu1_ms_M1[2], sqrt(wind_var1_M1[2, 2]))
lines(seq(min(xvals), max(xvals), length.out = 100), yvals, col = "#d9534f", lwd = 2)

yvals <- dnorm(seq(min(xvals), max(xvals), length.out = 100), mu_ms_M2[2], sqrt(wind_var_M2[2, 2]))
lines(seq(min(xvals), max(xvals), length.out = 100), yvals, col = "#93c44b", lwd = 2)

yvals <- dnorm(seq(min(xvals), max(xvals), length.out = 100), mu1_ms_M3[2], sqrt(wind_var_M3[2, 2]))
lines(seq(min(xvals), max(xvals), length.out = 100), yvals, col = "#3fb1e2", lwd = 2)

yvals <- dnorm(seq(min(xvals), max(xvals), length.out = 100), mu1_ms_M4[2], sqrt(wind_var_M4[2, 2]))
lines(seq(min(xvals), max(xvals), length.out = 100), yvals, col = "#916fcd", lwd = 2)

screen(13)
par(pty = 's')
par(mai=c(0.4, 0.4, 0.4, 0.8))

hist(unique(c(WIND_ARRAY[, 3, ])), freq = FALSE, xlab = "", ylab = "", main = NULL, cex.axis = 1.5, pch = 20, cex = 0.5, xlim = xvals)
mtext(bquote(paste(v[x], " (925 hPa) ")) , side = 2, line = 3, adj = 0.5, cex = 2)

yvals <- dnorm(seq(min(xvals), max(xvals), length.out = 100), mu2_ms_M1[1], sqrt(wind_var2_M1[1, 1]))
lines(seq(min(xvals), max(xvals), length.out = 100), yvals, col = "#d9534f", lwd = 2)

yvals <- dnorm(seq(min(xvals), max(xvals), length.out = 100), mu_ms_M2[1], sqrt(wind_var_M2[1, 1]))
lines(seq(min(xvals), max(xvals), length.out = 100), yvals, col = "#93c44b", lwd = 2)

yvals <- dnorm(seq(min(xvals), max(xvals), length.out = 100), mu2_ms_M3[1], sqrt(wind_var_M3[3, 3]))
lines(seq(min(xvals), max(xvals), length.out = 100), yvals, col = "#3fb1e2", lwd = 2)

yvals <- dnorm(seq(min(xvals), max(xvals), length.out = 100), mu2_ms_M4[1], sqrt(wind_var_M4[3, 3]))
lines(seq(min(xvals), max(xvals), length.out = 100), yvals, col = "#916fcd", lwd = 2)

screen(18)
par(pty = 's')
par(mai=c(0.4, 0.4, 0.4, 0.4))

hist(unique(c(WIND_ARRAY[, 4, ])), freq = FALSE, xlab = "", ylab = "", main = NULL, cex.axis = 1.5, pch = 20, cex = 0.5, xlim = xvals)
mtext(bquote(paste(v[y], " (925 hPa) ")) , side = 2, line = 3, adj = 0.5, cex = 2)

yvals <- dnorm(seq(min(xvals), max(xvals), length.out = 100), mu2_ms_M1[2], sqrt(wind_var2_M1[2, 2]))
lines(seq(min(xvals), max(xvals), length.out = 100), yvals, col = "#d9534f", lwd = 2)

yvals <- dnorm(seq(min(xvals), max(xvals), length.out = 100), mu_ms_M2[2], sqrt(wind_var_M2[2, 2]))
lines(seq(min(xvals), max(xvals), length.out = 100), yvals, col = "#93c44b", lwd = 2)

yvals <- dnorm(seq(min(xvals), max(xvals), length.out = 100), mu2_ms_M3[2], sqrt(wind_var_M3[4, 4]))
lines(seq(min(xvals), max(xvals), length.out = 100), yvals, col = "#3fb1e2", lwd = 2)

yvals <- dnorm(seq(min(xvals), max(xvals), length.out = 100), mu2_ms_M4[2], sqrt(wind_var_M4[4, 4]))
lines(seq(min(xvals), max(xvals), length.out = 100), yvals, col = "#916fcd", lwd = 2)


screen(4)
par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))

plot(unique(round(cbind(c(WIND_ARRAY[, 2, ]), c(WIND_ARRAY[, 1, ])), 0)), xlab = "", ylab = "", cex.axis = 1.5, pch = 20, cex = 0.5, col = "#808080", ylim = c(-max(WIND_ARRAY), max(WIND_ARRAY)), xlim = c(-max(WIND_ARRAY), max(WIND_ARRAY)))
abline(h = 0, v = 0, col = 1, lwd = 2, lty = 2)
mtext(bquote(paste(v[y], " (880 hPa) ")), side = 3, line = 0, adj = 0.5, cex = 2)


ellipse(mu = c(mu1_ms_M1[2], mu1_ms_M1[1]), sigma = wind_var1_M1[c(2,1), c(2,1)], alpha = .05, npoints = 250, col = "#d9534f", lwd = 2) 
points(matrix(c(mu1_ms_M1[2], mu1_ms_M1[1]), ncol = 2), pch = 8, lwd = 2, cex = 2, col = "#d9534f")

ellipse(mu = c(mu_ms_M2[2], mu_ms_M2[1]), sigma = wind_var_M2[c(2,1), c(2,1)], alpha = .05, npoints = 250, col = "#93c44b", lwd = 2) 
points(matrix(c(mu_ms_M2[2], mu_ms_M2[1]), ncol = 2), pch = 8, lwd = 2, cex = 2, col = "#93c44b")

ellipse(mu = c(mu1_ms_M3[2], mu1_ms_M3[1]), sigma = wind_var_M3[c(2, 1), c(2, 1)], alpha = .05, npoints = 250, col = "#3fb1e2", lwd = 2) 
points(matrix(c(mu1_ms_M3[2], mu1_ms_M3[1]), ncol = 2), pch = 8, lwd = 2, cex = 2, col = "#3fb1e2")

ellipse(mu = c(mu1_ms_M4[2], mu1_ms_M4[1]), sigma = wind_var_M4[c(2, 1), c(2, 1)], alpha = .05, npoints = 250, col = "#916fcd", lwd = 2) 
points(matrix(c(mu1_ms_M4[2], mu1_ms_M4[1]), ncol = 2), pch = 8, lwd = 2, cex = 2, col = "#916fcd")

screen(14)
par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))

plot(unique(round(cbind(c(WIND_ARRAY[, 4, ]), c(WIND_ARRAY[, 3, ])), 0)), xlab = "", ylab = "", cex.axis = 1.5, pch = 20, cex = 0.5, col = "#808080", ylim = c(-max(WIND_ARRAY), max(WIND_ARRAY)), xlim = c(-max(WIND_ARRAY), max(WIND_ARRAY)))
abline(h = 0, v = 0, col = 1, lwd = 2, lty = 2)

ellipse(mu = c(mu2_ms_M1[2], mu2_ms_M1[1]), sigma = wind_var2_M1[c(2,1), c(2,1)], alpha = .05, npoints = 250, col = "#d9534f", lwd = 2) 
points(matrix(c(mu2_ms_M1[2], mu2_ms_M1[1]), ncol = 2), pch = 8, lwd = 2, cex = 2, col = "#d9534f")

ellipse(mu = c(mu_ms_M2[2], mu_ms_M2[1]), sigma = wind_var_M2[c(2,1), c(2,1)], alpha = .05, npoints = 250, col = "#93c44b", lwd = 2) 
points(matrix(c(mu_ms_M2[2], mu_ms_M2[1]), ncol = 2), pch = 8, lwd = 1, cex = 2, col = "#93c44b")

ellipse(mu = c(mu2_ms_M3[2], mu2_ms_M3[1]), sigma = wind_var_M3[c(4, 3), c(4, 3)], alpha = .05, npoints = 250, col = "#3fb1e2", lwd = 2) 
points(matrix(c(mu2_ms_M3[2], mu2_ms_M3[1]), ncol = 2), pch = 8, lwd = 2, cex = 2, col = "#3fb1e2")

ellipse(mu = c(mu2_ms_M4[2], mu2_ms_M4[1]), sigma = wind_var_M4[c(4, 3), c(4, 3)], alpha = .05, npoints = 250, col = "#916fcd", lwd = 2) 
points(matrix(c(mu2_ms_M4[2], mu2_ms_M4[1]), ncol = 2), pch = 8, lwd = 2, cex = 2, col = "#916fcd")

screen(5)
par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))

plot(unique(round(cbind(c(WIND_ARRAY[, 3, ]), c(WIND_ARRAY[, 1, ])), 0)), xlab = "", ylab = "", cex.axis = 1.5, pch = 20, cex = 0.5, col = "#808080", ylim = c(-max(WIND_ARRAY), max(WIND_ARRAY)), xlim = c(-max(WIND_ARRAY), max(WIND_ARRAY)), yaxt = 'n', xaxt = 'n')
abline(h = 0, v = 0, col = 1, lwd = 2, lty = 2)
mtext(bquote(paste(v[x], " (925 hPa) ")), side = 3, line = 0, adj = 0.5, cex = 2)

ellipse(mu = c(mu2_ms_M3[1], mu1_ms_M3[1]), sigma = wind_var_M3[c(3,1), c(3,1)], alpha = .05, npoints = 250, col = "#3fb1e2", lwd = 2) 
points(matrix(c(mu2_ms_M3[1], mu1_ms_M3[1]), ncol = 2), pch = 8, lwd = 2, cex = 2, col = "#3fb1e2")
  
ellipse(mu = c(mu2_ms_M4[1], mu1_ms_M4[1]), sigma = wind_var_M4[c(3,1), c(3,1)], alpha = .05, npoints = 250, col = "#916fcd", lwd = 2) 
points(matrix(c(mu2_ms_M4[1], mu1_ms_M4[1]), ncol = 2), pch = 8, lwd = 2, cex = 2, col = "#916fcd")

screen(10)
par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))

plot(unique(round(cbind(c(WIND_ARRAY[, 4, ]), c(WIND_ARRAY[, 2, ])), 0)), xlab = "", ylab = "", cex.axis = 1.5, pch = 20, cex = 0.5, col = "#808080", ylim = c(-max(WIND_ARRAY), max(WIND_ARRAY)), xlim = c(-max(WIND_ARRAY), max(WIND_ARRAY)), yaxt = 'n', xaxt = 'n')
abline(h = 0, v = 0, col = 1, lwd = 2, lty = 2)

ellipse(mu = c(mu2_ms_M3[2], mu1_ms_M3[2]), sigma = wind_var_M3[c(4, 2), c(4, 2)], alpha = .05, npoints = 250, col = "#3fb1e2", lwd = 2) 
points(matrix(c(mu2_ms_M3[2], mu1_ms_M3[2]), ncol = 2), pch = 8, lwd = 2, cex = 2, col = "#3fb1e2")

ellipse(mu = c(mu2_ms_M4[2], mu1_ms_M4[2]), sigma = wind_var_M4[c(4, 2), c(4, 2)], alpha = .05, npoints = 250, col = "#916fcd", lwd = 2) 
points(matrix(c(mu2_ms_M4[2], mu1_ms_M4[2]), ncol = 2), pch = 8, lwd = 2, cex = 2, col = "#916fcd")

screen(6)
par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))

plot(unique(round(cbind(c(WIND_ARRAY[, 4, ]), c(WIND_ARRAY[, 1, ])), 0)), xlab = "", ylab = "", cex.axis = 1.5, pch = 20, cex = 0.5, col = "#808080", ylim = c(-max(WIND_ARRAY), max(WIND_ARRAY)), xlim = c(-max(WIND_ARRAY), max(WIND_ARRAY)), yaxt = 'n', xaxt = 'n')
abline(h = 0, v = 0, col = 1, lwd = 2, lty = 2)
mtext(bquote(paste(v[y], " (925 hPa) ")), side = 3, line = 0, adj = 0.5, cex = 2)

ellipse(mu = c(mu2_ms_M3[2], mu1_ms_M3[1]), sigma = wind_var_M3[c(4,1), c(4,1)], alpha = .05, npoints = 250, col = "#3fb1e2", lwd = 2) 
points(matrix(c(mu2_ms_M3[2], mu1_ms_M3[1]), ncol = 2), pch = 8, lwd = 2, cex = 2, col = "#3fb1e2")

ellipse(mu = c(mu2_ms_M4[2], mu1_ms_M4[1]), sigma = wind_var_M4[c(4, 1), c(4, 1)], alpha = .05, npoints = 250, col = "#916fcd", lwd = 2) 
points(matrix(c(mu2_ms_M4[2], mu1_ms_M4[1]), ncol = 2), pch = 8, lwd = 2, cex = 2, col = "#916fcd")

screen(9)
par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))

plot(unique(round(cbind(c(WIND_ARRAY[, 3, ]), c(WIND_ARRAY[, 2, ])), 0)), xlab = "", ylab = "", cex.axis = 1.5, pch = 20, cex = 0.5, col = "#808080", ylim = c(-max(WIND_ARRAY), max(WIND_ARRAY)), xlim = c(-max(WIND_ARRAY), max(WIND_ARRAY)))
abline(h = 0, v = 0, col = 1, lwd = 2, lty = 2)

ellipse(mu = c(mu2_ms_M3[1], mu1_ms_M3[2]), sigma = wind_var_M3[c(3,2), c(3,2)], alpha = .05, npoints = 250, col = "#3fb1e2", lwd = 2) 
points(matrix(c(mu2_ms_M3[1], mu1_ms_M3[2]), ncol = 2), pch = 8, lwd = 2, cex = 2, col = "#3fb1e2")

ellipse(mu = c(mu2_ms_M4[1], mu1_ms_M4[2]), sigma = wind_var_M4[c(3,2), c(3,2)], alpha = .05, npoints = 250, col = "#916fcd", lwd = 2) 
points(matrix(c(mu2_ms_M4[1], mu1_ms_M4[2]), ncol = 2), pch = 8, lwd = 2, cex = 2, col = "#916fcd")

screen(15)

plot(0,0, 'n', ylab = '', xlab = '', cex.axis = 2, yaxt = 'n', xaxt = 'n', bty="n")
legend(-1.2, 1.3, legend = c("M1", "M2", "M3", "M4"), lty = rep(1, 4), lwd = 4, col = c("#d9534f", "#93c44b", "#3fb1e2", "#916fcd"), cex = 2, bty="n")

close.screen( all=TRUE)
dev.off()

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

ellipse(mu = c(mu1_ms_M1[2], mu1_ms_M1[1]), sigma = wind_var1_M1[c(2,1), c(2,1)], alpha = .05, npoints = 250, col = "#d9534f", lwd = 2) 
points(matrix(c(mu1_ms_M1[2], mu1_ms_M1[1]), ncol = 2), pch = 8, lwd = 2, cex = 2, col = "#d9534f")

ellipse(mu = c(mu_ms_M2[2], mu_ms_M2[1]), sigma = wind_var_M2[c(2,1), c(2,1)], alpha = .05, npoints = 250, col = "#93c44b", lwd = 2) 
points(matrix(c(mu_ms_M2[2], mu_ms_M2[1]), ncol = 2), pch = 8, lwd = 2, cex = 2, col = "#93c44b")

ellipse(mu = c(mu1_ms_M3[2], mu1_ms_M3[1]), sigma = wind_var_M3[c(2, 1), c(2, 1)], alpha = .05, npoints = 250, col = "#3fb1e2", lwd = 2) 
points(matrix(c(mu1_ms_M3[2], mu1_ms_M3[1]), ncol = 2), pch = 8, lwd = 2, cex = 2, col = "#3fb1e2")

ellipse(mu = c(mu1_ms_M4[2], mu1_ms_M4[1]), sigma = wind_var_M4[c(2, 1), c(2, 1)], alpha = .05, npoints = 250, col = "#916fcd", lwd = 2) 
points(matrix(c(mu1_ms_M4[2], mu1_ms_M4[1]), ncol = 2), pch = 8, lwd = 2, cex = 2, col = "#916fcd")

screen(11)
par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))

plot(unique(round(cbind(c(WIND_ARRAY[, 4, ]), c(WIND_ARRAY[, 3, ])), 0)), xlab = "", ylab = "", cex.axis = 1.5, pch = 20, cex = 0.5, col = "#808080", ylim = c(-max(WIND_ARRAY), max(WIND_ARRAY)), xlim = c(-max(WIND_ARRAY), max(WIND_ARRAY)))
abline(h = 0, v = 0, col = 1, lwd = 2, lty = 2)
mtext(bquote(paste(v[x], " (925 hPa) ")) , side = 2, line = 3, adj = 0.5, cex = 2)

ellipse(mu = c(mu2_ms_M1[2], mu2_ms_M1[1]), sigma = wind_var2_M1[c(2,1), c(2,1)], alpha = .05, npoints = 250, col = "#d9534f", lwd = 2) 
points(matrix(c(mu2_ms_M1[2], mu2_ms_M1[1]), ncol = 2), pch = 8, lwd = 2, cex = 2, col = "#d9534f")

ellipse(mu = c(mu_ms_M2[2], mu_ms_M2[1]), sigma = wind_var_M2[c(2,1), c(2,1)], alpha = .05, npoints = 250, col = "#93c44b", lwd = 2) 
points(matrix(c(mu_ms_M2[2], mu_ms_M2[1]), ncol = 2), pch = 8, lwd = 1, cex = 2, col = "#93c44b")

ellipse(mu = c(mu2_ms_M3[2], mu2_ms_M3[1]), sigma = wind_var_M3[c(4, 3), c(4, 3)], alpha = .05, npoints = 250, col = "#3fb1e2", lwd = 2) 
points(matrix(c(mu2_ms_M3[2], mu2_ms_M3[1]), ncol = 2), pch = 8, lwd = 2, cex = 2, col = "#3fb1e2")

ellipse(mu = c(mu2_ms_M4[2], mu2_ms_M4[1]), sigma = wind_var_M4[c(4, 3), c(4, 3)], alpha = .05, npoints = 250, col = "#916fcd", lwd = 2) 
points(matrix(c(mu2_ms_M4[2], mu2_ms_M4[1]), ncol = 2), pch = 8, lwd = 2, cex = 2, col = "#916fcd")

screen(4)
par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))

plot(unique(round(cbind(c(WIND_ARRAY[, 3, ]), c(WIND_ARRAY[, 1, ])), 0)), xlab = "", ylab = "", cex.axis = 1.5, pch = 20, cex = 0.5, col = "#808080", ylim = c(-max(WIND_ARRAY), max(WIND_ARRAY)), xlim = c(-max(WIND_ARRAY), max(WIND_ARRAY)), yaxt = 'n', xaxt = 'n')
abline(h = 0, v = 0, col = 1, lwd = 2, lty = 2)
mtext(bquote(paste(v[x], " (925 hPa) ")), side = 3, line = 0, adj = 0.5, cex = 2)

ellipse(mu = c(mu2_ms_M3[1], mu1_ms_M3[1]), sigma = wind_var_M3[c(3,1), c(3,1)], alpha = .05, npoints = 250, col = "#3fb1e2", lwd = 2) 
points(matrix(c(mu2_ms_M3[1], mu1_ms_M3[1]), ncol = 2), pch = 8, lwd = 2, cex = 2, col = "#3fb1e2")
  
ellipse(mu = c(mu2_ms_M4[1], mu1_ms_M4[1]), sigma = wind_var_M4[c(3,1), c(3,1)], alpha = .05, npoints = 250, col = "#916fcd", lwd = 2) 
points(matrix(c(mu2_ms_M4[1], mu1_ms_M4[1]), ncol = 2), pch = 8, lwd = 2, cex = 2, col = "#916fcd")

screen(8)
par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))

plot(unique(round(cbind(c(WIND_ARRAY[, 4, ]), c(WIND_ARRAY[, 2, ])), 0)), xlab = "", ylab = "", cex.axis = 1.5, pch = 20, cex = 0.5, col = "#808080", ylim = c(-max(WIND_ARRAY), max(WIND_ARRAY)), xlim = c(-max(WIND_ARRAY), max(WIND_ARRAY)), yaxt = 'n', xaxt = 'n')
abline(h = 0, v = 0, col = 1, lwd = 2, lty = 2)

ellipse(mu = c(mu2_ms_M3[2], mu1_ms_M3[2]), sigma = wind_var_M3[c(4, 2), c(4, 2)], alpha = .05, npoints = 250, col = "#3fb1e2", lwd = 2) 
points(matrix(c(mu2_ms_M3[2], mu1_ms_M3[2]), ncol = 2), pch = 8, lwd = 2, cex = 2, col = "#3fb1e2")

ellipse(mu = c(mu2_ms_M4[2], mu1_ms_M4[2]), sigma = wind_var_M4[c(4, 2), c(4, 2)], alpha = .05, npoints = 250, col = "#916fcd", lwd = 2) 
points(matrix(c(mu2_ms_M4[2], mu1_ms_M4[2]), ncol = 2), pch = 8, lwd = 2, cex = 2, col = "#916fcd")

screen(5)
par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))

plot(unique(round(cbind(c(WIND_ARRAY[, 4, ]), c(WIND_ARRAY[, 1, ])), 0)), xlab = "", ylab = "", cex.axis = 1.5, pch = 20, cex = 0.5, col = "#808080", ylim = c(-max(WIND_ARRAY), max(WIND_ARRAY)), xlim = c(-max(WIND_ARRAY), max(WIND_ARRAY)), yaxt = 'n', xaxt = 'n')
abline(h = 0, v = 0, col = 1, lwd = 2, lty = 2)
mtext(bquote(paste(v[y], " (925 hPa) ")), side = 3, line = 0, adj = 0.5, cex = 2)

ellipse(mu = c(mu2_ms_M3[2], mu1_ms_M3[1]), sigma = wind_var_M3[c(4,1), c(4,1)], alpha = .05, npoints = 250, col = "#3fb1e2", lwd = 2) 
points(matrix(c(mu2_ms_M3[2], mu1_ms_M3[1]), ncol = 2), pch = 8, lwd = 2, cex = 2, col = "#3fb1e2")

ellipse(mu = c(mu2_ms_M4[2], mu1_ms_M4[1]), sigma = wind_var_M4[c(4, 1), c(4, 1)], alpha = .05, npoints = 250, col = "#916fcd", lwd = 2) 
points(matrix(c(mu2_ms_M4[2], mu1_ms_M4[1]), ncol = 2), pch = 8, lwd = 2, cex = 2, col = "#916fcd")

screen(7)
par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))

plot(unique(round(cbind(c(WIND_ARRAY[, 3, ]), c(WIND_ARRAY[, 2, ])), 0)), xlab = "", ylab = "", cex.axis = 1.5, pch = 20, cex = 0.5, col = "#808080", ylim = c(-max(WIND_ARRAY), max(WIND_ARRAY)), xlim = c(-max(WIND_ARRAY), max(WIND_ARRAY)))
abline(h = 0, v = 0, col = 1, lwd = 2, lty = 2)
mtext(bquote(paste(v[y], " (880 hPa) ")) , side = 2, line = 3, adj = 0.5, cex = 2)

ellipse(mu = c(mu2_ms_M3[1], mu1_ms_M3[2]), sigma = wind_var_M3[c(3,2), c(3,2)], alpha = .05, npoints = 250, col = "#3fb1e2", lwd = 2) 
points(matrix(c(mu2_ms_M3[1], mu1_ms_M3[2]), ncol = 2), pch = 8, lwd = 2, cex = 2, col = "#3fb1e2")

ellipse(mu = c(mu2_ms_M4[1], mu1_ms_M4[2]), sigma = wind_var_M4[c(3,2), c(3,2)], alpha = .05, npoints = 250, col = "#916fcd", lwd = 2) 
points(matrix(c(mu2_ms_M4[1], mu1_ms_M4[2]), ncol = 2), pch = 8, lwd = 2, cex = 2, col = "#916fcd")

screen(9)

plot(0,0, 'n', ylab = '', xlab = '', cex.axis = 2, yaxt = 'n', xaxt = 'n', bty="n")
legend(-0.7, 0.4, legend = c("M1", "M2", "M3", "M4"), lty = rep(1, 4), lwd = 3, col = c("#d9534f", "#93c44b", "#3fb1e2", "#916fcd"), cex = 2, bty="n")

close.screen( all=TRUE)
dev.off()

#########################################################     EXAMPLE OF KM/3HR TO M/S CONVERSION    #########################################################

mu1 <- c(-0.073211462, -0.005373963)
mu2 <- c(0.02514516, 0.01563370)

mu1_ms <- (mu1 * 100 / 3) * (0.2778)
mu2_ms <- (mu2 * 100 / 3) * (0.2778) 

wind_var_ms2 <- (wind_var * 10000 / 3^2) * 0.07716049


#########################################################     READING REAL DATA ESTIMATES    #########################################################

p <- read.table(paste(root, 'Results/multivariate_stationary_real_M1_A_FULL', sep = ''), header = FALSE, sep = " ") %>% as.matrix()

mu1 <- p[9:10]
mu2 <- p[11:12]

mu1_ms_M1 <- (mu1 * 100 / 3) * (0.2778)
mu2_ms_M1 <- (mu2 * 100 / 3) * (0.2778) 

wind_var_chol <- matrix(c(p[13], p[14], 0, p[15]), ncol = 2, byrow = T)
wind_var1_M1 <- t(wind_var_chol) %*% wind_var_chol

wind_var_chol <- matrix(c(p[16], p[17], 0, p[18]), ncol = 2, byrow = T)
wind_var2_M1 <- t(wind_var_chol) %*% wind_var_chol

p <- read.table(paste(root, 'Results/multivariate_stationary_real_M2_A_FULL', sep = ''), header = FALSE, sep = " ") %>% as.matrix()

mu <- p[10:11]

mu_ms_M2 <- (mu * 100 / 3) * (0.2778)

wind_var_chol <- matrix(c(p[12], p[13], 0, p[14]), ncol = 2, byrow = T)
wind_var_M2 <- t(wind_var_chol) %*% wind_var_chol

p <- read.table(paste(root, 'Results/multivariate_stationary_real_M3_A_FULL', sep = ''), header = FALSE, sep = " ") %>% as.matrix()

mu1 <- p[10:11]
mu2 <- p[12:13]
wind_var_chol <- matrix(c(p[14], p[15], p[16], p[17], 0, p[18], p[19], p[20], 0, 0, p[21], p[22], 0, 0, 0, p[23]), ncol = 4, byrow = T)
wind_var_M3 <- t(wind_var_chol) %*% wind_var_chol

mu1_ms_M3 <- (mu1 * 100 / 3) * (0.2778)
mu2_ms_M3 <- (mu2 * 100 / 3) * (0.2778) 

mu1_ms_M4 <- mu1_ms_M3 - 0.05
mu2_ms_M4 <- mu2_ms_M3 + 0.05
wind_var_M4 <- wind_var_M3 + 0.5







EMP <- unique(round(cbind(c(WIND_ARRAY[, 1, ]), c(WIND_ARRAY[, 2, ]), c(WIND_ARRAY[, 3, ]), c(WIND_ARRAY[, 4, ])), 0))

EMP_mean <- colMeans(EMP)

EMP_mean * 3 / ( 100 * 0.2778 )

test <- cov(unique(round(cbind(c(WIND_ARRAY[, 3, ]), c(WIND_ARRAY[, 1, ])), 0)))

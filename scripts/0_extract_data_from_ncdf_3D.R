
directory <- '/home/salvanmo/Desktop/'

root <- paste(directory, 'Spatio-Temporal-Cross-Covariance-Functions-under-the-Lagrangian-Framework/', sep = '')

source(file = paste(root, "Functions/load_packages.R",sep=''))

dname1 <- "BCPHOBIC"

saudi<- map("world", "Saudi", fill = TRUE)
IDs <- sapply(strsplit(saudi$names, ":"), function(x) x[1])
saudi <- map2SpatialPolygons(saudi, IDs=IDs, proj4string=CRS("+proj=longlat +datum=WGS84"))

for(yr in 1980:2019){

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

	for(mnth in 1:1){    #7,12
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
		ncname <- paste("/home/salvanmo/Downloads/MERRA/3d/MERRA2_", merra_ind, ".inst3_3d_aer_Nv.", yr, mo, "01.SUB.nc", sep='')

		ncin <- nc_open(ncname)

		DAT_array <- ncvar_get(ncin,dname1)
		u_array <- DAT_array[,, 63, ]
		v_array <- DAT_array[,, 72, ]

		# get longitude and latitude
		lon <- ncvar_get(ncin,"lon")
		lat <- ncvar_get(ncin,"lat")
		lon.lat <- expand.grid(lon, lat)

		nc_close(ncin)

		U <- u_array
		V <- v_array

		test1 <- data.frame(rep(lon.lat[,1], 8), rep(lon.lat[,2], 8), c(U), c(V))

		for(day in 2:mnth_end){
			cat('READING NETCDF DATA ===> year: ', yr, 'month: ', mnth, 'day: ', day, '\n')
			if(day > 9){
				ncname <- paste("/home/salvanmo/Downloads/MERRA/3d/MERRA2_", merra_ind, ".inst3_3d_aer_Nv.", yr, mo, day,".SUB.nc", sep='')
			}else{
				ncname <- paste("/home/salvanmo/Downloads/MERRA/3d/MERRA2_", merra_ind, ".inst3_3d_aer_Nv.", yr, mo, "0",day,".SUB.nc", sep='')
			}
			ncin <- nc_open(ncname)

			DAT_array <- ncvar_get(ncin,dname1)
			u_array <- DAT_array[,, 63, ]
			v_array <- DAT_array[,, 72, ]

			nc_close(ncin)

			U <- u_array
			V <- v_array

			test1 <- rbind(test1, data.frame(rep(lon.lat[,1], 8), rep(lon.lat[,2], 8),  c(U),  c(V)))
		}
		colnames(test1) <- c('lon', 'lat', 'Y1', 'Y2')
		spdf <- SpatialPointsDataFrame(coords = test1[, c("lon", "lat")], data = test1, proj4string = CRS("+proj=longlat +datum=WGS84"))
		saudi_data_orig <- data.frame(spdf[!is.na(over(spdf, as(saudi, "SpatialPolygons"))), ])

		data_temp <- matrix(saudi_data_orig[, 3], ncol = 550, byrow = T)
		#data_temp <- matrix(log(saudi_data_orig[, 3]), ncol = 550, byrow = T)
		data_temp[is.infinite(data_temp)] <- min(data_temp[!is.infinite(data_temp)])
		data_array1 <- rbind(data_array1, data_temp)

		data_temp <- matrix(saudi_data_orig[, 4], ncol = 550, byrow = T)
		#data_temp <- matrix(log(saudi_data_orig[, 4]), ncol = 550, byrow = T)
		data_temp[is.infinite(data_temp)] <- min(data_temp[!is.infinite(data_temp)])
		data_array2 <- rbind(data_array2, data_temp)
	}

	data_array3 <- saudi_data_orig[1:550, 1:2]
		
	write.table(data_array1, file = paste(root, "Data/ncdf/layer1_", yr, sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
	write.table(data_array2, file = paste(root, "Data/ncdf/layer2_", yr, sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)
	write.table(data_array3, file = paste(root, "Data/ncdf/LOCS-3D-dataset", sep = ''), sep = " ", row.names = FALSE, col.names = FALSE)

}

saudi<- map("world", c("Jordan", "Iraq", "Syria", "Lebanon", "Israel", "Kenya", "Eritrea", "Ethiopia", "South Sudan", "Sudan", "Egypt", "UAE", "Saudi", "Oman", "Yemen", "Somalia", "Djibouti", "Pakistan", "India", "Kuwait"), fill = TRUE)
IDs <- sapply(strsplit(saudi$names, ":"), function(x) x[1])
saudi <- map2SpatialPolygons(saudi, IDs=IDs, proj4string=CRS("+proj=longlat +datum=WGS84"))

TESTT <- test1

test1 <- TESTT[2 * 6862 + 1:6862, ]
colnames(test1) <- c('lon', 'lat', 'Z1', 'Z2')
spdf <- SpatialPointsDataFrame(coords = test1[, c("lon", "lat")], data = test1,
                               proj4string = CRS("+proj=longlat +datum=WGS84"))

saudi_data_orig <- data.frame(spdf[is.na(over(spdf, as(saudi, "SpatialPolygons"))), ])

IND <- which(saudi_data_orig[, 1] >= 43)
saudi_data_orig <- saudi_data_orig[IND,]
IND <- which(saudi_data_orig[, 2] <= 20)
saudi_data_orig <- saudi_data_orig[IND,]
IND <- which(saudi_data_orig[, 1] < 57 & saudi_data_orig[, 2] > 20)
saudi_data_orig <- saudi_data_orig[-IND,]
IND <- which(saudi_data_orig[, 1] > 70)
saudi_data_orig <- saudi_data_orig[-IND,]

saudi_data_orig[, 3] <- log(saudi_data_orig[, 3])
saudi_data_orig[, 4] <- log(saudi_data_orig[, 4])

ZZ <- log(VAL[1:550] / 39)
ZZ2 <- log(VAL[550 * 1 + 1:550] / 39)

        jpeg(file = paste(root, '/Figures/application_arabian_sea.jpg', sep = ''), width = 1600, height = 800)

split.screen( rbind(c(0.05,0.93,0.1,0.95), c(0.93,0.99,0.1,0.95)))
split.screen( figs = c( 1, 2 ), screen = 1 )

screen(3)
par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))
quilt.plot(saudi_data_orig[1:550, 1], saudi_data_orig[1:550, 2], ZZ, nx = 30, ny = 30, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 2)
map("worldHires", xlim = c(26.719, 85.078), ylim = c(5.625, 42.188), lwd = 0.75, add = T)
mtext('Longitude', side = 1, line = 4, adj = 0.5,  cex = 2.5, font = 2)
mtext('Latitude', side = 2, line = 4, adj = 0.5, cex = 2.5, font = 2)
mtext('U Component', side = 3, line = 1, adj = 0.5, cex = 3, font = 2)

screen(4)
par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))
quilt.plot(saudi_data_orig[1:550, 1], saudi_data_orig[1:550, 2], ZZ2, nx = 30, ny = 30, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 2, yaxt = 'n')
map("worldHires", xlim = c(26.719, 85.078), ylim = c(5.625, 42.188), lwd = 0.75, add = T)
mtext('Longitude', side = 1, line = 4, adj = 0.5, cex = 2.5, font = 2)
mtext('V Component', side = 3, line = 1, adj = 0.5, cex = 3, font = 2)

close.screen( all=TRUE)
dev.off()

jpeg(file = paste(root, '/Figures/application_arabian_sea.jpg', sep = ''), width = 1600, height = 800)

#plot(fit1$residuals)
#plot(data_array1[, 2])
plot(ZMAT[, 1])
dev.off()

demean <- demean2 <- matrix(, ncol = 550, nrow = 248)

for(ll in 1:550){
	fit1 <- arima(log(data_array1[, ll] + 1e-1), c(1, 0, 0))
	demean[, ll] <- fit1$residuals
	fit1 <- arima(log(data_array1[, ll] + 1e-1), c(1, 0, 0))
	demean2[, ll] <- fit1$residuals
}

ZZ <- demean[1, ]
ZZ2 <- demean[10, ]

ZZ <- data_array1[100, ]
ZZ2 <- data_array1[101, ]

ZZ <- ZZ2 <- log(data_array1[1, ]) - mean(log(data_array1[1, ])) - pred_tps

jpeg(file = paste(root, '/Figures/application_arabian_sea.jpg', sep = ''), width = 1600, height = 800)

split.screen( rbind(c(0.05,0.93,0.1,0.95), c(0.93,0.99,0.1,0.95)))
split.screen( figs = c( 1, 2 ), screen = 1 )

screen(3)
par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))
quilt.plot(saudi_data_orig[1:550, 1], saudi_data_orig[1:550, 2], ZZ, nx = 30, ny = 30, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 2)
map("worldHires", xlim = c(26.719, 85.078), ylim = c(5.625, 42.188), lwd = 0.75, add = T)
mtext('Longitude', side = 1, line = 4, adj = 0.5,  cex = 2.5, font = 2)
mtext('Latitude', side = 2, line = 4, adj = 0.5, cex = 2.5, font = 2)
mtext('U Component', side = 3, line = 1, adj = 0.5, cex = 3, font = 2)

screen(4)
par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))
quilt.plot(saudi_data_orig[1:550, 1], saudi_data_orig[1:550, 2], ZZ2, nx = 30, ny = 30, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 2, yaxt = 'n')
map("worldHires", xlim = c(26.719, 85.078), ylim = c(5.625, 42.188), lwd = 0.75, add = T)
mtext('Longitude', side = 1, line = 4, adj = 0.5, cex = 2.5, font = 2)
mtext('V Component', side = 3, line = 1, adj = 0.5, cex = 3, font = 2)

close.screen( all=TRUE)
dev.off()

X <- data.matrix(cbind(rep(1, 550), saudi_data_orig[550 * 11 + 1:550, 1:2], log(saudi_data_orig[550 * 11 + 1:550, 3])))

X <- data.matrix(cbind(rep(1, 550 * 10), saudi_data_orig[550 * 10 + 1:(550 * 10), 1:2], log(saudi_data_orig[550 * 9 + 1:(550 * 10), 3]), log(saudi_data_orig[550 * 10 + 1:(550 * 10), 3])))

start <- 200
TT <- 24

X <- data.matrix(cbind(rep(1, 550 * TT), saudi_data_orig[550 * start + 1:(550 * TT), 1:2], rep(1:TT, each = 550), log(saudi_data_orig[550 * start + 1:(550 * TT), 3])))

X1 <- rep(1, 550 * TT)
X2 <- (saudi_data_orig[550 * start + 1:(550 * TT), 1] - mean(saudi_data_orig[550 * start + 1:(550 * TT), 1])) / sd(saudi_data_orig[550 * start + 1:(550 * TT), 1])
X3 <- (saudi_data_orig[550 * start + 1:(550 * TT), 2] - mean(saudi_data_orig[550 * start + 1:(550 * TT), 2])) / sd(saudi_data_orig[550 * start + 1:(550 * TT), 2])
X4 <- (TEMP_vec[550 * start + 1:(550 * TT)] - mean(TEMP_vec[550 * start + 1:(550 * TT)])) / sd(TEMP_vec[550 * start + 1:(550 * TT)])
X5 <- (RH_vec[550 * start + 1:(550 * TT)] - mean(RH_vec[550 * start + 1:(550 * TT)])) / sd(RH_vec[550 * start + 1:(550 * TT)])
X6 <- (PS_vec[550 * start + 1:(550 * TT)] - mean(PS_vec[550 * start + 1:(550 * TT)])) / sd(PS_vec[550 * start + 1:(550 * TT)])
#X4 <- (TEMP[start + 1, ] - mean(TEMP[start + 1, ])) / sd(TEMP[start + 1, ])
#X5 <- (RH[start + 1, ] - mean(RH[start + 1, ])) / sd(RH[start + 1, ])
#X6 <- (PS[start + 1, ] - mean(PS[start + 1, ])) / sd(PS[start + 1, ])
X7 <- log(saudi_data_orig[550 * (start - 1) + 1:(550 * TT), 3])
X8 <- log(saudi_data_orig[550 * start + 1:(550 * TT), 3])

X <- cbind(X1, X2, X3, X4, X5, X6, X7, X8)

beta_hat <- solve(t(X[, 1:7]) %*% X[, 1:7]) %*% t(X[, 1:7]) %*% X[, 8]
mean_Z <- X[, 1:7] %*% beta_hat
err <- mean_Z - X[, 8]

ZZ <- ZZ2 <- err


ZMAT <- matrix(err, ncol = 550, byrow = T)
conso <- 2
ZZ <- colMeans(ZMAT[1:conso, ])
ZZ2 <- colMeans(ZMAT[conso + 1:conso, ])


ZZ <- ZZ2 <- colMeans(ZMAT)

##################################################################################


yr <- 2018

dname1 <- "PS"

#RH <- PS <- V <- U <- TEMP <- NULL

PS <- NULL

if(yr < 1992){
	merra_ind <- 100
}else if(yr >= 1992 & yr < 2001){
	merra_ind <- 200
}else if (yr >= 2001 & yr < 2011){
	merra_ind <- 300
}else{
	merra_ind <- 400
}

for(mnth in 1:1){    #7,12
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
	ncname <- paste("/home/salvanmo/Downloads/MERRA/3d/MERRA2_", merra_ind, ".inst3_3d_asm_Np.", yr, mo, "01.SUB.nc", sep='')

	ncin <- nc_open(ncname)

	DAT_array <- ncvar_get(ncin,dname1)
	u_array <- DAT_array#[,, 42, ]
	v_array <- DAT_array#[,, 42, ]

	# get longitude and latitude
	lon <- ncvar_get(ncin,"lon")
	lat <- ncvar_get(ncin,"lat")
	lon.lat <- expand.grid(lon, lat)

	nc_close(ncin)

	U <- u_array
	V <- v_array

	test1 <- data.frame(rep(lon.lat[,1], 8), rep(lon.lat[,2], 8), c(U), c(V))

	for(day in 2:mnth_end){
		cat('READING NETCDF DATA ===> year: ', yr, 'month: ', mnth, 'day: ', day, '\n')
		if(day > 9){
			ncname <- paste("/home/salvanmo/Downloads/MERRA/3d/MERRA2_", merra_ind, ".inst3_3d_asm_Np.", yr, mo, day,".SUB.nc", sep='')
		}else{
			ncname <- paste("/home/salvanmo/Downloads/MERRA/3d/MERRA2_", merra_ind, ".inst3_3d_asm_Np.", yr, mo, "0",day,".SUB.nc", sep='')
		}
		ncin <- nc_open(ncname)

		DAT_array <- ncvar_get(ncin,dname1)
		u_array <- DAT_array#[,, 42, ]
		v_array <- DAT_array#[,, 42, ]

		nc_close(ncin)

		U <- u_array
		V <- v_array

		test1 <- rbind(test1, data.frame(rep(lon.lat[,1], 8), rep(lon.lat[,2], 8),  c(U),  c(V)))
	}

	colnames(test1) <- c('lon', 'lat', 'Y1', 'Y2')
	spdf <- SpatialPointsDataFrame(coords = test1[, c("lon", "lat")], data = test1, proj4string = CRS("+proj=longlat +datum=WGS84"))
	saudi_data_orig <- data.frame(spdf[!is.na(over(spdf, as(saudi, "SpatialPolygons"))), ])

	data_temp <- matrix(saudi_data_orig[, 3], ncol = 550, byrow = T)
	PS <- rbind(PS, data_temp)
}

RH_vec <- PS_vec <- TEMP_vec <- NULL

for(tt in 1:248){
	RH_vec <- c(RH_vec, RH[tt, ])
	TEMP_vec <- c(TEMP_vec, TEMP[tt, ])
	PS_vec <- c(PS_vec, PS[tt, ])
}

#####################################################################


	ncname <- paste("/home/salvanmo/Downloads/MERRA2_400.inst3_3d_aer_Nv.20180101.SUB.nc", sep='')

	ncin <- nc_open(ncname)

	DAT_array <- ncvar_get(ncin,dname1)
	u_array <- DAT_array[,, 72, ]
	v_array <- DAT_array[,, 72, ]

	# get longitude and latitude
	lon <- ncvar_get(ncin,"lon")
	lat <- ncvar_get(ncin,"lat")
	lon.lat <- expand.grid(lon, lat)

	nc_close(ncin)

	U <- u_array
	V <- v_array

	test1 <- data.frame(rep(lon.lat[,1], 8), rep(lon.lat[,2], 8), c(U), c(V))
	for(day in 2:mnth_end){
		cat('READING NETCDF DATA ===> year: ', yr, 'month: ', mnth, 'day: ', day, '\n')
		if(day > 9){
			ncname <- paste("/home/salvanmo/Downloads/MERRA2_", merra_ind, ".inst3_3d_aer_Nv.", yr, mo, day,".SUB.nc", sep='')
		}else{
			ncname <- paste("/home/salvanmo/Downloads/MERRA2_", merra_ind, ".inst3_3d_aer_Nv.", yr, mo, "0",day,".SUB.nc", sep='')
		}
		ncin <- nc_open(ncname)

		DAT_array <- ncvar_get(ncin,dname1)
		u_array <- DAT_array[,, 72, ]
		v_array <- DAT_array[,, 72, ]

		nc_close(ncin)

		U <- u_array
		V <- v_array

		test1 <- rbind(test1, data.frame(rep(lon.lat[,1], 8), rep(lon.lat[,2], 8),  c(U),  c(V)))
	}

	colnames(test1) <- c('lon', 'lat', 'Y1', 'Y2')
	spdf <- SpatialPointsDataFrame(coords = test1[, c("lon", "lat")], data = test1, proj4string = CRS("+proj=longlat +datum=WGS84"))

	#saudi<- map("state", region = c("Iowa", "Wisconsin"), fill = TRUE)
	#saudi<- map("state", region = c("Missouri", "Iowa", "Illinois"), fill = TRUE)
	#saudi<- map("state", region = c("Missouri", "Iowa", "Illinois", "Indiana", "Wisconsin", "Minnesota"), fill = TRUE)
	#saudi<- map("state", region = c("Ohio", "West Virginia", "Kentucky", "Indiana"), fill = TRUE)
	#saudi<- map("state", region = c("Alabama", "Georgia", "South Carolina"), fill = TRUE)
	#saudi<- map("state", region = c("Colorado", "New Mexico", "Arizona", "Utah"), fill = TRUE)
	#saudi<- map("state", region = c("North Carolina", "Virginia", "West Virginia", "Kentucky", "Ohio"), fill = TRUE)
	saudi<- map("world", "France", fill = TRUE)
	IDs <- sapply(strsplit(saudi$names, ":"), function(x) x[1])
	saudi <- map2SpatialPolygons(saudi, IDs=IDs, proj4string=CRS("+proj=longlat +datum=WGS84"))

	saudi_data_orig <- data.frame(spdf[!is.na(over(spdf, as(saudi, "SpatialPolygons"))), ])

jpeg(file = paste(root, '/Figures/application_arabian_sea.jpg', sep = ''), width = 1600, height = 800)

split.screen( rbind(c(0.05,0.93,0.1,0.95), c(0.93,0.99,0.1,0.95)))
split.screen( figs = c( 1, 2 ), screen = 1 )

screen(3)
par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))
quilt.plot(saudi_data_orig[1:n, 1], saudi_data_orig[1:n, 2],  saudi_data_orig[1:n, 3], nx = 10, ny = 10, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 2)
map("worldHires", xlim = c(26.719, 85.078), ylim = c(5.625, 42.188), lwd = 0.75, add = T)
mtext('Longitude', side = 1, line = 4, adj = 0.5,  cex = 2.5, font = 2)
mtext('Latitude', side = 2, line = 4, adj = 0.5, cex = 2.5, font = 2)
mtext('U Component', side = 3, line = 1, adj = 0.5, cex = 3, font = 2)

screen(4)
par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))
quilt.plot(saudi_data_orig[1:n, 1], saudi_data_orig[1:n, 2],  saudi_data_orig[1:n, 3], nx = 30, ny = 30, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 2, yaxt = 'n')
map("worldHires", xlim = c(26.719, 85.078), ylim = c(5.625, 42.188), lwd = 0.75, add = T)
mtext('Longitude', side = 1, line = 4, adj = 0.5, cex = 2.5, font = 2)
mtext('V Component', side = 3, line = 1, adj = 0.5, cex = 3, font = 2)

close.screen( all=TRUE)
dev.off()

TT <- 8 * 31
n <- nrow(saudi_data_orig) / TT
X1 <- rep(1, n * (TT - 1))
X2 <- (saudi_data_orig[n + 1:(n *  (TT - 1)), 1] - mean(saudi_data_orig[n + 1:(n *  (TT - 1)), 1])) / sd(saudi_data_orig[n + 1:(n *  (TT - 1)), 1])
X3 <- (saudi_data_orig[n + 1:(n *  (TT - 1)), 2] - mean(saudi_data_orig[n + 1:(n *  (TT - 1)), 2])) / sd(saudi_data_orig[n + 1:(n *  (TT - 1)), 2])
X4 <- log(saudi_data_orig[1:(n *  (TT - 1)), 3])
inf_ind <- which(is.infinite(X4))
not_inf_ind <- which(!is.infinite(X4))

alternate_ind <- NULL
for(aa in 1:length(inf_ind)){
	alternate_ind <- c(alternate_ind, not_inf_ind[which.min((inf_ind[aa] - not_inf_ind)^2)])
}

X4[inf_ind] <- X4[alternate_ind]

X5 <- log(saudi_data_orig[n + 1:(n *  (TT - 1)), 3])
inf_ind <- which(is.infinite(X5))
not_inf_ind <- which(!is.infinite(X5))

alternate_ind <- NULL
for(aa in 1:length(inf_ind)){
	alternate_ind <- c(alternate_ind, not_inf_ind[which.min((inf_ind[aa] - not_inf_ind)^2)])
}

X5[inf_ind] <- X5[alternate_ind]


X <- cbind(X1, X2, X3, X4, X5)

beta_hat <- solve(t(X[, 1:4]) %*% X[, 1:4]) %*% t(X[, 1:4]) %*% X[, 5]
mean_Z <- X[, 1:4] %*% beta_hat
err <- mean_Z - X[, 5]

#err <- X[, 5]
ZZ <- err[1 * n + 1:n]
ZZ2 <- err[3 * n + 1:n]

ZZ <- res_mat[1, ]
ZZ2 <- res_mat2[1, ]

jpeg(file = paste(root, '/Figures/application_arabian_sea.jpg', sep = ''), width = 1600, height = 800)

split.screen( rbind(c(0.05,0.93,0.1,0.95), c(0.93,0.99,0.1,0.95)))
split.screen( figs = c( 1, 2 ), screen = 1 )

screen(3)
par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))
quilt.plot(saudi_data_orig[1:n, 1], saudi_data_orig[1:n, 2],  ZZ, nx = 20, ny = 20, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 2)
map("worldHires", xlim = c(26.719, 85.078), ylim = c(5.625, 42.188), lwd = 0.75, add = T)
mtext('Longitude', side = 1, line = 4, adj = 0.5,  cex = 2.5, font = 2)
mtext('Latitude', side = 2, line = 4, adj = 0.5, cex = 2.5, font = 2)
mtext('U Component', side = 3, line = 1, adj = 0.5, cex = 3, font = 2)

screen(4)
par(pty = 's')
par(mai=c(0.2,0.2,0.2,0.2))
quilt.plot(saudi_data_orig[1:n, 1], saudi_data_orig[1:n, 2],  ZZ2, nx = 20, ny = 20, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 2, yaxt = 'n')
map("worldHires", xlim = c(26.719, 85.078), ylim = c(5.625, 42.188), lwd = 0.75, add = T)
mtext('Longitude', side = 1, line = 4, adj = 0.5, cex = 2.5, font = 2)
mtext('V Component', side = 3, line = 1, adj = 0.5, cex = 3, font = 2)

close.screen( all=TRUE)
dev.off()

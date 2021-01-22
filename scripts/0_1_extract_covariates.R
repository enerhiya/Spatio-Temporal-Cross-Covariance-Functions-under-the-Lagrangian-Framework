
directory <- '/home/salvanmo/Desktop/'

root <- paste(directory, 'Spatio-Temporal-Cross-Covariance-Functions-under-the-Lagrangian-Framework/', sep = '')

source(file = paste(root, "Functions/load_packages.R",sep=''))
source(file = paste(root, "Functions/auxiliary_functions.R",sep=''))

saudi<- map("world", "Saudi", fill = TRUE)
IDs <- sapply(strsplit(saudi$names, ":"), function(x) x[1])
saudi <- map2SpatialPolygons(saudi, IDs=IDs, proj4string=CRS("+proj=longlat +datum=WGS84"))

for (yr in 1980:2019){

	data_array <- array(, dim = c(248, 550, 4, 2))

	for(VAR in 1:2){
		

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

		ncin <- nc_open(ncname)
		
		dname1 <- "RH"
		dname2 <- "T"
		dname3 <- "U"
		dname4 <- "V"

		u_array <- ncvar_get(ncin,dname1)
		v_array <- ncvar_get(ncin,dname2)
		a_array <- ncvar_get(ncin,dname3)
		b_array <- ncvar_get(ncin,dname4)
		
		# get longitude and latitude
		lon <- ncvar_get(ncin,"lon")
		lat <- ncvar_get(ncin,"lat")

		nc_close(ncin)

		if(VAR == 1)	lev <- 65	else	lev <- 68

		U <- u_array[,, lev, ]
		V <- v_array[,, lev, ]
		A <- a_array[,, lev, ]
		B <- b_array[,, lev, ]

		lon.lat <- expand.grid(lon,lat)
		lon_new <- matrix(lon.lat[, 1], ncol = length(lat))
		lat_new <- matrix(lon.lat[, 2], ncol = length(lat))

		test1 <- data.frame(rep(lon.lat[,1], 8), rep(lon.lat[,2], 8),  c(U),  c(V), c(A), c(B))

		for(day in 2:31){

			cat('READING NETCDF DATA ===> year: ', yr, 'month: ', mnth, 'day: ', day, '\n')
			if(day > 9){
				ncname <- paste("/home/salvanmo/Downloads/MERRA2_", merra_ind, ".inst3_3d_asm_Nv.", yr, mo, day,".SUB.nc", sep='')
			}else{
				ncname <- paste("/home/salvanmo/Downloads/MERRA2_", merra_ind, ".inst3_3d_asm_Nv.", yr, mo, "0",day,".SUB.nc", sep='')
			}
			ncin <- nc_open(ncname)

			u_array <- ncvar_get(ncin,dname1)
			v_array <- ncvar_get(ncin,dname2)

			nc_close(ncin)

			U <- u_array[,, lev, ]
			V <- v_array[,, lev, ]
			A <- a_array[,, lev, ]
			B <- b_array[,, lev, ]

			test1 <- rbind(test1, data.frame(rep(lon.lat[,1], 8), rep(lon.lat[,2], 8),  c(U),  c(V), c(A), c(B)))
		}

		colnames(test1) <- c('lon', 'lat', 'Y1', 'Y2')
		spdf <- SpatialPointsDataFrame(coords = test1[, c("lon", "lat")], data = test1, proj4string = CRS("+proj=longlat +datum=WGS84"))
		saudi_data_orig <- data.frame(spdf[!is.na(over(spdf, as(saudi, "SpatialPolygons"))), ])

		N <- nrow(saudi_data_orig)/(8 * mnth_end)

		data_array_temp <- array(, dim = c(248, 550, 4))
		for(tt in 1:4){
			data_temp <- matrix(saudi_data_orig[, tt + 2], ncol = N, byrow = T)
			data_array_temp[, , tt] <- data_temp
		}
		
		data_array[, , , VAR] <- data_array_temp

	}
	save(data_array, file = paste(root, "Data/ncdf/covariates_", yr, '.RData', sep = ''))

}



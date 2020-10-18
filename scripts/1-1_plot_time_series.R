#plot timeseries

#args <- commandArgs(TRUE)

#yr <- args[1]

ref_point <- 100

directory <- '/home/salvanmo/Desktop/'

root <- paste(directory, 'Spatio-Temporal-Cross-Covariance-Functions-under-the-Lagrangian-Framework/', sep = '')

source(file = paste(root, "Functions/load_packages.R",sep=''))

pdf(file = paste(root, 'Figures/timeseries-with-trend.pdf', sep = ''), width = 25, height = 20)

for(yr in 1980:1980){

	cat(yr, '\n')

	dat <- read.table(paste(root, 'Data/ncdf/DUSMASS25_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
	dat2 <- read.table(paste(root, 'Data/ncdf/BCSMASS_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()

	DAT <- dat[1:744, ]
	DAT2 <- dat2[1:744, ]

	#aggregate <- 6
	#TT <- floor(nrow(dat) / aggregate)
	#TT <- floor(744 / aggregate)

	#DAT <- DAT2 <- matrix(, nrow = TT, ncol = 550)

	#for(aa in 1:TT){
	#	DAT[aa, ] <- apply(dat[(aa - 1) * aggregate + 1:aggregate, ], 2, mean)
	#	DAT2[aa, ] <- apply(dat2[(aa - 1) * aggregate + 1:aggregate, ], 2, mean)
	#}

	dat_trend <- read.table(paste(root, 'Results/estimated_mean/DUSMASS25_trend_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
	dat2_trend <- read.table(paste(root, 'Results/estimated_mean/BCSMASS_trend_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()

	par(mfrow = c(2,1))

	par(mar=c(5, 5, 5, 5))

	plot(DAT[, ref_point], type = 'l', xaxt = 'n', xlab = '', ylab = 'log Dust Mass Concentration', cex.lab = 2)
	lines(dat_trend[, ref_point] + colMeans(DAT)[ref_point], lwd = 2, col = 'red')
	mtext(paste(yr, ' / Site: ', ref_point, sep = ''), side = 3, cex = 3, col = 'blue', font = 2, line = 1.5)	

	plot(DAT2[, ref_point], type = 'l', xlab = 'Hour', ylab = 'log Black Carbon Concentration', cex.lab = 2)
	lines(dat2_trend[, ref_point] + colMeans(DAT2)[ref_point], lwd = 2, col = 'red')

}

dev.off()


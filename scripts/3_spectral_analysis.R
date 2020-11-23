
directory <- '/home/salvanmo/Desktop/'

root <- paste(directory, 'Spatio-Temporal-Cross-Covariance-Functions-under-the-Lagrangian-Framework/', sep = '')

source(file = paste(root, "Functions/load_packages.R",sep=''))
source(file = paste(root, "Functions/auxiliary_functions.R",sep=''))

yr <- 1987

dat <- read.table(paste(root, 'Data/ncdf/layer1_residuals_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
dat2 <- read.table(paste(root, 'Data/ncdf/layer2_residuals_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
dat3 <- read.table(paste(root, 'Data/ncdf/LOCS-3D-dataset', sep = ''), header = FALSE, sep = " ") %>% as.matrix()

station_index <- c(37, 160, 413, 490)

pdf(file = paste(root, 'Figures/spacetime-spectra.pdf', sep = ''), width = 25, height = 10)
split.screen( rbind(c(0.05,0.96,0.1,0.95), c(0.98,0.99,0.1,0.95)))
split.screen( figs = c( 4, 2 ), screen = 1 )

count_screen <- 3
for(nn in 1:length(station_index)){
	
	screen(count_screen)
	count_screen <- count_screen + 1
	par(mai=c(.3,.3,.3,.3))

	raw.spec <- spec.pgram(ts.union(as.ts(dat[, station_index[nn]]), as.ts(dat2[, station_index[nn]])), spans = c(100, 100), plot = F)
	plot(raw.spec, plot.type = "coherency", main = '', ylab = '', xlab = '')
	mtext(paste('Location ', station_index[nn], ': (', round(dat3[station_index[nn], 1], 2), ', ', round(dat3[station_index[nn], 2], 2), ')', sep = ''), side = 2, line = 4, adj = 0.5,  cex = 1, font = 2, col = 'blue')
	mtext('Spectrum', side = 2, line = 2, adj = 0.5,  cex = 1, font = 2)
	if(nn == 1){
		mtext('Coherence', side = 3, line = 1.5, adj = 0.5,  cex = 1.5, font = 2)
	}
	if(nn == length(station_index)){
		mtext('Frequency', side = 1, line = 2, adj = 0.5,  cex = 1, font = 2)
	}				
	
	screen(count_screen)
	count_screen <- count_screen + 1
	
	par(mai=c(.3,.3,.3,.3))
	
	plot(raw.spec, plot.type = "phase", main = '', ylab = '', xlab = '')
	if(nn == 1){
		mtext('Phase', side = 3, line = 1.5, adj = 0.5,  cex = 1.5, font = 2)
	}
	if(nn == length(station_index)){
		mtext('Frequency', side = 1, line = 2, adj = 0.5,  cex = 1, font = 2)
	}				
}
close.screen( all=TRUE)
dev.off()


pdf(file = paste(root, 'Figures/spacetime-spectra.pdf', sep = ''), width = 10, height = 10)
spectrum(data.frame(dat2[, 1], dat2[, 2]))$coh
dev.off()


pdf(file = paste(root, 'Figures/spacetime-spectra.pdf', sep = ''), width = 10, height = 10)
coherence <- coh(dat2[, 1], dat2[, 10], f=1000, plot=TRUE)
dev.off()



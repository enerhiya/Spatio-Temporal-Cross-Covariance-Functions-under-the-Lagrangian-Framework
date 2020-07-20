
plot_func <- function(data, file_name, variable_name){
	zlim_range <- range(data[, 3] - mean(data[, 3]))
	zlim_range <- max(abs(ceiling(zlim_range)), abs(floor(zlim_range)))

	jpeg(file = paste(file_name, sep = ''), width = 800, height = 800)

	split.screen( rbind(c(0.05,0.93,0.1,0.95), c(0.92,0.99,0.1,0.95)))

	screen(1)
	quilt.plot(data[, 2], data[, 1], data[,3] - mean(data[,3]), zlim = c(-zlim_range, zlim_range), nx = 200, ny = 200, ylab = '', xlab = '', cex.lab = 4, add.legend = F, cex.axis = 2, ylim = c(25, 50), xlim = c(-90, -60))
	#US( add=TRUE, lwd=2)
	#map("worldHires", "usa", xlim=c(-140,-110), ylim=c(48,64))
	map("state", xlim = range(data[, 2]), ylim = range(data[, 1]), lwd = 0.75, add = T)
	mtext('Longitude', side = 1, line = 4, adj = 0.5,  cex = 2.5, font = 2)
	mtext('Latitude', side = 2, line = 4, adj = 0.5, cex = 2.5, font = 2)
	mtext(variable_name, side = 3, line = 1, adj = 0.5, cex = 3, font = 2)

	screen(2)
	x1 <- c(0.025,0.12,0.12,0.025)
	y1 <- c(0.25,0.25,0.7,0.7)
	legend.gradient2(cbind(x1,y1), title = "", limits = round(seq(-zlim_range, zlim_range,length.out = 5), 1), cex = 2)

	close.screen( all=TRUE)
	dev.off()
}

write_to_txt <- function(data, file_name){

	#comma-separated if you save the bivariate realizations data.
	#space-separated if you save the cross-covariance matrix.	

	if(ncol(data) > 2){
		write.table(data, file = file_name, sep = " ", row.names = FALSE, col.names = FALSE)
	}else{	
		write.table(data, file = file_name, sep = ",", row.names = FALSE, col.names = FALSE)
	}
}

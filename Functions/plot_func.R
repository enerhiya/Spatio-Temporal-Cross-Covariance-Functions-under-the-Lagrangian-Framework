directory <- '/home/salvanmo/Desktop/'

root <- paste(directory, 'Spatio-Temporal-Cross-Covariance-Functions-under-the-Lagrangian-Framework/', sep = '')

source(file = paste(root, "Functions/auxiliary_functions.R", sep = ''))

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


fig1 <- function(COVARIANCES, LOCS, file_name){

	#COVARIANCES is a list of matrices with dimension nrow(LOCS) * T x 3

  	n <- nrow(LOCS)
  	N <- sqrt(n)

  	pdf(file = file_name, width = 14, height = 14)
  	
	new_grid_x <- matrix(LOCS[, 1], N, N)
  	new_grid_y <- matrix(LOCS[, 1], N, N, byrow = T)
	
	split.screen( rbind(c(0.08,0.48,0.53,0.95), c(0.53,0.93,0.53,0.95), c(0.08,0.48,0.05,0.47), c(0.53,0.93,0.05,0.47), c(0.92,0.99,0.1,0.9)))

	split.screen( figs = c( 3, 3 ), screen = 1 )
	split.screen( figs = c( 3, 3 ), screen = 2 )
	split.screen( figs = c( 3, 3 ), screen = 3 )
	split.screen( figs = c( 3, 3 ), screen = 4 )

	for(model in 1:4){
		
		vals_temp <- COVARIANCES[[model]]
		zlim_range <- c(0, 1)
		
		for(variable in 1:3){
			for(tt in 1:3){
				
				screen(5 + (model - 1) * 9 + (variable - 1) * 3 + tt)	
				par(pty="s")
				par(mai=c(0.02, 0.02, 0.02, 0.02))
	
				poly.image(new_grid_x, new_grid_y, matrix(vals_temp[(tt - 1) * n + 1:n, variable], N, N), zlim = zlim_range, xlab = " ", ylab = " ", xaxt = 'n', yaxt = 'n')

				if(tt == 1 & (model == 1 | model == 3)){
					mtext(expression(s[y]), side = 2, line = 2, cex = 1)
                                        axis(2, at = seq(-0.5, 0.5, length.out = 6))
				}else if(variable == 3 & (model == 3 | model == 4)){
					mtext(expression(s[x]), side = 1, line = 2, cex = 1)
                                        axis(1, at = seq(-0.5, 0.5, length.out = 6))
				}
				if(variable == 3 & model == 3 & tt == 1){
                                        axis(1, at = seq(-0.5, 0.5, length.out = 6))
                                        mtext(expression(s[x]), side = 1, line = 2, cex = 1)
                                }

				if(variable == 1){
					if(tt == 1){
						mtext(expression(t[1]), side = 3, line = 0.1, adj = 0.5, cex = 1, font=2)
					}else if(tt == 2){
						mtext(expression(t[2]), side = 3, line = 0.1, adj = 0.5, cex = 1, font=2)
					}else{
						mtext(expression(t[3]), side = 3, line = 0.1, adj = 0.5, cex = 1, font=2)
					}
				}

				if(variable == 1 & tt == 1 & (model == 1 | model == 3)){
                                        mtext(expression(C[11]), side = 2, line = 3, cex = 1.5, font=2, col="#000080")
				}		
					 	
				if(variable == 2 & tt == 1 & (model == 1 | model == 3)){
                                        mtext(expression(C[22]), side = 2, line = 3, cex = 1.5, font=2, col="#000080")
				}		

				if(variable == 3 & tt == 1 & (model == 1 | model == 3)){
                                        mtext(expression(C[12]), side = 2, line = 3, cex = 1.5, font=2, col="#000080")
				}		
			}
		}
	}


  	screen(5)
  	x1 <- c(0.025,0.09,0.09,0.025) + 0.3
  	y1 <- c(0.25,0.25,0.75,0.75)
	legend.gradient3(cbind(x1,y1), title = "", limits = seq(0, 1,length.out = 3))

  	screen(1)
  	mtext(expression((a)), side = 3, line = 4.5, adj = 0.5, cex = 1.5, font=2, col="#000080")

  	screen(2)
  	mtext(expression((b)), side = 3, line = 4.5, adj = 0.5, cex = 1.5, font=2, col="#000080")

  	screen(3)
  	mtext(expression((c)), side = 3, line = 4.5, adj = 0.5, cex = 1.5, font=2, col="#000080")

  	screen(4)
  	mtext(expression((d)), side = 3, line = 4.5, adj = 0.5, cex = 1.5, font=2, col="#000080")

  	close.screen( all=TRUE)
  	dev.off()

}

fig2 <- function(REALIZATIONS, LOCS, file_name){

	#REALIZATIONS is a list of matrices with dimension nrow(LOCS) * T x 2

  	n <- nrow(LOCS)
  	N <- sqrt(n)

  	pdf(file = file_name, width = 14, height = 10)
  	
	new_grid_x <- matrix(LOCS[, 1], N, N)
  	new_grid_y <- matrix(LOCS[, 1], N, N, byrow = T)
  	
	split.screen( rbind(c(0.08,0.48,0.53,0.95), c(0.53,0.93,0.53,0.95), c(0.08,0.48,0.05,0.47), c(0.53,0.93,0.05,0.47), c(0.92,0.99,0.1,0.9)))

	split.screen( figs = c( 2, 3 ), screen = 1 )
	split.screen( figs = c( 2, 3 ), screen = 2 )
	split.screen( figs = c( 2, 3 ), screen = 3 )
	split.screen( figs = c( 2, 3 ), screen = 4 )

	for(model in 1:4){
		
		vals_temp <- REALIZATIONS[[model]]
		zlim_range <- range(vals_temp)
		
		for(variable in 1:2){
			for(tt in 1:3){
				
				screen(5 + (model - 1) * 6 + (variable - 1) * 3 + tt)	
				par(pty="s")
				par(mai=c(0.02, 0.02, 0.02, 0.02))
	
				poly.image(new_grid_x, new_grid_y, matrix(vals_temp[(tt - 1) * n + 1:n, variable], N, N), zlim = zlim_range, xlab = " ", ylab = " ", xaxt = 'n', yaxt = 'n')

				if(tt == 1 & (model == 1 | model == 3)){
					mtext(expression(s[y]), side = 2, line = 2, cex = 1)
                                        axis(2, at = seq(0, 1, length.out = 6))
				}else if(variable == 2 & (model == 3 | model == 4)){
					mtext(expression(s[x]), side = 1, line = 2, cex = 1)
                                        axis(1, at = seq(0, 1, length.out = 6))
				}
				if(variable == 2 & model == 3 & tt == 1){
                                        axis(1, at = seq(0, 1, length.out = 6))
                                        mtext(expression(s[x]), side = 1, line = 2, cex = 1)
                                }

				if(variable == 1){
					if(tt == 1){
						mtext(expression(t[1]), side = 3, line = 0.1, adj = 0.5, cex = 1, font=2)
					}else if(tt == 2){
						mtext(expression(t[2]), side = 3, line = 0.1, adj = 0.5, cex = 1, font=2)
					}else{
						mtext(expression(t[3]), side = 3, line = 0.1, adj = 0.5, cex = 1, font=2)
					}
				}

				if(variable == 1 & tt == 1 & (model == 1 | model == 3)){
                                        mtext(expression(Z[1]), side = 2, line = 3, cex = 1.5, font=2, col="#000080")
				}		
					 	
				if(variable == 2 & tt == 1 & (model == 1 | model == 3)){
                                        mtext(expression(Z[2]), side = 2, line = 3, cex = 1.5, font=2, col="#000080")
				}		

			}
		}
	}


  	screen(5)
  	x1 <- c(0.025,0.09,0.09,0.025) + 0.3
  	y1 <- c(0.25,0.25,0.75,0.75)
	legend.gradient2(cbind(x1,y1), title = "", limits = seq(-3, 3,length.out = 5))

  	screen(1)
  	mtext(expression((a)), side = 3, line = 4.5, adj = 0.5, cex = 1.5, font=2, col="#000080")

  	screen(2)
  	mtext(expression((b)), side = 3, line = 4.5, adj = 0.5, cex = 1.5, font=2, col="#000080")

  	screen(3)
  	mtext(expression((c)), side = 3, line = 4.5, adj = 0.5, cex = 1.5, font=2, col="#000080")

  	screen(4)
  	mtext(expression((d)), side = 3, line = 4.5, adj = 0.5, cex = 1.5, font=2, col="#000080")

  	close.screen( all=TRUE)
  	dev.off()

}


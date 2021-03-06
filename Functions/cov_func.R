#---------STATIONARY------------#
frozen_matern_cov_parapply <- function(wind, max_time_lag = 1, q = 2, LOCS = locs_demean, theta = space_params){

	###################												###################
	################### RETURNS a q * nrow(LOCS) * (max_time_lag + 1) x q * nrow(LOCS) * (max_time_lag + 1) matrix 	###################
	###################												###################	

	w <- wind
	
	locs <- coords <- LOCS
  
  	if (max_time_lag > 0){
		for (tt in 1:max_time_lag){
			temp_locs <- cbind(coords[, 1] - tt * w[1], coords[, 2] - tt * w[2])
			locs <- rbind(locs, temp_locs)
		}
  	} 
  
	nu <- theta[4:5]
	beta <- theta[3]
	sigma2 <- theta[1:2]
  
	dist0 <- parDist(x = locs, method = "euclidean") %>% as.matrix()	

	S <- matrix(NA,  q * nrow(dist0), q * nrow(dist0))
  
	for(i in 1:q){
		for(j in 1:i){

			temp <- (i - 1) * nrow(dist0) + 1:nrow(dist0)
	      		temp1 <- (j - 1) * nrow(dist0) + 1:nrow(dist0)
	      
	      		if(i == j){
		
				temp2 <- ifelse(dist0 != 0, sigma2[i] * (dist0 / beta)^nu[i] * besselK(dist0 / beta, nu[i]) / (2^(nu[i] - 1) * gamma(nu[i])), sigma2[i])
				S[temp, temp1] <- temp2
		
	      		}
	      
		      	if(i != j){
			
				nu1 <- nu[i]
				nu2 <- nu[j]
				nu3 <- (nu1 + nu2)/2
			
				rho <- theta[6] * (gamma(nu1 + 3/2) / gamma(nu1))^(1/2) * (gamma(nu2 + 3/2) / gamma(nu2))^(1/2) * gamma(nu3) / (gamma(nu3 + 3/2))
			
				temp3 <- (dist0 / beta)^nu3 * besselK(dist0 / beta, nu3)/(2^(nu3 - 1) * gamma(nu3)) * sqrt(sigma2[i] * sigma2[j]) * rho
				temp3[is.na(temp3)] <- sqrt(sigma2[i] * sigma2[j]) * rho
				S[temp, temp1] <- temp3
				S[temp1, temp] <- t(temp3)
		      }
	    	}
	  }
  
	return(c(S))
}
frozen_matern_cov <- function(theta, wind = NULL, max_time_lag = 0, q = 2, LOCS){

	###################												###################
	################### RETURNS a q * nrow(LOCS) * (max_time_lag + 1) x q * nrow(LOCS) * (max_time_lag + 1) matrix 	###################
	###################												###################	

	w <- wind
	
	locs <- coords <- LOCS
  
  	if (max_time_lag > 0){
		for (tt in 1:max_time_lag){
			temp_locs <- cbind(coords[, 1] - tt * w[1], coords[, 2] - tt * w[2])
			locs <- rbind(locs, temp_locs)
		}
  	} 
  
	nu <- theta[4:5]
	beta <- theta[3]
	sigma2 <- theta[1:2]
  
	dist0 <- parDist(x = locs, method = "euclidean") %>% as.matrix()	

	S <- matrix(NA,  q * nrow(dist0), q * nrow(dist0))
  
	for(i in 1:q){
		for(j in 1:i){

			temp <- (i - 1) * nrow(dist0) + 1:nrow(dist0)
	      		temp1 <- (j - 1) * nrow(dist0) + 1:nrow(dist0)
	      
	      		if(i == j){
		
				temp2 <- ifelse(dist0 != 0, sigma2[i] * (dist0 / beta)^nu[i] * besselK(dist0 / beta, nu[i]) / (2^(nu[i] - 1) * gamma(nu[i])), sigma2[i])
				S[temp, temp1] <- temp2
		
	      		}
	      
		      	if(i != j){
			
				nu1 <- nu[i]
				nu2 <- nu[j]
				nu3 <- (nu1 + nu2)/2
			
				rho <- theta[6] * (gamma(nu1 + 3/2) / gamma(nu1))^(1/2) * (gamma(nu2 + 3/2) / gamma(nu2))^(1/2) * gamma(nu3) / (gamma(nu3 + 3/2))
			
				temp3 <- (dist0 / beta)^nu3 * besselK(dist0 / beta, nu3)/(2^(nu3 - 1) * gamma(nu3)) * sqrt(sigma2[i] * sigma2[j]) * rho
				temp3[is.na(temp3)] <- sqrt(sigma2[i] * sigma2[j]) * rho
				S[temp, temp1] <- temp3
				S[temp1, temp] <- t(temp3)
		      }
	    	}
	  }
  
	return(S)
}

frozen_matern_cov_rep_I <- function(theta, wind, max_time_lag = 0, q = 2, LOCS){

	###################												###################
	################### RETURNS a q * nrow(LOCS) * (max_time_lag + 1) x q * nrow(LOCS) * (max_time_lag + 1) matrix 	###################
	###################												###################	

	w <- wind
	
	locs <- coords <- LOCS
  
  	if (max_time_lag > 0){
		for (tt in 1:max_time_lag){
			temp_locs <- cbind(coords[, 1] - tt * w[1], coords[, 2] - tt * w[2])
			locs <- rbind(locs, temp_locs)
		}
  	} 
  
	nu <- theta[4:5]
	beta <- theta[3]
	sigma2 <- theta[1:2]
  
	dist0 <- parDist(x = locs, method = "euclidean") %>% as.matrix()	

	S <- matrix(NA,  q * nrow(dist0), q * nrow(dist0))
  
	for(i in 1:q){
		for(j in 1:i){

			temp <- seq(i, nrow(dist0) * q, by = q)
	      		temp1 <- seq(j, nrow(dist0) * q, by = q)
	      
	      		if(i == j){
		
				temp2 <- ifelse(dist0 != 0, sigma2[i] * (dist0 / beta)^nu[i] * besselK(dist0 / beta, nu[i]) / (2^(nu[i] - 1) * gamma(nu[i])), sigma2[i])
				S[temp, temp1] <- temp2
		
	      		}
	      
		      	if(i != j){
			
				nu1 <- nu[i]
				nu2 <- nu[j]
				nu3 <- (nu1 + nu2)/2
			
				rho <- theta[6] * (gamma(nu1 + 3/2) / gamma(nu1))^(1/2) * (gamma(nu2 + 3/2) / gamma(nu2))^(1/2) * gamma(nu3) / (gamma(nu3 + 3/2))
			
				temp3 <- (dist0 / beta)^nu3 * besselK(dist0 / beta, nu3)/(2^(nu3 - 1) * gamma(nu3)) * sqrt(sigma2[i] * sigma2[j]) * rho
				temp3[is.na(temp3)] <- sqrt(sigma2[i] * sigma2[j]) * rho
				S[temp, temp1] <- temp3
				S[temp1, temp] <- t(temp3)
		      }
	    	}
	  }
  
	return(S)
}

frozen_matern_cov_for_heatmap <- function(theta, wind, q = 2){

	###################															###################
	################### RETURNS a q * nrow(LOCS) * (max_time_lag + 1) x 3 matrix of covariance values with respect to location (0,0) 	###################
	###################															###################	

	N <- 51
        n <- N^2
        TT <- 3
        grid_x <- seq(from = -0.5, to = 0.5, length.out = N)
	sim_grid_locations <- expand.grid(grid_x, grid_x) %>% as.matrix()

	w <- wind
	
	nu <- theta[4:5]
	beta <- theta[3]
	sigma2 <- theta[1:2]

	nu1 <- nu[1]
	nu2 <- nu[2]
	nu3 <- (nu1 + nu2)/2
	rho <- theta[6] * (gamma(nu1 + 3/2) / gamma(nu1))^(1/2) * (gamma(nu2 + 3/2) / gamma(nu2))^(1/2) * gamma(nu3) / (gamma(nu3 + 3/2))

	S <- matrix(NA,  n * TT, 3)

	for(i in 1:3){
		for(tt in 0:(TT - 1)){

			for(l in 1:nrow(sim_grid_locations)){

				temp_locs <- rbind(cbind(0, 0), cbind(sim_grid_locations[l, 1] - tt * w[1], sim_grid_locations[l, 2] - tt * w[2]))
				dist0 <- dist(temp_locs) %>% as.numeric()
				
				if(i < 3) temp2 <- ifelse(dist0 != 0, sigma2[i] * (dist0 / beta)^nu[i] * besselK(dist0 / beta, nu[i]) / (2^(nu[i] - 1) * gamma(nu[i])), sigma2[i])
				else temp2 <- ifelse(dist0 != 0, (dist0 / beta)^nu3 * besselK(dist0 / beta, nu3)/(2^(nu3 - 1) * gamma(nu3)) * sqrt(sigma2[1] * sigma2[2]) * rho, sqrt(sigma2[1] * sigma2[2]) * rho)
				S[tt * n + l, i] <- temp2 
			}
		}
	}

	return(S)
}



#------------------------------- END ----------------------------#



matern_cov_soph <- function(theta, wind, max_time_lag, q, new_locations = locations, meters = T, nug_eff = F, kap = F){
  
  w <- wind
  
  loc1 <- coords1 <- new_locations
  
  if (max_time_lag == 0){
    loc1 <- loc1
  } else {
    for (tt in 1:max_time_lag){
      temploc <- matrix(, ncol=2, nrow=nrow(coords1))
      for(rr in 1:nrow(coords1)){
        temploc[rr,] <- c(coords1[rr,1] - tt*w[1], coords1[rr,2] - tt*w[2])
      }
      loc1 <- rbind(loc1, temploc)
    }
  }
  
  loc2 <- coords2 <- cbind(new_locations[,1] - kap[1], new_locations[,2] - kap[2])
  
  if (max_time_lag == 0){
    loc2 <- loc2
  } else {
    for (tt in 1:max_time_lag){
      temploc <- matrix(, ncol=2, nrow=nrow(coords2))
      for(rr in 1:nrow(coords2)){
        temploc[rr,] <- c(coords2[rr,1] - tt*w[1], coords2[rr,2] - tt*w[2])
      }
      loc2 <- rbind(loc2, temploc)
    }
  }
  loc <- rbind(loc1, loc1, loc2)
  
  if(meters == T){
    dist0 <- spDists(loc, longlat=F)/1000
  }else{
    dist0 <- spDists(loc, longlat=F)
  }
  
  nu <- theta[1:2]
  beta <- theta[3]
  var <- theta[4:5]
  rho <- theta[6]
  
  if(nug_eff == T){
    nug <- theta[7:8]
  }else{
    nug <- c(0, 0)
  }
  
  S=matrix(NA,  q*dim(dist0)[1], q*dim(dist0)[1])
  
  for(i in 1:q){
    for(j in 1:i){
      
      temp=(i-1)*dim(dist0)[1]+1:dim(dist0)[1]
      temp1=(j-1)*dim(dist0)[1]+1:dim(dist0)[1]
      
      if(i==j){
        
        temp2=ifelse(dist0!=0,var[i]*(dist0/beta)^nu[i] * besselK(dist0/beta, nu[i])/(2^(nu[i]-1)*gamma(nu[i])),var[i]+nug[i])
        S[temp,temp1]=temp2
        
      }
      
      if(i != j){
        
        nu1 <- nu[i]
        nu2 <- nu[j]
        nu3 <- (nu1 + nu2)/2
        
        #rho=Beta[i,j]*(gamma(nu1+3/2)/gamma(nu1))^(1/2) * (gamma(nu2+3/2)/gamma(nu2))^(1/2)*gamma(nu3)/(gamma(nu3+3/2))
        
        temp3 <- (dist0/beta)^nu3 * besselK(dist0/beta,nu3)/(2^(nu3-1)*gamma(nu3))*sqrt(var[i] * var[j])*rho
        temp3[is.na(temp3)] <- sqrt(var[i] * var[j])*rho
        S[temp,temp1] <- temp3
        S[temp1,temp] <- t(temp3)
      }
    }
  }
  
  S1 <- rbind(cbind(S[1:nrow(loc1), 1:nrow(loc1)], S[1:nrow(loc1), (nrow(loc1)*3 + 1):(nrow(loc1)*4)]),
              cbind(S[(nrow(loc1)*3 + 1):(nrow(loc1)*4), 1:nrow(loc1)], S[(nrow(loc1)*3 + 1):(nrow(loc1)*4), (nrow(loc1)*3 + 1):(nrow(loc1)*4)]))
  return(S1)
}

matern_cov_old <- function(theta, wind, max_time_lag, q, new_locations = locations, meters = T, nug_eff, kap = matrix(c(704400+ 100, 205700 + 100, 100, 100), ncol = 2, nrow = 2, byrow=T)){
  
  w <- wind
  
  loc1 <- coords1 <- cbind(new_locations[,1] + kap[1,1] - kap[2,1], new_locations[,2] + kap[1,2] - kap[2,2])
  
  if (max_time_lag == 0){
    loc1 <- loc1
  } else {
    for (tt in 1:max_time_lag){
      temploc <- matrix(, ncol=2, nrow=nrow(coords1))
      for(rr in 1:nrow(coords1)){
        temploc[rr,] <- c(coords1[rr,1] - tt*w[1], coords1[rr,2] - tt*w[2])
      }
      loc1 <- rbind(loc1, temploc)
    }
  }
  
  loc2 <- coords2 <- cbind(new_locations[,1] - kap[1,1] + kap[2,1], new_locations[,2] - kap[1,2] + kap[2,2])
  
  if (max_time_lag == 0){
    loc2 <- loc2
  } else {
    for (tt in 1:max_time_lag){
      temploc <- matrix(, ncol=2, nrow=nrow(coords2))
      for(rr in 1:nrow(coords2)){
        temploc[rr,] <- c(coords2[rr,1] - tt*w[1], coords2[rr,2] - tt*w[2])
      }
      loc2 <- rbind(loc2, temploc)
    }
  }
  loc <- rbind(loc1, loc2)
  
  if(meters == T){
    dist0 <- spDists(loc, longlat=F)/1000
  }else{
    dist0 <- spDists(loc, longlat=F)
  }
  
  nu <- theta[1:2]
  beta <- theta[3]
  var <- theta[4:5]
  rho <- theta[6]
  
  if(nug_eff == T){
    nug <- theta[7:8]
  }else{
    nug <- c(0, 0)
  }
  
  S=matrix(NA,  q*dim(dist0)[1], q*dim(dist0)[1])
  
  for(i in 1:q){
    for(j in 1:i){
      
      temp=(i-1)*dim(dist0)[1]+1:dim(dist0)[1]
      temp1=(j-1)*dim(dist0)[1]+1:dim(dist0)[1]
      
      if(i==j){
        
        temp2=ifelse(dist0!=0,var[i]*(dist0/beta)^nu[i] * besselK(dist0/beta, nu[i])/(2^(nu[i]-1)*gamma(nu[i])),var[i]+nug[i])
        S[temp,temp1]=temp2
        
      }
      
      if(i != j){
        
        nu1 <- nu[i]
        nu2 <- nu[j]
        nu3 <- (nu1 + nu2)/2
        
        #rho=Beta[i,j]*(gamma(nu1+3/2)/gamma(nu1))^(1/2) * (gamma(nu2+3/2)/gamma(nu2))^(1/2)*gamma(nu3)/(gamma(nu3+3/2))
        
        temp3 <- (dist0/beta)^nu3 * besselK(dist0/beta,nu3)/(2^(nu3-1)*gamma(nu3))*sqrt(var[i] * var[j])*rho
        temp3[is.na(temp3)] <- sqrt(var[i] * var[j])*rho
        S[temp,temp1] <- temp3
        S[temp1,temp] <- t(temp3)
      }
    }
  }
  
  S1 <- rbind(cbind(S[1:nrow(loc1), 1:nrow(loc1)], S[1:nrow(loc1), (nrow(loc1)*3 + 1):(nrow(loc1)*4)]),
              cbind(S[(nrow(loc1)*3 + 1):(nrow(loc1)*4), 1:nrow(loc1)], S[(nrow(loc1)*3 + 1):(nrow(loc1)*4), (nrow(loc1)*3 + 1):(nrow(loc1)*4)]))
  return(S1)
}

matern_random_cov <- function(theta, wind, wind_var, max_time_lag, q, new_locations, meters = T, nug_eff, kap){
  
  nu <- theta[1:2]
  beta <- theta[3]
  rho <- theta[4]
  var <- theta[5:6]
  if(nug_eff == T){
    nug <- theta[7:8]
  }else{
    nug <- c(0, 0)
  }
  
  if(meters == T){
    w <- wind/1000
    Sigma <- matrix(c(wind_var[1:2], wind_var[2:3]), ncol = 2)/1000
    loc <- coords <- new_locations/1000
    loc2 <- coords2 <- cbind(new_locations[,1] - kap[1], new_locations[,2] - kap[2])/1000
  }else{
    w <- wind
    Sigma <- matrix(c(wind_var[1:2], wind_var[2:3]), ncol=2)
    loc <- coords <- new_locations
    loc2 <- coords2 <- cbind(new_locations[,1] - kap[1], new_locations[,2] - kap[2])
  }
  
  SS <- list()
  
  S <- matrix(NA,  q*nrow(coords)*(max_time_lag + 1), q*nrow(coords)*(max_time_lag + 1))
  
  for(i in 1:q){
    for(j in 1:i){
      
      temp2 <- (i-1)*(nrow(coords)*(max_time_lag + 1)) + 1:(nrow(coords)*(max_time_lag + 1))
      temp1 <- (j-1)*(nrow(coords)*(max_time_lag + 1)) + 1:(nrow(coords)*(max_time_lag + 1))
      
      if(i == j){
        
        for(tt in 0:max_time_lag){
          temploc <- matrix(, ncol=nrow(coords), nrow=nrow(coords))
          for(rr in 1:nrow(coords)){
            for(ss in 1:nrow(coords)){
              cat(tt,rr,ss,'\n')
              h <- c(coords[rr,1]-coords[ss,1],coords[rr,2]-coords[ss,2])
              emp_cov1 <- c(h[1], h[2], tt)
              Int.func <- function(c, hvec){   
                y.fun  <- function(y) y^(nu[i])*exp(-y)*dmvn(X = hvec[1:2], mu = hvec[3]*w, sigma = (hvec[3]^2*Sigma + beta^2*2*y*diag(2)))
                sapply(c, y.fun)
              }
              lai <- function(xxxx) integrate(Int.func, lower = 0, upper = Inf, hvec = xxxx, abs.tol = 1e-18, rel.tol = 1e-18)$val
              temp <- lai(emp_cov1) 
              
              temploc[rr,ss] <- ifelse(tt == 0 & h[1] == 0 & h[2] == 0, var[i], var[i]*4*pi*temp*beta^2/gamma(nu[i]))
            }
          }
          SS[[tt + 1]] <- temploc
        }
        S2 <- toeplitz_mat(SS)
        S[temp2,temp1] <- S2 
        
      }
      
      if(i != j){
        
        nu1 <- nu[i]
        nu2 <- nu[j]
        nu3 <- (nu1 + nu2)/2
        
        #rho=rot*(gamma(nu1+3/2)/gamma(nu1))^(1/2) * (gamma(nu2+3/2)/gamma(nu2))^(1/2)*gamma(nu3)/(gamma(nu3+3/2))
        
        for(tt in 0:max_time_lag){
          temploc <- matrix(, ncol = nrow(coords), nrow = nrow(coords))
          for(rr in 1:nrow(coords)){
            for(ss in 1:nrow(coords)){
              cat(tt, rr, ss,'\n')
              
              h <- c(coords2[rr,1] - coords2[ss,1], coords2[rr,2] - coords2[ss,2])
              emp_cov1 <- c(h[1], h[2], tt)
              Int.func <- function(c, hvec){   
                y.fun  <- function(y) y^(nu3)*exp(-y)*dmvn(X = hvec[1:2], mu = hvec[3]*w, sigma = (hvec[3]^2*Sigma + beta^2*2*y*diag(2)))
                sapply(c, y.fun)
              }
              lai <- function(xxxx) integrate(Int.func, lower = 0, upper = Inf, hvec = xxxx, abs.tol = 1e-18, rel.tol = 1e-18)$val
              temp <- lai(emp_cov1) 
              
              temploc[rr, ss] <- ifelse(tt == 0 & h[1] == 0 & h[2] == 0, sqrt(var[i] * var[j])*rho, sqrt(var[i] * var[j])*rho*4*pi*temp*beta^2/gamma(nu3))
            }
          }
          SS[[tt + 1]] <- temploc
        }
        S2 <- toeplitz_mat(SS)
        S[temp2, temp1] <- S2
        S[temp1, temp2] <- t(S2)
      }
    }
  }
  return(S)
}

lmc_cov <- function(theta, wind, max_time_lag, q, new_locations = locations, meters = T, nug_eff){
  
  nu <- theta[1:2]
  beta <- theta[3:4]
  var <- theta[5:6]
  
  if(nug_eff == T){
    nug <- theta[7:8]
  }else{
    nug <- c(0, 0)
  }
  
  alpha <- matrix(c(theta[7], theta[8], theta[9], theta[10]), ncol=2, byrow=T)
  
  S <- list()
  for(i in 1:q){
    
    if(meters == T){
      w <- matrix(wind, ncol = 2, byrow = T)/1000
      loc <- coords <- locations/1000
    }else{
      w <- matrix(wind, ncol = 2, byrow = T)
      loc <- coords <- locations
    }
    
    if (max_time_lag == 0){
      loc <- coords
    } else {
      for (tt in 1:max_time_lag){
        temploc <- matrix(, ncol=2, nrow=nrow(coords))
        for(rr in 1:nrow(coords)){
          temploc[rr,] <- c(coords[rr,1] - tt*w[i,1], coords[rr,2] - tt*w[i,2])
        }
        loc <- rbind(loc, temploc)
      }
    }
    dist0 <- spDists(loc, longlat = F)
  
    SS <- ifelse(dist0 != 0, var[i]*(dist0/beta[i])^nu[i] * besselK(dist0/beta[i], nu[i])/(2^(nu[i] - 1)*gamma(nu[i])), var[i] + nug[i])
    
    S[[i]] <- SS
  }
  S1 <- rbind(cbind(alpha[1,1]^2*S[[1]] + alpha[1,2]^2*S[[2]], alpha[1,1]*alpha[2,1]*S[[1]] + alpha[1,2]*alpha[2,2]*S[[2]]),
              cbind(alpha[1,1]*alpha[2,1]*S[[1]] + alpha[1,2]*alpha[2,2]*S[[2]], alpha[2,1]^2*S[[1]] + alpha[2,2]^2*S[[2]]))
  
  return(S1)
}

lmc_random_cov <- function(theta, wind, wind_var, max_time_lag, q, new_locations, meters = T, nug_eff){
  
  nu <- theta[1:2]
  beta <- theta[3:4]
  var <- theta[5:6]
  
  if(nug_eff == T){
    nug <- theta[7:8]
  }else{
    nug <- c(0, 0)
  }
  
  alpha <- matrix(c(theta[7], theta[8], theta[9], theta[10]), ncol=2, byrow=T)
  
  if(meters == T){
    w <- matrix(wind, ncol = 2, byrow = T)/1000
    loc <- coords <- locations/1000
    sigma <- wind_var
    sigma[[1]] <- sigma[[1]]/1000
    sigma[[2]] <- sigma[[2]]/1000
  }else{
    w <- matrix(wind, ncol = 2, byrow = T)
    loc <- coords <- locations
    sigma <- wind_var
  }
  
  S <- list()
  temploc <- denom <- list()
  
  for(i in 1:q){
    
    temploc[[1]] <- loc <- spDists(coords, longlat = F)
    denom[[1]] <- matrix(1, ncol = ncol(temploc[[1]]), nrow = nrow(temploc[[1]]))
    
    if (max_time_lag == 0){
      loc <- spDists(coords, longlat = F)
    } else {
      for (tt in 1:max_time_lag){
        
        temploc.temp <- matrix(, ncol=nrow(coords), nrow=nrow(coords))
        for(rr in 1:nrow(coords)){
          for(ss in 1:nrow(coords)){
            temploc.temp[rr,ss] <- sqrt((coords[rr,] - coords[ss,] - tt*w[i,])%*%solve(diag(2) + sigma[[i]])%*%matrix((coords[rr,] - coords[ss,] - tt*w[i,]), ncol=1))
          }
        }
        
        temploc[[tt + 1]] <- temploc.temp
        denom[[tt + 1]] <- sqrt(det(diag(2) + tt^2*sigma[[i]]))
      }
    }
    dist0 <- toeplitz_mat(temploc)
    denom.fin <- toeplitz_mat(denom)
    
    SS <- ifelse(dist0 != 0, var[i]*(dist0/beta[i])^nu[i] * besselK(dist0/beta[i], nu[i])/(2^(nu[i]-1)*gamma(nu[i])), var[i] + nug[i])
    
    S[[i]] <- SS/denom.fin
  }
  S1 <- rbind(cbind(alpha[1,1]^2*S[[1]] + alpha[1,2]^2*S[[2]], alpha[1,1]*alpha[2,1]*S[[1]] + alpha[1,2]*alpha[2,2]*S[[2]]),
              cbind(alpha[1,1]*alpha[2,1]*S[[1]] + alpha[1,2]*alpha[2,2]*S[[2]], alpha[2,1]^2*S[[1]] + alpha[2,2]^2*S[[2]]))
  
  return(S1)
}

matern_allard <- function(theta, max_time_lag, q, new_locations = locations, meters = T, nug_eff){
  
   loc <- coords <- new_locations
   tloc <- list()
  
  if (max_time_lag == 0){
    loc <- loc
    tloc.temp <- matrix(0, ncol = nrow(coords), nrow = nrow(coords))
    tloc[[1]] <- tloc.temp
  } else {
    tloc.temp <- matrix(0, ncol = nrow(coords), nrow = nrow(coords))
    tloc[[1]] <- tloc.temp
    for (tt in 1:max_time_lag){
      loc <- rbind(loc, coords)
      tloc.temp <- matrix(tt, ncol = nrow(coords), nrow = nrow(coords))
      tloc[[tt + 1]] <- tloc.temp
    }
  }
  
  if(meters == T){
    dist1 <- spDists(loc, longlat = F)/1000
  }else{
    dist1 <- spDists(loc, longlat = F)
  }
  
  nu <- theta[1:2]
  beta <- theta[3]
  var <- theta[4:5]
  rho <- theta[6]
  
  if(nug_eff == T){
    nug <- theta[7:8]
  }else{
    nug <- c(0, 0)
    alpha <- theta[7]
    a <- 1
    b <- theta[8]
  }
  
  dist2 <- toeplitz_mat(tloc)
  dist0 <- dist1/(alpha*dist2^(2*a) + 1)^(b/2)
  
  S=matrix(NA,  q*dim(dist0)[1], q*dim(dist0)[1])
  
  for(i in 1:q){
    for(j in 1:i){
      
      temp=(i-1)*dim(dist0)[1]+1:dim(dist0)[1]
      temp1=(j-1)*dim(dist0)[1]+1:dim(dist0)[1]
      
      if(i == j){
        
        temp2=ifelse(dist0 != 0, var[i]*(dist0/beta)^nu[i] * besselK(dist0/beta,nu[i])/(2^(nu[i]-1)*gamma(nu[i]))/((alpha*(dist2)^(2*a)+1)),(var[i]+nug[i])/(alpha*(dist2)^(2*a)+1))
        #temp2=ifelse(dist0!=0,var[i]*(dist0/beta[i])^nu[i] * besselK(dist0/beta[i],nu[i])/(2^(nu[i]-1)*gamma(nu[i])),var[i]+nug[i])
        
        S[temp,temp1]=temp2
        
      }
      
      if(i != j){
        
        nu1 <- nu[i]
        nu2 <- nu[j]
        nu3 <- (nu1+nu2)/2
        
        #rho=Beta[i,j]*(gamma(nu1+3/2)/gamma(nu1))^(1/2) * (gamma(nu2+3/2)/gamma(nu2))^(1/2)*gamma(nu3)/(gamma(nu3+3/2))
        
        lai=ifelse(dist0 != 0, (dist0/beta)^nu3 * besselK(dist0/beta, nu3)/(2^(nu3 - 1)*gamma(nu3))*sqrt(var[i] * var[j])*rho/((alpha*(dist2)^(2*a) + 1)), sqrt(var[i] * var[j])*rho/((alpha*(dist2)^(2*a)+1)))
        S[temp,temp1] <- lai
        S[temp1,temp] <- t(lai)
        
      }
    }
  }
  return(S)
}

cov_lagrangian <- function(theta, wind, max_time_lag, q, new_locations = locations, meters = T, nug_eff){
  
  w <- wind
  
  loc <- coords <- new_locations
  
  if (max_time_lag == 0){
    loc <- loc
  } else {
    for (tt in 1:max_time_lag){
      temploc <- matrix(, ncol=2, nrow=nrow(coords))
      for(rr in 1:nrow(coords1)){
        temploc[rr,] <- c(coords[rr,1] - tt*w[1], coords[rr,2] - tt*w[2])
      }
      loc <- rbind(loc, temploc)
    }
  }
  
  if(meters == T){
    dist0 <- spDists(loc, longlat = F)/1000
  }else{
    dist0 <- spDists(loc, longlat = F)
  }
  
  mu <- theta[1:2]
  nu <- theta[3]
  scale <- theta[4]
  beta <- c(1,1)
  
  S <- matrix(NA,  q*dim(dist0)[1], q*dim(dist0)[1])
  
  for(i in 1:q){
    for(j in 1:i){
      
      temp <- (i - 1)*dim(dist0)[1] + 1:dim(dist0)[1]
      temp1 <- (j - 1)*dim(dist0)[1] + 1:dim(dist0)[1]
      
      if(i == j){
        
        temp2 <- pmax((1 - dist0/scale), 0)^(nu + mu[i])
        S[temp, temp1] <- temp2
        
      }
      
      if(i != j){
        
        mu1 <- mu[i]
        mu2 <- mu[j]
        mu3 <- (mu1 + mu2)/2
        
        beta3 <- (gamma(1 + mu3)/gamma(1 + nu + mu3))*sqrt((gamma(1 + nu + mu1)*gamma(1 + nu + mu2))/(gamma(1 + mu1)*gamma(1 + mu2)))
        
        lai <- beta3*pmax((1 - dist0/scale), 0)^(nu + mu3)
        S[temp, temp1] <- lai
        S[temp1, temp] <- t(lai)
        
      }
    }
  }
  return(S)
}


#---------NONSTATIONARY---------#

matern_cov_regular_grid <-function(theta,wind,time){
  
  w <- wind
  
  t=time
  q=2
  
  # create a spatial autocorrelation signature
  # coordinate list
  
  loc <- coords <- sim_grid_locations
  
  n <- nrow(sim_grid_locations)
  
  if (t==1){
    loc <- loc
  } else {
    for (tt in 1:(t-1)){
      temploc <- matrix(,ncol=2,nrow=nrow(coords))
      for(rr in 1:nrow(coords)){
        temploc[rr,] <- c(coords[rr,1]-tt*w[1],coords[rr,2]-tt*w[2])
      }
      loc <- rbind(loc, temploc)
    }
  }
  
  locations <- loc
  
  theta2 <- function (n,beta0,beta1,beta2,beta3,beta4) {
    theta3 <- beta0 + beta1*(locations[,1] - .5) + beta2*(locations[,2]-.5) + 
      beta3*(locations[,1] - .5)^2 + beta4*(locations[,2] - .5)^2
    theta3 <- matrix(theta3,nrow=nrow(locations),ncol=1)
    return(theta3)
  }
  
  #log.lam1.1<-theta2(n,-3,1,1,-6,-7)
  #log.lam1.2<-theta2(n,-5,1,1,6,-4)
  #logit.phi.1<-theta2(n,0,1,-2,0,1)
  
  #log.lam2.1<-theta2(n,-1.65,0.5,0.5,0,0)
  #log.lam2.2<-theta2(n,-2.8,-1,2,0,-7)
  #logit.phi.2<-theta2(n,-3,-1,2,0,-1)
  
  #log.lam2.1<-theta2(n,-3,-1,-1,-6,-7)
  #log.lam2.2<-theta2(n,-5,-1,-1,6,-4)
  #logit.phi.2<-theta2(n,0,-1,-2,0,1)
  
  log.lam1.1<-theta2(n,-3,1,1,-6,-7)
  log.lam1.2<-theta2(n,-5,1,1,2,-12)
  logit.phi.1<-theta2(n,0,1,-2,0,1)
  
  log.lam2.1<-theta2(n,-3,-1,-1,-6,-7)
  log.lam2.2<-theta2(n,-5,-1,-1,2,-12)
  logit.phi.2<-theta2(n,0,-1,-2,0,1)
  
  KERNEL_LIST <- list()
  
  kernel.local <- array(0, dim = c(2, 2, nrow(locations)))
  for(i in 1:nrow(locations)){
    lam1 <- exp(log.lam1.1[i,])
    lam2 <- exp(log.lam1.2[i,])
    phi <- (pi/2)*exp(logit.phi.1[i,])/(1+exp(logit.phi.1[i,]))
    Pmat <- matrix(c(cos(phi), -sin(phi), sin(phi), cos(phi)), nrow = 2, byrow = T)
    Dmat <- diag(c(lam1, lam2))
    Sigma <- Pmat %*% Dmat %*% t(Pmat)
    kernel.local[, ,i] <-  Sigma
  }
  
  KERNEL_LIST[[1]] <- kernel.local
  
  for(i in 1:nrow(locations)){
    lam1 <- exp(log.lam2.1[i,])
    lam2 <- exp(log.lam2.2[i,])
    phi <- (pi/2)*exp(logit.phi.2[i,])/(1+exp(logit.phi.2[i,]))
    Pmat <- matrix(c(cos(phi), -sin(phi), sin(phi), cos(phi)), nrow = 2, byrow = T)
    Dmat <- diag(c(lam1, lam2))
    Sigma <- Pmat %*% Dmat %*% t(Pmat)
    kernel.local[, ,i] <-  Sigma
  }
  
  KERNEL_LIST[[2]] <- kernel.local
  
  ##Calculate Matern form Nonstationary Covariance function 
  FIN_Sigma.mat <- list()
  dist0 <- list()
  
  for(KK in 1:2){
    Sigma.mat <- matrix(rep(NA, (n*t)^2), nrow = n*t)
    Q.mat <- matrix(rep(NA, (n*t)^2), nrow = n*t)
    Inv_ij <- matrix(rep(NA,4),2,2)
    
    for (i in 1:nrow(locations)) {
      #Sigma.mat[i, i] <- 1
      #Q.mat[i, i] <- 0
      Kernel_i <- KERNEL_LIST[[KK]][, , i]
      det_i <- Kernel_i[1,1] * Kernel_i[2,2] - Kernel_i[1,2] * Kernel_i[2,1]
      for (j in 1:nrow(locations)) {
        Kernel_j <- KERNEL_LIST[[KK]][, , j]
        det_j <- Kernel_j[1,1] * Kernel_j[2,2] - Kernel_j[1,2] * Kernel_j[2,1]
        Kernel_ij <- 0.5 * (Kernel_i + Kernel_j)
        Inv_ij[1,1] <- Kernel_ij[2,2] 
        Inv_ij[2,2] <- Kernel_ij[1,1] 
        Inv_ij[2,1] <- - Kernel_ij[2,1] 
        Inv_ij[1,2] <- - Kernel_ij[1,2] 
        det_ij <- Kernel_ij[1,1] * Kernel_ij[2,2] - Kernel_ij[1,2] * Kernel_ij[2,1]
        x <- c(locations[i,1] - locations[j,1], locations[i,2] - locations[j,2])
        Sigma.mat[i, j] <- sqrt(sqrt(det_i * det_j)/det_ij)
        Q.mat[i, j] <- sqrt(t(x) %*% Inv_ij %*% x/det_ij)
        #Sigma.mat[j, i] <- Sigma.mat[i, j]
        #Q.mat[j, i] <- Q.mat[i, j]
      }
    }
    FIN_Sigma.mat[[KK]] <- Sigma.mat
    dist0[[KK]] <- Q.mat
  }
  
  for (i in 1:nrow(locations)) {
    #Sigma.mat[i, i] <- 1
    #Q.mat[i, i] <- 0
    Kernel_i <- KERNEL_LIST[[1]][, , i]
    det_i <- Kernel_i[1,1] * Kernel_i[2,2] - Kernel_i[1,2] * Kernel_i[2,1]
    for (j in 1:nrow(locations)) {
      Kernel_j <- KERNEL_LIST[[2]][, , j]
      det_j <- Kernel_j[1,1] * Kernel_j[2,2] - Kernel_j[1,2] * Kernel_j[2,1]
      Kernel_ij <- 0.5 * (Kernel_i + Kernel_j)
      Inv_ij[1,1] <- Kernel_ij[2,2] 
      Inv_ij[2,2] <- Kernel_ij[1,1] 
      Inv_ij[2,1] <- - Kernel_ij[2,1] 
      Inv_ij[1,2] <- - Kernel_ij[1,2] 
      det_ij <- Kernel_ij[1,1] * Kernel_ij[2,2] - Kernel_ij[1,2] * Kernel_ij[2,1]
      x <- c(locations[i,1] - locations[j,1], locations[i,2] - locations[j,2])
      Sigma.mat[i, j] <- sqrt(sqrt(det_i * det_j)/det_ij)
      Q.mat[i, j] <- sqrt(t(x) %*% Inv_ij %*% x/det_ij)
      #Sigma.mat[j, i] <- Sigma.mat[i, j]
      #Q.mat[j, i] <- Q.mat[i, j]
    }
  }
  FIN_Sigma.mat[[3]] <- Sigma.mat
  dist0[[3]] <- Q.mat
  
  for (i in 1:nrow(locations)) {
    #Sigma.mat[i, i] <- 1
    #Q.mat[i, i] <- 0
    Kernel_i <- KERNEL_LIST[[2]][, , i]
    det_i <- Kernel_i[1,1] * Kernel_i[2,2] - Kernel_i[1,2] * Kernel_i[2,1]
    for (j in 1:nrow(locations)) {
      Kernel_j <- KERNEL_LIST[[1]][, , j]
      det_j <- Kernel_j[1,1] * Kernel_j[2,2] - Kernel_j[1,2] * Kernel_j[2,1]
      Kernel_ij <- 0.5 * (Kernel_i + Kernel_j)
      Inv_ij[1,1] <- Kernel_ij[2,2] 
      Inv_ij[2,2] <- Kernel_ij[1,1] 
      Inv_ij[2,1] <- - Kernel_ij[2,1] 
      Inv_ij[1,2] <- - Kernel_ij[1,2] 
      det_ij <- Kernel_ij[1,1] * Kernel_ij[2,2] - Kernel_ij[1,2] * Kernel_ij[2,1]
      x <- c(locations[i,1] - locations[j,1], locations[i,2] - locations[j,2])
      Sigma.mat[i, j] <- sqrt(sqrt(det_i * det_j)/det_ij)
      Q.mat[i, j] <- sqrt(t(x) %*% Inv_ij %*% x/det_ij)
      #Sigma.mat[j, i] <- Sigma.mat[i, j]
      #Q.mat[j, i] <- Q.mat[i, j]
    }
  }
  FIN_Sigma.mat[[4]] <- Sigma.mat
  dist0[[4]] <- Q.mat
  
  nu=theta[1:2]
  beta=theta[3]
  rot=theta[4]
  Beta=matrix(0,q,q)
  diag(Beta)=1
  Beta[2,1]=rot
  Beta[1,2]=rot
  
  S=matrix(NA,  q*dim(dist0[[1]])[1], q*dim(dist0[[1]])[1])
  
  for(i in 1:q){
    for(j in 1:q){
      
      temp=(i-1)*dim(dist0[[1]])[1]+1:dim(dist0[[1]])[1]
      temp1=(j-1)*dim(dist0[[1]])[1]+1:dim(dist0[[1]])[1]
      
      if(i==j){
        
        temp2=ifelse(dist0[[i]]!=0,FIN_Sigma.mat[[i]]*(dist0[[i]]/beta)^nu[i] * besselK(dist0[[i]]/beta,nu[i])/(2^(nu[i]-1)*gamma(nu[i])),FIN_Sigma.mat[[i]])
        S[temp,temp1]=temp2
        
      }
      
      if(i!=j & i<j){
        
        nu1=nu[i]
        nu2=nu[j]
        nu3=(nu[i]+nu[j])/2
        
        rho=Beta[i,j]*(gamma(nu1+3/2)/gamma(nu1))^(1/2) * (gamma(nu2+3/2)/gamma(nu2))^(1/2)*gamma(nu3)/(gamma(nu3+3/2))
        #rho <- Beta[i,j]*gamma(nu3)/sqrt(gamma(nu1)*gamma(nu2))
        
        lai=ifelse(dist0[[3]]!=0 ,(dist0[[3]]/beta)^nu3 * besselK(dist0[[3]]/beta,nu3)/(2^(nu3-1)*gamma(nu3))*sqrt(FIN_Sigma.mat[[i]] * FIN_Sigma.mat[[j]])*rho,sqrt(FIN_Sigma.mat[[i]] * FIN_Sigma.mat[[j]])*rho)
        S[temp,temp1]=lai
        #S[temp1,temp]=t(lai)
      }
      
      if(i!=j & i>j){
        
        nu1=nu[i]
        nu2=nu[j]
        nu3=(nu[i]+nu[j])/2
        
        rho=Beta[i,j]*(gamma(nu1+3/2)/gamma(nu1))^(1/2) * (gamma(nu2+3/2)/gamma(nu2))^(1/2)*gamma(nu3)/(gamma(nu3+3/2))
        #rho <- Beta[i,j]*gamma(nu3)/sqrt(gamma(nu1)*gamma(nu2))
        
        lai=ifelse(dist0[[4]]!=0 ,(dist0[[4]]/beta)^nu3 * besselK(dist0[[4]]/beta,nu3)/(2^(nu3-1)*gamma(nu3))*sqrt(FIN_Sigma.mat[[i]] * FIN_Sigma.mat[[j]])*rho,sqrt(FIN_Sigma.mat[[i]] * FIN_Sigma.mat[[j]])*rho)
        S[temp,temp1]=lai
        #S[temp1,temp]=t(lai)
      }
    }
  }
  
  return(S)
}

matern_cov_regular_grid_v4 <-function(theta,wind,time){
  
  w <- wind
  
  t=time
  q=2
  
  # create a spatial autocorrelation signature
  # coordinate list
  
  loc <- coords <- sim_grid_locations
  
  n <- nrow(sim_grid_locations)
  
  if (t==1){
    loc <- loc
  } else {
    for (tt in 1:(t-1)){
      temploc <- matrix(,ncol=2,nrow=nrow(coords))
      for(rr in 1:nrow(coords)){
        temploc[rr,] <- c(coords[rr,1]-tt*w[1],coords[rr,2]-tt*w[2])
      }
      loc <- rbind(loc, temploc)
    }
  }
  
  locations <- loc
  
  theta2 <- function (n,beta0,beta1,beta2,beta3,beta4) {
    theta3 <- beta0 + beta1*(locations[,1] - .5) + beta2*(locations[,2]-.5) + 
      beta3*(locations[,1] - .5)^2 + beta4*(locations[,2] - .5)^2
    theta3 <- matrix(theta3,nrow=nrow(locations),ncol=1)
    return(theta3)
  }
  
  #log.lam1.1<-theta2(n,-3,1,1,-6,-7)
  #log.lam1.2<-theta2(n,-5,1,1,6,-4)
  #logit.phi.1<-theta2(n,0,1,-2,0,1)
  
  #log.lam2.1<-theta2(n,-1.65,0.5,0.5,0,0)
  #log.lam2.2<-theta2(n,-2.8,-1,2,0,-7)
  #logit.phi.2<-theta2(n,-3,-1,2,0,-1)
  
  #log.lam2.1<-theta2(n,-3,-1,-1,-6,-7)
  #log.lam2.2<-theta2(n,-5,-1,-1,6,-4)
  #logit.phi.2<-theta2(n,0,-1,-2,0,1)
  
  log.lam1.1<-theta2(n,-3,1,1,-6,-7)
  log.lam1.2<-theta2(n,-5,1,1,2,-12)
  logit.phi.1<-theta2(n,0,1,-2,0,1)
  
  log.lam2.1<-theta2(n,-3,-1,-1,-6,-7)
  log.lam2.2<-theta2(n,-5,-1,-1,2,-12)
  logit.phi.2<-theta2(n,0,-1,-2,0,1)
  
  KERNEL_LIST <- list()
  
  kernel.local <- array(0, dim = c(2, 2, nrow(locations)))
  for(i in 1:nrow(locations)){
    lam1 <- exp(log.lam1.1[i,])
    lam2 <- exp(log.lam1.2[i,])
    phi <- (pi/2)*exp(logit.phi.1[i,])/(1+exp(logit.phi.1[i,]))
    Pmat <- matrix(c(cos(phi), -sin(phi), sin(phi), cos(phi)), nrow = 2, byrow = T)
    Dmat <- diag(c(lam1, lam2))
    Sigma <- Pmat %*% Dmat %*% t(Pmat)
    kernel.local[, ,i] <-  Sigma
  }
  
  KERNEL_LIST[[1]] <- kernel.local
  
  for(i in 1:nrow(locations)){
    lam1 <- exp(log.lam2.1[i,])
    lam2 <- exp(log.lam2.2[i,])
    phi <- (pi/2)*exp(logit.phi.2[i,])/(1+exp(logit.phi.2[i,]))
    Pmat <- matrix(c(cos(phi), -sin(phi), sin(phi), cos(phi)), nrow = 2, byrow = T)
    Dmat <- diag(c(lam1, lam2))
    Sigma <- Pmat %*% Dmat %*% t(Pmat)
    kernel.local[, ,i] <-  Sigma
  }
  
  KERNEL_LIST[[2]] <- kernel.local
  
  ##Calculate Matern form Nonstationary Covariance function 
  FIN_Sigma.mat <- list()
  dist0 <- list()
  
  for(KK in 1:2){
    Sigma.mat <- matrix(rep(NA, (n*t)^2), nrow = n*t)
    Q.mat <- matrix(rep(NA, (n*t)^2), nrow = n*t)
    Inv_ij <- matrix(rep(NA,4),2,2)
    
    for (i in 1:nrow(locations)) {
      #Sigma.mat[i, i] <- 1
      #Q.mat[i, i] <- 0
      Kernel_i <- KERNEL_LIST[[KK]][, , i]
      det_i <- Kernel_i[1,1] * Kernel_i[2,2] - Kernel_i[1,2] * Kernel_i[2,1]
      for (j in 1:nrow(locations)) {
        Kernel_j <- KERNEL_LIST[[KK]][, , j]
        det_j <- Kernel_j[1,1] * Kernel_j[2,2] - Kernel_j[1,2] * Kernel_j[2,1]
        Kernel_ij <- 0.5 * (Kernel_i + Kernel_j)
        Inv_ij[1,1] <- Kernel_ij[2,2] 
        Inv_ij[2,2] <- Kernel_ij[1,1] 
        Inv_ij[2,1] <- - Kernel_ij[2,1] 
        Inv_ij[1,2] <- - Kernel_ij[1,2] 
        det_ij <- Kernel_ij[1,1] * Kernel_ij[2,2] - Kernel_ij[1,2] * Kernel_ij[2,1]
        x <- c(locations[i,1] - locations[j,1], locations[i,2] - locations[j,2])
        Sigma.mat[i, j] <- sqrt(sqrt(det_i * det_j)/det_ij)
        Q.mat[i, j] <- sqrt(t(x) %*% Inv_ij %*% x/det_ij)
        #Sigma.mat[j, i] <- Sigma.mat[i, j]
        #Q.mat[j, i] <- Q.mat[i, j]
      }
    }
    FIN_Sigma.mat[[KK]] <- Sigma.mat
    dist0[[KK]] <- Q.mat
  }
  
  for (i in 1:nrow(locations)) {
    #Sigma.mat[i, i] <- 1
    #Q.mat[i, i] <- 0
    Kernel_i <- KERNEL_LIST[[1]][, , i]
    det_i <- Kernel_i[1,1] * Kernel_i[2,2] - Kernel_i[1,2] * Kernel_i[2,1]
    for (j in 1:nrow(locations)) {
      Kernel_j <- KERNEL_LIST[[2]][, , j]
      det_j <- Kernel_j[1,1] * Kernel_j[2,2] - Kernel_j[1,2] * Kernel_j[2,1]
      Kernel_ij <- 0.5 * (Kernel_i + Kernel_j)
      Inv_ij[1,1] <- Kernel_ij[2,2] 
      Inv_ij[2,2] <- Kernel_ij[1,1] 
      Inv_ij[2,1] <- - Kernel_ij[2,1] 
      Inv_ij[1,2] <- - Kernel_ij[1,2] 
      det_ij <- Kernel_ij[1,1] * Kernel_ij[2,2] - Kernel_ij[1,2] * Kernel_ij[2,1]
      x <- c(locations[i,1] - locations[j,1], locations[i,2] - locations[j,2])
      Sigma.mat[i, j] <- sqrt(sqrt(det_i * det_j)/det_ij)
      Q.mat[i, j] <- sqrt(t(x) %*% Inv_ij %*% x/det_ij)
      #Sigma.mat[j, i] <- Sigma.mat[i, j]
      #Q.mat[j, i] <- Q.mat[i, j]
    }
  }
  FIN_Sigma.mat[[3]] <- Sigma.mat
  dist0[[3]] <- Q.mat
  
  for (i in 1:nrow(locations)) {
    #Sigma.mat[i, i] <- 1
    #Q.mat[i, i] <- 0
    Kernel_i <- KERNEL_LIST[[2]][, , i]
    det_i <- Kernel_i[1,1] * Kernel_i[2,2] - Kernel_i[1,2] * Kernel_i[2,1]
    for (j in 1:nrow(locations)) {
      Kernel_j <- KERNEL_LIST[[1]][, , j]
      det_j <- Kernel_j[1,1] * Kernel_j[2,2] - Kernel_j[1,2] * Kernel_j[2,1]
      Kernel_ij <- 0.5 * (Kernel_i + Kernel_j)
      Inv_ij[1,1] <- Kernel_ij[2,2] 
      Inv_ij[2,2] <- Kernel_ij[1,1] 
      Inv_ij[2,1] <- - Kernel_ij[2,1] 
      Inv_ij[1,2] <- - Kernel_ij[1,2] 
      det_ij <- Kernel_ij[1,1] * Kernel_ij[2,2] - Kernel_ij[1,2] * Kernel_ij[2,1]
      x <- c(locations[i,1] - locations[j,1], locations[i,2] - locations[j,2])
      Sigma.mat[i, j] <- sqrt(sqrt(det_i * det_j)/det_ij)
      Q.mat[i, j] <- sqrt(t(x) %*% Inv_ij %*% x/det_ij)
      #Sigma.mat[j, i] <- Sigma.mat[i, j]
      #Q.mat[j, i] <- Q.mat[i, j]
    }
  }
  FIN_Sigma.mat[[4]] <- Sigma.mat
  dist0[[4]] <- Q.mat
  
  nu=theta[1:2]
  beta=theta[3]
  rot=theta[4]
  Beta=matrix(0,q,q)
  diag(Beta)=1
  Beta[2,1]=rot
  Beta[1,2]=rot
  
  S=matrix(NA,  q*dim(dist0[[1]])[1], q*dim(dist0[[1]])[1])
  
  for(i in 1:q){
    for(j in 1:q){
      
      temp=(i-1)*dim(dist0[[1]])[1]+1:dim(dist0[[1]])[1]
      temp1=(j-1)*dim(dist0[[1]])[1]+1:dim(dist0[[1]])[1]
      
      if(i==j){
        
        temp2=ifelse(dist0[[i]]!=0,FIN_Sigma.mat[[i]]*(dist0[[i]]/beta)^nu[i] * besselK(dist0[[i]]/beta,nu[i])/(2^(nu[i]-1)*gamma(nu[i])),FIN_Sigma.mat[[i]])
        S[temp,temp1]=temp2
        
      }
      
      if(i!=j & i<j){
        
        nu1=nu[i]
        nu2=nu[j]
        nu3=(nu[i]+nu[j])/2
        
        rho=Beta[i,j]*(gamma(nu1+3/2)/gamma(nu1))^(1/2) * (gamma(nu2+3/2)/gamma(nu2))^(1/2)*gamma(nu3)/(gamma(nu3+3/2))
        #rho <- Beta[i,j]*gamma(nu3)/sqrt(gamma(nu1)*gamma(nu2))
        
        lai=ifelse(dist0[[3]]!=0 ,(dist0[[3]]/beta)^nu3 * besselK(dist0[[3]]/beta,nu3)/(2^(nu3-1)*gamma(nu3))*sqrt(FIN_Sigma.mat[[i]] * FIN_Sigma.mat[[j]])*rho,sqrt(FIN_Sigma.mat[[i]] * FIN_Sigma.mat[[j]])*rho)
        S[temp,temp1]=lai
        #S[temp1,temp]=t(lai)
      }
      
      if(i!=j & i>j){
        
        nu1=nu[i]
        nu2=nu[j]
        nu3=(nu[i]+nu[j])/2
        
        rho=Beta[i,j]*(gamma(nu1+3/2)/gamma(nu1))^(1/2) * (gamma(nu2+3/2)/gamma(nu2))^(1/2)*gamma(nu3)/(gamma(nu3+3/2))
        #rho <- Beta[i,j]*gamma(nu3)/sqrt(gamma(nu1)*gamma(nu2))
        
        lai=ifelse(dist0[[4]]!=0 ,(dist0[[4]]/beta)^nu3 * besselK(dist0[[4]]/beta,nu3)/(2^(nu3-1)*gamma(nu3))*sqrt(FIN_Sigma.mat[[i]] * FIN_Sigma.mat[[j]])*rho,sqrt(FIN_Sigma.mat[[i]] * FIN_Sigma.mat[[j]])*rho)
        S[temp,temp1]=lai
        #S[temp1,temp]=t(lai)
      }
    }
  }
  
  return(S)
}
#--------------------------------------------------------------------------#

matern_cov_regular_grid_v2_for_estimation_sim_step1 <-function(theta,Q.mat1,Q.mat2,Q.mat3){
  q=2
  dist0 <- list()
  
  dist0[[1]] <- Q.mat1
  dist0[[2]] <- Q.mat2
  dist0[[3]] <- Q.mat3
  
  nu=theta[1:2]
  beta=theta[3]
  rot=theta[4]
  Beta=matrix(0,q,q)
  diag(Beta)=1
  Beta[2,1]=rot
  Beta[1,2]=rot
  
  var=theta[5:6]
  
  S=matrix(NA,  q*dim(dist0[[1]])[1], q*dim(dist0[[1]])[1])
  
  for(i in 1:q){
    for(j in 1:i){
      
      temp=(i-1)*dim(dist0[[1]])[1]+1:dim(dist0[[1]])[1]
      temp1=(j-1)*dim(dist0[[1]])[1]+1:dim(dist0[[1]])[1]
      
      if(i==j){
        
        temp2=ifelse(dist0[[i]]!=0,(dist0[[i]]/beta)^nu[i] * besselK(dist0[[i]]/beta,nu[i])/(2^(nu[i]-1)*gamma(nu[i])),var[i])
        #diag(temp2)=var[i]+nug[i]
        S[temp,temp1]=temp2
        
      }
      
      if(i !=j){
        
        nu1=nu[i]
        nu2=nu[j]
        nu3=(nu[i]+nu[j])/2
        
        rho=Beta[i,j]*(gamma(nu1+3/2)/gamma(nu1))^(1/2) * (gamma(nu2+3/2)/gamma(nu2))^(1/2)*gamma(nu3)/(gamma(nu3+3/2))
        
        lai=ifelse(dist0[[3]]!=0 ,(dist0[[3]]/beta)^nu3 * besselK(dist0[[3]]/beta,nu3)/(2^(nu3-1)*gamma(nu3))*sqrt(var[i]*var[j])*rho,sqrt(var[i]*var[j])*rho)
        S[temp,temp1]=lai
        S[temp1,temp]=t(lai)
      }
    }
  }
  
  return(S)
}

matern_cov_regular_grid_for_estimation_sim_v2 <-function(theta,Q.mat){
  q=2
  dist0 <- Q.mat
  
  nu=theta[1:2]
  beta=theta[3]
  rot=theta[4]
  Beta=matrix(0,q,q)
  diag(Beta)=1
  Beta[2,1]=rot
  Beta[1,2]=rot
  
  var=theta[5:6]
  
  S=matrix(NA,  q*dim(dist0)[1], q*dim(dist0)[1])
  
  for(i in 1:q){
    for(j in 1:i){
      
      temp=(i-1)*dim(dist0)[1]+1:dim(dist0)[1]
      temp1=(j-1)*dim(dist0)[1]+1:dim(dist0)[1]
      
      if(i==j){
        
        temp2=ifelse(dist0!=0,(dist0/beta)^nu[i] * besselK(dist0/beta,nu[i])/(2^(nu[i]-1)*gamma(nu[i])),var[i])
        #diag(temp2)=var[i]+nug[i]
        S[temp,temp1]=temp2
        
      }
      
      if(i !=j){
        
        nu1=nu[i]
        nu2=nu[j]
        nu3=(nu[i]+nu[j])/2
        
        rho=Beta[i,j]*(gamma(nu1+3/2)/gamma(nu1))^(1/2) * (gamma(nu2+3/2)/gamma(nu2))^(1/2)*gamma(nu3)/(gamma(nu3+3/2))
        
        lai=ifelse(dist0!=0 ,(dist0/beta)^nu3 * besselK(dist0/beta,nu3)/(2^(nu3-1)*gamma(nu3))*sqrt(var[i]*var[j])*rho,sqrt(var[i]*var[j])*rho)
        S[temp,temp1]=lai
        S[temp1,temp]=t(lai)
      }
    }
  }
  
  return(S)
}

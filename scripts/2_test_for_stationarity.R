#test_for_stationarity

directory <- '/home/salvanmo/Desktop/'

root <- paste(directory, 'Spatio-Temporal-Cross-Covariance-Functions-under-the-Lagrangian-Framework/', sep = '')

source(file = paste(root, "Functions/load_packages.R",sep=''))

compute_splag <- function(LOCS){

	N <- nrow(LOCS)
	xlags <- ylags <- matrix(, ncol = N, nrow = N)

	for(i in 1:N){
		xlags[i, ] <- LOCS[i, 1] - LOCS[, 1]
		ylags[i, ] <- LOCS[i, 2] - LOCS[, 2]
	}

	splag <- list()

	for(ll in 1:nrow(lag_targ)){
		locs_sub <- NULL

		for(i in 1:N){
			locs_temp <- which(xlags[i, ] == lag_targ[ll, 1] & ylags[i, ] == lag_targ[ll, 2])
			if(length(locs_temp) != 0) locs_sub <- rbind(locs_sub, cbind(rep(i, length(locs_temp)), locs_temp))
		}
		splag[[ll]] <- locs_sub
	}

	return(splag)

}

get_A<-function(Z,lag,tlag,index)
{
  nTime = dim(Z)[1]
  A1 = A2 = 0
  Sn = dim(lag)[1]
  for(s in 1:Sn)
  {
    for(t in index)
    {
      A1 = A1 + Z[t,lag[s,1]]*Z[t+tlag,lag[s,2]]
      A2 = A2 + Z[t+nTime/2,lag[s,1]]*Z[t+tlag+nTime/2,lag[s,2]]
    }
  }
  A1 = A1/Sn/length(index)
  A2 = A2/Sn/length(index)
  return(list(A1,A2))
}

get_G<-function(Z,splag,ulag,index)
{
  m = length(splag)
  G = integer(2*m)
  
  for(slag in 1:length(splag))
  {
    tmp = get_A(Z,splag[[slag]],ulag[slag],index);
    G[slag] = tmp[[1]]
    G[slag + m] = tmp[[2]]
  }
  
  return(G)
}

get_A_space<-function(Z1, Z2,lag1, lag2, index)
{
  nTime = dim(Z)[1]
  A1 = A2 = 0
  Sn1 = dim(lag1)[1]
  Sn2 = dim(lag2)[1]
  for(s in 1:Sn1)
  {
    for(t in index)
    {
      A1 = A1 + Z1[t, lag1[s,1]] * Z1[t, lag1[s,2]]
    }
  }
  for(s in 1:Sn2)
  {
    for(t in index)
    {
      A2 = A2 + Z2[t, lag2[s,1]] * Z2[t, lag2[s,2]]
    }
  }
  A1 = A1/Sn1/length(index)
  A2 = A2/Sn2/length(index)
  return(list(A1,A2))
}

get_G_space<-function(Z1, Z2,splag1, splag2, index)
{
  m = length(splag1)
  G = integer(2*m)
  
  for(slag in 1:length(splag1))
  {
    tmp = get_A_space(Z1, Z2,splag1[[slag]],splag2[[slag]],index);
    G[slag] = tmp[[1]]
    G[slag + m] = tmp[[2]]
  }
  
  return(G)
}

stationary.test_space<-function(Z1, Z2, splag1, splag2)
{
  m = length(splag1)
  
  nTime = dim(Z1)[1]
  index = 1:nTime 
  
  G = get_G_space(Z1, Z2, splag1, splag2, index) * sqrt(length(index))
  
  G_ = NULL

  sample_size <- 50

  #for(i in 1:(floor((nTime/2-max(ulag))/10)))
  for(i in 1:(floor(nTime/sample_size)))
  {
    index = ((i-1)*sample_size+1):(i*sample_size)
    
    G_ = cbind(G_,get_G_space(Z1, Z2, splag1, splag2,index)*sqrt(length(index)))
    
    #index = index + max(ulag)

    #G_ = cbind(G_,get_G_space(Z1, Z2, splag1, splag2,index)*sqrt(length(index)))
  }
  
  Sigma = cov(t(G_))
  
  X=cbind(diag(m),-diag(m))
  
  ts = t(X%*%G)%*%solve(X%*%Sigma%*%t(X),X%*%G)
  return(1-pchisq(ts,df=m))

}

stationary.test<-function(Z,splag,ulag)
{
  m = length(splag)
  
  nTime = dim(Z)[1]
  index = 1:(nTime/2 - max(ulag)) 
  
  G = get_G(Z,splag,ulag,index) * sqrt(length(index))
  
  G_ = NULL
  for(i in 1:(floor((nTime/2-max(ulag))/10)))
  {
    index = ((i-1)*10+1):(i*10)
    
    G_ = cbind(G_,get_G(Z,splag,ulag,index)*sqrt(length(index)))
    
    index = index + max(ulag)
    
    G_ = cbind(G_,get_G(Z,splag,ulag,index)*sqrt(length(index)))
  }
  
  Sigma = cov(t(G_))
  
  X=cbind(diag(m),-diag(m))
  
  ts = t(X%*%G)%*%solve(X%*%Sigma%*%t(X),X%*%G)
  return(1-pchisq(ts,df=m))

}

########################

lag_targ <- matrix(c(c(0.625, 0), c(1.250, 0)), ncol = 2, byrow = T)
ulag <- rep(0, nrow(lag_targ))
lag_targ <- matrix(c(c(0, 0.5), c(0, 1)), ncol = 2, byrow = T)
ulag <- rep(0, nrow(lag_targ))
lag_targ <- matrix(c(c(0.625, 0), c(0, 0.5), c(0.625, 0.5)), ncol = 2, byrow = T)
ulag <- rep(0, nrow(lag_targ))
lag_targ <- matrix(c(c(0.625, 0), c(0, 0.5), c(0.625, 0.5)), ncol = 2, byrow = T)
ulag <- rep(1, nrow(lag_targ))

PVALS <- matrix(, ncol = 2, nrow = 40)
for(yr in 1980:2019){
	cat(yr, '\n')
	dat <- read.table(paste(root, 'Data/ncdf/layer1_residuals_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
	dat2 <- read.table(paste(root, 'Data/ncdf/layer2_residuals_' , yr, sep = ''), header = FALSE, sep = " ") %>% as.matrix()
	dat3 <- read.table(paste(root, 'Data/ncdf/LOCS-3D-dataset', sep = ''), header = FALSE, sep = " ") %>% as.matrix()

	Z <- dat
	Z2 <- dat2

	for(zz in 1:nrow(Z)){
		Z[zz, ] <- Z[zz, ] - mean(Z[zz, ])
		Z2[zz, ] <- Z2[zz, ] - mean(Z2[zz, ])
	}

	splag <- compute_splag(dat3)

	pvalues1 <- stationary.test(Z, splag, ulag)
	pvalues2 <- stationary.test(Z2, splag, ulag)

	PVALS[yr - 1979, 1] <- pvalues1
	PVALS[yr - 1979, 2] <- pvalues2
}

ind1 <- which(dat3[, 1] <= 44)

ind1 <- which(dat3[, 2] <= 23.5)

locs1 <- dat3[ind1, ]
locs2 <- dat3[-ind1, ]


ZA <- Z[, ind1]
ZB <- Z[, -ind1]
Z2A <- Z2[, ind1]
Z2B <- Z2[, -ind1]
for(zz in 1:nrow(Z)){
	ZA[zz, ] <- ZA[zz, ] - mean(ZA[zz, ])
	ZB[zz, ] <- ZB[zz, ] - mean(ZB[zz, ])
	Z2A[zz, ] <- Z2A[zz, ] - mean(Z2A[zz, ])
	Z2B[zz, ] <- Z2B[zz, ] - mean(Z2B[zz, ])
}
stationary.test_space(ZA, ZB, splag1, splag2)
stationary.test_space(Z2A, Z2B, splag1, splag2)

splag1 <- compute_splag(locs1)
splag2 <- compute_splag(locs2)
stationary.test(ZA, splag1, ulag)
stationary.test(ZB, splag2, ulag)


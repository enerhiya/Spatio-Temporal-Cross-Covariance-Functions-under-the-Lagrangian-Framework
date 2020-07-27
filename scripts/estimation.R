directory <- '/home/salvanmo/Desktop/'

root <- paste(directory, 'Spatio-Temporal-Cross-Covariance-Functions-under-the-Lagrangian-Framework/', sep = '')

source(file = paste(root, "Functions/load_packages.R", sep = ''))
source(file = paste(root, "Functions/auxiliary_functions.R", sep = ''))
source(file = paste(root, "Functions/cov_func.R", sep = ''))
source(file = paste(root, "Functions/plot_func.R", sep = ''))

load_default_data = F

if(load_default_data){

	cov1 <- read.table(paste(root, "Data/frozen_matern_cov_for_heatmap", sep = ''), header = FALSE, sep = " ") %>% as.matrix() 
	r1 <- read.table(paste(root, "Data/frozen_matern_realizations", sep = ''), header = FALSE, sep = ",") %>% as.matrix() 
	locs <- read.table(paste(root, "Data/simulation_locs", sep = ''), header = FALSE, sep = " ") %>% as.matrix() 
	TT <- max(unique(locs[, 3]))
	n <- nrow(locs)/TT	


}else{

	N <- 20
	n <- N^2
	t <- 10
	grid_x <- seq(from = 0, to = 1, length.out = N)
	grid_y <- seq(from = 0, to = 1, length.out = N)
	sim_grid_locations <- expand.grid(grid_x, grid_y) %>% as.matrix()

	wind <- c(0.1001, 0.1001)
	theta <- c(1, 1, 0.23, 0.5, 1, 0.8)

	cov1 <- frozen_matern_cov_rep_I(theta, wind, max_time_lag = t - 1, LOCS = sim_grid_locations)

	set.seed(1235)
	r1 <- rmvn(1, rep(0, ncol(cov1)), cov1, ncores = 25)

	ind_var1 <- seq(1, ncol(cov1), 2)
	ind_var2 <- seq(2, ncol(cov1), 2)
	
	jpeg(file = paste(root, 'Figures/frozen_matern_realizations.jpg', sep = ''), width = 1600, height = 800)

        par(mfrow = c(2, 2))

        for(tt in 1:2){

                image.plot(matrix(r1[1, ind_var1[(tt - 1) * n + 1:n]], ncol = N, nrow = N))
                image.plot(matrix(r1[1, ind_var2[(tt - 1) * n + 1:n]], ncol = N, nrow = N))

        }

        dev.off()


        locs_rand_sample <- sim_grid_locations
        Z_rand_sample <- r1
        #Z_rand_sample <- matrix((Z_rand_sample - mean(Z_rand_sample)) / sd(Z_rand_sample), ncol = 1)
        #Z_rand_sample <- matrix(Z_rand_sample - mean(Z_rand_sample), ncol = 1)
	
        NEGLOGLIK <- function(p){
	  	theta <- c(exp(p[1:5]), p[6])
		wind <- p[7:8]		

          	print(c(theta, wind))

          	Sigma <- frozen_matern_cov(theta, wind, max_time_lag = t - 1, LOCS = sim_grid_locations)
          	Sigma.inv <- solve(Sigma)
          	det.Sigma <- determinant(Sigma,log=T)$mod[1]
          	out <- 0.5 * det.Sigma + 0.5 * t(Z_rand_sample) %*% Sigma.inv %*% Z_rand_sample
          	return(out)
        }

	init <- c(log(sd(Z_rand_sample[1:(n * t),])), log(sd(Z_rand_sample[n * t + 1:(n * t),])), log(0.1), log(0.5), log(1), cor(Z_rand_sample[1:(n * t),], Z_rand_sample[n * t + 1:(n * t),]), 0.05, 0.05)

        fit1 <- optim(par = init, fn = NEGLOGLIK, control = list(trace = 5, maxit = 10000)) #

        for(rep in 1:100){
          fit1 <- optim(par = fit1$par, fn = NEGLOGLIK, control = list(trace = 5))
        }

	cut_off_t <- 3

        FULL_VECCHIA_NEGLOGLIK <- function(p){
	  	theta <- c(exp(p[1:5]), p[6])
		wind <- p[7:8]		

          	print(c(theta, wind))

          	Sigma <- frozen_matern_cov_rep_I(theta, wind, max_time_lag = cut_off_t - 1, LOCS = sim_grid_locations)
	
		R11 <- Sigma[1:(n * q * (cut_off_t - 1)), 1:(n * q * (cut_off_t - 1))]
		R12 <- Sigma[(n * q * (cut_off_t - 1)) + 1:(n * q), 1:(n * q * (cut_off_t - 1))]
		#R12 <- Sigma[1:(n * q * (cut_off_t - 1)), (n * q * (cut_off_t - 1)) + 1:(n * q)]
		R22 <- Sigma[(n * q * (cut_off_t - 1)) + 1:(n * q), (n * q * (cut_off_t - 1)) + 1:(n * q)]

		Z_full <-  matrix(Z_rand_sample[1:(n * q * cut_off_t)] - mean(Z_rand_sample[1:(n * q * cut_off_t)]), ncol = 1)

		Sigma.inv <- solve(Sigma)
                det.Sigma <- determinant(Sigma,log=T)$mod[1]
                out <- 0.5 * det.Sigma + 0.5 * t(Z_full) %*% Sigma.inv %*% Z_full

		for(tt in 1:(t - cut_off_t)){
			Z_cond <-  matrix(Z_rand_sample[n * (tt - 1) + 1:(n * q * (cut_off_t - 1))] - mean(Z_rand_sample[n * (tt - 1) + 1:(n * q * (cut_off_t - 1))]), ncol = 1)
			Z_forward <- matrix(Z_rand_sample[n * (tt - 1) + (n * q * (cut_off_t - 1)) + 1:(n * q)] - mean(Z_rand_sample[n * (tt - 1) + (n * q * (cut_off_t - 1)) + 1:(n * q)]), ncol = 1)

			C_t <- R12 %*% solve(R11) 
			V_t <- R22 - C_t %*% t(R12) 
		
			err <- Z_forward - C_t %*% Z_cond
			#sum(err^2)

			cholmat <- tryCatch(t(chol(V_t)) , error = function(a) numeric(0) )
    			if( length(cholmat) == 0 ){
        			# sometimes the covmat is not numerically positive definite.
        			# we probably need a better solution for this.
        			return(-9999999)
        			print("One of the Choleskys failed")
    			}

			z       <- forwardsolve(cholmat, err)
		    	# sum of log conditional variances
		    	logsig  <- 2 * sum(log(diag(cholmat)))
		    	# simply add them up to get conditional loglikelihood
		    	out  <- out + length(z)/2 * log(2 * pi) + 1/2 * logsig + 1/2 * sum(z^2)
						
			#Sigma.inv <- solve(V_t)
			#det.Sigma <- determinant(V_t,log=T)$mod[1]
			#out <- out + 0.5 * det.Sigma + 0.5 * t(err) %*% Sigma.inv %*% err
		}
          	return(out)
        }


	init <- c(log(1), log(1), log(0.3), log(0.5), log(1), cor(Z_rand_sample[ind_var1], Z_rand_sample[ind_var2]), 0.05, 0.05)
        fit1 <- optim(par = init, fn = FULL_VECCHIA_NEGLOGLIK, control = list(trace = 5, maxit = 1000)) #

        for(rep in 1:100){
          fit1 <- optim(par = fit1$par, fn = FULL_VECCHIA_NEGLOGLIK, control = list(trace = 5))
        }
######## --------------------------------------------------------- #########



library(psych)
######## ------  CONDITIONAL LOGLIKELIHOOD  ------ #######

	ZZ <- Z_rand_sample[1:(n * 3 * 2)] %*% t(Z_rand_sample[1:(n * 3 * 2)])
	
	for(tt in 1:17){

		ZZ <- ZZ + Z_rand_sample[n * tt * 2 + 1:(n * 3 * 2)] %*% t(Z_rand_sample[n * tt * 2 + 1:(n * 3 * 2)])

	}

	G <- chol(ZZ)

        COND_NEGLOGLIK <- function(p){
	  	theta <- c(exp(p[1:5]), p[6])
		wind <- p[7:8]		

          	print(c(theta, wind))

		#t = 2, i.e., use two previous timesteps
		ttt <- 2

          	Sigma <- frozen_matern_cov_rep_I(theta, wind, max_time_lag = ttt, LOCS = sim_grid_locations)
		#r1 <- rmvn(1, rep(0, ncol(Sigma)), Sigma, ncores = 25)

		n <- nrow(sim_grid_locations)

		conditional_ind <- 1:(n * ttt * 2)
		R11_b <- Sigma[conditional_ind, conditional_ind]
		
		unobserved_ind <- (n * ttt * 2) + 1:(n * 2)
		R12_b <- Sigma[conditional_ind, unobserved_ind]
		R22 <- Sigma[unobserved_ind, unobserved_ind]

		#pred <- C_t %*% matrix(r1[1, conditional_ind], ncol = 1)	
		#err <- pred - matrix(r1[1, unobserved_ind], ncol = 1)	
		
		R11_f <- Sigma[n * 2 + conditional_ind, n * 2 + conditional_ind]
		R12_f <- Sigma[1:(n * 2), n * 2 + conditional_ind]

		C_b <- t(R12_b) %*% solve(R11_b)
		C_f <- R12_f %*% solve(R11_f)
		
		a_f <- R22 - C_b %*% R12_b
		a_b <- R22 - C_f %*% t(R12_f)

		A <- B <- A_temp <- B_temp <- diag(2 * n) 

		R <- Sigma[1:(n * 2), 1:(n * 2)]		

		beta <-	t(B) %*% R

		k_f <- solve(a_f) %*% beta 	
		k_b <- solve(a_b) %*% t(beta) 	
		
		ZERO_MAT <-  matrix(0, ncol = n * 2, nrow = n * 2)

		for(tt in 1:2){
		
			A <- rbind(A_temp, ZERO_MAT) + rbind(ZERO_MAT, B_temp) %*% k_f
			B <- rbind(ZERO_MAT, B_temp) + rbind(A_temp, ZERO_MAT) %*% k_b

			A_temp <- A
			B_temp <- B	
			a_f <- a_f + beta %*% k_b	
			a_b <- a_b + t(beta) %*% k_f	
		}

		V_t <- a_f
		C_t <- A

		H <- tryCatch(t(chol(V_t)) , error = function(a) numeric(0) )
		if( length(H) == 0 ){
			# sometimes the covmat is not numerically positive definite.
			# we probably need a better solution for this.
			return(-9999999)
			print("One of the Choleskys failed")
		}

		det_V <- 2 * sum(log(diag(H)))

		third_term <- tr(B %*% solve(V_t)  %*% t(B) %*% ZZ) 

		loglik <- n * (10 - 2) * log(2 * pi) / 2 + (10 - 2) / 2 * det_V + 0.5 * 1

          	return(loglik)
        }
        
	fit2 <- optim(par = init, fn = COND_NEGLOGLIK, control = list(trace = 5, maxit = 100)) #
	fit2 <- optim(par = c(rep(log(1), 5), 0.1, 0.2, 0.4), fn = COND_NEGLOGLIK, control = list(trace = 5, maxit = 10000)) #

mvnMultiCondLik <- function( covparms, covfun, y, locs, likinds, returnz = FALSE ){

    # computes the ordered conditional multivariate normal density for the
    # data in the likinds position of y

    # for example, if y has length 4, and likinds = c(2,4),
    # this returns the log density for the second observation
    # given the first plus the log density for the fourth
    # observation given the first three

    leny    <- length(y)
    distmat <- rdist(locs,locs)
    covmat  <- covfun( locs, covparms, returnD1 = FALSE )
    cholmat <- tryCatch(t(chol(covmat)) , error = function(a) numeric(0) )
    if( length(cholmat) == 0 ){
        # sometimes the covmat is not numerically positive definite.
        # we probably need a better solution for this.
        return(-9999999)
        print("One of the Choleskys failed")
    }

    # decorrelating transform
    z       <- forwardsolve(cholmat,y)
    # sum of log conditional variances
    logsig  <- 2*sum(log(diag(cholmat)[likinds]))
    # simply add them up to get conditional loglikelihood
    loglik  <- -length(likinds)/2*log(2*pi) - 1/2*logsig - 1/2*sum(z[likinds]^2)
}
}

        FULL_NEGLOGLIK <- function(p){
	  	theta <- c(exp(p[1:5]), p[6])
		wind <- p[7:8]		

          	print(c(theta, wind))

          	Sigma <- frozen_matern_cov(theta, wind, max_time_lag = t - 1, LOCS = sim_grid_locations)
		cholmat <- tryCatch(t(chol(Sigma)) , error = function(a) numeric(0) )
		if( length(cholmat) == 0 ){
			# sometimes the covmat is not numerically positive definite.
			# we probably need a better solution for this.
			return(-9999999)
			print("One of the Choleskys failed")
		}

    		z <- forwardsolve(cholmat, Z_rand_sample)
    		
		# sum of log conditional variances

		loglik <- 0

		for(cc in 1:9){

			lastk <- n

			inds2 <- n * cc + 1:n
			logsig <- 2*sum(log(diag(cholmat)[inds2]))
			loglik  <- loglik + lastk/2*log(2*pi) + 1/2*logsig + 1/2*sum(z[inds2]^2)
		}

    		
    		# simply add them up to get conditional loglikelihood

          	return(loglik)
        }
#2000 iterations
#p <- c(-0.32351046, -0.30254643, -1.64192653, -0.73872861,  0.01612647,  0.80238952,  0.09987040,  0.09982339)
#0.72360439 0.73893418 0.19360669 0.47772090 1.01625721 0.80238952 0.09987040 0.09982339

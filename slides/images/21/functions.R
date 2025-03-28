# functions to be used to tune hyperparameters for HSGP under matern 32 kernel
m_QE <- function(c,rho,S) ceiling(3.42 * c / (rho/S))
rho_QE <- function(c,m,S) round(S * 3.42 * c / m, 3) 

c_vs_rho_QE <- function(rho,S){
  c =  4.5*rho/S
  if(c < 1.2)
    c = 1.2
  c
}

diagnostic <- function(rho,rho_hat) rho_hat + 0.01 >= rho

# helper functions to run HSGP 

runHSGP <- function(rs,mod,i,seed=1,testshort=F,testfull=F){
  
  if (testshort){
    t_iter <- 100; t_burn <- 50
    refresh <- max(10, t_iter/10)
  } else {
    t_iter <- rs$niter; t_burn <- rs$nburn
    
    if (testfull){
      refresh <- max(10, t_iter/10)
    } else {
      refresh <- 0
    }
  }

  HSGP_data <- list(N=dim(data)[1],q=dim(rs$u_new)[1],p=dim(rs$X)[2],Y=rs$Y,
               n=dim(rs$u)[1],X=rs$X,Ids=rs$Ids,u=rs$u,u_new=rs$u_new,L=rs$c[[i]]*rs$S,M=rs$m[[i]])
 
  rs$fit[[i]] <- sampling(mod, data = HSGP_data, init=0.5, chains = rs$nchain, 
                          warmup = t_burn,iter = t_iter, core=rs$nchain, 
                          thin = rs$nthin, seed = seed, include=T, refresh=refresh,
                          pars=c("beta","alpha","sigma","tau","rho","theta_new","tau_adj")) 

  return(rs)
}

updateHSGP <- function(rs,i){
  S <- rs$S
  d <- 2
  
  # we assume each GP has only one lengthscale parameter, for all the dimensions
  # but HSGP setup would require one lengthscale for each dimension during configuration
  # because the length of each dimension could be different
  
  ci <- mi <- rhoij <- numeric(d)
  check <- F
  
  if (i==1){
    rhoi <- rs$rho0
    
    for (j in 1:d){
      ci[j] <- c_vs_rho_QE(rho=rhoi, S=S[j])
      mi[j] <- m_QE(c=ci[j], rho=rhoi, S=S[j]) 
    }
    
  } else {
    
    rhoi_hat <- round(summary(rs$fit[[i-1]], pars ="rho")$summary[,1], 2)
    diagi <- diagnostic(rs$rho[[i-1]], rhoi_hat)
    
    if (i==2){
      
      if (diagi){
        mi <- rs$m[[i-1]] + 2
        for (j in 1:d){
          ci[j] <- c_vs_rho_QE(rho=rhoi_hat, S=S[j])
          rhoij[j] <- rho_QE(ci[j], m=mi[j], S=S[j])
        }
        rhoi <- max(rhoij)
      } else {
        rhoi <- rhoi_hat
        for (j in 1:d){
          ci[j] <- c_vs_rho_QE(rho=rhoi, S=S[j])
          mi[j] <- m_QE(c=ci[j], rho=rhoi, S=S[j])
        }
      }
      
    } else {
      
      if (diagi & !rs$diag[[i-2]]){
        mi <- rs$m[[i-1]] + 2
        for (j in 1:d){
          ci[j] <- c_vs_rho_QE(rho=rhoi_hat, S=S[j])
          rhoij[j] <- rho_QE(ci[j], m=mi[j], S=S[j])
        }
        rhoi <- max(rhoij)
      } else if (diagi & rs$diag[[i-2]]){
        check <- T
        rhoi <- rs$rho[[i-1]]
        ci <- rs$c[[i-1]]
        mi <- rs$m[[i-1]]
      } else if(!diagi){
        rhoi <- rhoi_hat
        for (j in 1:d){
          ci[j] <- c_vs_rho_QE(rho=rhoi, S=S[j])
          # for the first 5 runs before things stabilize, we don't apply the flooring
          if (i>5){
            mi[j] <- max(m_QE(c=ci[j], rho=rhoi, S=S[j]), rs$m[[i-1]][j])
          } else {
            mi[j] <- m_QE(c=ci[j], rho=rhoi, S=S[j])
          }
        }       
      }
      
    }
    
    rs$rho_hat[[i-1]] <- rhoi_hat
    rs$diag[[i-1]] <- diagi
  }
  
  if (i==1){
    rs$fit <- list()
    rs$rho <- list() # lengthscale
    rs$c <- list()  # boundary factor
    rs$m <- list()  # number of basis functions
    rs$rho_hat <- list()  # HSGP lengthscale estimate
    rs$diag <- list()
    rs$check <- list()
  }
  rs$check[[i]] <- check
  rs$rho[[i]] <- rhoi
  rs$c[[i]] <- ci
  rs$m[[i]] <- mi
  
  return(rs)
}


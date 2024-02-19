
# Loading required libraries:
library("MASS")
library("stochvol")
library("dbarts")

############## Aux. Functions  #####################

# lag variables
mlag <- function(X,lag){
  p <- lag
  X <- as.matrix(X)
  Traw <- nrow(X)
  N <- ncol(X)
  Xlag <- matrix(NA,Traw,p*N)
  for (ii in 1:p){
    Xlag[(p+1):Traw,(N*(ii-1)+1):(N*ii)]=X[(p+1-ii):(Traw-ii),(1:N)]
  }
  return(Xlag)
}

get.hs <- function(bdraw,lambda.hs,nu.hs,tau.hs,zeta.hs){
  k <- length(bdraw)
  if (is.na(tau.hs)){
    tau.hs <- 1   
  }else{
    tau.hs <- invgamma::rinvgamma(1,shape=(k+1)/2,rate=1/zeta.hs+sum(bdraw^2/lambda.hs)/2) 
  }
  
  lambda.hs <- invgamma::rinvgamma(k,shape=1,rate=1/nu.hs+bdraw^2/(2*tau.hs))
  
  nu.hs <- invgamma::rinvgamma(k,shape=1,rate=1+1/lambda.hs)
  zeta.hs <- invgamma::rinvgamma(1,shape=1,rate=1+1/tau.hs)
  
  ret <- list("psi"=(lambda.hs*tau.hs),"lambda"=lambda.hs,"tau"=tau.hs,"nu"=nu.hs,"zeta"=zeta.hs)
  return(ret)
}


############ Main ###################

VARBART = function(Yraw, p, fhorz, nburn = 5000, nsave = 5000, thinfac = 1, sv = TRUE){
  
  ############# Creating data matrices ##################
  
  nthin = round(thinfac * nsave)
  ntot = nburn + nsave
  thin.set = floor(seq(nburn+1,ntot,length.out=nthin))
  
  iter.update = 250
  in.thin  = 0
  
  Ymu = apply(Yraw, 2, mean,na.rm=T)
  Ysd = apply(Yraw, 2, sd,na.rm=T)
  Yraw = apply(Yraw, 2, function(x){(x-mean(x,na.rm=T))/sd(x,na.rm=T)})
  
  # create design matrices Y/X
  
  X = cbind(mlag(Yraw,p))[(p+1):nrow(Yraw),]
  
  Y = Yraw[(p+1):nrow(Yraw), ] # leave original input matrix unchanged
  
  colnames(X) = paste(rep(colnames(Yraw), p),
                      sort(rep(paste(".t-", 1:p, sep = ""),ncol(Yraw))), sep = "" )
  
  
  N = ncol(Yraw)
  TT = nrow(Y)
  K  = ncol(X)
  
  ############### Prior Parameters ##################
  
  # The model is the Following :
  # Yt = F(X) + e A0
  # Yt: TxN, where N is the number of End.
  # X: TxK dimmension, K = Np 
  
  ## Initializing the coefficients and Sd
  PAI.OLS = solve(crossprod(X)) %*% crossprod(X,Y) 
  Sig.OLS = crossprod(Y-X%*%PAI.OLS)/TT
  
  # Now the time varying varyances, this will be used for the forecasting!
  Sig.t = array(0, dim= c(TT,N,N))
  for(tt in 1:TT) Sig.t[tt,,] = Sig.OLS
  
  #Initialization of Covariance Objects:
  
  th.A0  = matrix(1,N,N)
  eta = matrix(NA, TT, N)
  H = matrix(0, TT, N)
  d.A0 = diag(N)
  
  # Stochastic Volatility:
  sv_priors = list()
  sv_draw = list()
  h_latent = list()
  
  if(sv){
    for(nn in 1:N){
      sv_draw[[nn]] = list(mu = 0, phi = 0.99, sigma = 0.01, nu = Inf, rho = 0, beta = NA, latent0 = 0)
      h_latent[[nn]] = rep(0,TT)
      sv_priors[[nn]] = stochvol::specify_priors(
        mu     = sv_normal(mean =0, sd = 10),
        phi    = sv_beta(shape1 =5 , shape2 = 1.5),
        sigma2 = sv_gamma(shape = 0.5, rate = 10),
        nu     = sv_infinity(),
        rho    = sv_constant(0))
    }
  }else{
    for(mm in 1:N){
      sv_draw[[mm]] = list(mu = 0, phi = 0.99, sigma = 0.01, nu = Inf, rho = 0, beta = NA, latent0 = 0)
      h_latent[[mm]] = rep(0,TT)
      sv_priors[[mm]] <- specify_priors(
        mu = sv_constant(0),
        phi = sv_constant(1-1e-12),
        sigma2 = sv_constant(1e-12),
        nu = sv_infinity(),
        rho = sv_constant(0)
      )
    }
    
  }
  
  
  # Bart Initialization 
  
  cgm.level=0.95
  cgm.exp=2
  sd.mu=2
  num.trees=250
  prior.sig = c(nrow(Y)/2, 0.75)
  
  control = dbartsControl(verbose = FALSE, keepTrainingFits = TRUE,
                          useQuantiles = FALSE,
                          keepTrees = FALSE, n.samples = ntot,
                          n.cuts = 100L, n.burn = nburn,
                          n.trees = num.trees, n.chains = 1,
                          n.threads = 1, n.thin = 1L, 
                          printEvery = 1,
                          printCutoffs = 0L, rngKind = "default",
                          rngNormalKind = "default",
                          updateState = FALSE)
  
  bart.sampler.ls = list()
  svdraw.ls = list()
  
  for(jj in 1:N){
    sigma.init = 1
    bart.sampler.ls[[jj]] = dbarts(Y[,jj]~X, control = control,
                                   tree.prior = cgm(cgm.exp, cgm.level),
                                   node.prior = normal(sd.mu),
                                   n.samples = nsave, weights=rep(1,TT),
                                   sigma=sqrt(Sig.OLS[jj,jj]),
                                   resid.prior = chisq(prior.sig[[1]],
                                                       prior.sig[[2]]))
  }
  bart.sampler.run = list()
  sigma.mat = matrix(NA, N, 1)
  count.mat = matrix(0, K, N)
  row.names(count.mat) = colnames(X)
  colnames(count.mat) = colnames(Yraw)
 
  # Init. for the Horseshoe:
  lambda.A0 = 1
  nu.A0 = 1
  tau.A0 =1
  zeta.A0 =1
  prior.cov = rep(1, N*(N-1)/2)


  sv_parms_mat <- matrix(NA, N, 4)
  
  Y.fit.BART = Y*0
  
  # Storage objects
  Y.store = array(NA, dim=c(nthin, TT, N))
  H_store = array(NA, dim=c(nthin, TT,N))
  
  forecast.store = array(NA, dim = c(nthin,fhorz,N ))
  Hfcst_store =  array(NA, dim= c(nthin,fhorz,N))
  
  # Starting the Gibbs Sampling!
  
  pb = txtProgressBar(min = 0, max = ntot, style = 3)
  start = Sys.time()
  
  for(nrep in 1:ntot){
    for(nn in 1:N){
      
      if(nn >1){
        eta_nn = eta[,1:(nn -1), drop = FALSE]
        A0_nn  = d.A0[nn,1:(nn-1)]
        bart.sampler.ls[[nn]]$setResponse(Y[,nn] - eta_nn%*%A0_nn)
      }
      
      rep.nn = bart.sampler.ls[[nn]]$run(0L,1L)
      bart.sampler.run[[nn]] = rep.nn
      sigma.mat[nn,] = rep.nn$sigma
      if(any(is.na(rep.nn$train))){
        break
      }
      Y.fit.BART[,nn] = rep.nn$train
      eta[,nn] = Y[,nn] - rep.nn$train
      count.mat[,nn] = rep.nn$varcount
      
      if(nn >1){
        norm.nn = as.numeric(exp(-.5*h_latent[[nn]]) * 1/sigma.mat[nn,])
        u_nn    = eta[,1:(nn-1), drop =F]*norm.nn
        eta.nn  = eta[,nn]*norm.nn
        if (nn == 2){
          v0.inv = 1/th.A0[nn,1]
        }else{
          v0.inv = diag(1/th.A0[nn,1:(nn-1)])
        } 
        V.cov = solve(crossprod(u_nn) + v0.inv)
        mu.cov = V.cov %*% crossprod(u_nn, eta.nn)
        mu.cov.draw = mu.cov + t(chol(V.cov)) %*% rnorm(ncol(V.cov)) 
        d.A0[nn,1:(nn-1)] = mu.cov.draw
      }
    }
    
    shocks = eta %*%t(solve(d.A0))
    if(sv){
      for (mm in 1:N){
        svdraw_mm =  svsample_general_cpp(shocks[,mm]/sigma.mat[mm], 
                                          startpara = sv_draw[[mm]], startlatent = h_latent[[mm]],
                                          priorspec = sv_priors[[mm]])
        sv_draw[[mm]][c("mu", "phi", "sigma")] = as.list(svdraw_mm$para[, c("mu", "phi", "sigma")])
        h_latent[[mm]] = svdraw_mm$latent
        sv_parms_mat[mm, ] = c(svdraw_mm$para[, c("mu", "phi", "sigma")], svdraw_mm$latent[TT])
        normalizer = as.numeric(exp(-.5*svdraw_mm$latent))
        weights.new = as.numeric(exp(-svdraw_mm$latent))
        dat = dbartsData(formula = Y[,mm]~X,weights=weights.new)
        bart.sampler.ls[[mm]]$setData(dat)
        H[,mm] = log(sigma.mat[mm]^2) + svdraw_mm$latent
      }
    }else{
      H[,mm] <- log(sigma.mat[mm]^2)
    }
    
    H[H<log(1e-6)] = log(1e-6)
    
    for(tt in 1:TT){
      aux.s = exp(H[tt,])
      s.t   = t(d.A0)%*%crossprod(diag(aux.s), d.A0)
      Sig.t[tt,,] = s.t 
    }
    # Updating the Shrinkage priors
    hs_draw = get.hs(bdraw=d.A0[lower.tri(d.A0)],lambda.hs=lambda.A0,nu.hs=nu.A0,tau.hs=tau.A0,zeta.hs=zeta.A0)
    lambda.A0 = hs_draw$lambda
    nu.A0  = hs_draw$nu
    tau.A0 = hs_draw$tau
    zeta.A0   = hs_draw$zeta
    prior.cov = hs_draw$psi
    
    th.A0[lower.tri(th.A0)] = prior.cov
    th.A0[th.A0>10] = 10
    th.A0[th.A0<1e-8] = 1e-8
    
    #th.A[th.A < 1e-12] <- 1e-12

    ################## Creating and Updating the Forecast ####################  
    
    if(nrep %in% thin.set){
      
      in.thin = in.thin +1
      H_store[in.thin,,] = H
      
      if(fhorz>0){
        Yfc = matrix(NA, fhorz, N)
        Hfc = matrix(NA,fhorz,N)
        
        X.hat = c( Y[TT,], X[TT,1:(N*(p-1))] )
        if(sv == TRUE){
          HT = H[TT,] - log(as.numeric(sigma.mat)^2)  
        }else{
          HT = H[TT,]
        }
        Sig.T = Sig.t[TT,,]
        tree.pred = matrix(0,N)
        for(hh in 1:fhorz){
          HT  = log(as.numeric(sigma.mat)^2) + (sv_parms_mat[, 1] + sv_parms_mat[ , 2] * (HT - sv_parms_mat[,1]) + sv_parms_mat[ , 3]*rnorm(N))
          Hfc[hh,] = exp(HT)
          for(j in 1:N){
            tree.pred[j] = bart.sampler.ls[[j]]$predict(X.hat)
          }
          #HT  = log(as.numeric(sigma.mat)^2) + (sv_parms_mat[, 1] + sv_parms_mat[ , 2] * (HT - sv_parms_mat[,1]) + sv_parms_mat[ , 3]*rnorm(N))
          #Sig.T = diag(exp(HT))
          if(sv == T){
            Sig.T = t(d.A0)%*%crossprod(diag(exp(HT)), d.A0)
          }
          Y.tp1 = as.numeric(tree.pred) + t(chol(Sig.T)) %*%rnorm(N)
          X.hat = c(Y.tp1, X.hat[1:(N*(p-1))])
          Yfc[hh,] = Y.tp1
          
        }
        
        forecast.store[in.thin,,] <- (Yfc*t(matrix(Ysd,N,fhorz)))+t(matrix(Ymu,N,fhorz))
        Hfcst_store[in.thin,,] <-  Hfc*t(matrix(Ysd,N,fhorz))
        
      }
      
    }
    
    setTxtProgressBar(pb, nrep)
    if (nrep %% iter.update==0){
      end =  Sys.time()
      message(paste0("\n Average time for single draw over last ",iter.update," draws ", round(as.numeric(end-start)/iter.update, digits=4), " seconds, currently at draw ", nrep))
      start = Sys.time() 
      
      par(mfrow = c(3,3))
      for(i in 1:N){
        ts.plot(cbind(Y.fit.BART[,i], Y[,i]), col=c("red","black") ,xlab = "", main = colnames(Yraw)[i])
      }
    }
  }
  dimnames(forecast.store) = list(paste0("mcmc",1:nthin),paste0("t+", 1:fhorz),colnames(Y))
  dimnames(Hfcst_store) =  dimnames(forecast.store)
  return_obj = list("fcst"=forecast.store,"Hfcst"=Hfcst_store, "Y" = Y, "Ysd" = Ysd, "Ymu" = Ymu, "tree.splits" = count.mat, "Y.fit.BART"= Y.fit.BART)
  return(return_obj)
}

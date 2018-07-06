

require(tclust)

m_t_clust = function(data,alpha.fixed,K,modelNames,iter.max=50,restr.factor=50,alpha.init=.20,restr.fact.init=50,n.start.init=100,initial.values = NULL){
  
  #Common values for each model
  
  
  
  n         = nrow(data)
  p         = ncol(data)
  
  #parameters to be consdered in the while loop
  
  loglik    = c()
  epsilon   = 1e-10
  iter      = 2
  
  war       = 0
  #init=tkmeans(data,k=K,alpha=0.20)
  
  if(is.null(initial.values)){
    init=tclust(data,k=K,alpha=alpha.init,restr.fact = restr.fact.init,nstart = n.start.init)
    
    init$cov  = array(dim=c(p,p,K))
    init$pro  = table(init$cluster)
    
    for(j in 1:K){
      init$cov[,,j]=cov(data[init$cluster==j,])
    }
    
    #fills the vector containing all the volumes
    sigmasq   =c()
    
    for(j in 1:K){
      sigmasq[j] = (det(init$cov[,,j]))^(1/p)
    }
    
    init$sigmasq = sigmasq
    #fills the matrix with the eigen values 
    shape =matrix(nrow=p,ncol=K)
    
    for(k in 1:K){
      shape[,k] = eigen(init$cov[,,k]/sigmasq[k])$values
    }
    init$shape = shape
    #fills the array with the matrix of the eigenvectors
    orientation = array(dim=c(p,p,K))
    
    for(k in 1:K){
      orientation[,,k] = eigen(init$cov[,,k])$vectors
    }
    
    init$orientation = orientation
    #fills the array with the Cholesky decomposition
    
    cholsigma = array(,dim=c(p,p,K))
    
    for(k in 1:K){
      cholsigma[,,k] = chol(init$cov[,,k])
    }
    init$cholsigma = cholsigma
    
  }else{init=initial.values
  sigmasq   = init$sigmasq
  shape     = init$shape
  cholsigma = init$cholsigma
  orientation=init$orientation
  
  }
  
  if(modelNames == "EII"){
    
    parameters     = list(mean=init$centers,
                          variance=list(sigmasq=det(cov(data[init$cluster!=0,]))^(1/p),
                                        Sigma=cov(data[init$cluster!=0,])/(det(cov(data[init$cluster!=0,]))^(1/p)),
                                        pro=init$weights))
    
    #E-step
    assig = unmap(init$clus[init$clus!=0])   
    #M-step
    mstep          = try(mstepEII(data[init$clus!=0,],z=assig,warn=T))
    parameters     = mstep$parameters
    
    
    #E-step
    discriminants  = d_function(data,parameters,"EII",K)$discriminants  
    trimming       = trim(data,K,D=discriminants,alpha.fixed)
    class          = trimming$classification
    cont           = sum(log(trimming$contribution[trimming$class!=0]))
    loglik[1]      = cont
    
    #Performs trimming on the observations having lowest contribution
    
    mstep          = try(mstepEII(data[class!=0,],z=unmap(class[class!=0]),warn=T))
    parameters     = mstep$parameters
    discriminants =  d_function(data,parameters,"EII",K)$discriminants  
    trimming       = trim(data,K,D=discriminants,alpha.fixed)
    class          = trimming$classification
    cont           = sum(log(trimming$contribution[class!=0]))
    loglik[2]      = cont
    
    while(  (loglik[iter]-loglik[iter-1]) >epsilon){#fix epsilon and run
      
      #Performs EM algorithm
      mstep          = try(mstepEII(data[class!=0,],z=unmap(class[class!=0]),warn=T))
      parameters     = mstep$parameters
      
      #Update the trimming
      D        = d_function(data,parameters=parameters,modelNames = "EEI",K)
      trimming = trim(data,K,D=D$discriminants,alpha.fixed)
      class    = trimming$classification
      cont     = sum(log(trimming$contribution[class!=0]))
      
      loglik[iter+1]   = cont
      if(is.na(loglik[iter+1]))(loglik[iter+1]=loglik[1])
      iter=iter+1
      if(iter==iter.max) break      
    }
    
    
  }  
  
  if(modelNames == "VII"){
    
    parameters     = list(mean=init$centers,variance=list(sigmasq=sigmasq,
                                                          sigma=init$cov ), pro=init$weights)
    
    #E-step
    assig = unmap(init$clus[init$clus!=0])   
    #M-step
    mstep          = try(mstepVII(data[init$clus!=0,],z=assig,warn=T))
    parameters     = mstep$parameters
    
    
    #E-step
    discriminants  = d_function(data,parameters,"VII",K)$discriminants  
    trimming       = trim(data,K,D=discriminants,alpha.fixed)
    class          = trimming$classification
    cont           = sum(log(trimming$contribution[trimming$class!=0]))
    loglik[1]      = cont
    
    #Performs trimming on the observations having lowest contribution
    
    mstep          = try(mstepVII(data[class!=0,],z=unmap(class[class!=0]),warn=T))
    parameters     = mstep$parameters
    discriminants =  d_function(data,parameters,"VII",K)$discriminants  
    trimming       = trim(data,K,D=discriminants,alpha.fixed)
    class          = trimming$classification
    cont           = sum(log(trimming$contribution[class!=0]))
    loglik[2]      = cont
    
    while(  (loglik[iter]-loglik[iter-1]) >epsilon){#fix epsilon and run
      
      #Performs EM algorithm
      mstep          = try(mstepVII(data[class!=0,],z=unmap(class[class!=0]),warn=T))
      parameters     = mstep$parameters
      
      #Update the trimming
      D        = d_function(data,parameters=parameters,modelNames = "VII",K)
      trimming = trim(data,K,D=D$discriminants,alpha.fixed)
      class    = trimming$classification
      cont     = sum(log(trimming$contribution[class!=0]))
      
      loglik[iter+1]   = cont
      if(is.na(loglik[iter+1]))(loglik[iter+1]=loglik[1])
      iter=iter+1
      if(iter==iter.max) break      
    }
    
    
    
  }
  
  if(modelNames == "EEI"){
    
    parameters    = list(mean=init$centers,variance=list(scale=sigmasq[which.max(sigmasq)],
                                                         Sigma=init$cov[,,which.max(sigmasq)], 
                                                         shape=shape[,which.max(sigmasq)] ),
                         pro=init$weights)
    #E-step
    assig = unmap(init$clus[init$clus!=0])   
    #M-step
    mstep          = try(mstepEEI(data[init$clus!=0,],z=assig,warn=T))
    parameters     = mstep$parameters
    
    
    #E-step
    discriminants  = d_function(data,parameters,"EEI",K)$discriminants  
    trimming       = trim(data,K,D=discriminants,alpha.fixed)
    class          = trimming$classification
    cont           = sum(log(trimming$contribution[trimming$class!=0]))
    loglik[1]      = cont
    
    #Performs trimming on the observations having lowest contribution
    
    mstep          = try(mstepEEI(data[class!=0,],z=unmap(class[class!=0]),warn=T))
    parameters     = mstep$parameters
    discriminants =  d_function(data,parameters,"EEI",K)$discriminants  
    trimming       = trim(data,K,D=discriminants,alpha.fixed)
    class          = trimming$classification
    cont           = sum(log(trimming$contribution[class!=0]))
    loglik[2]      = cont
    
    while(  (loglik[iter]-loglik[iter-1]) >epsilon){#fix epsilon and run
      
      #Performs EM algorithm
      mstep          = try(mstepEEI(data[class!=0,],z=unmap(class[class!=0]),warn=T))
      parameters     = mstep$parameters
      
      #Update the trimming
      D        = d_function(data,parameters=parameters,modelNames = "EEI",K)
      trimming = trim(data,K,D=D$discriminants,alpha.fixed)
      class    = trimming$classification
      cont     = sum(log(trimming$contribution[class!=0]))
      
      loglik[iter+1]   = cont
      if(is.na(loglik[iter+1]))(loglik[iter+1]=loglik[1])
      iter=iter+1
      if(iter==iter.max) break      
    }
    
  }
  
  if(modelNames == "VEI"){
    
    parameters    = list(mean=init$centers,variance=list(scale=sigmasq,
                                                         shape=shape[,which.max(sigmasq)],
                                                         sigma=init$cov),
                         pro=init$weights)
    #E-step
    assig = unmap(init$clus[init$clus!=0])   
    #M-step
    mstep          = try(mstepVEI(data[init$clus!=0,],z=assig,warn=T))
    parameters     = mstep$parameters
    
    
    #E-step
    discriminants  = d_function(data,parameters,"VEI",K)$discriminants  
    trimming       = trim(data,K,D=discriminants,alpha.fixed)
    class          = trimming$classification
    cont           = sum(log(trimming$contribution[trimming$class!=0]))
    loglik[1]      = cont
    
    #Performs trimming on the observations having lowest contribution
    
    mstep          = try(mstepVEI(data[class!=0,],z=unmap(class[class!=0]),warn=T))
    parameters     = mstep$parameters
    discriminants =  d_function(data,parameters,"VEI",K)$discriminants  
    trimming       = trim(data,K,D=discriminants,alpha.fixed)
    class          = trimming$classification
    cont           = sum(log(trimming$contribution[class!=0]))
    loglik[2]      = cont
    
    while(  (loglik[iter]-loglik[iter-1]) >epsilon){#fix epsilon and run
      
      #Performs EM algorithm
      mstep          = try(mstepVEI(data[class!=0,],z=unmap(class[class!=0]),warn=T))
      parameters     = mstep$parameters
      
      #Update the trimming
      D        = d_function(data,parameters=parameters,modelNames = "VEI",K)
      trimming = trim(data,K,D=D$discriminants,alpha.fixed)
      class    = trimming$classification
      cont     = sum(log(trimming$contribution[class!=0]))
      
      loglik[iter+1]   = cont
      if(is.na(loglik[iter+1]))(loglik[iter+1]=loglik[1])
      iter=iter+1
      if(iter==iter.max) break      
    }
    
  }
  
  if(modelNames == "EVI"){
    
    parameters    = list(mean=init$centers,variance=list(scale=sigmasq[which.max(sigmasq)],
                                                         shape=shape,
                                                         sigma=init$cov),
                         pro=init$weights )
    #E-step
    assig = unmap(init$clus[init$clus!=0])   
    #M-step
    mstep = try(mstepEVI(data[init$clus!=0,],z=assig,warn=T))
    
    #E-step
    discriminants  = d_function(data,parameters,"EVI",K)$discriminants  
    trimming       = trim(data,K,D=discriminants,alpha.fixed)
    class          = trimming$classification
    cont           = sum(log(trimming$contribution[trimming$class!=0]))
    loglik[1]      = cont
    
    #Performs trimming on the observations having lowest contribution
    
    mstep          = try(mstepEVI(data[class!=0,],z=unmap(class[class!=0]),warn=T))
    parameters     = mstep$parameters
    discriminants  = d_function(data,parameters,"EVI",K)$discriminants  
    trimming       = trim(data,K,D=discriminants,alpha.fixed)
    class          = trimming$classification
    cont           = sum(log(trimming$contribution[class!=0]))
    loglik[2]      = cont
    
    while(  (loglik[iter]-loglik[iter-1]) >epsilon){#fix epsilon and run
      
      #Performs EM algorithm
      mstep          = try(mstepEVI(data[class!=0,],z=unmap(class[class!=0]),warn=T))
      parameters     = mstep$parameters
      
      #Update the trimming
      D        = d_function(data,parameters=parameters,modelNames = "EVI",K)
      trimming = trim(data,K,D=D$discriminants,alpha.fixed)
      class    = trimming$classification
      cont     = sum(log(trimming$contribution[class!=0]))
      
      loglik[iter+1]   = cont
      if(is.na(loglik[iter+1]))(loglik[iter+1]=loglik[1])
      iter=iter+1
      if(iter==iter.max) break      
    }
    
    
  }  
  
  if(modelNames == "VVI"){
    
    parameters     = list(mean=init$centers,
                          variance=list(scale=sigmasq,
                                        sigma=init$cov,shape=shape),
                          pro=init$weights )
    
    #E-step
    assig = unmap(init$clus[init$clus!=0])   
    #M-step
    mstep          = try(mstepVVI(data[init$clus!=0,],z=assig,warn=T))
    parameters     = mstep$parameters
    
    #ER
    ni.ini     = c(unlist(table(init$clus[init$clus!=0])))    
    autovalues = matrix(nrow=p,ncol=K) 
    autovectors= array(dim=c(p,p,K))
    for(k in 1:K){
      autovalues[,k]   = eigen(mstep$parameters$variance$sigma[,,k])$values
      autovectors[,,k] = eigen(mstep$parameters$variance$sigma[,,k])$vectors
    }
    M_n = max(autovalues)
    m_n = min(autovalues)
    
    if(M_n/m_n >restr.factor){
      war = 1
      autovalues = .restr2_eigenv(autovalues,ni.ini,restr.factor,zero.tol = 1e-16)
      for(k in 1:K){
        mstep$parameters$variance$sigma[,,k] = autovectors[,,k]%*%(diag(autovalues[,k],nrow=p,ncol=p))%*%t(autovectors[,,k])
        mstep$parameters$variance$scale[k]   = det(mstep$parameters$variance$sigma[,,k])^(1/p)
        mstep$parameters$variance$shape[,k]  = eigen(mstep$parameters$variance$sigma[,,k]/mstep$parameters$variance$scale[k])$values
      }
      warning("solution contrained with the eigenvalue ratio")
    }
    
    
    #E-step
    discriminants  = d_function(data,parameters,"VVI",K)$discriminants  
    trimming       = trim(data,K,D=discriminants,alpha.fixed)
    class          = trimming$classification
    cont           = sum(log(trimming$contribution[trimming$class!=0]))
    loglik[1]      = cont
    
    #Performs trimming on the observations having lowest contribution
    
    mstep          = try(mstepVVI(data[class!=0,],z=unmap(class[class!=0]),warn=T))
    
    #ER
    ni.ini     = c(unlist(table(class[class!=0])))    
    autovalues = matrix(nrow=p,ncol=K) 
    autovectors= array(dim=c(p,p,K))
    for(k in 1:K){
      autovalues[,k]   = eigen(mstep$parameters$variance$sigma[,,k])$values
      autovectors[,,k] = eigen(mstep$parameters$variance$sigma[,,k])$vectors
    }
    M_n = max(autovalues)
    m_n = min(autovalues)
    
    if(M_n/m_n >restr.factor){
      war = 1
      autovalues = .restr2_eigenv(autovalues,ni.ini,restr.factor,zero.tol = 1e-16)
      for(k in 1:K){
        mstep$parameters$variance$sigma[,,k] = autovectors[,,k]%*%(diag(autovalues[,k],nrow=p,ncol=p))%*%t(autovectors[,,k])
        mstep$parameters$variance$scale[k]   = det(mstep$parameters$variance$sigma[,,k])^(1/p)
        mstep$parameters$variance$shape[,k]  = eigen(mstep$parameters$variance$sigma[,,k]/mstep$parameters$variance$scale[k])$values
      }
      warning("solution contrained with the eigenvalue ratio")
    }
    
    parameters     = mstep$parameters
    discriminants  = d_function(data,parameters,"VVI",K)$discriminants  
    trimming       = trim(data,K,D=discriminants,alpha.fixed)
    class          = trimming$classification
    cont           = sum(log(trimming$contribution[class!=0]))
    loglik[2]      = cont
    
    while(  (loglik[iter]-loglik[iter-1]) >epsilon){#fix epsilon and run
      
      #Performs EM algorithm
      mstep          = try(mstepVVI(data[class!=0,],z=unmap(class[class!=0]),warn=T))
      ni.ini     = c(unlist(table(class[class!=0])))    
      autovalues = matrix(nrow=p,ncol=K) 
      autovectors= array(dim=c(p,p,K))
      for(k in 1:K){
        autovalues[,k]   = eigen(mstep$parameters$variance$sigma[,,k])$values
        autovectors[,,k] = eigen(mstep$parameters$variance$sigma[,,k])$vectors
      }
      M_n = max(autovalues)
      m_n = min(autovalues)
      
      if(M_n/m_n >restr.factor){
        war = 1
        autovalues = .restr2_eigenv(autovalues,ni.ini,restr.factor,zero.tol = 1e-16)
        for(k in 1:K){
          mstep$parameters$variance$sigma[,,k] = autovectors[,,k]%*%(diag(autovalues[,k],nrow=p,ncol=p))%*%t(autovectors[,,k])
          mstep$parameters$variance$scale[k]   = det(mstep$parameters$variance$sigma[,,k])^(1/p)
          mstep$parameters$variance$shape[,k]  = eigen(mstep$parameters$variance$sigma[,,k]/mstep$parameters$variance$scale[k])$values
        }
        warning("solution contrained with the eigenvalue ratio")
      }
      
      parameters     = mstep$parameters      
      #Update the trimming
      D        = d_function(data,parameters=parameters,modelNames = "VVI",K)
      trimming = trim(data,K,D=D$discriminants,alpha.fixed)
      class    = trimming$classification
      cont     = sum(log(trimming$contribution[class!=0]))
      
      loglik[iter+1]   = cont
      if(is.na(loglik[iter+1]))(loglik[iter+1]=loglik[1])
      iter=iter+1
      if(iter==iter.max) break      
    }
    
    
  }
  
  if(modelNames == "EEE"){
    
    parameters    = list(mean=init$centers,variance=list(sigmasq=sigmasq,
                                                         Sigma=init$cov[,,which.max(sigmasq)],
                                                         cholSigma=cholsigma[,,which.max(sigmasq)]),
                         pro=init$weights )
    
    #E-step
    assig = unmap(init$clus[init$clus!=0])   
    #M-step
    mstep          = try(mstepEEE(data[init$clus!=0,],z=assig,warn=T))
    parameters     = mstep$parameters
    
    #E-step
    discriminants  = d_function(data,parameters,"EEE",K)$discriminants  
    trimming       = trim(data,K,D=discriminants,alpha.fixed)
    class          = trimming$classification
    cont           = sum(log(trimming$contribution[trimming$class!=0]))
    loglik[1]      = cont
    
    #Performs trimming on the observations having lowest contribution
    
    mstep          = try(mstepEEE(data[class!=0,],z=unmap(class[class!=0]),warn=T))
    parameters     = mstep$parameters
    discriminants  = d_function(data,parameters,"EEE",K)$discriminants  
    trimming       = trim(data,K,D=discriminants,alpha.fixed)
    class          = trimming$classification
    cont           = sum(log(trimming$contribution[class!=0]))
    loglik[2]      = cont
    
    while(  (loglik[iter]-loglik[iter-1]) >epsilon){#fix epsilon and run
      
      #Performs EM algorithm
      mstep          = try(mstepEEE(data[class!=0,],z=unmap(class[class!=0]),warn=T))
      parameters     = mstep$parameters
      
      #Update the trimming
      D        = d_function(data,parameters=parameters,modelNames = "EEE",K)
      trimming = trim(data,K,D=D$discriminants,alpha.fixed)
      class    = trimming$classification
      cont     = sum(log(trimming$contribution[class!=0]))
      
      loglik[iter+1]   = cont
      if(is.na(loglik[iter+1]))(loglik[iter+1]=loglik[1])
      iter=iter+1
      if(iter==iter.max) break      
    }
    
    
  }
  
  if(modelNames == "EVE"){
    
    parameters    = list(mean=init$centers,
                         variance=list(scale=sigmasq,cholsigma=cholsigma,
                                       sigma=init$cov,shape=shape,orientation=orientation),
                         pro=init$weights )
    
    #E-step
    assig = unmap(init$clus[init$clus!=0])   
    #M-step
    mstep          = try(mstepEVE(data[init$clus!=0,],z=assig,warn=T))
    parameters     = mstep$parameters
    
    #ER
    ni.ini     = c(unlist(table(init$clus[init$clus!=0])))    
    autovalues = matrix(nrow=p,ncol=K) 
    autovectors= array(dim=c(p,p,K))
    
    for(k in 1:K){
      autovalues[,k]   = eigen(mstep$parameters$variance$sigma[,,k])$values
      autovectors[,,k] = eigen(mstep$parameters$variance$sigma[,,k])$vectors
    }
    
    M_n = max(autovalues)
    m_n = min(autovalues)
    
    if(M_n/m_n >restr.factor){
      war = 1
      autovalues = .restr2_eigenv(autovalues,ni.ini,restr.factor,zero.tol = 1e-16)
      for(k in 1:K){
        mstep$parameters$variance$sigma[,,k]     = autovectors[,,k]%*%(diag(autovalues[,k],nrow=p,ncol=p))%*%t(autovectors[,,k])
        #???mstep$parameters$variance$cholsigma[,,k] = chol(mstep$parameters$variance$sigma[,,k])
      }
      warning("solution contrained with the eigenvalue ratio")
    }
    
    #E-step
    parameters     = mstep$parameters
    discriminants  = d_function(data,parameters,"EVE",K)$discriminants  
    trimming       = trim(data,K,D=discriminants,alpha.fixed)
    class          = trimming$classification
    cont           = sum(log(trimming$contribution[trimming$class!=0]))
    loglik[1]      = cont
    
    #Performs trimming on the observations having lowest contribution
    
    mstep          = try(mstepEVE(data[class!=0,],z=unmap(class[class!=0]),warn=T))
    #ER
    ni.ini     = c(unlist(table(class[class!=0])))    
    autovalues = matrix(nrow=p,ncol=K) 
    autovectors= array(dim=c(p,p,K))
    
    for(k in 1:K){
      autovalues[,k]   = eigen(mstep$parameters$variance$sigma[,,k])$values
      autovectors[,,k] = eigen(mstep$parameters$variance$sigma[,,k])$vectors
    }
    
    M_n = max(autovalues)
    m_n = min(autovalues)
    
    if(M_n/m_n >restr.factor){
      war = 1
      autovalues = .restr2_eigenv(autovalues,ni.ini,restr.factor,zero.tol = 1e-16)
      for(k in 1:K){
        mstep$parameters$variance$sigma[,,k] = autovectors[,,k]%*%(diag(autovalues[,k],nrow=p,ncol=p))%*%t(autovectors[,,k])
        #   mstep$parameters$variance$cholsigma[,,k] =chol( mstep$parameters$variance$sigma[,,k])
      }
      warning("solution contrained with the eigenvalue ratio")
    }
    
    parameters     = mstep$parameters
    discriminants  = d_function(data,parameters,"EVE",K)$discriminants  
    trimming       = trim(data,K,D=discriminants,alpha.fixed)
    class          = trimming$classification
    cont           = sum(log(trimming$contribution[class!=0]))
    loglik[2]      = cont
    
    while(  (loglik[iter]-loglik[iter-1]) >epsilon){#fix epsilon and run
      
      #Performs EM algorithm
      mstep          = try(mstepEVE(data[class!=0,],z=unmap(class[class!=0]),warn=T))
      ni.ini     = c(unlist(table(class[class!=0])))    
      autovalues = matrix(nrow=p,ncol=K) 
      autovectors= array(dim=c(p,p,K))
      for(k in 1:K){
        autovalues[,k]   = eigen(mstep$parameters$variance$sigma[,,k])$values
        autovectors[,,k] = eigen(mstep$parameters$variance$sigma[,,k])$vectors
      }
      M_n = max(autovalues)
      m_n = min(autovalues)
      
      if(M_n/m_n >restr.factor){
        war = 1
        autovalues = .restr2_eigenv(autovalues,ni.ini,restr.factor,zero.tol = 1e-16)
        for(k in 1:K){
          mstep$parameters$variance$sigma[,,k]     = autovectors[,,k]%*%(diag(autovalues[,k],nrow=p,ncol=p))%*%t(autovectors[,,k])
          #     mstep$parameters$variance$cholsigma[,,k] = chol(mstep$parameters$variance$sigma[,,k])    
        }
        warning("solution contrained with the eigenvalue ratio")
      }
      
      parameters     = mstep$parameters      
      #Update the trimming
      D        = d_function(data,parameters=parameters,modelNames = "EVE",K)
      trimming = trim(data,K,D=D$discriminants,alpha.fixed)
      class    = trimming$classification
      cont     = sum(log(trimming$contribution[class!=0]))
      
      loglik[iter+1]   = cont
      if(is.na(loglik[iter+1]))(loglik[iter+1]=loglik[1])
      iter=iter+1
      if(iter==iter.max) break      
      
    }
    
  }
  
  if(modelNames == "VEE"){
    
    parameters    = list(mean=init$centers,
                         variance=list(scale=sigmasq[which.max(sigmasq)],
                                       sigma=init$cov,shape=shape,orientation=orientation),
                         pro=init$weights )
    #E-step
    assig = unmap(init$clus[init$clus!=0])   
    #M-step
    mstep          = try(mstepVEE(data[init$clus!=0,],z=assig,warn=T))
    parameters     = mstep$parameters
    
    #E-step
    discriminants  = d_function(data,parameters,"VEE",K)$discriminants  
    trimming       = trim(data,K,D=discriminants,alpha.fixed)
    class          = trimming$classification
    cont           = sum(log(trimming$contribution[trimming$class!=0]))
    loglik[1]      = cont
    
    #Performs trimming on the observations having lowest contribution
    
    mstep          = try(mstepVEE(data[class!=0,],z=unmap(class[class!=0]),warn=T))
    parameters     = mstep$parameters
    discriminants  = d_function(data,parameters,"VEE",K)$discriminants  
    trimming       = trim(data,K,D=discriminants,alpha.fixed)
    class          = trimming$classification
    cont           = sum(log(trimming$contribution[class!=0]))
    loglik[2]      = cont
    
    while(  (loglik[iter]-loglik[iter-1]) >epsilon){#fix epsilon and run
      
      #Performs EM algorithm
      mstep          = try(mstepVEE(data[class!=0,],z=unmap(class[class!=0]),warn=T))
      parameters     = mstep$parameters
      
      #Update the trimming
      D        = d_function(data,parameters=parameters,modelNames = "VEE",K)
      trimming = trim(data,K,D=D$discriminants,alpha.fixed)
      class    = trimming$classification
      cont     = sum(log(trimming$contribution[class!=0]))
      
      loglik[iter+1]   = cont
      if(is.na(loglik[iter+1]))(loglik[iter+1]=loglik[1])
      iter=iter+1
      if(iter==iter.max) break      
    }# debug
    
    
  }
  
  if(modelNames == "VVE"){
    
    parameters    = list(mean=init$centers,
                         variance=list(scale=sigmasq,cholsigma=cholsigma,
                                       sigma=init$cov,shape=shape,orientation=orientation),
                         pro=init$weights )
    
    #E-step
    assig = unmap(init$clus[init$clus!=0])   
    #M-step
    mstep          = try(mstepVVE(data[init$clus!=0,],z=assig,warn=T))
    parameters     = mstep$parameters
    
    #ER
    ni.ini     = c(unlist(table(init$clus[init$clus!=0])))    
    autovalues = matrix(nrow=p,ncol=K) 
    autovectors= array(dim=c(p,p,K))
    for(k in 1:K){
      autovalues[,k]   = eigen(mstep$parameters$variance$sigma[,,k])$values
      autovectors[,,k] = eigen(mstep$parameters$variance$sigma[,,k])$vectors
    }
    M_n = max(autovalues)
    m_n = min(autovalues)
    
    if(M_n/m_n >restr.factor){
      war = 1
      autovalues = .restr2_eigenv(autovalues,ni.ini,restr.factor,zero.tol = 1e-16)
      for(k in 1:K){
        mstep$parameters$variance$sigma[,,k]     = autovectors[,,k]%*%(diag(autovalues[,k],nrow=p,ncol=p))%*%t(autovectors[,,k])
        #mstep$parameters$variance$cholsigma[,,k] = chol(mstep$parameters$variance$sigma[,,k])
      }
      warning("solution contrained with the eigenvalue ratio")
    }
    
    #E-step
    parameters     = mstep$parameters
    discriminants  = d_function(data,parameters,"VVE",K)$discriminants  
    trimming       = trim(data,K,D=discriminants,alpha.fixed)
    class          = trimming$classification
    cont           = sum(log(trimming$contribution[trimming$class!=0]))
    loglik[1]      = cont
    
    #Performs trimming on the observations having lowest contribution
    
    mstep          = try(mstepVVE(data[class!=0,],z=unmap(class[class!=0]),warn=T))
    #ER
    ni.ini     = c(unlist(table(class[class!=0])))    
    autovalues = matrix(nrow=p,ncol=K) 
    autovectors= array(dim=c(p,p,K))
    
    for(k in 1:K){
      autovalues[,k]   = eigen(mstep$parameters$variance$sigma[,,k])$values
      autovectors[,,k] = eigen(mstep$parameters$variance$sigma[,,k])$vectors
    }
    
    M_n = max(autovalues)
    m_n = min(autovalues)
    
    if(M_n/m_n >restr.factor){
      war = 1
      autovalues = .restr2_eigenv(autovalues,ni.ini,restr.factor,zero.tol = 1e-16)
      for(k in 1:K){
        mstep$parameters$variance$sigma[,,k] = autovectors[,,k]%*%(diag(autovalues[,k],nrow=p,ncol=p))%*%t(autovectors[,,k])
        #mstep$parameters$variance$cholsigma[,,k] =chol( mstep$parameters$variance$sigma[,,k])
      }
      warning("solution contrained with the eigenvalue ratio")
    }
    
    parameters     = mstep$parameters
    discriminants  = d_function(data,parameters,"VVE",K)$discriminants  
    trimming       = trim(data,K,D=discriminants,alpha.fixed)
    class          = trimming$classification
    cont           = sum(log(trimming$contribution[class!=0]))
    loglik[2]      = cont
    
    while(  (loglik[iter]-loglik[iter-1]) >epsilon){#fix epsilon and run
      
      #Performs EM algorithm
      mstep          = try(mstepVVE(data[class!=0,],z=unmap(class[class!=0]),warn=T))
      ni.ini     = c(unlist(table(class[class!=0])))    
      autovalues = matrix(nrow=p,ncol=K) 
      autovectors= array(dim=c(p,p,K))
      for(k in 1:K){
        autovalues[,k]   = eigen(mstep$parameters$variance$sigma[,,k])$values
        autovectors[,,k] = eigen(mstep$parameters$variance$sigma[,,k])$vectors
      }
      M_n = max(autovalues)
      m_n = min(autovalues)
      
      if(M_n/m_n >restr.factor){
        war = 1
        autovalues = .restr2_eigenv(autovalues,ni.ini,restr.factor,zero.tol = 1e-16)
        for(k in 1:K){
          mstep$parameters$variance$sigma[,,k]     = autovectors[,,k]%*%(diag(autovalues[,k],nrow=p,ncol=p))%*%t(autovectors[,,k])
          #mstep$parameters$variance$cholsigma[,,k] = chol(mstep$parameters$variance$sigma[,,k])    
        }
        warning("solution contrained with the eigenvalue ratio")
      }
      
      parameters     = mstep$parameters      
      #Update the trimming
      D        = d_function(data,parameters=parameters,modelNames = "VVE",K)
      trimming = trim(data,K,D=D$discriminants,alpha.fixed)
      class    = trimming$classification
      cont     = sum(log(trimming$contribution[class!=0]))
      
      loglik[iter+1]   = cont
      if(is.na(loglik[iter+1]))(loglik[iter+1]=loglik[1])
      iter=iter+1
      if(iter==iter.max) break      
      
    }
    
  }
  
  if(modelNames == "EEV"){
    
    parameters    = list(mean=init$centers,
                         variance=list(scale=sigmasq[which.max(sigmasq)],
                                       sigma=init$cov,shape=shape,orientation=orientation),
                         pro=init$weights )
    #E-step
    assig = unmap(init$clus[init$clus!=0])   
    #M-step
    mstep          = try(mstepEEV(data[init$clus!=0,],z=assig,warn=T))
    parameters     = mstep$parameters
    
    #E-step
    discriminants  = d_function(data,parameters,"EEV",K)$discriminants  
    trimming       = trim(data,K,D=discriminants,alpha.fixed)
    class          = trimming$classification
    cont           = sum(log(trimming$contribution[trimming$class!=0]))
    loglik[1]      = cont
    
    #Performs trimming on the observations having lowest contribution
    
    mstep          = try(mstepEEV(data[class!=0,],z=unmap(class[class!=0]),warn=T))
    parameters     = mstep$parameters
    discriminants  = d_function(data,parameters,"EEV",K)$discriminants  
    trimming       = trim(data,K,D=discriminants,alpha.fixed)
    class          = trimming$classification
    cont           = sum(log(trimming$contribution[class!=0]))
    loglik[2]      = cont
    
    while(  (loglik[iter]-loglik[iter-1]) >epsilon){#fix epsilon and run
      
      #Performs EM algorithm
      mstep          = try(mstepEEV(data[class!=0,],z=unmap(class[class!=0]),warn=T))
      parameters     = mstep$parameters
      
      #Update the trimming
      D        = d_function(data,parameters=parameters,modelNames = "EEV",K)
      trimming = trim(data,K,D=D$discriminants,alpha.fixed)
      class    = trimming$classification
      cont     = sum(log(trimming$contribution[class!=0]))
      
      loglik[iter+1]   = cont
      if(is.na(loglik[iter+1]))(loglik[iter+1]=loglik[1])
      iter=iter+1
      if(iter==iter.max) break      
    }
    
    
    
  }
  
  if(modelNames == "VEV"){
    
    parameters    = list(mean=init$centers,
                         variance=list(scale=sigmasq,sigma=init$cov,shape=shape,orientation=orientation),
                         pro=init$weights )
    
    #E-step
    assig = unmap(init$clus[init$clus!=0])   
    #M-step
    mstep          = try(mstepVEV(data[init$clus!=0,],z=assig,warn=T))
    parameters     = mstep$parameters
    
    #E-step
    discriminants  = d_function(data,parameters,"VEV",K)$discriminants  
    trimming       = trim(data,K,D=discriminants,alpha.fixed)
    class          = trimming$classification
    cont           = sum(log(trimming$contribution[trimming$class!=0]))
    loglik[1]      = cont
    
    #Performs trimming on the observations having lowest contribution
    
    mstep          = try(mstepVEV(data[class!=0,],z=unmap(class[class!=0]),warn=T))
    parameters     = mstep$parameters
    discriminants  = d_function(data,parameters,"VEV",K)$discriminants  
    trimming       = trim(data,K,D=discriminants,alpha.fixed)
    class          = trimming$classification
    cont           = sum(log(trimming$contribution[class!=0]))
    loglik[2]      = cont
    
    while(  (loglik[iter]-loglik[iter-1]) >epsilon){#fix epsilon and run
      
      #Performs EM algorithm
      mstep          = try(mstepVEV(data[class!=0,],z=unmap(class[class!=0]),warn=T))
      parameters     = mstep$parameters
      
      #Update the trimming
      D        = d_function(data,parameters=parameters,modelNames = "VEV",K)
      trimming = trim(data,K,D=D$discriminants,alpha.fixed)
      class    = trimming$classification
      cont     = sum(log(trimming$contribution[class!=0]))
      
      loglik[iter+1]   = cont
      if(is.na(loglik[iter+1]))(loglik[iter+1]=loglik[1])
      iter=iter+1
      if(iter==iter.max) break      
    }
    
  }
  
  if(modelNames == "EVV"){
    
    parameters    = list(mean=init$centers,
                         variance=list(scale=sigmasq,cholsigma=cholsigma,
                                       sigma=init$cov,shape=shape,orientation=orientation),
                         pro=init$weights )
    
    #E-step
    assig = unmap(init$clus[init$clus!=0])   
    #M-step
    mstep          = try(mstepEVV(data[init$clus!=0,],z=assig,warn=T))
    parameters     = mstep$parameters
    
    #ER
    ni.ini     = c(unlist(table(init$clus[init$clus!=0])))    
    autovalues = matrix(nrow=p,ncol=K) 
    autovectors= array(dim=c(p,p,K))
    for(k in 1:K){
      autovalues[,k]   = eigen(mstep$parameters$variance$sigma[,,k])$values
      autovectors[,,k] = eigen(mstep$parameters$variance$sigma[,,k])$vectors
    }
    M_n = max(autovalues)
    m_n = min(autovalues)
    
    if(M_n/m_n >restr.factor){
      war = 1
      autovalues = .restr2_eigenv(autovalues,ni.ini,restr.factor,zero.tol = 1e-16)
      for(k in 1:K){
        mstep$parameters$variance$sigma[,,k]     = autovectors[,,k]%*%(diag(autovalues[,k],nrow=p,ncol=p))%*%t(autovectors[,,k])
        #mstep$parameters$variance$cholsigma[,,k] = chol(mstep$parameters$variance$sigma[,,k])
      }
      warning("solution contrained with the eigenvalue ratio")
    }
    
    #E-step
    parameters     = mstep$parameters
    discriminants  = d_function(data,parameters,"EVV",K)$discriminants  
    trimming       = trim(data,K,D=discriminants,alpha.fixed)
    class          = trimming$classification
    cont           = sum(log(trimming$contribution[trimming$class!=0]))
    loglik[1]      = cont
    
    #Performs trimming on the observations having lowest contribution
    
    mstep          = try(mstepEVV(data[class!=0,],z=unmap(class[class!=0]),warn=T))
    #ER
    ni.ini     = c(unlist(table(class[class!=0])))    
    autovalues = matrix(nrow=p,ncol=K) 
    autovectors= array(dim=c(p,p,K))
    
    for(k in 1:K){
      autovalues[,k]   = eigen(mstep$parameters$variance$sigma[,,k])$values
      autovectors[,,k] = eigen(mstep$parameters$variance$sigma[,,k])$vectors
    }
    
    M_n = max(autovalues)
    m_n = min(autovalues)
    
    if(M_n/m_n >restr.factor){
      war = 1
      autovalues = .restr2_eigenv(autovalues,ni.ini,restr.factor,zero.tol = 1e-16)
      for(k in 1:K){
        mstep$parameters$variance$sigma[,,k] = autovectors[,,k]%*%(diag(autovalues[,k],nrow=p,ncol=p))%*%t(autovectors[,,k])
        #mstep$parameters$variance$cholsigma[,,k] =chol( mstep$parameters$variance$sigma[,,k])
      }
      warning("solution contrained with the eigenvalue ratio")
    }
    
    parameters     = mstep$parameters
    discriminants  = d_function(data,parameters,"EVV",K)$discriminants  
    trimming       = trim(data,K,D=discriminants,alpha.fixed)
    class          = trimming$classification
    cont           = sum(log(trimming$contribution[class!=0]))
    loglik[2]      = cont
    
    while(  (loglik[iter]-loglik[iter-1]) >epsilon){#fix epsilon and run
      
      #Performs EM algorithm
      mstep          = try(mstepEVV(data[class!=0,],z=unmap(class[class!=0]),warn=T))
      ni.ini     = c(unlist(table(class[class!=0])))    
      autovalues = matrix(nrow=p,ncol=K) 
      autovectors= array(dim=c(p,p,K))
      for(k in 1:K){
        autovalues[,k]   = eigen(mstep$parameters$variance$sigma[,,k])$values
        autovectors[,,k] = eigen(mstep$parameters$variance$sigma[,,k])$vectors
      }
      M_n = max(autovalues)
      m_n = min(autovalues)
      
      if(M_n/m_n >restr.factor){
        war = 1
        autovalues = .restr2_eigenv(autovalues,ni.ini,restr.factor,zero.tol = 1e-16)
        for(k in 1:K){
          mstep$parameters$variance$sigma[,,k]     = autovectors[,,k]%*%(diag(autovalues[,k],nrow=p,ncol=p))%*%t(autovectors[,,k])
          #mstep$parameters$variance$cholsigma[,,k] = chol(mstep$parameters$variance$sigma[,,k])    
        }
        warning("solution contrained with the eigenvalue ratio")
      }
      
      parameters     = mstep$parameters      
      #Update the trimming
      D        = d_function(data,parameters=parameters,modelNames = "EVV",K)
      trimming = trim(data,K,D=D$discriminants,alpha.fixed)
      class    = trimming$classification
      cont     = sum(log(trimming$contribution[class!=0]))
      
      loglik[iter+1]   = cont
      if(is.na(loglik[iter+1]))(loglik[iter+1]=loglik[1])
      iter=iter+1
      if(iter==iter.max) break      
      
    }
    
  }
  
  if(modelNames == "VVV"){
    
    parameters    = list(mean=init$centers,
                         variance=list(scale=sigmasq,cholsigma=cholsigma,
                                       sigma=init$cov,shape=shape,orientation=orientation),
                         pro=init$weights )
    
    #E-step
    assig = unmap(init$clus[init$clus!=0])   
    #M-step
    mstep          = try(mstepVVV(data[init$clus!=0,],z=assig,warn=T))
    parameters     = mstep$parameters
    
    #ER
    ni.ini     = c(unlist(table(init$clus[init$clus!=0])))    
    autovalues = matrix(nrow=p,ncol=K) 
    autovectors= array(dim=c(p,p,K))
    for(k in 1:K){
      autovalues[,k]   = eigen(mstep$parameters$variance$sigma[,,k])$values
      autovectors[,,k] = eigen(mstep$parameters$variance$sigma[,,k])$vectors
    }
    M_n = max(autovalues)
    m_n = min(autovalues)
    
    if(M_n/m_n >restr.factor){
      war = 1
      autovalues = .restr2_eigenv(autovalues,ni.ini,restr.factor,zero.tol = 1e-16)
      for(k in 1:K){
        mstep$parameters$variance$sigma[,,k]     = autovectors[,,k]%*%(diag(autovalues[,k],nrow=p,ncol=p))%*%t(autovectors[,,k])
        mstep$parameters$variance$cholsigma[,,k] = chol(mstep$parameters$variance$sigma[,,k])
      }
      warning("solution contrained with the eigenvalue ratio")
    }
    
    #E-step
    parameters     = mstep$parameters
    discriminants  = d_function(data,parameters,"VVV",K)$discriminants  
    trimming       = trim(data,K,D=discriminants,alpha.fixed)
    class          = trimming$classification
    cont           = sum(log(trimming$contribution[trimming$class!=0]))
    loglik[1]      = cont
    
    #Performs trimming on the observations having lowest contribution
    
    mstep          = try(mstepVVV(data[class!=0,],z=unmap(class[class!=0]),warn=T))
    #ER
    ni.ini     = c(unlist(table(class[class!=0])))    
    autovalues = matrix(nrow=p,ncol=K) 
    autovectors= array(dim=c(p,p,K))
    
    for(k in 1:K){
      autovalues[,k]   = eigen(mstep$parameters$variance$sigma[,,k])$values
      autovectors[,,k] = eigen(mstep$parameters$variance$sigma[,,k])$vectors
    }
    
    M_n = max(autovalues)
    m_n = min(autovalues)
    
    if(M_n/m_n >restr.factor){
      war = 1
      autovalues = .restr2_eigenv(autovalues,ni.ini,restr.factor,zero.tol = 1e-16)
      for(k in 1:K){
        mstep$parameters$variance$sigma[,,k] = autovectors[,,k]%*%(diag(autovalues[,k],nrow=p,ncol=p))%*%t(autovectors[,,k])
        mstep$parameters$variance$cholsigma[,,k] =chol( mstep$parameters$variance$sigma[,,k])
      }
      warning("solution contrained with the eigenvalue ratio")
    }
    
    parameters     = mstep$parameters
    discriminants  = d_function(data,parameters,"VVV",K)$discriminants  
    trimming       = trim(data,K,D=discriminants,alpha.fixed)
    class          = trimming$classification
    cont           = sum(log(trimming$contribution[class!=0]))
    loglik[2]      = cont
    
    while(  (loglik[iter]-loglik[iter-1]) >epsilon){#fix epsilon and run
      
      #Performs EM algorithm
      mstep          = try(mstepVVV(data[class!=0,],z=unmap(class[class!=0]),warn=T))
      ni.ini     = c(unlist(table(class[class!=0])))    
      autovalues = matrix(nrow=p,ncol=K) 
      autovectors= array(dim=c(p,p,K))
      for(k in 1:K){
        autovalues[,k]   = eigen(mstep$parameters$variance$sigma[,,k])$values
        autovectors[,,k] = eigen(mstep$parameters$variance$sigma[,,k])$vectors
      }
      M_n = max(autovalues)
      m_n = min(autovalues)
      
      if(M_n/m_n >restr.factor){
        war = 1
        autovalues = .restr2_eigenv(autovalues,ni.ini,restr.factor,zero.tol = 1e-16)
        for(k in 1:K){
          mstep$parameters$variance$sigma[,,k]     = autovectors[,,k]%*%(diag(autovalues[,k],nrow=p,ncol=p))%*%t(autovectors[,,k])
          mstep$parameters$variance$cholsigma[,,k] = chol(mstep$parameters$variance$sigma[,,k])    
        }
        warning("solution contrained with the eigenvalue ratio")
      }
      
      parameters     = mstep$parameters      
      #Update the trimming
      D        = d_function(data,parameters=parameters,modelNames = "VVV",K)
      trimming = trim(data,K,D=D$discriminants,alpha.fixed)
      class    = trimming$classification
      cont     = sum(log(trimming$contribution[class!=0]))
      
      loglik[iter+1]   = cont
      if(is.na(loglik[iter+1]))(loglik[iter+1]=loglik[1])
      iter=iter+1
      if(iter==iter.max) break      
      
    }
  }
  
  return(list(data=data,classification=trimming$classification,init=init,
              loglikelihood=loglik,contribution=trimming$contribution,
              parameters=mstep$parameters,war=war))
  
}


# These are the internal functions involved by the main function. The user
# is not supposed to directely deal with them.

is.NullOb <- function(x) is.null(x) | all(sapply(x, is.null))

rmNullObs <- function(x) {
  x <- Filter(Negate(is.NullOb), x)
  lapply(x, function(x) if (is.list(x)) rmNullObs(x) else x)
}

d_function = function(data,parameters,modelNames,K){
  
  n  = nrow(data)
  p  = ncol(data)
  D  = matrix(nrow=n,ncol=K)
  
  
  if(modelNames=="EII" | modelNames == "EEI" | modelNames=="EEE"){
    
    for(j in 1:K){
      D[,j] = c(dmvnorm(data,mean=parameters$mean[,j],
                        sigma=parameters$variance$Sigma)*parameters$pro[j])
    }
  }else{
    for(j in 1:K){
      D[,j] = c(dmvnorm(data,mean=parameters$mean[,j],
                        sigma=parameters$variance$sigma[,,j])*parameters$pro[j])
    }
  }
  
  return(list(discriminants=D))
  
}


trim = function(data,K,D,alpha.fixed){
  
  n=nrow(data)
  p=ncol(data)
  
  D=matrix(unlist(D),nrow=nrow(data),ncol=K)
  classification = apply(D,1,which.max)
  cont = apply(D,1,max)
  cont_ord = c(rank(cont))
  classification[which(cont_ord<=ceiling(n*(alpha.fixed)))]=0
  
  return(list(data=data,contribution=cont,classification=classification))
}

#function to apply the ER constraint

.restr2_eigenv <- function (autovalues, ni.ini, restr.fact, zero.tol)
{
  ev <- autovalues
  
  ###### function parameters:
  ###### ev: matrix containin eigenvalues  								#n proper naming - changed autovalue to ev (eigenvalue)
  ###### ni.ini: current sample size of the clusters						#n proper naming - changed ni.ini to cSize (cluster size)
  ###### factor: the factor parameter in tclust program 
  ###### init: 1 if we are applying restrictions during the smart inicialization   
  ######       0 if we are applying restrictions during the C steps execution
  
  ###### some inicializations
  
  if (!is.matrix (ev))												#n	checking for malformed ev - argument (e.g. a vector of determinants instead of a matrix)
    if (is.atomic (ev))												#n
      ev <- t (as.numeric (ev))									#n
    else															#n
      ev <- as.matrix (ev)										#n
    
    stopifnot (ncol (ev) == length (ni.ini))							#n	check wether the matrix autovalues and ni.ini have the right dimension.
    
    d <- t (ev)
    
    p <- nrow (ev) 
    K <- ncol (ev)	
    
    n <- sum(ni.ini)
    
    nis <- matrix(data=ni.ini,nrow=K,ncol=p)
    
    #m	MOVED: this block has been moved up a bit. see "old position of "block A" for it's old occurrence
    idx.nis.gr.0 <- nis > zero.tol										#n	as nis € R we have to be carefull when checking against 0
    used.ev <- ni.ini > zero.tol										#n
    ev.nz <- ev[,used.ev]												#n	non-zero eigenvalues
    #m
    #o	if ((max (d[nis > 0]) <= zero.tol))									#m
    #i	if ((max (d[idx.nis.gr.0]) <= zero.tol))							#n	"idx.nis.gr.0" instead of (nis > 0)
    if ((max (ev.nz) <= zero.tol))										#n	simplify syntax
      return (matrix (0, nrow = p, ncol = K))							#m
    #m
    ###### we check if the  eigenvalues verify the restrictions			#m
    #m
    #o	if (max (d[nis > 0]) / min (d[nis > 0]) <= restr.fact)				#m
    #i	if (max (d[idx.nis.gr.0]) / min (d[idx.nis.gr.0]) <= restr.fact)	#n	"idx.nis.gr.0" instead of (nis > 0)
    if (max (ev.nz) / min (ev.nz) <= restr.fact)						#n	simplify syntax
      
    {																	#m
      #o		d[!idx.nis.gr.0] <- mean (d[idx.nis.gr.0])						#n	"idx.nis.gr.0" instead of (nis > 0)
      ev[,!used.ev] <- mean (ev.nz)									#n	simplify syntax
      return (ev)														#m
    }																	#m
    
    ###### d_ is the ordered set of values in which the restriction objective function change the definition
    ###### points in d_ correspond to  the frontiers for the intervals in which this objective function has the same definition
    ###### ed is a set with the middle points of these intervals 
    
    #o	d_ <- sort (c (d, d / restr.fact))
    d_ <- sort (c (ev, ev / restr.fact))								#n	using ev instead of d
    dim <- length (d_)													##2do: continue here cleaning up .restr2_eigenv
    d_1 <- d_
    d_1[dim+1] <- d_[dim] * 2
    d_2 <- c (0, d_)
    ed <- (d_1 + d_2) / 2
    dim <- dim + 1;
    
    ##o
    ##o	old position of "block A"
    ##o
    
    ###### the only relevant eigenvalues are those belong to a clusters with sample size greater than 0.
    ###### eigenvalues corresponding to a clusters with 0 individuals has no influence in the objective function
    ###### if all the eigenvalues are 0 during the smart initialization we assign to all the eigenvalues the value 1  
    
    ###### we build the sol array 
    ###### sol[1],sol[2],.... this array contains the critical values of the interval functions which defines the m objective function  
    ###### we use the centers of the interval to get a definition for the function in each interval
    ###### this set with the critical values (in the array sol) contains the optimum m value  
    
    t <- s <- r <- array(0, c(K, dim))
    sol <- sal <- array(0, c(dim))
    
    for (mp_ in 1:dim) 
    {
      for (i in 1:K)
      {  
        r[i,mp_] <- sum ((d[i,] < ed[mp_])) + sum((d[i,] > ed[mp_]*restr.fact))
        s[i,mp_] <- sum (d[i,]*(d[i,] < ed[mp_]))
        t[i,mp_] <- sum (d[i,]*(d[i,] > ed[mp_] * restr.fact))
      }
      
      sol[mp_] <- sum (ni.ini / n * (s[,mp_] + t[,mp_] / restr.fact)) / (sum(ni.ini /n * (r[, mp_])))
      
      e <-	sol[mp_] * (d < sol[mp_]) +
        d * (d >= sol[mp_]) * (d <= restr.fact * sol[mp_]) + 
        (restr.fact*sol[mp_]) * (d > restr.fact * sol[mp_])
      o <- -1/2 * nis / n * (log(e) + d / e)
      
      sal[mp_] <- sum(o)
    }
    
    ###### m is the optimum value for the eigenvalues procedure 
    #o	eo <- which.max (c (sal))						## remove c ()
    #	m <- sol[eo]
    m <- sol[which.max (sal)]						#n
    ###### based on the m value we get the restricted eigenvalues
    
    t (m * (d < m) + d * (d >= m) * (d <= restr.fact * m) + (restr.fact * m) * (d > restr.fact * m))	#
}



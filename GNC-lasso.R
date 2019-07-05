library(Matrix)
library(irlba)
library(MASS)
library(glasso)




gnclasso <- function(X,A,lambda,alpha,tol=1e-4,verbose=FALSE,W.init=NULL,Theta.init=NULL,M.init=NULL,L.U=NULL,L.tau=NULL,low.rank=NULL){
  A <- as.matrix(A)
  D <- diag(colSums(A))
  L <- D-A
  n <- nrow(A)
  p <- ncol(X)
  if(ncol(A)!=n) stop("Invalid adjacency matrix!")
  if(is.null(W.init)){ 
    S <- t(X)%*%X/n
    diag(S) <- diag(S)+lambda
    W.init <- S
  }
  if(is.null(Theta.init)){
    Theta.init <- solve(W.init)
  }	
  if(is.null(L.U)){
    if(is.null(low.rank)){
      L.eigen <- eigen(L)
      L.U <- L.eigen$vectors
      L.tau <- L.eigen$values + 0.01*min(L.eigen$values)
    }else{
      L.eigen <-	partial_eigen(L, n = low.rank, symmetric = TRUE)
      L.tau <- L.eigen$values + 0.01*min(L.eigen$values)
      L.U <- L.eigen$vectors
    }
  }
  if(is.null(M.init)){
    M.init <- matrix(0,nrow=n,ncol=p)
    tmp <- t(L.U)%*%X
    tmp <- tmp*(1/(1+alpha*L.tau))
    M.init <- L.U%*%tmp
  }
  M <- M.init
  X.centered <- X-M
  S <- t(X.centered)%*%X.centered/n
  result.g <- glasso(s=S,rho=lambda,trace=verbose)
  result <- list()
  edge.df <- (sum(abs(result.g$wi)>1e-8)-p)/2
  if(verbose) print(paste(edge.df,"edges"))
  cohesion.df <- p*sum(1/(1+alpha*L.tau))
  df <- cohesion.df + edge.df
  result$df <- df
  loglike = determinant(result.g$wi, logarithm = TRUE)$modulus[1] - sum(diag(result.g$wi%*%S))
  BIC = loglike - df*log(n)
  AIC = loglike - df*2
  GCV = norm(X.centered,"F")^2/(p*(n-cohesion.df/p)^2)
  result$loglike = loglike
  result$bic = BIC
  result$aic = AIC
  result$gcv = GCV
  result$cohesion.df <- cohesion.df
  result$Theta = as(result.g$wi,"dgCMatrix")
  result$M <- M
  result$num.edge <- edge.df
  result$df <- df
  result$Sigma <- result.g$w
  return(result)
}




gnclassopath <- function(X,A,lambda.seq,alpha,tol=1e-4,verbose=FALSE,M.init=NULL,L.U=NULL,L.tau=NULL,low.rank=NULL){
  A <- as.matrix(A)
  D <- diag(colSums(A))
  L <- D-A
  n <- nrow(A)
  p <- ncol(X)
  if(ncol(A)!=n) stop("Invalid adjacency matrix!")
  ptm <- proc.time()
  if(is.null(L.U)){
    if(is.null(low.rank)){
      L.eigen <- eigen(L)
      L.U <- L.eigen$vectors
      L.tau <- L.eigen$values + 0.01*min(L.eigen$values)
    }else{
      L.eigen <-	partial_eigen(L, n = low.rank, symmetric = TRUE)
      L.tau <- L.eigen$values + 0.01*min(L.eigen$values)
      L.U <- L.eigen$vectors
    }
  }
  if(is.null(M.init)){
    M.init <- matrix(0,nrow=n,ncol=p)
    tmp <- t(L.U)%*%X
    tmp <- tmp*(1/(1+alpha*L.tau))
    M.init <- L.U%*%tmp
  }
  M <- M.init
  m.time <- proc.time() - ptm
  X.centered <- X-M
  S <- t(X.centered)%*%X.centered/n
  lambda.seq <- sort(lambda.seq)
  sigma.time <- system.time(result.g <- glassopath(s=S,rholist=lambda.seq,trace=verbose))
  result.list <- list()
  cohesion.df <- p*sum(1/(1+alpha*L.tau))
  result.list$models <- list()
  for(dd in 1:length(lambda.seq)){
    
    result <- list()
    edge.df <- (sum(abs(result.g$wi[,,dd])>1e-8)-p)/2
    if(verbose) print(paste(edge.df,"edges"))
    df <- cohesion.df + edge.df
    result$df <- df
    loglike = determinant(result.g$wi[,,dd], logarithm = TRUE)$modulus[1] - sum(diag(result.g$wi[,,dd]%*%S))
    BIC = loglike - df*log(n)
    AIC = loglike - df*2
    GCV = n*norm(X.centered,"F")^2/((n-cohesion.df/p)^2)
    result$loglike = loglike
    result$bic = BIC
    result$aic = AIC
    result$gcv = GCV
    result$cohesion.df <- cohesion.df
    result$Theta = as(result.g$wi[,,dd],"dgCMatrix")
    result$M <- M
    result$num.edge <- edge.df
    result$df <- df
    result$W <- result.g$w[,,dd]
    result.list$models[[dd]] <- result
  }
  result.list$gcv <- result$gcv
  result.list$cohesion.df <- result$cohesion.df
  result.list$df.seq <- unlist(lapply(result.list$models,function(x) return(x$df)))
  result.list$bic.seq <- unlist(lapply(result.list$models,function(x) return(x$bic)))
  result.list$aic.seq <- unlist(lapply(result.list$models,function(x) return(x$aic)))
  result.list$num.edge.seq <- unlist(lapply(result.list$models,function(x) return(x$num.edge)))
  result.list$lambda.seq <- sort(lambda.seq,decreasing=TRUE)
  result.list$m.time <- m.time
  result.list$sigma.time <- sigma.time
  return(result.list)
}





gnc.M.semi <- function(X,A,alpha,L.U=NULL,L.tau=NULL,low.rank=NULL,train.index){
  A <- as.matrix(A)
  D <- diag(colSums(A))
  L <- D-A
  n <- nrow(A)
  p <- ncol(X)
  if(ncol(A)!=n) stop("Invalid adjacency matrix!")
  
  if(is.null(L.U)){
    if(is.null(low.rank)){
      L.eigen <- eigen(L)
      L.U <- L.eigen$vectors
      L.tau <- L.eigen$values + 0.01*min(L.eigen$values)
    }else{
      L.eigen <-	partial_eigen(L, n = low.rank, symmetric = TRUE)
      L.tau <- L.eigen$values + 0.01*min(L.eigen$values)
      L.U <- L.eigen$vectors
    }
  }
  X.train <- X[train.index,]
  X.test <- X[-train.index,]
  test.index <- setdiff(1:n,train.index)
  L.U.train <- L.U[train.index,]
  B <- solve(t(L.U.train)%*%L.U.train+alpha*diag(L.tau),t(L.U.train)%*%X.train)
  M <- L.U.train%*%B
  M.test <- L.U[test.index,]%*%B
  M.all <- matrix(0,n,p)
  M.all[train.index,] <- M
  M.all[test.index,] <- M.test
  
  err <- norm(M.test-X.test,"F")^2/(length(test.index)*p)
  return(list(M.train=M,M.test=M.test,err=err,M.all=M.all))
}

gnc.M.cv <- function(X,A,alpha.seq,V,proportion=0.1,L.U=NULL,L.tau=NULL,low.rank=NULL){
  n <- nrow(A)
  p <- ncol(X)
  K <- V
  m <- floor(n*proportion)
  A <- as.matrix(A)
  D <- diag(colSums(A))
  L <- D-A  
  if(ncol(A)!=n) stop("Invalid adjacency matrix!")
  
  if(is.null(L.U)){
    if(is.null(low.rank)){
      L.eigen <- eigen(L)
      L.U <- L.eigen$vectors
      L.tau <- L.eigen$values + 0.00001*min(L.eigen$values)
    }else{
      L.eigen <-	partial_eigen(L, n = low.rank, symmetric = TRUE)
      L.tau <- L.eigen$values + 0.00001*min(L.eigen$values)
      L.U <- L.eigen$vectors
    }
  }
  cv.err.mat <- matrix(0,K,length(alpha.seq))
  for(k in 1:K){
    print(paste("CV iteration: ",k))
    train.index <- sample(n,m)
    for(i in 1:length(alpha.seq)){
      cv.semi <- gnc.M.semi(X,A,alpha.seq[i],L.U=L.U,L.tau=L.tau,low.rank=low.rank,train.index=train.index)
      cv.err.mat[k,i] <- cv.semi$err
    }
  }
  return(list(cv.err.mat=cv.err.mat,opt.index=which.min(colMeans(cv.err.mat))))
}


glasso.cv <- function(X,M,lambda.seq,V,proportion=0.1){
  n <- nrow(X)
  p <- ncol(X)
  K <- V
  m <- floor(n*proportion)
  lambda.seq <- sort(lambda.seq,decreasing=FALSE)
  cv.err.mat <- matrix(0,K,length(lambda.seq))
  X.centered <- X-M
  cv.loglike.mat <- matrix(0,K,length(lambda.seq))
  for(k in 1:K){
    print(paste("CV iteration: ",k))
    train.index <- sample(n,m)
    X.train <- X.centered[train.index,]
    test.index <- setdiff(1:n,train.index)
    X.test <- X.centered[test.index,]
    S.train <- t(X.train)%*%X.train/length(train.index)
    S.test <- t(X.test)%*%X.test/length(test.index)
    t0.path <- glassopath(s=S.train,rholist=lambda.seq,trace=FALSE)
    for(dd in 1:length(lambda.seq)){
      loglike <- determinant(t0.path$wi[,,dd], logarithm = TRUE)$modulus[1] - sum(diag(t0.path$wi[,,dd]%*%S.test)) 
      cv.loglike.mat[k,dd] <- loglike
    }
  }
  return(list(cv.loglike.mat=cv.loglike.mat,opt.index=which.max(colMeans(cv.loglike.mat))))
}



gnc.M <- function(X,A,alpha,L.U=NULL,L.tau=NULL,low.rank=NULL){
  A <- as.matrix(A)
  D <- diag(colSums(A))
  L <- D-A
  n <- nrow(A)
  p <- ncol(X)
  if(ncol(A)!=n) stop("Invalid adjacency matrix!")
  
  if(is.null(L.U)){
    if(is.null(low.rank)){
      L.eigen <- eigen(L)
      L.U <- L.eigen$vectors
      L.tau <- L.eigen$values + 0.00001*min(L.eigen$values)
    }else{
      L.eigen <-	partial_eigen(L, n = low.rank, symmetric = TRUE)
      L.tau <- L.eigen$values + 0.00001*min(L.eigen$values)
      L.U <- L.eigen$vectors
    }
  }
  
  M.init <- matrix(0,nrow=n,ncol=p)
  tmp <- t(L.U)%*%X
  tmp <- tmp*(1/(1+alpha*L.tau))
  M.init <- L.U%*%tmp
  
  M <- M.init
  return(M)
} 




gnc.M.gcv <- function(X,A,alpha,L.U=NULL,L.tau=NULL,low.rank=NULL){
  A <- as.matrix(A)
  D <- diag(colSums(A))
  L <- D-A
  n <- nrow(A)
  p <- ncol(X)
  if(ncol(A)!=n) stop("Invalid adjacency matrix!")
  
  if(is.null(L.U)){
    if(is.null(low.rank)){
      L.eigen <- eigen(L)
      L.U <- L.eigen$vectors
      L.tau <- L.eigen$values + 0.00001*min(L.eigen$values)
    }else{
      L.eigen <-	partial_eigen(L, n = low.rank, symmetric = TRUE)
      L.tau <- L.eigen$values + 0.00001*min(L.eigen$values)
      L.U <- L.eigen$vectors
    }
  }

  AA <- t(L.U)*(1/(1+alpha*L.tau))
  AA <- L.U%*%AA

  M <- AA%*%X
  Err <- norm(X-M,"F")^2/(n*p)
  denom <- (1-sum(diag(AA))/n)^2
  gcv <- Err/denom
  
  
  return(list(M=M,gcv=gcv))
} 





gnc.M.gcv.seq <- function(X,A,alpha.seq,L.U=NULL,L.tau=NULL,low.rank=NULL){
  A <- as.matrix(A)
  D <- diag(colSums(A))
  L <- D-A
  n <- nrow(A)
  p <- ncol(X)
  if(ncol(A)!=n) stop("Invalid adjacency matrix!")
  
  if(is.null(L.U)){
    if(is.null(low.rank)){
      L.eigen <- eigen(L)
      L.U <- L.eigen$vectors
      L.tau <- L.eigen$values + 0.00001*min(L.eigen$values)
    }else{
      L.eigen <-	partial_eigen(L, n = low.rank, symmetric = TRUE)
      L.tau <- L.eigen$values + 0.00001*min(L.eigen$values)
      L.U <- L.eigen$vectors
    }
  }
  gcv.seq <- rep(0,length(alpha.seq))
  for(i in 1:length(alpha.seq)){
      alpha <- alpha.seq[i]
  AA <- t(L.U)*(1/(1+alpha*L.tau))
  AA <- L.U%*%AA

  M <- AA%*%X
  Err <- norm(X-M,"F")^2/(n*p)
  denom <- (1-sum(diag(AA))/n)^2
  gcv <- Err/denom
      gcv.seq[i] <- gcv
  }
  
  return(gcv.seq)
} 






gnc.M.sup <- function(X,A,alpha,L.U=NULL,L.tau=NULL,low.rank=NULL,train.index){
  n <- nrow(A)
  test.index <- setdiff(1:n,train.index)
  n.train <- length(train.index)
  X <- X[c(train.index,test.index),]
  X.train <- X[1:n.train,]
  X.test <- X[(n.train+1):n,]
  A <- A[c(train.index,test.index),c(train.index,test.index)]
  A.train <- as.matrix(A[1:n.train,1:n.train])
  A.full <- as.matrix(A)
  D.train <- diag(colSums(A.train))
  D.full <- diag(colSums(A.full))
  L.train <- D.train-A.train
  L.full <- D.full-A.full
  p <- ncol(X)
  if(ncol(A)!=n) stop("Invalid adjacency matrix!")
  
  if(is.null(L.U)){
    if(is.null(low.rank)){
      L.eigen <- eigen(L.train)
      L.U <- L.eigen$vectors
      L.tau <- L.eigen$values + 0.01*min(L.eigen$values)
    }else{
      L.eigen <-	partial_eigen(L, n = low.rank, symmetric = TRUE)
      L.tau <- L.eigen$values + 0.01*min(L.eigen$values)
      L.U <- L.eigen$vectors
    }
  }
  
  M.train <- matrix(0,nrow=n.train,ncol=p)
  tmp <- t(L.U)%*%X.train
  tmp <- tmp*(1/(1+alpha*L.tau))
  M.train <- L.U%*%tmp
  
  M.test <- matrix(0,nrow=length(test.index),ncol=p)
  L11 <- L.full[(n.train+1):n,(n.train+1):n]
  L12 <- L.full[(n.train+1):n,1:n.train]
  M.test <- -solve(L11,L12%*%M.train)
  M.all <- matrix(0,n,p)
  M.all[train.index,] <- M.train
  M.all[test.index,] <- M.test
  
  err <- norm(M.test-X.test,"F")^2/(length(test.index)*p)
  return(list(M.train=M,M.test=M.test,err=err,M.all=M.all))
} 




gnc.M.cv.sup <- function(X,A,alpha.seq,V,proportion=0.1,L.U=NULL,L.tau=NULL,low.rank=NULL){
  n <- nrow(A)
  p <- ncol(X)
  K <- V
  m <- floor(n*proportion)
  A <- as.matrix(A)
  D <- diag(colSums(A))
  L <- D-A  
  if(ncol(A)!=n) stop("Invalid adjacency matrix!")
  
  cv.err.mat <- matrix(0,K,length(alpha.seq))
  for(k in 1:K){
    print(paste("CV iteration: ",k))
    train.index <- sample(n,m)
    for(i in 1:length(alpha.seq)){
      cv.semi <- gnc.M.sup(X,A,alpha.seq[i],train.index=train.index)
      cv.err.mat[k,i] <- cv.semi$err
    }
  }
  return(list(cv.err.mat=cv.err.mat,opt.index=which.min(colMeans(cv.err.mat))))
}

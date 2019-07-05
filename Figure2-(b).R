library(netcoh)
### running the proposed two-stage estimation does not need the package. But iterative version needs the packaage.
### you should use the version of netcoh included in the folder (version 0.35)
library(pROC)
library(glasso)
library(huge)
library(doParallel)
registerDoParallel(cores=25)

source("GNC-lasso.R")
source("GridGenFunc.R")

m <- 20
set.seed(100)
m <- 20
GridData <- GridNetDataGen(m=m)


n <- m^2
p <- 500

library(lattice)


set.seed(999)
Adj <- GridData$A
library(igraph)
g <- graph.adjacency(Adj,mode="undirected")

D <- diag(colSums(Adj))

L <- D-Adj

library(irlba)

eig <- eigen(L,symmetric=TRUE)

m1 <- eig$vectors[,400]
m2 <- eig$vectors[,399]



mm <- sqrt(0.95)*m1 + sqrt(0.05)*m2



M <- matrix(0,n,p)
for(j in 1:p){
  #M[,j] <- GridData$value/max(GridData$value)#*sample(c(-1,1),size=1)
    M[,j] <- mm*sqrt(n)*0.5#mm/(max(mm)-min(mm))
}



#### unequal variance setting



library(RSpectra)
library(irlba)


gnc0.err <- gnc1.err <- gnc3.err <- rep(0,100)
gnc0.AUC <- glasso.AUC <- gnc.AUC <- rep(0,100)
gnc0.fp <- gnc0.tp <- glasso.fp <- glasso.tp <- gnc1.fp <- gnc1.tp <- gnc3.fp <- gnc3.tp <- NULL
set.seed(100)

n <- 400
set.seed(999)
sim.data <- huge.generator(n=n*100,d=p,prob=0.01,v=0.05)



result <- foreach(ii = 1:100)%dopar%{
  #write(ii,file="log.txt")
  X <- 0.5*sim.data$data[(((ii-1)*n+1):(ii*n)),]
  
  new.X <- X+M
  new.X <- scale(new.X,center=FALSE,scale=FALSE)
  mean.X <- colMeans(new.X)
  sd.X <- apply(new.X,2,sd)
  scaled.M <- t((t(M)-mean.X)/sd.X)
  X.std <- scale(new.X,center=TRUE,scale=TRUE)
  new.S <- t(X.std)%*%X.std/n
  lambda.seq=exp(seq(log(0.3),log(0.105),length.out=50))
  #system.time(t0.path <- rncglassopath(X=X.std,A=Adj,lambdaseq=lambda.seq,alpha=0.05*n/log(n),safe=TRUE,unequal.scale = FALSE,verbose=TRUE,max.iter=1))
  #t0.path <- gnclassopath(X=X.std,A=Adj,lambda.seq=lambda.seq,alpha=0.1*n/(log(n)))

  d <- mean(colSums(Adj))
  rho.seq <- exp(seq(log(0.001),log(0.8),length.out=20))
  alpha.seq <- rho.seq*n/d
  system.time(test.M.cv <- gnc.M.cv(X.std,Adj,alpha.seq,V=10))
  alpha.cv <- alpha.seq[test.M.cv$opt.index]
  #gcv.seq <- gnc.M.gcv.seq(X.std,Adj,alpha.seq)
  #alpha.cv <- alpha.seq[which.min(gcv.seq)]
  
  lambda.seq3=exp(seq(log(0.2),log(0.065),length.out=50))
  system.time(t3.path <- gnclassopath(X=X.std,A=Adj,lambda.seq=lambda.seq3,alpha=alpha.cv))
  t3.path$num.edge

 lambda.seq2=exp(seq(log(0.195),log(0.065),length.out=50))
  system.time(t2.path <- gnclassopath(X=X.std,A=Adj,lambda.seq=lambda.seq2,alpha=0.06*n/d))
  t2.path$num.edge
 
   lambda.seq1=exp(seq(log(0.195),log(0.06),length.out=50))
  system.time(t1.path <- rncglassopath(X=X.std,A=Adj,lambdaseq=lambda.seq1,alpha=0.078*n/log(n),safe=FALSE,unequal.scale = FALSE,verbose=TRUE,max.iter=50))
  t1.path$num.edge

  t0.err <- t1.err <- t2.err <- t3.err <- list()
  t0.graph <- t1.graph <- t2.graph <- t3.graph <- list()

  
  for(k in 1:length(t1.path$models)){
    A1 <- matrix(0,p,p)
    A1[abs(as.matrix(t1.path$models[[k]]$Theta))>1e-7] <- 1
    diag(A1) <- 0
    t1.graph[[k]] <- A1
    t1.err[[k]] <- max(abs(t1.path$models[[k]]$M-scaled.M))
    }
  png("tmp.png")
  t1.roc <- huge.roc(t1.graph, sim.data$theta, verbose = TRUE)
  dev.off()

  for(k in 1:length(t2.path$models)){
    A2 <- matrix(0,p,p)
    A2[abs(as.matrix(t2.path$models[[k]]$Theta))>1e-7] <- 1
    diag(A1) <- 0
    t2.graph[[k]] <- A2
    t2.err[[k]] <- max(abs(t2.path$models[[k]]$M-scaled.M))
    }
  png("tmp.png")
  t2.roc <- huge.roc(t2.graph, sim.data$theta, verbose = TRUE)
  dev.off()

  for(k in 1:length(t3.path$models)){
    A3 <- matrix(0,p,p)
    A3[abs(as.matrix(t3.path$models[[k]]$Theta))>1e-7] <- 1
    diag(A3) <- 0
    t3.graph[[k]] <- A3
    t3.err[[k]] <- max(abs(t3.path$models[[k]]$M-scaled.M))
    }
  png("tmp.png")
  t3.roc <- huge.roc(t3.graph, sim.data$theta, verbose = TRUE)
  dev.off()
  #diag(t0.path$models[[40]]$W)
  #diag(t1.path$models[[40]]$W)

  #plot(t2.roc$fp,t2.roc$tp,type="b")
  
  
  #gnc0.AUC[ii] <- t0.roc$AUC
  #gnc0.fp <- t0.roc$fp
  #gnc0.tp <- t0.roc$tp
  
  
  gnc3.fp <- t3.roc$fp
  gnc3.tp <- t3.roc$tp
  gnc2.tp <- t2.roc$tp
  gnc2.fp <- t2.roc$fp
  gnc1.fp <- t1.roc$fp
  gnc1.tp <- t1.roc$tp
  tmp <- list(gnc3.fp=gnc3.fp,gnc3.tp=gnc3.tp,gnc1.fp=gnc1.fp,gnc1.tp=gnc1.tp,gnc2.tp=gnc2.tp,gnc2.fp=gnc2.fp)
}
# 
# lw0 <- loess(gnc0.tp~gnc0.fp,degree=1,span=0.7)
# lw3 <- loess(gnc3.tp~gnc3.fp,degree=1,span=0.8)
# lw1 <- loess(gnc1.tp~gnc1.fp,degree=1,span=0.7)
# lw2 <- loess(glasso.tp~glasso.fp,degree=1,span=0.7)
# plot.x <- seq(0,0.002,length.out=1000)
# plot.y0 <- predict(lw0,plot.x)
# plot.y1 <- predict(lw1,plot.x)
# plot.y2 <- predict(lw2,plot.x)
# plot.y3 <- predict(lw3,plot.x)



save(result,file="GridSimulationResults-HighDim_LowSNR_CompareOracle-005.Rda")


# 
# pdf("GridSimulation_WithTwoStage.pdf")
# plot(plot.x,plot.y0,type="l",col="blue",lwd=2,xlab="False Positive",ylab="True Positive",xlim=c(0,0.002),ylim=c(0,1),lty=3)
# lines(plot.x,plot.y1,col="red",lty=1,lwd=2)
# lines(plot.x,plot.y3,col="purple",lty=4,lwd=2)
# lines(plot.x,plot.y2,col="green",lwd=2,lty=2)
# abline(a=0,b=1,lty=2)
# legend("bottomright",c("oracle two-stage gnc-lasso","oracle iterative gnc-lasso","two-stage gnc-lasso","glasso"),lty=c(3,1,4,2),col=c("blue","red","purple","green"),lwd=c(2,2,2,2))
# dev.off()


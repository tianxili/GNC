library(igraph)

library(doParallel)
registerDoParallel(cores=24)
source("GNC-lasso.R")

load("DataForGNC-Plot-Combined.Rda")


CCScale <- function(X,max.iter=20){
    old.X <- X
    diff.seq <- NULL
    for(k in 1:max.iter){
        X1 <- scale(old.X,TRUE,TRUE)
        X <- t(scale(t(X1),TRUE,TRUE))
        diff.seq <- c(diff.seq,norm(X-old.X,"F"))
        old.X <- X
    }
    print(diff.seq)
    return(old.X)
}




x <- X[,1:300]
summary(x[,1:5])


n <- nrow(x)
p <- ncol(x)
dim(x)

x.std <- CCScale(x)

x.reg.std <- scale(x,TRUE,TRUE)
S <- t(x.reg.std)%*%x.reg.std/n

norm(A-t(A)) 
d <- mean(colSums(A))
d

rho.seq <- exp(seq(log(0.001),log(0.02),length.out=30))
alpha.seq <- rho.seq*n/d
set.seed(500)

test.M.cv <- gnc.M.cv(x.std,A,alpha.seq,V=10)
test.M.gcv <- gnc.M.gcv.seq(x.std,A,alpha.seq)
colMeans(test.M.cv$cv.err.mat)
test.M.gcv

test.M.cv$opt.index
alpha.seq




alpha.cv <- alpha.seq[test.M.cv$opt.index]
#alpha.cv <- alpha.seq[which.min(test.M.gcv)]

M.cv <- gnc.M(x.std,A,alpha.cv)

x.centered <- x.std-M.cv

S.cv <- t(x.centered)%*%x.centered/n

lambda.seq <- exp(seq(log(0.15),log(0.3),length.out=100))
#lambda.seq <- exp(seq(log(0.0009),log(0.002),length.out=100))

t1.path <- glassopath(S.cv,rholist=lambda.seq)
t1.num.edge <- rep(0,dim(t1.path$wi)[3])
for(k in 1:dim(t1.path$wi)[3]){
	t1.num.edge[k] <- (sum(abs(t1.path$wi[,,k])>1e-7)-p)/2
}
t1.num.edge





A1 <- matrix(0,p,p)
A1[abs(t1.path$wi[,,49])>1e-7] <- 1
diag(A1) <- 0
g1 <- graph.adjacency(A1,mode="undirected")
V(g1)$label <- colnames(x)
Theta <- t1.path$wi[,,32]
rownames(Theta) <- colnames(Theta) <- colnames(x)

rholist <- exp(seq(log(0.35),log(0.6),length.out=100))
t2.path <- glassopath(S,rholist=rholist)
t2.num.edge <- rep(0,dim(t2.path$wi)[3])
for(k in 1:dim(t2.path$wi)[3]){
	t2.num.edge[k] <- (sum(abs(t2.path$wi[,,k])>1e-7)-p)/2
}
t2.num.edge



A2 <- matrix(0,p,p)
#A2[abs(t2.path$wi[,,Theta.cv2$opt.index])>1e-7] <- 1
A2[abs(t2.path$wi[,,85])>1e-7] <- 1
#A2[abs(t2.path$wi[,,50])>1e-7] <- 1
diag(A2) <- 0
g2 <- graph.adjacency(A2,mode="undirected")
V(g2)$label <- colnames(x)




A0 <- A1+A2
niso <- which(colSums(A0)>0)

A0 <- A0[niso,niso]
A1 <- A1[niso,niso]
A2 <- A2[niso,niso]
labels <- colnames(x[,niso])
g1 <- graph.adjacency(A1,mode="undirected")
V(g1)$label <- labels
g2 <- graph.adjacency(A2,mode="undirected")
V(g2)$label <- labels

g0 <- graph.adjacency(A0,mode="undirected")
lo <- layout.auto(g0)


pdf("CombinedData-NewScaling-Glasso-25.pdf",height=10,width=10)
V(g2)$color <- "orange"
plot(g2,edge.width=3,vertex.size=160,vertex.label.cex=0.75,layout=lo*3,rescale=FALSE,xlim=c(30,55),ylim=c(0,55),vertex.frame.color="orange")
dev.off()

pdf("CombinedData-NewScaling-GNC-25.pdf",height=10,width=10)
V(g1)$color <- "orange"
plot(g1,edge.width=3,vertex.size=160,vertex.label.cex=0.75,layout=lo*3,rescale=FALSE,xlim=c(30,55),ylim=c(0,55),vertex.frame.color="orange")
dev.off()





#gnc <- t1.path5$models[[42]]
M <- M.cv

dim(M)
#M <- scale(M,TRUE,TRUE)

colnames(M) <- colnames(x.std)
rownames(M) <- rownames(X)

M <- M[,niso]


library(MASS)

mydata <- t(M)
d <- dist(mydata)+0.001 # euclidean distances between the rows

fit <- cmdscale(d,eig=TRUE, k=3) # k is the number of dim
fit # view results

# plot solution

pdf("TermMDS.pdf")
set.seed(800)
x <- fit$points[,1]+rnorm(nrow(mydata))*0.35
y <- fit$points[,2]+rnorm(nrow(mydata))*0.35
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
  main="Term	MDS",	type="n",xlim=c(-24,12))
text(x, y, labels = row.names(mydata), cex=.7)
dev.off()

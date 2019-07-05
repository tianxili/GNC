
grid.to.graph <- function(i,j,m){
  return((j-1)*m + i)
}
graph.to.grid <- function(v,m){
  i <- v%%m
  j <- (v-i)/m + 1
  if(i==0) i<- 24
  return(c(i,j))
}


Grid2Mat <- function(i,j,m){
  return(i+(j-1)*m)
}
GridNetDataGen <- function(m){
  n <- m^2
  Adj <- matrix(0,n,n)
  value <- rep(0,n)
  x <- 1:m
  value.mat <- outer(x,x,function(a,b) (a/m)*(b/m))
  
  for(i in 1:m){
    for(j in 1:m){
      v <- grid.to.graph(i,j,m)
      value[v] <- value.mat[i,j]
      if(i>1){
        v.upper <- grid.to.graph(i-1,j,m)
        Adj[v,v.upper] <- 1
      }
      if(i<m){
        v.lower <- grid.to.graph(i+1,j,m)
        Adj[v,v.lower] <- 1
      }
      if(j>1){
        v.left <- grid.to.graph(i,j-1,m)
        Adj[v,v.left] <- 1
      }
      if(j<m){
        v.right <- grid.to.graph(i,j+1,m)
        Adj[v,v.right] <- 1
      }
    }
  }
  return(list(A=Adj,value=value,value.mat=value.mat))
}
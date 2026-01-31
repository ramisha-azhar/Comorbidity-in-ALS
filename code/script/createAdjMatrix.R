createAdjMatrix <- function(network){
  
  source <- network$source
  target <- network$target
  
  nodi <- unique(c(source,target))
  
  N <- length(nodi)
  
  adj <- matrix(0,N,N)
  colnames(adj) <- nodi
  row.names(adj) <- nodi
  
  for(i in 1:N){
    ind <- which(source %in% nodi[i])
    M <- length(ind)
    if ( M > 0 ){
      for (k in 1:M){
        j <- which(nodi %in% target[ind[k]])
        adj[i,j] = 1
        adj[j,i] =  adj[i,j] 
      }
    }
  }
  return(adj)
}


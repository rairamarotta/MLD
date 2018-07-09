normal_var_con <- function(y, m0, C0, V, W, F, G, D){
  
  #Auxilio Desconto
  if(is.null(dim(G)) == FALSE){
    matrixaux <- G == 0 
    if(G[1,1] == 1 & G[1,2] == 1 & G[2,2] == 1 ) {matrixaux[2] <- 0}
    tira <- which(matrixaux == 1)
    mantem <- which(matrixaux == 0)}
  
  #Definindo os objetos
  n <-nrow(F) ; r <-1
  T <- length(y)
  mt <- matrix(0,nrow=n,ncol=T+1)
  Ct <- array(rep(diag(n),T+1),dim=c(n,n,T+1))
  Rt <- array(rep(0,T+1),dim=c(n,n,T+1))
  Pt <- array(rep(0,T+1),dim=c(n,n,T+1))
  Wt <- array(rep(0,T+1),dim=c(n,n,T+1))
  ft <- matrix(0,nrow=T+1,ncol=r)
  at <- matrix(0,nrow=n,ncol=T+1)
  Qt <- matrix(0,nrow=T+1,ncol=r)
  et <- matrix(0,nrow=T+1,ncol=r)
  At <- matrix(0,nrow=n,ncol=T+1)
  
  ## Passo t=1
  
  # Priori em t=1
  if(is.null(D) == TRUE){
    at[,1] <- G%*%m0
    Rt[,,1] <- G%*%C0%*%(t(G)) + W} else{
      if(is.null(dim(D)) == TRUE){
        at[,1] <- G%*%m0
        Pt[,,1] <- G%*%C0%*%(t(G))
        Wt[,,1] <- (D^2 -1)*Pt[,,1]
        Rt[,,1] <- D%*%Pt[,,1]%*%t(D)
      } else{
        at[,1] <- G%*%m0
        Pt[,,1] <- G%*%C0%*%(t(G))
        aux1 <- (D^2)*Pt[,,1]
        aux2 <- Pt[,,1]
        aux.W <- D^2- matrix(1,ncol=n, nrow = n)
        
        Wt[,,1] <-  Pt[,,1]*aux.W
        for(i in 1: length(mantem)){
          aux2[mantem[i]] <- 0}
        Rt[,,1] <- aux1 + aux2}}
  
  # Previsăo 1 passo-a-frente
  
  ft[1,] <- t(F[,1])%*%at[,1]
  Qt[1,] <- t(F[,1])%*%Rt[,,1]%*%F[,1] + V
  
  # Posteriori em t = 1
  
  At[,1] <- Rt[,,1]%*%F[,1]*(1/Qt[1,])
  et[1,] <- y[1]-ft[1,]
  
  mt[,1] <- at[,1]+At[,1]*et[1,]
  Ct[,,1] <-(Rt[,,1]-At[,1]%*%t(At[,1])*Qt[1,])
  
  
  for(t in 2:(T+1)){			
    
    #Priori em t=1
    if(is.null(D) == TRUE){
      at[,t] <- G%*%mt[,t-1]
      Rt[,,t] <- G%*%Ct[,,t-1]%*%(t(G)) + W} else{
        if(is.null(dim(D)) == TRUE){
          at[,t] <- G%*%mt[,t-1]
          Pt[,,t] <- G%*%Ct[,,t-1]%*%(t(G))
          Wt[,,t] <- (D^2 -1)*Pt[,,t]
          Rt[,,t] <- D%*%Pt[,,t]%*%t(D)
        } else{
          at[,t] <- G%*%mt[,t-1]
          Pt[,,t] <- G%*%Ct[,,t-1]%*%(t(G))
          aux1 <- (D^2)*Pt[,,t]
          aux2 <- Pt[,,t]
          Wt[,,t] <-  Pt[,,t]*aux.W
          for(i in 1: length(mantem)){
            aux2[mantem[i]] <- 0}
          Rt[,,t] <- aux1 + aux2}}
    
    # Previsăo 1 passo-a-frente
    ft[t,] <- t(F[,t])%*%at[,t]
    Qt[t,] <- t(F[,t])%*%Rt[,,t]%*%F[,t] + V
    
    # Posteriori em t
    At[,t] <- Rt[,,t]%*%F[,t]*(1/Qt[t,])
    et[t,] <- y[t]-ft[t,]
    mt[,t] <- at[,t]+At[,t]*et[t,]
    Ct[,,t] <- (Rt[,,t]-At[,t]%*%t(At[,t])*Qt[t,])
  }
  
  if(is.null(W) == TRUE){
    result <- list(mt,Ct,ft,Qt, et, Pt, Wt)
    names(result) <- c("mt", "Ct", "ft", "Qt", "et", "Pt", "Wt")} else{
      result <- list(mt,Ct,ft,Qt, et)  
      names(result) <- c("mt", "Ct", "ft", "Qt", "et")
    }
  return(result)
}



# Você pode escolher entrar com W ou D, assim como na funçao de estimaçao
# Observe que se você não especificar W, mas escolher o método aditivo
# Wt será calculado a cada passo t e a estimação será da forma Rt = Pt + Wt.
# Caso o método multiplicativo seja escolhido, Rt = Pt/delta

# calculo base pode ser "aditivo" ou "multiplicativo
funcao.previsao2 <- function(k.passos, mt,Ct, G.prev, F.prev, D, W, calculo.base = "aditivo"){
  
  r <- 1
  at.prev <- matrix(0,nrow=n,ncol=k+1)
  Rt.prev <- array(rep(0,k+1),dim=c(n,n,k+1))
  Pt.prev <- array(rep(0,k+1),dim=c(n,n,k+1))
  Wt.prev <- array(rep(0,k+1),dim=c(n,n,k+1))
  ft.k <-    matrix(0,nrow=k+1,ncol=r)
  Qt.k <-    matrix(0,nrow=k+1,ncol=r)
  
  #Auxilio Desconto
  if(is.null(dim(G.prev)) == FALSE){
    matrixaux <- G.prev == 0 
    if(G.prev[1,1] == 1 & G.prev[1,2] == 1 & G.prev[2,2] == 1 ) {matrixaux[2] <- 0}
    tira <- which(matrixaux == 1)
    mantem <- which(matrixaux == 0)}
  
  # Passo 0
  at.prev[,1] <- mt[, T]
  Rt.prev[,,1] <- Ct[,, T]
  Pt.prev[,,1] <-  Rt.prev[,,1]
  
  # Passo 1
  at.prev[,2] <- G.prev%*%at.prev[,1]
  if(is.null(D) == TRUE){
    Rt.prev[,,2] <- G.prev%*%Rt.prev[,,1]%*%t(G.prev) + W} else{
      if(is.null(dim(D)) == TRUE){
        
        if(calculo.base == "multiplicativo"){
          Rt.prev[,,2]<- D%*%G.prev%*%Rt.prev[,,1]%*%t(G.prev)%*%t(D)} else{
            Pt[,,2] <- G.prev%*%Rt.prev[,,1]%*%t(G.prev)
            aux.W <- D^2 - 1
            aux1 <- (D^2)*Pt[,,2]
            aux2 <- Pt[,,2]
            Wt[,,2] <-  Pt[,,2]*aux.W
            Rt.prev[,,2] <-  Pt[,,2] +  Wt[,,2]
          }
      } else {
        
        if(calculo.base == "multiplicativo"){
          Pt.prev[,,2] <- G.prev%*%Rt.prev[,,1]%*%t(G.prev)
          aux.W <- D^2- matrix(1,ncol=n, nrow = n)
          aux1 <- (D^2)*Pt.prev[,,2]
          aux2 <- Pt.prev[,,2]
          Wt.prev[,,2] <-  Pt.prev[,,2]*aux.W
          for(i in 1: length(mantem)){
            aux2[mantem[i]] <- 0}
          Rt.prev[,,2] <- aux1 + aux2} else  {
            Pt.prev[,,2] <- G.prev%*%Rt.prev[,,1]%*%t(G.prev)
            aux.W <- D^2 - 1
            aux1 <- (D^2)*Pt.prev[,,2]
            aux2 <- Pt.prev[,,2]
            Wt.prev[,,2] <-  Pt.prev[,,2]*aux.W
            Rt.prev[,,2] <-  Pt.prev[,,2] +  Wt.prev[,,2]
          }
        
      }
      
      
    }
  ft.k[2,] <-  t(F.prev[,1])%*%at.prev[,2] 
  Qt.k[2,] <-  t(F.prev[,1])%*%Rt.prev[,,2]%*%F.prev[,1] + V 
  
  for(i in 3:(k+1)){
    at.prev[,i] <- G.prev%*%at.prev[,i-1]
    if(is.null(D) == TRUE){
      Rt.prev[,,i] <- G.prev%*%Rt.prev[,,i-1]%*%t(G.prev) + W} else{
        if(is.null(dim(D)) == TRUE){
          
          if(calculo.base == "multiplicativo"){
            Rt.prev[,,i]<- D%*%G.prev%*%Rt.prev[,,i-1]%*%t(G.prev)%*%t(D)} else{
              Pt.prev[,,i] <- G.prev%*%Rt.prev[,,i-1]%*%t(G.prev)
              aux.W <- D^2 - 1
              aux1 <- (D^2)*Pt.prev[,,i]
              aux2 <- Pt[,,i]
              Wt.prev[,,i] <-  Pt[,,2]*aux.W
              Rt.prev[,,i] <-  Pt.prev[,,i] +  Wt.prev[,,i]
            }
        } else {
          
          if(calculo.base == "multiplicativo"){
            Pt.prev[,,i] <- G.prev%*%Rt.prev[,,i-1]%*%t(G.prev)
            aux.W <- D^2- matrix(1,ncol=n, nrow = n)
            aux1 <- (D^2)*Pt.prev[,,i]
            aux2 <- Pt.prev[,,i]
            Wt.prev[,,i] <-  Pt.prev[,,i]*aux.W
            for(j in 1: length(mantem)){
              aux2[mantem[j]] <- 0}
            Rt.prev[,,i] <- aux1 + aux2} else {
              Pt.prev[,,i] <- G.prev%*%Rt.prev[,,i-1]%*%t(G.prev)
              aux.W <- D^2 - 1
              aux1 <- (D^2)*Pt.prev[,,i]
              aux2 <- Pt.prev[,,i]
              Wt.prev[,,i] <-  Pt.prev[,,i]*aux.W
              Rt.prev[,,i] <-  Pt.prev[,,i] +  Wt.prev[,,i]
            }
          
        }
        
        
      }
    ft.k[i,] <-  t(F.prev[,i-1])%*%at.prev[,i] 
    Qt.k[i,] <-  t(F.prev[,i-1])%*%Rt.prev[,,i]%*%F.prev[,i-1] + V 
  }
  
  aux.prev <- list(at.prev, Rt.prev, ft.k, Qt.k)
  names(aux.prev) <- c("at.prev", "Rt.prev", "ft.k", "Qt.k")
  return(aux.prev)
}
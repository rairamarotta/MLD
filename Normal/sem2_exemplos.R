# Códigos da Semana 2 -------
# Autora: Raíra Marotta
# e-mail: raira@dme.ufrj.br
# Estes códigos encontram-se em https://github.com/rairamarotta


# W/V ------------

gera.serie <- function(T, nivel.mu,Wt,Vt, plot = TRUE){
  set.seed(12345)
  mu <- NULL; mu[1] <- nivel.mu
  for(i in 2:T){
    mu[i] <- mu[i-1] + rnorm(1,0, sqrt(Wt))}
  yt <- rnorm(T, mu, sqrt(Vt))
  
  if(plot == TRUE){
    plot(yt, type = "l", ylim=c(min(yt) - 3, max(yt) + 3), ylab = " ", xlab = " ")
    lines(mu, lty = 2, col = 4)
    legend("topleft", lty = c(1,2), col =c(1,4), legend = c(expression(y[t]), expression(mu[t])), bty = "n")}
  return(data.frame(mu,yt))
}

T <- 100; nivel.mu <- 25; Wt = 0.05; Vt = 1
dados <- gera.serie(T, nivel.mu,Wt,Vt)

T <- 100; nivel.mu <- 25; Wt = 0.5; Vt = 1
dados2 <- gera.serie(T, nivel.mu,Wt,Vt)


# Modelo Normal ------
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
# Modelo de Primeira ordem ---------

library(dlm)
library(readxl)
candy <- read_excel("Migon/candy.xlsx")
y <- candy$sales

m0 <- 0; C0<- 100; V <-0.1; W<- 0.01

# Definindo a matriz F
F1 <- matrix(1,nrow=1,ncol=(length(y)+1))
F <- F1
F

resultados <- normal_var_con(y, m0, C0, V, W, F, G = 1, D = NULL)

# Estimação do nível

LS <- qnorm(0.975,resultados$mt[1,],sqrt(resultados$Ct[1,1,]))
LI <- qnorm(0.025,resultados$mt[1,],sqrt(resultados$Ct[1,1,]))
LS <- LS[1:length(LS)-1]
LI <- LI[1:length(LI)-1]

plot(y,pch = 20, ylab = " ", xlab = " ", ylim= c(min(LI) - 1, max(LS)+1))
lines(resultados$mt[1,], col = 2)
lines(LS, col = 4, lty = 2)
lines(LI, col = 4, lty = 2)


# Previsão a 1 passo
LS <- dropFirst(qnorm(0.975, resultados$ft, sqrt(resultados$Qt)))
LI <- dropFirst(qnorm(0.025, resultados$ft, sqrt(resultados$Qt)))

plot(y,pch = 20, ylab = " ", xlab = " ", ylim= c(0,15))
lines(dropFirst(resultados$ft), col = 2)
lines(LS, col = 4, lty = 2)
lines(LI, col = 4, lty = 2)
legend("topleft", legend = c(expression(y[t]),expression(f[t])),
       bty = "n", pch = c(20, NA), lty = c(NA,1), col = c(1,2),
       cex = 1.2)


# Modelo de Segunda  ordem ---------

# Definindo a matriz F
n <- 2
F1 <- matrix(c(1,0),nrow=n,ncol=(length(y)+1))
F <- F1
F

# Definindo a matriz G
G <- matrix(c(1,0,1,1),2,2)
G

# Definindo as quantidades iniciais
m0 <- c(0,0)
C0 <- diag(100,n,n) 
W = diag(c(0.02, 0.01))

resultados <- normal_var_con(y, m0, C0, V, W, F, G , D = NULL)

# Estimação do nível

LS <- qnorm(0.975,resultados$mt[1,],sqrt(resultados$Ct[1,1,]))
LI <- qnorm(0.025,resultados$mt[1,],sqrt(resultados$Ct[1,1,]))
LS <- LS[1:length(LS)-1]
LI <- LI[1:length(LI)-1]

plot(y,pch = 20, ylab = " ", xlab = " ", ylim= c(min(LI) - 1, max(LS)+1))
lines(resultados$mt[1,], col = 2)
lines(LS, col = 4, lty = 2)
lines(LI, col = 4, lty = 2)

# Previsão a 1 passo
LS <- dropFirst(qnorm(0.975, resultados$ft, sqrt(resultados$Qt)))
LI <- dropFirst(qnorm(0.025, resultados$ft, sqrt(resultados$Qt)))

plot(y,pch = 20, ylab = " ", xlab = " ", ylim= c(0,15))
lines(dropFirst(resultados$ft), col = 2)
lines(LS, col = 4, lty = 2)
lines(LI, col = 4, lty = 2)
legend("topleft", legend = c(expression(y[t]),expression(f[t])),
       bty = "n", pch = c(20, NA), lty = c(NA,1), col = c(1,2),
       cex = 1.2)

# Fator de crescimento
LS <- qnorm(0.975,resultados$mt[2,],sqrt(resultados$Ct[2,2,]))
LI <- qnorm(0.025,resultados$mt[2,],sqrt(resultados$Ct[2,2,]))
LS <- LS[1:length(LS)-1]
LI <- LI[1:length(LI)-1]

plot(resultados$mt[2,], type = "l", col =2, ylim= c(-2,2))
lines(LS, col = 4, lty = 2)
lines(LI, col = 4, lty = 2)


# Descontos -------------

# Desconto 1
m0 <- 0; C0<- 100; V <-0.1; W = NULL
F <- matrix(1,nrow=1,ncol=(length(y)+1))
G <- 1
D <- 1
resultados <- normal_var_con(y, m0, C0, V, W, F, G , D)

# Previsão a 1 passo
plot(y,pch = 20, ylab = " ", xlab = " ", ylim= c(0,15))
lines(dropFirst(resultados$ft), col = 2)
legend("topleft", legend = c(expression(y[t]),expression(f[t])),
       bty = "n", pch = c(20, NA), lty = c(NA,1), col = c(1,2),
       cex = 1.2)


# Desconto 0.9
m0 <- 0; C0<- 100; V <-0.1; W = NULL
F <- matrix(1,nrow=1,ncol=(length(y)+1))
G <- 1
D <- 1/sqrt(0.9)
resultados <- normal_var_con(y, m0, C0, V, W, F, G , D)

# Previsão a 1 passo

plot(y,pch = 20, ylab = " ", xlab = " ", ylim= c(5,13))
lines(dropFirst(resultados$ft), col = 2)
legend("topleft", legend = c(expression(y[t]),expression(f[t])),
       bty = "n", pch = c(20, NA), lty = c(NA,1), col = c(1,2),
       cex = 1.2)

# Desconto 0.8
m0 <- 0; C0<- 100; V <-0.1; W = NULL
F <- matrix(1,nrow=1,ncol=(length(y)+1))
G <- 1
D <- 1/sqrt(0.8)
resultados <- normal_var_con(y, m0, C0, V, W, F, G , D)

# Previsão a 1 passo

plot(y,pch = 20, ylab = " ", xlab = " ", ylim= c(5,13))
lines(dropFirst(resultados$ft), col = 2)
legend("topleft", legend = c(expression(y[t]),expression(f[t])),
       bty = "n", pch = c(20, NA), lty = c(NA,1), col = c(1,2),
       cex = 1.2)


# Regressão múltipla ---------

# Figura comparativa
par(mfrow=c(2,1), mar = c(4,2,2,1))
candydata <- ts(candy,start=c(1976,1), freq = 12)
plot(candydata[,3], type = "l", main = "Sales", xlab = " ")
plot(candydata[,2], type = "l", main = "Price", xlab = " ")

#X-Y
par(mfrow=c(1,1), mar = c(4,4,2,1))
plot(y,candy$price, xlab = "Price", ylab = "Sales", pch = 20)

# Modelo
n <- 2
m0  <- c(0,0); C0<- diag(100,2); V <-0.1; W = NULL
T <- length(y)

#Definindo F
F <- matrix(c(rep(1, length(y)),candy$price),ncol=T, nrow=n, byrow=T)
F <- cbind(F, c(1,1)) #gambiarra no código (pq antes estava indo até t+1)

#Definindo G
G <- diag(1,2)

#Definindo D
D1 <- 1/sqrt(0.9)
D2 <- 1/sqrt(0.99)
D<- bdiag(D1,D2)

# Resultados
resultados <- normal_var_con(y,m0, C0, V, W, F, G , D)

# Previsao
LS <- dropFirst(qnorm(0.975,resultados$ft ,sqrt(resultados$Qt)))
LI <- dropFirst(qnorm(0.025,resultados$ft ,sqrt(resultados$Qt)))

plot(y,pch = 20, ylab = " ", xlab = " ", ylim= c(3,15))
lines(dropFirst(resultados$ft), col = 2)
lines(LS, col = 4, lty = 2)
lines(LI, col = 4, lty = 2)
legend("topleft", legend = c(expression(y[t]),expression(f[t])),
       bty = "n", pch = c(20, NA), lty = c(NA,1), col = c(1,2),
       cex = 1.2)

# Efeito do Preço

LS <- qnorm(0.975,resultados$mt[2,],sqrt(resultados$Ct[2,2,]))
LI <- qnorm(0.025,resultados$mt[2,],sqrt(resultados$Ct[2,2,]))
LS <- LS[1:length(LS)-1]
LI <- LI[1:length(LI)-1]

par(mfrow=c(2,1), mar = c(4,2,2,1))
plot(resultados$mt[2,], type = "l", col =2, ylim=c(-5,3),
     xlab = expression(paste("E[",beta[t],"|D"[t], "]")))
lines(LS, col = 4, lty = 2)
lines(LI, col = 4, lty = 2)

# Efeito do preço x Preço
plot(na.omit(resultados$mt[2,])*candy$price, type = "l", col =2, ylim=c(-5,3),
     xlab = expression(paste("E[",beta[t],"|D"[t], "] . x"[t])))
lines((candy$price)*LS, col = 4, lty = 2)
lines((candy$price)*LI, col = 4, lty = 2)

# Previsao k passos  ------

# Você pode escolher entrar com W ou D, assim como na funçao de estimaçao
# Observe que se você não especificar W, mas escolher o método aditivo
# Wt será calculado a cada passo t e a estimação será da forma Rt = Pt + Wt.
# Caso o método multiplicativo seja escolhido, Rt = Pt/delta

# calculo base pode ser "aditivo" ou "multiplicativo
funcao.previsao <- function(k.passos, mt,Ct, G.prev, F.prev, D, W, calculo.base = "aditivo"){
  
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
k <- 12
W = diag(0.01,2)
G.prev <- matrix(c(1,0,0,0), ncol = n)
F.prev <- matrix(rep(c(1,0), k), ncol = k, nrow = n)

previsao <- funcao.previsao(k.passos = k,
                            mt = resultados$mt,
                            Ct = resultados$Ct, G.prev, F.prev, D, W = NULL,
                            calculo.base = "multiplicativo")

LS <- c(y,dropFirst(qnorm(0.975,previsao$ft.k,sqrt(previsao$Qt.k))))
LI <- c(y,dropFirst(qnorm(0.025,previsao$ft.k,sqrt(previsao$Qt.k))))
y.prev <- c(y, rep(NA, k))
par(mfrow=c(1,1))
plot(y.prev, pch = 20, ylim=c(6,13))
points(c(rep(NA, length(y)),previsao$ft.k), pch = 20, col = 2)
lines(LS, lty = 2, col = 2)
lines(LI, lty = 2, col = 2)


# Modelo com 1 par de harmonicos ------


# Definindo as quantidades
n <- 4
w <- (2*pi)/(12)
G1 <- diag(1,2)
G2 <- matrix(c(cos(w),-sin(w),sin(w),cos(w)),2,2)
G <- dlm::bdiag(G1, G2)

F <- matrix(c(rep(1, length(y)),candy$price, 
              rep(1, length(y)),  rep(0, length(y))),ncol=T, nrow=n, byrow=T)
F <- cbind(F, c(1,1,1,1)) #gambiarra no código (pq antes estava indo até t+1)

#Desconto da tendência
D1 <- 1/sqrt(0.97)

#Desconto da covariável
D2 <- 1/sqrt(0.99)

#Desconto da sazonalidade
D3 <- matrix(1/sqrt(0.99),2,2)
D <- dlm::bdiag(D1,D2, D3)

F.prev <- F
F.prev[2,] <- 0 #Não temos mais dados da covariável, por isso, zero
F.prev <- F.prev[,1:12]
G.prev <- G

m0 <- rep(0,n)
C0 <- diag(100,n)

#Resultados
resultados <- normal_var_con(y,m0, C0, V, W, F, G , D)

#Previsão
LS <- dropFirst(qnorm(0.975,resultados$ft ,sqrt(resultados$Qt)))
LI <- dropFirst(qnorm(0.025,resultados$ft ,sqrt(resultados$Qt)))

plot(y,pch = 20, ylab = " ", xlab = " ", ylim= c(3,15))
lines(dropFirst(resultados$ft), col = 2)
lines(LS, col = 4, lty = 2)
lines(LI, col = 4, lty = 2)
legend("topleft", legend = c(expression(y[t]),expression(f[t])),
       bty = "n", pch = c(20, NA), lty = c(NA,1), col = c(1,2),
       cex = 1.2)

# Efeito do harmonico

LS <- qnorm(0.975,resultados$mt[3,],sqrt(resultados$Ct[3,3,]))
LI <- qnorm(0.025,resultados$mt[3,],sqrt(resultados$Ct[3,3,]))
LS <- LS[1:length(LS)-1]
LI <- LI[1:length(LI)-1]

plot(resultados$mt[3,], type = "l", col =2, ylim= c(-2,2), xlab = " ", ylab = " ")
lines(LS, col = 4, lty = 2)
lines(LI, col = 4, lty = 2)

# Previsao k passos

prev.har <- funcao.previsao(k.passos = k,
                            mt = resultados$mt,
                            Ct = resultados$Ct, G.prev, F.prev, D, W = NULL,
                            calculo.base = "multiplicativo")

LS <- c(y,dropFirst(qnorm(0.95,prev.har$ft.k,sqrt(prev.har$Qt.k))))
LI <- c(y,dropFirst(qnorm(0.05,prev.har$ft.k,sqrt(prev.har$Qt.k))))
y.prev <- c(y, rep(NA, k))
plot(y.prev, pch = 20, ylim=c(6,13), xlab = " ", ylab = " ")
points(c(rep(NA, length(y)),prev.har$ft.k), pch = 20, col = 2)
lines(LS, lty = 2, col = 2)
lines(LI, lty = 2, col = 2)



# Modelo com 2 harmonicos ------


# Definindo as quantidades
n <- 6
w <- (2*pi)/(12)
G2 <- matrix(c(cos(w),-sin(w),sin(w),cos(w)),2,2)
G3 <- matrix(c(cos(2*w),-sin(2*w),sin(2*w),cos(2*w)),2,2)
G <- dlm::bdiag(G1, G2, G3)


F <- matrix(c(rep(1, length(y)),candy$price, 
              rep(1, length(y)),  rep(0, length(y)),
              rep(1, length(y)),  rep(0, length(y))),ncol=T, nrow=n, byrow=T)
F <- cbind(F, c(1,1,1,1,1,1)) #gambiarra no código (pq antes estava indo até t+1)

#Desconto da tendência
D1 <- 1/sqrt(0.97)

#Desconto da covariável
D2 <- 1/sqrt(0.99)

#Desconto da sazonalidade
D3 <- matrix(1/sqrt(0.99),2,2)
D4 <- matrix(1/sqrt(0.99),2,2)

#Desconto final
D <- dlm::bdiag(D1,D2,D3,D4)

#Previsão
F.prev <- F
F.prev[2,] <- 0 #Não temos mais dados da covariável, por isso, zero
F.prev <- F.prev[,1:12]
G.prev <- G

m0 <- rep(0,n)
C0 <- diag(100,n)

#Resultados
resultados <- normal_var_con(y,m0, C0, V, W, F, G , D)

#Previsão
LS <- dropFirst(qnorm(0.975,resultados$ft ,sqrt(resultados$Qt)))
LI <- dropFirst(qnorm(0.025,resultados$ft ,sqrt(resultados$Qt)))

plot(y,pch = 20, ylab = " ", xlab = " ", ylim= c(3,15))
lines(dropFirst(resultados$ft), col = 2)
lines(LS, col = 4, lty = 2)
lines(LI, col = 4, lty = 2)
legend("topleft", legend = c(expression(y[t]),expression(f[t])),
       bty = "n", pch = c(20, NA), lty = c(NA,1), col = c(1,2),
       cex = 1.2)

# Efeito do harmonico

LS <- qnorm(0.975,resultados$mt[3,],sqrt(resultados$Ct[3,3,]))
LI <- qnorm(0.025,resultados$mt[3,],sqrt(resultados$Ct[3,3,]))
LS <- LS[1:length(LS)-1]
LI <- LI[1:length(LI)-1]

plot(resultados$mt[3,], type = "l", col =2, ylim= c(-2,2), xlab = " ", ylab = " ")
lines(LS, col = 4, lty = 2)
lines(LI, col = 4, lty = 2)

# Previsao k passos

prev.har <- funcao.previsao(k.passos = k,
                            mt = resultados$mt,
                            Ct = resultados$Ct, G.prev, F.prev, D, W = NULL,
                            calculo.base = "multiplicativo")

LS <- c(y,dropFirst(qnorm(0.95,prev.har$ft.k,sqrt(prev.har$Qt.k))))
LI <- c(y,dropFirst(qnorm(0.05,prev.har$ft.k,sqrt(prev.har$Qt.k))))
y.prev <- c(y, rep(NA, k))
plot(y.prev, pch = 20, ylim=c(6,13), xlab = " ", ylab = " ")
points(c(rep(NA, length(y)),prev.har$ft.k), pch = 20, col = 2)
lines(LS, lty = 2, col = 2)
lines(LI, lty = 2, col = 2)


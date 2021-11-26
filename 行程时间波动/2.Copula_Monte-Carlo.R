library(fitdistrplus)
library(Hmisc)
library(fBasics)
library(copula)
library(riskSimul)
library(ghyp)


speedlog <- read.csv('C:\\Users\\g\\Desktop\\博论书写\\论文3\\TTV\\data\\case2_ttv_data.csv')
speedlog <- as.matrix(speedlog)

# ReturnCopula calculates the portfolio return for N replications of
# multinormal Z and gamma input Y. (Z: D*N, Y: length of N)

ReturnCopula<-function(Z,  Y,	portfobj){
  D <- dim(Z)[1]
  N <- dim(Z)[2]                                   
  nu <- portfobj[[4]]
  L <- portfobj[[5]]
  c <- portfobj[[6]]
  w <- portfobj[[7]] # analysis of parameters
  T <- L%*%Z/matrix(sqrt(Y/nu),D,N,byrow=TRUE)         # generating multi-T, %*%:multiplies two matrices
  marg <- matrix(0,D,N)
  for(i in 1:D){
    if(portfobj[[1]][i]=="t"){
      tpmg <- portfobj[[2]][i,]
      marg[i,] <- qt(pt(T[i,],nu),tpmg[3])*tpmg[2]+tpmg[1]  # generate t marginals of t-copula
      print(c(i,"t","is done!"))
      
    }else if(portfobj[[1]][i]=="GH"){
      gpmg <- portfobj[[3]][i,]
      gen <- pinvd.new(distr=udghyp(lambda=gpmg[1], 
                                    alpha=gpmg[2], 
                                    beta=gpmg[3], 
                                    delta=gpmg[4], 
                                    mu=gpmg[5]))
      marg[i,] <- uq(gen,pt(T[i,],nu))                     # generate GH marginals of t-copula
      print(c(i,"GH","is done!"))
    }
  }
  as.vector(t(w)%*%exp(c*marg))                             # weighted sum of log-scaled-marginals (portfolio return)
}


NVCopula<-function(N,portfobj, clock=F){		
  if(clock) start<-proc.time()                       # timer of the algorithm
  D<-length(portfobj[[5]])
  nu<-portfobj[[3]]                                                    
  if(portfobj[[1]]=="GH"){
    gen<-list()
    pmg <- portfobj[[2]]
    gen<-list()
    for(d in 1:D){
      gen[[d]]<-pinvd.new(distr=udghyp(lambda=pmg[d,1], alpha=pmg[d,2], beta=pmg[d,3], delta=pmg[d,4], mu=pmg[d,5]))	
    }
    portfobj[[2]]<-gen
  }
  Z<-matrix(rnorm(D*N),D,N)                         # multi normal input
  Y<-rgamma(N,shape=nu/2,scale=2)                   # chi-squared input
  return<-ReturnCopula(Z,Y,portfobj)                # copula return
  rr <- sort(return)
  LVAR <- rr[length(return)*0.05-1]
  GVAR <- rr[length(return)*0.95]
  ETTR <- mean(rr[(N*0.05):(N*0.95)])
  
  res<-rbind(log(LVAR),log(GVAR),log(ETTR))
  if(clock){
    finish<-proc.time()
    list(res,finish[3]-start[3])                    # return results with timing
  }else{
    list(res)                                       # return results without timing
  }
}

## calculate the LVAR,GVAR and ETTR with t marginals of T-copula 
preday <- function(dd){
  pmg<- matrix(NA,ncol=3,nrow=5)  
  colnames(pmg) <- c("mu","sigma","nu")
  rowi = 1
  for (i in names(speedlog)) {
    x1 = speedlog[i]
    x1 = x1[0:dd,]
    ft <- fitdistr(x1,'t',df=3)
    pmg[rowi,] <- c(ft$estimate[1],ft$estimate[2],5)
    rowi = rowi+1
  }
  
  M1 <- as.matrix(speedlog[1:dd,1:5])
  #u1 <- speedlog
  u1 <- pobs(M1) # convert to rank,and then normalization
  param_len <- ncol(u1)*(ncol(u1)-1)/2
  df <- 6
  tc.mpl <- fitCopula(tCopula(dim=ncol(u1), dispstr="un"),u1, method="ml", estimate.variance=TRUE,start=c(rep(0,param_len),df)) ## start = (rho[1:3], df)
  #summary(tc.mpl)
  
  R <- matrix(
    c(rep(1,25)),ncol = 5)
  
  # rewrite the up-triangular matrix
  param_ind <- 1
  for (i in 1:4) { # row
    for (j in (i+1):5) { # col
      R[i,j] <- tc.mpl@copula@parameters[param_ind]
      param_ind <- param_ind+1
    }
  }
  # rewrite the low-triangular matrix
  param_ind <- 1
  for (i in 1:4) { # col
    for (j in (i+1):5) {  # row
      R[j,i] <- tc.mpl@copula@parameters[param_ind]
      param_ind <- param_ind+1
    }
  }
  nu = tc.mpl@estimate[11]
  portfo <- new.portfobj(nu=nu,R=R,typemg="t",parmg=pmg,c=rep(1,5),w=c(0.195, 0.098, 0.168, 0.238, 0.301)) 
  
  nv_res2 <- NVCopula(N=10^4,portfobj=portfo)
  return(nv_res2)
}

# 閹笛嗩攽妫板嫭绁撮崙鑺ユ殶
res <- matrix(c(rep(1,62)), ncol = 3) #閸掓繂顫愰崠鏍波閺嬫粎鐓╅梼?
res_row <- 0
for (d in 44:(nrow(speedlog)-1)){
  resu = preday(d)
  res_row = res_row+1
  res[res_row,1]=resu[[1]][1]
  res[res_row,2]=resu[[1]][2]
  res[res_row,3]=resu[[1]][3]
}


### *_*_*_*_*_*_*_*_*_*_*_*_ Test sample *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_###
dd=62
t_para = 3
N_para = 2
GH_para = 5
link_nums = 12

## parameters for T marginal
t_pmg<- matrix(NA,ncol=t_para,nrow=link_nums)  
colnames(t_pmg) <- c("mu","sigma","nu")
rowi = 1
for (i in 1:link_nums) {
  x1 = speedlog[0:dd,i]
  ft <- fit.tuv(x1,nu=5)
  #ft <- fitdistr(x1,'t',df=8)
  t_pmg[rowi,] <- c(ft@mu,ft@sigma[1,1],-2*ft@lambda)
  rowi = rowi+1
}

## parameters for Normal margin
n_pmg<- matrix(NA,ncol=N_para,nrow=link_nums)  
colnames(n_pmg) <- c("mu","sigma")
rowi = 1
for (i in 1:link_nums) {
  x1 = speedlog[0:dd,i]
  fn <- fit.gaussuv(x1, na.rm = T, save.data = T)
  #fn <- fitdist(x1,'norm')
  n_pmg[rowi,] <- c(fn@mu, fn@sigma[1,1])
  rowi = rowi+1
}

## parameters for GH margin
g_pmg<- matrix(NA,ncol=GH_para,nrow=link_nums)  
colnames(g_pmg) <- c("lambda","alpha","beta","delta","mu")
rowi = 1
for (i in 1:link_nums) {
  x1 = speedlog[0:dd,i]
  fg <- fit.ghypuv(x1, lambda = 1, alpha.bar = 0.5, mu = mean(x1), sigma = mad(x1), gamma = 0.5)
  #print('fg is done')
  g_pmg[rowi,] <- c(fg@lambda, fg@alpha.bar, fg@gamma, fg@sigma, fg@mu)
  rowi = rowi+1
}

## Correlation matrix based copula
M1 <- speedlog[0:dd,]
#u1 <- speedlog
u1 <- pobs(M1) # convert to rank,and then normalization
param_len <- ncol(u1)*(ncol(u1)-1)/2
df <- 20
tc.mpl <- fitCopula(tCopula(dim=ncol(u1), dispstr="un"),
                    u1, 
                    method="ml", 
                    estimate.variance=TRUE,
                    start=c(rep(0,param_len),df)) ## start = (rho[1:3], df)
summary(tc.mpl)

# ralation matrix 
R <- matrix(c(rep(1,link_nums^2)),ncol = link_nums)

# rewrite the up-triangular matrix
param_ind <- 1
for (i in 1:(link_nums-1)) { # row
  for (j in (i+1):link_nums) { # col
    R[i,j] <- tc.mpl@copula@parameters[param_ind]
    param_ind <- param_ind+1
  }
}
# rewrite the low-triangular matrix
param_ind <- 1
for (i in 1:(link_nums-1)) { # col
  for (j in (i+1):link_nums) {  # row
    R[j,i] <- tc.mpl@copula@parameters[param_ind]
    param_ind <- param_ind+1
  }
}


## Parameter portfolio 
nu = tail(tc.mpl@estimate,1)  # the last one param
m1 = c("t","t","t","t","t","Norm","t","t","Norm","t","t","GH","GH","t","t","t")
m2 = c("t","t","t","t","t","t","t","t","t","t","t","t")
m2 = c("Norm","Norm","Norm","Norm","Norm","Norm","Norm","Norm","Norm","Norm","Norm","Norm")

portfo <- list(typemg=m2,
               parmg1=t_pmg,
               parmg2=g_pmg,
               parmg3=n_pmg,
               nu=nu,
               R=t(chol(R)), #Cholesky decomposition
               c=rep(1,link_nums),
               w=rep(1/link_nums,link_nums))


## execute simulation
MC_simul <- function(N,portfobj, clock=F){		
  if(clock) start<-proc.time()                                 # timer of the algorithm
  D <- length(portfobj[[7]])
  nu <- portfobj[[5]]  
  L <- portfobj[[6]]
  c <- portfobj[[7]]
  w <- portfobj[[8]] # analysis of parameters
  Z <- matrix(rnorm(D*N),D,N)                                  # random generator for multi normal distribution
  Y<-rgamma(N,shape=nu/2,scale = 2)                             # random generator for gamma distribution
  #Y <- rchisq(N,df=nu,ncp = 0)   
  # random generator for chi-squared distribution
  T <- L%*%Z/matrix(sqrt(Y/nu),D,N,byrow=TRUE)         # generating multi-T, %*%:multiplies two matrices
  marg <- matrix(0,D,N)
  for(i in 1:D){
    if(portfobj[[1]][i]=="t"){
      tpmg <- portfobj[[2]][i,]
      marg[i,] <- qt(pt(T[i,],nu),tpmg[3])*tpmg[2]+tpmg[1]  # generate t marginals of t-copula
      print(c(i,"t","is done!"))
      
    }else if(portfobj[[1]][i]=="GH"){
      gpmg <- portfobj[[3]][i,]
      gen <- pinvd.new(distr=udghyp(lambda=gpmg[1], 
                                    alpha=gpmg[2], 
                                    beta=gpmg[3], 
                                    delta=gpmg[4], 
                                    mu=gpmg[5]))
      marg[i,] <- uq(gen,pt(T[i,],nu))                     # generate GH marginals of t-copula
      print(c(i,"GH","is done!"))
    }else if(portfobj[[1]][i]=="Norm"){
      npmg <- portfobj[[4]][i,]
      marg[i,] <- qnorm(pnorm(Z[i,]),mean=npmg[1], sd=npmg[2])      # generate Norm marginals of t-copula
      print(c(i,"Norm","is done!"))
    }
  }
  as.vector(t(w)%*%exp(c*marg))                         # weighted sum of log-scaled-marginals (portfolio return)
}

simul <- MC_simul(N=10^4,portfobj=portfo)
hist(simul[simul<5],breaks=10)
length(simul[simul<40])
simul <- as.data.frame(simul)
write.csv(simul,"C:\\Users\\g\\Desktop\\博论书写\\论文3\\TTV\\data\\return2.csv",row.names = FALSE)
#nv_res2 <- NVCopula(N=10^4,portfobj=portfo)
#nv_res2





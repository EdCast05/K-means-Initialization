##*****************************************
##*
##* @file: ELpp.R
##* 
##* 
##* Perform EL++. Return centers than 
##* can be used for initialization for k-means.
##* 
##* Author:
##* 
##* Eduar A. Castañeda Molina, Student
##* Israel A. Almodóvar-Rivera, PhD
##* University of Puerto Rico at Mayagüez
##*
##* 
##*
##*****************************************
##*
elppv2 <- function(X,K){  ## Data X, K desired centers
  library(emplik)
  n <- nrow(X)
  idx <- sample(1:n,size=1) ## First centroid choosen with Prob = 1/n
  c1 <- X[idx,]
  M1 <- t(c1)              ## call it mu_1
  Prob.Step <- matrix(0,nrow=n,ncol=K-1) ## store probability values
  while((nrow(M1) < K)){
    if(nrow(M1)==1){
      psi <- sapply(1:n,function(i) apply(t(X[i,]-c1),1,norm,type="2"))
      Mpsi <- max(psi)     ## Maximun of distances between data and mu_1
      EL <- el.test(psi,mu=Mpsi)
      pr <- EL$wts      ## Pr(choose) weights of el.test
      idx2 <- sample(1:n,size=1,prob = pr) ##next center is chosen with Prob = pr
      M1 <- rbind(M1,X[idx2,])  
      Prob.Step[,1] <- pr
    } else{
      m <- nrow(M1)        ## new number of means
      psi <- sapply(1:n,function(i) min(sapply(1:m,function(k) apply(t(X[i,]-M1[k,]),1,norm,type="2"))))
      Mpsi <- max(psi)
      EL <- el.test(psi,mu=Mpsi)
      pr <- as.vector(EL$wts)         ## weights of el.test
      idx2 <- sample(1:n,size=1,prob = pr)
      M1 <- rbind(M1,X[idx2,])
      Prob.Step[,m] <- pr
    }
  }
  list(Means = M1, Prob = Prob.Step)
}













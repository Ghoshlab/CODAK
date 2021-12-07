
############### Functions for cell type abundance comparison ################


library("MASS")
library("coda.base")
library("lme4")


#### Function for the KDC test ####


dcor.comp <- function(propmat, xvec, num.perm = 1e4, 
                      Xdist = "categorical", result = "test"){
  
  # propmat: n x num.cat matrix of proportions (n = sample size, num.cat = number of cell types)
  # xvec: an n x 1 vector containing the values of the predictor
  # num.perm: number of permutations
  # Xdist: specifying the kernel used for the predictor. 
  #       "categorical" - Hamming distance kernel, "LK" - Linear kernel, "Gaussian" - Gaussian kernel
  # result: "test" - provides test results, "cor" - provides value of dcor
  
  d <- dist(propmat, method = 'aitchison')
  
  alldist <- as.numeric(d)
  K <- exp(-as.matrix(d)/median(alldist))
  
  if (Xdist == "categorical"){
    d2 <- outer(xvec, xvec, 
                FUN = function(x,y){1 - as.numeric(x == y)})
    K2 <- exp(-d2)
  }
  
  if (Xdist == "Gaussian"){
    d2 <- outer(xvec, xvec, 
                FUN = function(x,y){abs(x-y)})
    K2 <- exp(-d2/median(as.numeric(d2)))
  }
  
  if (Xdist == "LK"){
    K2 <- outer(xvec, xvec, 
                FUN = function(x,y){x*y})
  }
  
  m <- nrow(K)
  H <- diag(m) - (1/m) * matrix(1, nrow=m, ncol=m)
  
  if ("cor" %in% result){
    corobs <- cor(as.numeric(H %*% K %*% H), 
                  as.numeric(H %*% K2 %*% H))
  }else{
    corobs <- NA
  }
  
  if ("test" %in% result){
    kdcobs <- (1/(m^2))*sum(K * H %*% K2 %*% H)
    #Do permutation
    kdcstore <- numeric(num.perm)
    for (i in 1:num.perm){
      xvec.curr <- sample(xvec)
      if (Xdist == "categorical"){
        d2 <- outer(xvec.curr, xvec.curr, 
                    FUN = function(x,y){1 - as.numeric(x == y)})
        K2 <- exp(-d2)
      }
      
      if (Xdist == "euclidean"){
        d2 <- outer(xvec.curr, xvec.curr, 
                    FUN = function(x,y){abs(x-y)})
        K2 <- exp(-d2/median(as.numeric(d2)))
      }
      
      if (Xdist == "LK"){
        K2 <- outer(xvec.curr, xvec.curr, 
                    FUN = function(x,y){x*y})
      }
      
      kdcstore[i] <- (1/(m^2))*sum(K * H %*% K2 %*% H)
    }
    pval <- (sum(kdcstore[-1] >= kdcobs)+1)/(num.perm)
  }else{
    kdcobs <- pval <- NA
  }
  
  
  return(list(kdcstat = kdcobs, pval = pval, dcor = corobs))
}

#kdcstat: KDC statistics, pval: P-value for the test, dcor: value of dcor







######### Function to calculate LOO dcor #########

LOO.comp <- function(propmat, xvec, Xdist = "categorical"){
  
  # function arguments as before
  
  out <- dcor.comp(propmat, xvec, Xdist = Xdist, result = "cor")
  
  dcor.loo <- numeric(ncol(propmat))
  for (i in 1:ncol(propmat)){
    propmat.loo <- propmat[, -i]
    propmat.loo <- propmat.loo/apply(propmat.loo, 1, sum)
    dcor.loo[i] <- dcor.comp(propmat.loo, xvec, Xdist = Xdist, result = "cor")$dcor
  }
  
  return(list(dcor = out$dcor, loo = dcor.loo))
}






######## Functions for weighted dcor ########

## Function to calculate weighted Aitchison distance
aitchison.w <- function(propmat, wt){
  
  # wt: weights, other arguments as before
  y <- t(t(propmat)/wt)
  sw <- sum(wt)
  gp <- apply(y, 1, function(x){
    exp(sum(wt*log(x))/sw)
  })
  
  fn <- function(i, j){
    sqrt(sum(wt*(log(y[i,]/gp[i])-log(y[j,]/gp[j]))^2))
  }
  dmat <- outer(1:nrow(y), 1:nrow(y), Vectorize(fn))
  return(dmat)
}




## This function calculates dcor for a given weight, to be passed into the optimization
calc.dcor.w <- function(g, wt, propmat, xvec, 
                        Xdist = "categorical", output = "cor"){
  
  #g: gamma, other arguments as before
  
  wt <- wt^g/sqrt(sum(wt^(2*g)))
  d <- aitchison.w(propmat, wt)
  alldist <- as.numeric(d)
  K <- exp(-as.matrix(d)/median(alldist))
  
  if (Xdist == "categorical"){
    d2 <- outer(xvec, xvec, 
                FUN = function(x,y){1 - as.numeric(x == y)})
    K2 <- exp(-d2)
  }
  
  if (Xdist == "euclidean"){
    d2 <- outer(xvec, xvec, 
                FUN = function(x,y){abs(x-y)})
    K2 <- exp(-d2/median(as.numeric(d2)))
  }
  
  if (Xdist == "LK"){
    K2 <- outer(xvec, xvec, 
                FUN = function(x,y){x*y})
  }
  
  m <- nrow(K)
  H <- diag(m) - (1/m) * matrix(1, nrow=m, ncol=m)
  
  if (output == "cor"){
    corobs <- cor(as.numeric(H %*% K %*% H), 
                  as.numeric(H %*% K2 %*% H))
    return(corobs)
  }
  
  if (output == "cov"){
    kdcstat <- (1/(m^2))*sum(K * H %*% K2 %*% H)
    return(kdcstat)
  }
}




#### Function to find optimal weights ####
find.opt.wt <- function(wt.in, propmat, xvec, Xdist = "categorical"){
  
  #wt.in: Initial weights, other arguments as before
  
  out <- optim(par = 1, fn = calc.dcor.w,
               wt = wt.in,
               propmat = propmat, xvec = xvec, Xdist = Xdist,
               control=list(fnscale=-1),
               method = "Brent", lower = 0, upper = 100)
  return(out$par)
}





#### Optimized weighted dcor ####
dcor.comp.w <- function(propmat, xvec, num.perm = 1e4, 
                        Xdist = "categorical", result = "test",
                        wt.in = NULL){
  
  # function arguments as before
  
  if (is.null(wt.in)){
    wt.in <- numeric(ncol(propmat))
    for (r in 1:ncol(propmat)){
      wt.in[r] <- calc.dcor.w(g = 1, rep(1, 2), 
                              cbind(propmat[,r],1-propmat[,r]), xvec, Xdist)
    }
  }
  
  g <- find.opt.wt(wt.in, propmat, xvec, Xdist)
  
  if ("cor" %in% result){
    corobs <- calc.dcor.w(g, wt.in, propmat, xvec, Xdist, output = "cor")
  }else{
    corobs <- NA
  }
  
  if ("test" %in% result){
    kdcobs <- calc.dcor.w(g, wt.in, propmat, xvec, Xdist, output = "cov")
    #Do permutation
    kdcstore <- numeric(num.perm)
    for (i in 1:num.perm){
      if(floor(i/100) == i/100){print(i)}
      xvec.curr <- sample(xvec)
      g.curr <- find.opt.wt(wt.in, propmat, xvec.curr, Xdist)
      
      kdcstore[i] <- calc.dcor.w(g.curr, wt.in, propmat, xvec.curr, Xdist, output = "cov")
    }
    pval <- (sum(kdcstore[-1] >= kdcobs)+1)/(num.perm)
  }else{
    kdcobs <- pval <- NA
  }
  
  wt.out <- wt.in^g/sqrt(sum(wt.in^(2*g)))
  
  return(list(kdcstat = kdcobs, pval = pval, dcor = corobs, wt.out = wt.out))
}







######### Functions for covariate adjustments when using dcor test #########

clr <- function(x){
  clrs <- log(x/exp(mean(log(x))))
  return(clrs)
}

invclr <- function(clrs){
  x <- exp(clrs)
  return(x/sum(x))
}

alr <- function(x, denomind = NA){
  if(is.na(denomind)){
    denomind <- length(x)
  }
  denom <- x[denomind]
  alrs <- log(x/denom)[-denomind]
  return(alrs)
}

invalr <- function(alrs, denomind = NA){
  if(is.na(denomind)){
    denomind <- length(alrs)+1
  }
  x <- numeric(length(alrs)+1)
  x[-denomind] <- exp(alrs)
  x[denomind] <- 1
  return(x/sum(x))
}







#### Function for computing dcor while adjusting for categorical covariate ####
# This method uses the stratified kernel method only
# See simulation codes for the alr/clr approach

dcor.comp.cov <- function(propmat, xvec, zmat = NA, num.perm = 1e4, 
                          Xdist = "categorical",
                          covadj = "SK",
                          result = "test"){
  
  # zmat: an n x 1 vector having the covariate values
  # covadj: Must be "SK"
  # other arguments as before
  
  d <- dist(propmat, method = 'aitchison')
  
  alldist <- as.numeric(d)
  K <- exp(-as.matrix(d)/median(alldist))
  
  if (Xdist == "categorical"){
    d2 <- outer(xvec, xvec, 
                FUN = function(x,y){1 - as.numeric(x == y)})
    K2 <- exp(-d2)
  }
  
  if (Xdist == "Gaussian"){
    d2 <- outer(xvec, xvec, 
                FUN = function(x,y){abs(x-y)})
    K2 <- exp(-d2/median(as.numeric(d2)))
  }
  
  if (Xdist == "LK"){
    K2 <- outer(xvec, xvec, 
                FUN = function(x,y){x*y})
  }
  
  m <- nrow(K)
  H <- diag(m) - (1/m) * matrix(1, nrow=m, ncol=m)
  Kc <- H %*% K %*% H
  K2c <- H %*% K2 %*% H
  
  
  #Finding the rows having same z
  n <- nrow(zmat)
  issame <- function(i,j,data) {sum(data[i,] - data[j,]) == 0}
  issamevec <- Vectorize(issame, vectorize.args=list("i","j"))
  
  keep <- outer(1:n, 1:n, issamevec, data = zmat)
  drop <- keep == F
  
  if (covadj == "SK"){
    Kc[drop] <- 0
    K2c[drop] <- 0
    
    xx <- as.numeric(Kc[row(Kc) != col(Kc)])
    yy <- as.numeric(K2c[row(K2c) != col(K2c)])
    
    if ("cor" %in% result){
      corobs <- cor(xx, yy)
    }else{
      corobs <- NA
    }
  }

  
  if ("test" %in% result){
    m <- length(xx)
    kdcobs <- cov(xx, yy)*((m-1)/m)
    
    #Do permutation
    kdcstore <- numeric(num.perm)
    for (i in 1:num.perm){
      xvec.curr <- sample(xvec)
      if (Xdist == "categorical"){
        d2 <- outer(xvec.curr, xvec.curr, 
                    FUN = function(x,y){1 - as.numeric(x == y)})
        K2 <- exp(-d2)
      }
      
      if (Xdist == "euclidean"){
        d2 <- outer(xvec.curr, xvec.curr, 
                    FUN = function(x,y){abs(x-y)})
        K2 <- exp(-d2/median(as.numeric(d2)))
      }
      
      if (Xdist == "LK"){
        K2 <- outer(xvec.curr, xvec.curr, 
                    FUN = function(x,y){x*y})
      }
      
      if (covadj == "SK"){
        K2c <- H %*% K2 %*% H
        
        K2c[drop] <- 0
        yy <- as.numeric(K2c[row(K2c) != col(K2c)])
      }
      
      kdcstore[i] <- cov(xx, yy)*((m-1)/m)
    }
    pval <- (sum(kdcstore[-1] >= kdcobs)+1)/(num.perm)
  }else{
    kdcobs <- pval <- NA
  }
  
  
  return(list(kdcstat = kdcobs, pval = pval, dcor = corobs))
}

###########################################################################







###### Doing the distance correlation test for the new data ######

#hold is an expression containing a which expression with 
  #the values of Group to hold
#xexp is an expression containing the name of the predictor

cellfreq.test.K <- function(propmat, infotable, hold, xexp, num.perm = 1e4){
  propmat.curr <- propmat[eval(hold), ]
  infotable.curr <- infotable[eval(hold), ]
  xvec <- eval(xexp)
  out <- dcor.comp(propmat.curr, xvec, num.perm, result = c("cor", "test"))
  return(out)
}

#Function with similar arguments, to be used to obtain the LOO dcor values
cellfreq.K.LOO <- function(propmat, infotable, hold, xexp){
  propmat.curr <- propmat[eval(hold), ]
  infotable.curr <- infotable[eval(hold), ]
  xvec <- eval(xexp)
  out <- dcor.comp(propmat.curr, xvec, result = "cor")
  
  dcor.loo <- numeric(ncol(propmat.curr))
  for (i in 1:ncol(propmat.curr)){
    propmat.loo <- propmat.curr[, -i]
    propmat.loo <- propmat.loo/apply(propmat.loo, 1, sum)
    #set.seed(myseed)
    dcor.loo[i] <- dcor.comp(propmat.loo, xvec, result = "cor")$dcor
  }
  
  dcor.loo <- c(out$dcor, dcor.loo)
  names(dcor.loo) <- c("Full", colnames(propmat.curr))
  return(dcor.loo)
}



#############################################################################




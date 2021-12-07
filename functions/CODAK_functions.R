
############### Functions for cell type abundance comparison ################


library("MASS")
library("coda.base")
library("lme4")
library("vegan")
library("MiRKAT")


#### Function for the KDC test ####


dcor.comp <- function(propmat, xvec, num.perm = 1e4, 
                      Xdist = "categorical", compdist = "aitchison",
                      gaussian = F, result = "test"){
  
  # propmat: n x num.cat matrix of proportions (n = sample size, num.cat = number of cell types)
  # xvec: an n x 1 vector containing the values of the predictor
  # num.perm: number of permutations
  # Xdist: specifying the kernel used for the predictor. 
  #       "categorical" - Hamming distance kernel, "LK" - Linear kernel, "gaussian" - Gaussian kernel
  # compdist: specifying the kernel used for the compositions. 
  #       "aitchison" - Aitchison distance,"euclidean" - Euclidean distance, 
  #       "BCD" - Bray Curtis distance
  # gaussian: T/F, whether to put the distance into a Gaussian Kernel (if FALSE, use D2K like mirkat)
  #         Don't use gaussian = T for compdist = "BCD"
  # result: a vector, "test" - provides test results, "cor" - provides value of dcor
  
  if (compdist == "BCD"){
    d <- vegdist(propmat, method="bray")
  }else{
    d <- dist(propmat, method = compdist)
  }
  
  if (gaussian == F){
    K <- D2K(as.matrix(d))
  }else{
    if (compdist == "BCD"){
      print("Gaussian kernel not recommended for BCD, using D2K instead")
      K <- D2K(as.matrix(d))
    }else{
      alldist <- as.numeric(d)
      K <- exp(-as.matrix(d^2)/median(alldist)^2)
    }
  }
  
  if (Xdist == "categorical"){
    d2 <- outer(xvec, xvec, 
                FUN = function(x,y){1 - as.numeric(x == y)})
    K2 <- D2K(as.matrix(d2))
  }
  
  if (Xdist == "gaussian"){
    d2 <- outer(xvec, xvec, 
                FUN = function(x,y){(x-y)^2})
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
    #Note that this is same as (1/(m^2))*sum(diag(K %*% H %*% K2 %*% H))
    
    #Do permutation
    kdcstore <- numeric(num.perm)
    
    for (i in 1:num.perm){
      ind <- sample(m)
      K.curr <- K[ind, ind]
      
      kdcstore[i] <- (1/(m^2))*sum(K.curr * H %*% K2 %*% H)
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







#### Function for computing dcor while adjusting for covariates ####

dcor.comp.cov <- function(propmat, xvec, zmat = NA, num.perm = 1e4, 
                          Xdist = "categorical",
                          covadj = "SK",
                          result = "test",
                          permutation = "shuffleRes"){
  
  # zmat: an n x 1 vector having the covariate values
  # covadj: Must be "SK", "alr" or "clr"
  # permutation: Options are "shuffleY", "shuffleRes" and "FL" (Freedman-Lane)
  # other arguments as before
  
  
  
  if (covadj == "SK"){
    
    xvec <- xvec[order(zmat)]
    propmat <- propmat[order(zmat), ]
    zmat <- as.matrix(zmat[order(zmat), ])
    tabz <- table(zmat)
    
    d <- dist(propmat, method = 'aitchison') #we assume compdist = "aitchison"
    
    K <- D2K(as.matrix(d)) #we assume gaussian = F
    
    if (Xdist == "categorical"){
      d2 <- outer(xvec, xvec, 
                  FUN = function(x,y){1 - as.numeric(x == y)})
      K2 <- exp(-d2)
    }
    
    if (Xdist == "gaussian"){
      d2 <- outer(xvec, xvec, 
                  FUN = function(x,y){(x-y)^2})
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
        
        if (length(tabz) == 2){
          ind <- c(sample(tabz[1]), tabz[1] + sample(tabz[2])) # works only for 2 grps
        }else{
          ind <- sample(tabz[1])
          for (l in 2:length(tabz)){
            ind <- c(ind, sum(tabz[1:(l-1)] + sample(tabz[l])))
          }
        }
        
        
        K.curr <- K[ind, ind]
        
        if (covadj == "SK"){
          Kc <- H %*% K.curr %*% H
          
          Kc[drop] <- 0
          xx <- as.numeric(Kc[row(Kc) != col(Kc)])
        }
        
        kdcstore[i] <- cov(xx, yy)*((m-1)/m)
      }
      pval <- (sum(kdcstore[-1] >= kdcobs)+1)/(num.perm)
    }else{
      kdcobs <- pval <- NA
    }
  }
  
  
  
  
  if (covadj == "alr" | covadj == "clr"){
    z <- zmat
    
    m <- length(z)
    zm <- cbind(rep(1, length(z)), z)
    Pz <- zm%*%solve(t(zm)%*%zm)%*%t(zm)
    resx <- as.numeric((diag(m) - Pz)%*%x)
    
    if (covadj == "alr"){
      transprop <- t(apply(propmat, 1, alr))
      
      resp <- (diag(m) - Pz)%*%transprop
      
      respropmat <- t(apply(resp, 1, invalr))
    }else{
      transprop <- t(apply(propmat, 1, clr))
      
      resp <- (diag(m) - Pz)%*%transprop
      
      respropmat <- t(apply(resp, 1, invclr))
    }

    
    if (permutation == "shuffleRes"){
      out <- dcor.comp(respropmat, resx, num.perm, 
                       Xdist = "LK", result = result)
      pval <- out$pval
      corobs <- out$dcor
      kdcobs <- out$kdcstat
    }
    
    if (permutation == "shuffleY"){
      outobs <- dcor.comp(respropmat, resx, num.perm = 1, 
                          Xdist = "LK", result = c("test", "cor"))
      corobs <- outobs$dcor
      kdcobs <- outobs$kdcstat
      
      corstore <- numeric(num.perm)
      
      for (i in 1:num.perm){
        propmat.curr <- propmat[sample(m),]
        
        if (covadj == "alr"){
          transprop1 <- t(apply(propmat.curr, 1, alr))
          
          resp1 <- (diag(m) - Pz)%*%transprop1
          
          respropmat <- t(apply(resp1, 1, invalr))
        }else{
          transprop1 <- t(apply(propmat.curr, 1, clr))
          
          resp1 <- (diag(m) - Pz)%*%transprop1
          
          respropmat <- t(apply(resp1, 1, invclr))
        }
        
        corstore[i] <- dcor.comp(respropmat, resx, num.perm = 1, 
                                 Xdist = "LK", result = c("cor"))$dcor
      }
      
      pval <- (sum(corstore[-1] >= corobs)+1)/(num.perm)
    }
    
    if (permutation == "FL"){
      outobs <- dcor.comp(respropmat, resx, num.perm = 1, 
                          Xdist = "LK", result = c("test", "cor"))
      corobs <- outobs$dcor
      kdcobs <- outobs$kdcstat
      
      corstore <- numeric(num.perm)
      
      for (i in 1:num.perm){
        resp.curr <- resp[sample(m),]
        
        if (covadj == "alr"){
          transprop1 <- resp.curr #+ Pz%*%transprop #This is not needed, provides same result
          
          resp1 <- (diag(m) - Pz)%*%transprop1
          
          respropmat <- t(apply(resp1, 1, invalr))
        }else{
          transprop1 <- resp.curr #+ Pz%*%transprop
          
          resp1 <- (diag(m) - Pz)%*%transprop1
          
          respropmat <- t(apply(resp1, 1, invclr))
        }
        
        corstore[i] <- dcor.comp(respropmat, resx, num.perm = 1, 
                                 Xdist = "LK", result = c("cor"))$dcor
      }
      
      pval <- (sum(corstore[-1] >= corobs)+1)/(num.perm)
    }
  }
  
  
  return(list(kdcstat = kdcobs, pval = pval, dcor = corobs))
}



###########################################################################







###### Doing the distance correlation test for the real data ######

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
    
    dcor.loo[i] <- dcor.comp(propmat.loo, xvec, result = "cor")$dcor
  }
  
  dcor.loo <- c(out$dcor, dcor.loo)
  names(dcor.loo) <- c("Full", colnames(propmat.curr))
  return(dcor.loo)
}



#############################################################################




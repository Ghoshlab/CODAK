

###### Functions for simulating data for cell type abundance comparison #######

library("MASS")
library("coda.base")
library("lme4")
library("SummarizedExperiment")
library("diffcyt")





##### Function for adding effects to cell type compositions #####

add.effects <- function(prop, maxeffect, ind.maxeff, 
                        n.nochange = 0, ind.nochange = NULL){
  
  # prop: proportions for the reference group (for x = 0 when x is continuous)
  # maxeffect: maximum effectsize in the prop scale
  # ind.maxeff: index of the celltype for which maximum effectsize to be added
  # n.nochange: how many celltypes don't have different abundance
  # ind.nochange: indices of the nochange cell types, a vector of length n.nochange
  
  n <- length(prop)
  addeff <- rep(NA, n)
  sum.nochange <- 0
  
  if (n.nochange > 0){
    addeff[ind.nochange] <- 0
    sum.nochange <-sum(prop[ind.nochange], na.rm = T)
  }
  
  addeff[ind.maxeff] <- maxeffect
  
  n2 <- n - 1 - n.nochange
  vec <- runif(n2, -abs(maxeffect) + maxeffect/n2, abs(maxeffect) - maxeffect/n2)
  temp <- which(sign(vec) == sign(sum(vec)))
  vec[temp] <- vec[temp] - sum(vec)/length(temp)
  vec <- vec - maxeffect/n2
  
  adj <- is.na(addeff)
  addeff[adj] <- vec
  
  prop.new <- prop + addeff
  
  if (sum(prop.new <= 0) > 0){
    temp <- which(prop.new <= 0)
    buffer <- 0.0001*length(temp)
    sumtoadd <- -sum(prop.new[temp])+ buffer
    temp2 <- which(vec > 0)
    vec[temp2] <- vec[temp2] - sumtoadd/length(temp2)
    
    addeff[adj] <- vec
    prop.new <- prop + addeff
    
    prop.new[temp] <- buffer/length(temp)
  }
  
  return(prop.new)
}








############### Function to obtain effect sizes from one given proportion vector using multilogit ###################

#Getting the effect size betas in the mlogit scale starting from the prop scale

getES.mlogit <- function(prop, ref = "last"){
  
  # ref: reference cell type
  
  n <-length(prop)
  if (ref == "last"){ref <- n}
  
  beta <- numeric(n)
  beta[ref] <- 0
  
  beta[-ref] <- log(prop[-ref]/prop[ref])
  
  return(beta)
}




############### Function to obtain proportions from effect sizes using multilogit ###################

getProps.mlogit <- function(beta, betaX = NULL, X = NULL){
  
  # beta: vector of coefficients for categories
  # betaX: matrix of coefficients for the predictors:
  #       rows are categories, columns are different predictors, 
  # X is a column vector of the predictor values, num.pred x 1
  # Note that although we use one predictor of interest, 
  #       a covariate is essentially another predictor
  
  n <- length(beta) #number of categories
  
  if(is.null(betaX) == 1){
    betaX <- matrix(0, nrow = n,ncol = 1)
    X <- 0
  }
  
  ref <- which(beta == 0)
  
  prop <- numeric(n)
  
  prop[ref] <- 1
  prop[-ref] <- exp(as.numeric(beta[-ref] + as.matrix(betaX[-ref,])%*%X))*prop[ref]
  
  prop <- prop/sum(prop)
  return(prop)
}






######## Function to simulate proportions based on effect sizes #########

simProps <- function(beta, betaX, Xmat, nvec, sigb = NULL){
  
  # rows of Xmat denote sample units, columns denote different X variables
  # nvec: a vector containing total number of cells per file
  # sigb: covariance matrix of the random effects (num.cat x num.cat)
  
  num.cat <- length(beta)
  y <- matrix(NA, nrow = nrow(Xmat), ncol = num.cat)
  
  if(is.null(sigb)){
    for (i in 1:nrow(Xmat)){
      probs <- getProps.mlogit(beta, betaX, as.matrix(Xmat[i, ]))
      y[i, ] <- t(rmultinom(n = 1, size = nvec[i], prob = probs))
    }
  }else{
    bmat <- matrix(0, nrow = nrow(Xmat), ncol = num.cat)
    ref <- which(beta == 0)
    bmat[, ref] <- rep(0, nrow(Xmat))
    bmat[, -ref] <- mvrnorm(nrow(Xmat), mu = rep(0, num.cat-1), sigb)
    for (i in 1:nrow(Xmat)){
      probs <- getProps.mlogit(beta + bmat[i, ], betaX, as.matrix(Xmat[i, ]))
      y[i, ] <- t(rmultinom(n = 1, size = nvec[i], prob = probs))
    }
  }
  
  return(y)
}

##########################################################################











########## Simulations Study codes ############


#### Fixed effects mlogit, not practical

simulation_study_1 <- function(beta, betaX, Xmat, nvec, num.sim, 
                               num.perm = 10000, printafter = 50,
                               Xdist = "categorical",
                               storeY = FALSE){
  #No random effect at all
  #rows of Xmat denote sample units, columns denote different X variables
  #storeY: stores the proportion matrix if TRUE
  
  pval.kernel <- pval.logistic <- numeric(num.sim)
  
  if (storeY == T){
    ylist <- list()
  }
  
  for (i in 1:num.sim){
    if (floor(i/printafter) == i/printafter){print(i)}
    
    y <- simProps(beta, betaX, Xmat, nvec)
    
    ng <- nrow(y)/2
    num.cat <- ncol(y)
    
    propmat <- y/rep(nvec, num.cat)
    x <- as.vector(Xmat)
    
    ptemp <- rep(NA, num.cat) 
    for (k in 1:num.cat){ 
      y.bin <- y[, k]
      fit <- glm(y.bin/nvec ~ x, weights = nvec, 
                 family = binomial(link = "logit"))
      ptemp[k] <- summary(fit)$coefficients[2, 4]
    }
    
    pval.logistic[i] <- min(p.adjust(ptemp, method = "bonferroni"))
    
    out <- dcor.comp(propmat, x, num.perm, Xdist = Xdist)
    pval.kernel[i] <- out$pval
    
    if (storeY == T){
      ylist[[i]] <- y
    }
  }
  
  if (storeY == T){
    retlist <-  list(logistic = pval.logistic, kernel = pval.kernel, ylist = ylist)
  }else{
    retlist <-  list(logistic = pval.logistic, kernel = pval.kernel)
  }
  
  return(retlist)
}

###############################################################################






# Mixed effects mlogit: OLRE added, realistic data

simulation_study_2 <- function(beta, betaX, Xmat, nvec, sigb, 
                               num.sim, num.perm = 10000, printafter = 100, 
                               Xdist = "categorical",
                               whichtests = c("Kernel", "LogisticM"),
                               adjmethod = "bonferroni",
                               lmertest = c("Wald", "LRT"), lmer.nump = 100, 
                               diffcyt.method = c("edgeR", "voom"),
                               storeY = FALSE,
                               standardize = TRUE){
  
  # rows of Xmat denote sample units, columns denote different X variables
  # sigb is the covariance matrix of the random effects (num.cat x num.cat)
  # whichtests: which tests to conduct
  # adjmethod: method to use for multiple testing adjustment
  # lmertest: Which logistic mixed models to use
  # lmer.nump: number of permutations for the LRT-permutation method
  # diffcyt.method: which diffcyt methods to use
  # standardize: standardizes continuous predictor values if TRUE
  
  pval.logistic <- pval.kernel <- pval.logisticM <- pval.diffcyt <- NA
  
  if ("Kernel" %in% whichtests){
    pval.kernel <- numeric(num.sim)
  }
  
  if ("Logistic" %in% whichtests){
    pval.logistic <- numeric(num.sim)
  }
  
  if ("LogisticM" %in% whichtests){
    pval.logisticM <- matrix(NA, nrow = num.sim, ncol = length(lmertest))
    colnames(pval.logisticM) <- lmertest
  }
  
  if ("diffcyt" %in% whichtests){
    pval.diffcyt <- matrix(NA, nrow = num.sim, ncol = length(diffcyt.method))
    colnames(pval.diffcyt) <- diffcyt.method
  }
  
  if (storeY == T){
    ylist <- list()
  }else{
    ylist = NA
  }
  
  
  
  #### Simulation begins ####
  
  for (i in 1:num.sim){
    if (floor(i/printafter) == i/printafter){print(i)}
    
    y <- simProps(beta, betaX, Xmat, nvec, sigb)
    
    ng <- nrow(y)/2
    num.cat <- ncol(y)
    
    propmat <- y/rep(nvec, num.cat)
    
    #x is the predictor of interest. Works for only one predictor.
    x <- as.vector(Xmat)
    
    
    #### Scaling for continuous predictors ####
    
    if (standardize == T & Xdist %in% c("Gaussian", "LK")){
      x <- (x - mean(x))/sd(x)
    }
    
    
    #### dcor test ####
    
    if ("Kernel" %in% whichtests){
      out <- dcor.comp(propmat, x, num.perm, Xdist = Xdist)
      pval.kernel[i] <- out$pval
    }  
    
    
    
    #### diffcyt tests ####
    if ("diffcyt" %in% whichtests){
      for (k in 1:length(diffcyt.method)){
        curr.method <- diffcyt.method[k]
        temp <- diffcyt.test(y, x, method = curr.method)
        pval.diffcyt[i, k] <- min(p.adjust(temp, method = adjmethod))
      }
    }
    
    
    
    
    #### Logistic or Logistic Mixed models ####
    
    if ("LogisticM" %in% whichtests | "Logistic" %in% whichtests){
      ptemp <- rep(NA, num.cat) 
      ptemp.mixed <- matrix(NA, nrow = num.cat, ncol = length(lmertest))
      
      for (k in 1:num.cat){ 
        y.bin <- y[, k]
        fit <- glm(y.bin/nvec ~ x, weights = nvec, 
                   family = binomial(link = "logit"))
        
        if ("Logistic" %in% whichtests){
          ptemp[k] <- summary(fit)$coefficients[2, 4]
        }  
        
        
        if ("LogisticM" %in% whichtests){
          id <- 1:length(x)
          fitmixed <- try(glmer(y.bin/nvec ~ x + (1|id), weights = nvec, 
                            family = binomial(link = "logit")))
          if (class(fitmixed) == "try-error"){
            ptemp.mixed[k, colind]<- NA
          }else{
            if ("Wald" %in% lmertest){
              colind <- which(lmertest == "Wald")
              ptemp.mixed[k, colind] <- summary(fitmixed)$coefficients[2, 4]
            }
            
            if ("LRT" %in% lmertest){
              colind <- which(lmertest == "LRT")
              temptry <- try(drop1(fitmixed, test="Chisq")$"Pr(Chi)"[2])
              if (class(fitmixed) == "try-error"){
                ptemp.mixed[k, colind]<- NA
              }else{
                ptemp.mixed[k, colind] <- temptry
              }
            }
            
            if ("Permutation" %in% lmertest){
              colind <- which(lmertest == "Permutation")
              lrstat <- drop1(fitmixed, test="Chisq")$LRT[2]
              lrstat.perm <- numeric(lmer.nump)
              for (p in 1:lmer.nump){
                xperm <- sample(x)
                fitmixed <- glmer(y.bin/nvec ~ xperm + (1|id), weights = nvec, 
                                  family = binomial(link = "logit"))
                lrstat.perm[p] <- drop1(fitmixed, test="Chisq")$LRT[2]
              }
              ptemp.mixed[k, colind] <- (sum(lrstat.perm >= lrstat)+1)/(lmer.nump+1)
            }#perm test loop ends
          }#try loop ends
        }#logisticM loop ends  
      }#categories loop ends
      
      if ("Logistic" %in% whichtests){
        pval.logistic[i] <- min(p.adjust(ptemp, method = adjmethod))
      } 
      
      if ("LogisticM" %in% whichtests){
        pval.logisticM[i, ] <- apply(ptemp.mixed, 2, function(w){
          return(min(p.adjust(w, method = adjmethod)))
        })
      } 
      
    }#logistic loop ends
    
    
    
    #### Storing the simulated data ####
    
    if (storeY == T){
      ylist[[i]] <- y
    }
  }#sim loop ends  
  
  
  retlist <-  list(logistic = pval.logistic, kernel = pval.kernel, 
                   logistic.mixed = pval.logisticM,
                   diffcyt = pval.diffcyt,
                   ylist = ylist)
  return(retlist)
}

##############################################################################








# Mixed effects mlogit: OLRE added, realistic data + Covariate

simulation_study_3 <- function(beta, betaX, Xmat, nvec, sigb, 
                               betaZ, Zmat, 
                               num.sim, num.perm = 10000, printafter = 100, 
                               Xdist = "categorical", covadjs = c("SK", "alr"),
                               whichtests = c("Kernel", "LogisticM"),
                               adjmethod = "bonferroni",
                               lmertest = c("Wald", "LRT"), lmer.nump = 100, 
                               diffcyt.method = c("edgeR", "voom"),
                               storeY = FALSE){
  
  #rows of Xmat denote sample units, columns denote different X variables
  #sigb is the covariance matrix of the random effects (num.cat x num.cat)
  
  pval.logistic <- pval.kernel <- pval.logisticM <- pval.diffcyt <- NA
  
  if ("Kernel" %in% whichtests){
    pval.kernel <- matrix(NA, nrow = num.sim, ncol = length(covadjs))
    colnames(pval.kernel) <- covadjs
  }
  
  
  if ("Logistic" %in% whichtests){
    pval.logistic <- numeric(num.sim)
  }
  
  if ("LogisticM" %in% whichtests){
    pval.logisticM <- matrix(NA, nrow = num.sim, ncol = length(lmertest))
    colnames(pval.logisticM) <- lmertest
  }
  
  if ("diffcyt" %in% whichtests){
    pval.diffcyt <- matrix(NA, nrow = num.sim, ncol = length(diffcyt.method))
    colnames(pval.diffcyt) <- diffcyt.method
  }
  
  if (storeY == T){
    ylist <- list()
  }else{
    ylist = NA
  }
  
  
  
  #### Simulation begins ####
  
  for (i in 1:num.sim){
    if (floor(i/printafter) == i/printafter){print(i)}
    
    y <- simProps(beta, cbind(betaX, betaZ), cbind(Xmat, Zmat), nvec, sigb)
    
    ng <- nrow(y)/2
    num.cat <- ncol(y)
    
    propmat <- y/rep(nvec, num.cat)
    
    #x is the predictor of interest. 
    #z is the covariate, works for only one covariate.
    x <- as.vector(Xmat)
    z <- as.vector(Zmat)
    
    #### Scaling for continuous predictors ####
    
    if (Xdist == "Gaussian" | Xdist == "LK"){
      x <- (x - mean(x))/sd(x)
    }
    
    
    #### dcor test ####
    
    if ("Kernel" %in% whichtests){
      for (j in 1:length(covadjs)){
        covadj.curr <- covadjs[j]
        
        if (covadj.curr %in% c("SK", "SK2", "SK3")){
          out <- try(dcor.comp.cov(propmat, x, Zmat, num.perm, 
                               Xdist = Xdist, covadj = covadj.curr), silent = T)
          if (class(out) == "try-error"){
            pval.kernel[i, j] <-NA
          }else{
            pval.kernel[i, j] <- out$pval
          }
        }
        
        if (covadj.curr == "alr"){
          alrs <- t(apply(propmat, 1, alr))
          fit1 <- try(lm(alrs ~ z), silent = T)
          if (class(fit1)[1] == "try-error"){
            pval.kernel[i, j] <-NA
          }else{
            fit2 <- lm(x ~ z)
            respropmat <- t(apply(residuals(fit1), 1, invalr))
            resx <- residuals(fit2)
            out <- dcor.comp(respropmat, resx, num.perm, 
                             Xdist = "LK")
            pval.kernel[i, j] <- out$pval
          }
        }
        
        
        if (covadj.curr == "clr"){
          clrs <- t(apply(propmat, 1, clr))
          fit1 <- try(lm(clrs ~ z), silent = T)
          if (class(fit1)[1] == "try-error"){
            pval.kernel[i, j] <-NA
          }else{
            fit2 <- lm(x ~ z)
            respropmat <- t(apply(residuals(fit1), 1, invclr))
            resx <- residuals(fit2)
            out <- dcor.comp(respropmat, resx, num.perm, 
                             Xdist = "LK")
            pval.kernel[i, j] <- out$pval
          }
        }
      }#End covadj loop
    }#End kernel loop
    
    
    #### diffcyt tests ####
    if ("diffcyt" %in% whichtests){
      for (k in 1:length(diffcyt.method)){
        curr.method <- diffcyt.method[k]
        temp <- diffcyt.test.cov(y, x, z, method = curr.method)
        pval.diffcyt[i, k] <- min(p.adjust(temp, method = adjmethod))
      }
    }
    
    
    
    
    
    #### Logistic or Logistic Mixed models ####
    
    if ("LogisticM" %in% whichtests | "Logistic" %in% whichtests){
      ptemp <- rep(NA, num.cat) 
      ptemp.mixed <- matrix(NA, nrow = num.cat, ncol = length(lmertest))
      
      for (k in 1:num.cat){ 
        y.bin <- y[, k]
        fit <- glm(y.bin/nvec ~ x + z, weights = nvec, 
                   family = binomial(link = "logit"))
        
        if ("Logistic" %in% whichtests){
          ptemp[k] <- summary(fit)$coefficients[2, 4]
        }  
        
        
        if ("LogisticM" %in% whichtests){
          id <- 1:length(x)
          fitmixed <- try(glmer(y.bin/nvec ~ x + z + (1|id), weights = nvec, 
                                family = binomial(link = "logit")))
          if (class(fitmixed) == "try-error"){
            ptemp.mixed[k, colind]<- NA
          }else{
            if ("Wald" %in% lmertest){
              colind <- which(lmertest == "Wald")
              ptemp.mixed[k, colind] <- summary(fitmixed)$coefficients[2, 4]
            }
            
            if ("LRT" %in% lmertest){
              colind <- which(lmertest == "LRT")
              temptry <- try(drop1(fitmixed, test="Chisq")$"Pr(Chi)"[2])
              if (class(temptry) == "try-error"){
                ptemp.mixed[k, colind]<- NA
              }else{
                ptemp.mixed[k, colind] <- temptry
              }
            }
            
            if ("Permutation" %in% lmertest){
              colind <- which(lmertest == "Permutation")
              lrstat <- drop1(fitmixed, test="Chisq")$LRT[2]
              lrstat.perm <- numeric(lmer.nump)
              for (p in 1:lmer.nump){
                xperm <- sample(x)
                fitmixed <- glmer(y.bin/nvec ~ xperm + z + (1|id), weights = nvec, 
                                  family = binomial(link = "logit"))
                lrstat.perm[p] <- drop1(fitmixed, test="Chisq")$LRT[2]
              }
              ptemp.mixed[k, colind] <- (sum(lrstat.perm >= lrstat)+1)/(lmer.nump+1)
            }#perm test loop ends
          }#try loop ends
        }#logisticM loop ends  
      }#categories loop ends
      
      if ("Logistic" %in% whichtests){
        pval.logistic[i] <- min(p.adjust(ptemp, method = adjmethod))
      } 
      
      if ("LogisticM" %in% whichtests){
        pval.logisticM[i, ] <- apply(ptemp.mixed, 2, function(w){
          return(min(p.adjust(w, method = adjmethod)))
        })
      } 
      
    }#logistic loop ends
    
    
    
    #### Storing the simulated data ####
    
    if (storeY == T){
      ylist[[i]] <- y
    }
  }#sim loop ends  
  
  
  retlist <-  list(logistic = pval.logistic, kernel = pval.kernel, 
                   logistic.mixed = pval.logisticM,
                   diffcyt = pval.diffcyt,
                   ylist = ylist)
  return(retlist)
}


############################################################################







# Simulation 3 function + potential dependence in the data 

simulation_study_4 <- function(beta, betaX, Xmat, nvec, sigb, 
                               betaZ, Zmat,
                               subjind, sig.subj,
                               num.sim, num.perm = 10000, printafter = 100, 
                               Xdist = "categorical", covadjs = c("SK", "alr"),
                               whichtests = c("Kernel", "LogisticM"),
                               adjmethod = "bonferroni",
                               lmertest = c("Wald", "LRT"), lmer.nump = 100, 
                               diffcyt.method = c("edgeR", "voom"),
                               storeY = FALSE){
  
  #rows of Xmat denote sample units, columns denote different X variables
  #sigb is the covariance matrix of the random effects (num.cat x num.cat)
  #subjind is the vector of subjects, sig.subj is the subject specific random
      # effect variance. 
  
  pval.logistic <- pval.kernel <- pval.logisticM <- pval.diffcyt <- pval.chisq <- NA
  
  if ("Kernel" %in% whichtests){
    pval.kernel <- matrix(NA, nrow = num.sim, ncol = length(covadjs))
    colnames(pval.kernel) <- covadjs
  }
  
  if ("Logistic" %in% whichtests){
    pval.logistic <- numeric(num.sim)
  }
  
  if ("LogisticM" %in% whichtests){
    pval.logisticM <- matrix(NA, nrow = num.sim, ncol = length(lmertest))
    colnames(pval.logisticM) <- lmertest
  }
  
  if ("diffcyt" %in% whichtests){
    pval.diffcyt <- matrix(NA, nrow = num.sim, ncol = length(diffcyt.method))
    colnames(pval.diffcyt) <- diffcyt.method
  }
  
  if (storeY == T){
    ylist <- list()
  }else{
    ylist = NA
  }
  
  
  
  #### Simulation begins ####
  
  for (i in 1:num.sim){
    if (floor(i/printafter) == i/printafter){print(i)}
    
    
    y <- simProps.d(beta, cbind(betaX, betaZ), cbind(Xmat, Zmat), nvec, sigb,
                  subjind, sig.subj)
    
    ng <- nrow(y)/2
    num.cat <- ncol(y)
    
    propmat <- y/rep(nvec, num.cat)
    
    #x is the predictor of interest. 
    #z is the covariate, works for only one covariate.
    x <- as.vector(Xmat)
    z <- as.vector(Zmat)
    
    #### Scaling for continuous predictors ####
    
    if (Xdist == "Gaussian" | Xdist == "LK"){
      x <- (x - mean(x))/sd(x)
    }
    
    
    #### dcor test ####
    
    if ("Kernel" %in% whichtests){
      for (j in 1:length(covadjs)){
        covadj.curr <- covadjs[j]
        
        if (covadj.curr == "SK"){
          out <- try(dcor.comp.cov(propmat, x, Zmat, num.perm, 
                                   Xdist = Xdist, covadj = covadj.curr), silent = T)
          if (class(out) == "try-error"){
            pval.kernel[i, j] <-NA
          }else{
            pval.kernel[i, j] <- out$pval
          }
        }
        
        if (covadj.curr == "alr"){
          alrs <- t(apply(propmat, 1, alr))
          fit1 <- try(lm(alrs ~ z))
          if (class(fit1)[1] == "try-error"){
            pval.kernel[i, j] <-NA
          }else{
            fit2 <- lm(x ~ z)
            respropmat <- t(apply(residuals(fit1), 1, invalr))
            resx <- residuals(fit2)
            out <- dcor.comp(respropmat, resx, num.perm, 
                             Xdist = "LK")
            pval.kernel[i, j] <- out$pval
          }
        }
        
        
        if (covadj.curr == "clr"){
          clrs <- t(apply(propmat, 1, clr))
          fit1 <- try(lm(clrs ~ z))
          if (class(fit1)[1] == "try-error"){
            pval.kernel[i, j] <-NA
          }else{
            fit2 <- lm(x ~ z)
            respropmat <- t(apply(residuals(fit1), 1, invclr))
            resx <- residuals(fit2)
            out <- dcor.comp(respropmat, resx, num.perm, 
                             Xdist = "LK")
            pval.kernel[i, j] <- out$pval
          }
        }
      }#End covadj loop
    }#End kernel loop
    
    
    
    #### diffcyt tests ####
    if ("diffcyt" %in% whichtests){
      for (k in 1:length(diffcyt.method)){
        curr.method <- diffcyt.method[k]
        temp <- try(diffcyt.test.cov(y, x, z, method = curr.method), silent = T)
        if (class(out) == "try-error"){
          pval.diffcyt[i, k] <- NA
        }else{  
          pval.diffcyt[i, k] <- min(p.adjust(temp, method = adjmethod))
        }
      }
    }
    
    
    
    #### Logistic or Logistic Mixed models ####
    
    if ("LogisticM" %in% whichtests | "Logistic" %in% whichtests){
      ptemp <- rep(NA, num.cat) 
      ptemp.mixed <- matrix(NA, nrow = num.cat, ncol = length(lmertest))
      
      for (k in 1:num.cat){ 
        y.bin <- y[, k]
        
        if ("Logistic" %in% whichtests){
          fit <- glm(y.bin/nvec ~ x + z, weights = nvec, 
                     family = binomial(link = "logit"))
          ptemp[k] <- summary(fit)$coefficients[2, 4]
        }  
        
        
        if ("LogisticM" %in% whichtests){
          id <- subjind
          fitmixed <- suppressMessages(try(glmer(y.bin/nvec ~ x + z + (1|id) + (1|id:z), 
                                weights = nvec, 
                                family = binomial(link = "logit")), silent = T))
          if (class(fitmixed) == "try-error"){
            ptemp.mixed[k, colind]<- NA
          }else{
            if ("Wald" %in% lmertest){
              colind <- which(lmertest == "Wald")
              ptemp.mixed[k, colind] <- summary(fitmixed)$coefficients[2, 4]
            }
            
            if ("LRT" %in% lmertest){
              colind <- which(lmertest == "LRT")
              temptry <- suppressMessages(try(drop1(fitmixed, test="Chisq")$"Pr(Chi)"[2], silent = T))
              if (class(temptry) == "try-error"){
                ptemp.mixed[k, colind]<- NA
              }else{
                ptemp.mixed[k, colind] <- temptry
              }
            }
            
            if ("Permutation" %in% lmertest){
              colind <- which(lmertest == "Permutation")
              lrstat <- drop1(fitmixed, test="Chisq")$LRT[2]
              lrstat.perm <- numeric(lmer.nump)
              for (p in 1:lmer.nump){
                xperm <- sample(x)
                fitmixed <- glmer(y.bin/nvec ~ xperm + z + (1|id), weights = nvec, 
                                  family = binomial(link = "logit"))
                lrstat.perm[p] <- drop1(fitmixed, test="Chisq")$LRT[2]
              }
              ptemp.mixed[k, colind] <- (sum(lrstat.perm >= lrstat)+1)/(lmer.nump+1)
            }#perm test loop ends
          }#try loop ends
        }#logisticM loop ends  
      }#categories loop ends
      
      if ("Logistic" %in% whichtests){
        pval.logistic[i] <- min(p.adjust(ptemp, method = adjmethod))
      } 
      
      if ("LogisticM" %in% whichtests){
        pval.logisticM[i, ] <- apply(ptemp.mixed, 2, function(w){
          return(min(p.adjust(w, method = adjmethod)))
        })
      } 
      
    }#logistic loop ends
    
    
    
    #### Storing the simulated data ####
    
    if (storeY == T){
      ylist[[i]] <- y
    }
  }#sim loop ends  
  
  
  retlist <-  list(logistic = pval.logistic, kernel = pval.kernel, 
                   chisq = pval.chisq, logistic.mixed = pval.logisticM,
                   diffcyt = pval.diffcyt,
                   ylist = ylist)
  return(retlist)
}





# simprops function for the dependent case

simProps.d <- function(beta, betaX, Xmat, nvec, sigb = NULL,
                       subjind = NULL, sig.subj = NULL){
  #rows of Xmat denote sample units, columns denote different X variables
  #sigb is the covariance matrix of the random effects (num.cat x num.cat)
  
  num.cat <- length(beta)
  y <- matrix(NA, nrow = nrow(Xmat), ncol = num.cat)
  
  if(is.null(sigb)){
    for (i in 1:nrow(Xmat)){
      probs <- getProps.mlogit(beta, betaX, as.matrix(Xmat[i, ]))
      y[i, ] <- t(rmultinom(n = 1, size = nvec[i], prob = probs))
    }
  }else{
    bmat <- matrix(0, nrow = nrow(Xmat), ncol = num.cat)
    ref <- which(beta == 0)
    bmat[, ref] <- rep(0, nrow(Xmat))
    bmat[, -ref] <- mvrnorm(nrow(Xmat), mu = rep(0, num.cat-1), sigb)
    
    cmat <- matrix(0, nrow = length(unique(subjind)), ncol = num.cat)
    ref <- which(beta == 0)
    cmat[, ref] <- rep(0, length(unique(subjind)))
    cmat[, -ref] <- mvrnorm(length(unique(subjind)), mu = rep(0, num.cat-1),
                            diag(num.cat-1)*sig.subj^2)
    
    for (i in 1:nrow(Xmat)){
      probs <- getProps.mlogit(beta + bmat[i, ] + cmat[subjind[i],], 
                               betaX, as.matrix(Xmat[i, ]))
      y[i, ] <- t(rmultinom(n = 1, size = nvec[i], prob = probs))
    }
  }
  
  return(y)
}

################################################################################






##### Wrapper to run the diffcyt methods #####

diffcyt.test <- function(ymat, x, method = "GLMM"){
  sid <- factor(paste0("sample", 1:length(x)))
  gid <- factor(paste0("group", (1+x)))
  colData <- DataFrame(sample_id=sid,
                       group_id=gid)
  se <- SummarizedExperiment(assays=list(counts=t(ymat)), 
                             colData = colData)
  cid <- factor(1:ncol(ymat))
  nc <- apply(ymat, 2, sum)
  elementMD <- DataFrame(cluster_id=cid,
                         n_cells=nc)
  elementMetadata(se) <- elementMD
  
  experiment_info <- data.frame(
    sample_id = sid, 
    group_id = gid, 
    stringsAsFactors = FALSE
  )
  
  if (method == "GLMM"){
    formula <- createFormula(experiment_info, 
                             cols_fixed = "group_id", 
                             cols_random = "sample_id")
    contrast <- createContrast(c(0, 1))
    out <- testDA_GLMM(se, formula, contrast)
    return(as.data.frame(rowData(out))$p_val)
  }
  
  if (method == "edgeR"){
    design <- createDesignMatrix(experiment_info, 
                                 cols_design = "group_id")
    contrast <- createContrast(c(0, 1))
    out <- testDA_edgeR(se, design, contrast)
    return(as.data.frame(rowData(out))$p_val)
  }
  
  if (method == "voom"){
    design <- createDesignMatrix(experiment_info, 
                                 cols_design = "group_id")
    contrast <- createContrast(c(0, 1))
    out <- testDA_voom(se, design, contrast)
    return(as.data.frame(rowData(out))$p_val)
  }
}


# diffcyt with covariates
diffcyt.test.cov <- function(ymat, x, z, method = "GLMM"){
  sid <- factor(paste0("sample", 1:length(x)))
  gid <- factor(paste0("group", (1+x)))
  colData <- DataFrame(sample_id=sid,
                       group_id=gid)
  se <- SummarizedExperiment(assays=list(counts=t(ymat)), 
                             colData = colData)
  cid <- factor(1:ncol(ymat))
  nc <- apply(ymat, 2, sum)
  elementMD <- DataFrame(cluster_id=cid,
                         n_cells=nc)
  elementMetadata(se) <- elementMD
  
  experiment_info <- data.frame(
    sample_id = sid, 
    group_id = gid,
    z = z,
    stringsAsFactors = FALSE
  )
  
  if (method == "GLMM"){
    formula <- createFormula(experiment_info, 
                             cols_fixed = c("group_id", "z"), 
                             cols_random = "sample_id")
    contrast <- createContrast(c(0, 1, 0))
    out <- testDA_GLMM(se, formula, contrast)
    return(as.data.frame(rowData(out))$p_val)
  }
  
  if (method == "edgeR"){
    design <- createDesignMatrix(experiment_info, 
                                 cols_design = c("group_id", "z"))
    contrast <- createContrast(c(0, 1, 0))
    out <- testDA_edgeR(se, design, contrast)
    return(as.data.frame(rowData(out))$p_val)
  }
  
  if (method == "voom"){
    design <- createDesignMatrix(experiment_info, 
                                 cols_design = c("group_id", "z"))
    contrast <- createContrast(c(0, 1, 0))
    out <- testDA_voom(se, design, contrast)
    return(as.data.frame(rowData(out))$p_val)
  }
}


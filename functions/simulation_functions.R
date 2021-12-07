

###### Functions for simulating data for cell type abundance comparison #######

library("MASS")
library("coda.base")
library("lme4")
library("SummarizedExperiment")
library("diffcyt")
library("energy")
library("vegan")
library("MiRKAT")




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
                               Xdist = "categorical", gaussian = F,
                               whichtests = c("Kernel", "LogisticM"),
                               adjmethod = "bonferroni",
                               lmertest = c("Wald", "LRT"), lmer.nump = 100, 
                               diffcyt.method = c("edgeR", "voom"),
                               mirkat.method = c("Opt"),
                               storeY = FALSE,
                               standardize = TRUE){
  
  # rows of Xmat denote sample units, columns denote different X variables
  # sigb is the covariance matrix of the random effects (num.cat x num.cat)
  # whichtests: which tests to conduct
  # adjmethod: method to use for multiple testing adjustment
  # lmertest: Which logistic mixed models to use
  # lmer.nump: number of permutations for the LRT-permutation method
  # diffcyt.method: which diffcyt methods to use
  # mirkat.method: which mirkat methods to use
  # standardize: standardizes continuous predictor values if TRUE
  
  pval.logistic <- pval.kernel <- pval.kernel.bc <- pval.euclidean <- pval.chisq <- 
    pval.logisticM <- pval.diffcyt <- pval.pm <- pval.dcor <- pval.mirkat <- NA
  
  if ("Kernel" %in% whichtests){ #CODAK with AD
    pval.kernel <- numeric(num.sim)
  }
  
  if ("Kernel.BC" %in% whichtests){ #CODAK with Bray-Curtis
    pval.kernel.bc <- numeric(num.sim)
  }
  
  if ("Euclidean" %in% whichtests){ #CODAK with ED, not recommended
    pval.euclidean <- numeric(num.sim)
  }
  
  if ("dcor" %in% whichtests){ #Original dcor, not recommended
    pval.dcor <- numeric(num.sim)
  }
  
  if ("permanova" %in% whichtests){ #permanova with AD
    pval.pm <- numeric(num.sim)
  }  
  
  if ("Logistic" %in% whichtests){ #Logistic model without OLRE, not practical
    pval.logistic <- numeric(num.sim)
  }
  
  if ("LogisticM" %in% whichtests){ #Logistic GLMM
    pval.logisticM <- matrix(NA, nrow = num.sim, ncol = length(lmertest))
    colnames(pval.logisticM) <- lmertest
  }
  
  if ("diffcyt" %in% whichtests){ #diffcyt
    pval.diffcyt <- matrix(NA, nrow = num.sim, ncol = length(diffcyt.method))
    colnames(pval.diffcyt) <- diffcyt.method
  }
  
  if ("mirkat" %in% whichtests){ #mirkat with AD, BCD or Opt
    ncolmir <- length(mirkat.method) + 2*as.numeric("Opt" %in% mirkat.method)
    pval.mirkat <- matrix(NA, nrow = num.sim, ncol = ncolmir)
    
    if ("Opt" %in% mirkat.method){
      colnamesmir <- c("AD", "BCD", mirkat.method)
      optcols <- c(1, 2, 2 + which(mirkat.method == "Opt"))
    }else{
      colnamesmir <- mirkat.method
    }
    colnames(pval.mirkat) <- colnamesmir
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
    
    if (standardize == T & Xdist %in% c("gaussian", "LK")){
      x <- (x - mean(x))/sd(x)
    }
    
    
    #### CODAK test ####
    
    if ("Kernel" %in% whichtests){
      out <- try(dcor.comp(propmat, x, num.perm, Xdist = Xdist, gaussian = gaussian,
                           result = c("test", "cor")), silent = T) #remove cor
      if (class(out) == "try-error"){
        pval.kernel[i] <- NA
      }else{
        pval.kernel[i] <- out$pval
      }
    }  
    
    if ("Kernel.BC" %in% whichtests){
      out <- try(dcor.comp(propmat, x, num.perm, Xdist = Xdist, 
                           compdist = "BCD", result = c("test")), silent = T) 
      if (class(out) == "try-error"){
        pval.kernel.bc[i] <- NA
      }else{
        pval.kernel.bc[i] <- out$pval
      }
    } 
    
    if ("Euclidean" %in% whichtests){
      out <- dcor.comp(propmat, x, num.perm, Xdist = Xdist, gaussian = gaussian,
                       compdist = "euclidean", result = c("test")) 
      pval.euclidean[i] <- out$pval
    }  
    
    
    #PERMANOVA
    if ("permanova" %in% whichtests){
      d <- dist(propmat, method = 'aitchison')
      
      pmfit <- try(adonis(d ~ x, permutations = num.perm), silent = T)
      if (class(pmfit) == "try-error"){
        pval.pm[i] <- NA
      }else{
        pval.pm[i] <- pmfit$aov.tab$"Pr(>F)"[1]
      }
    }
    
    
    
    #Original dcor
    if ("dcor" %in% whichtests){
      pval.dcor[i] <- dcov.test(propmat, x, R = num.perm)$p.value
    }
    
    
    #mirkat
    if ("mirkat" %in% whichtests){
      typ = "C"
      if (Xdist == "categorical"){typ = "D"}
      
      for (k in 1:length(mirkat.method)){
        curr.method <- mirkat.method[k]
        temp <- try(mirkat.test(y, x, method = curr.method, typ = typ, num.perm = num.perm), silent = T)
        
        if (class(temp) == "try-error"){
          temp <- data.frame(p_values = rep(NA, 1 + as.numeric("Opt" %in% mirkat.method)), 
                             omnibus_p = NA)
        }
        
        if (curr.method == "Opt"){
          pval.mirkat[i, optcols] <- c(temp$p_values, temp$omnibus_p[1])
        }else{
          pval.mirkat[i, (1:ncol(pval.mirkat))[-optcols][k]] <- temp$p_values
        }
      }
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
                   kernelbc = pval.kernel.bc, euclidean = pval.euclidean,
                   permanova = pval.pm, dcor = pval.dcor,
                   chisq = pval.chisq, logistic.mixed = pval.logisticM,
                   diffcyt = pval.diffcyt, mirkat = pval.mirkat,
                   ylist = ylist)
  return(retlist)
}

##############################################################################








# Mixed effects mlogit: OLRE added, realistic data + Covariate

simulation_study_3 <- function(beta, betaX, Xmat, nvec, sigb, 
                               betaZ = NULL, Zmat, 
                               num.sim, num.perm = 10000, printafter = 100, 
                               Xdist = "categorical", covadjs = c("SK", "alr"),
                               covpermute = c("", "shuffleRes"),
                               whichtests = c("Kernel", "LogisticM"),
                               adjmethod = "bonferroni",
                               lmertest = c("Wald", "LRT"), lmer.nump = 100, 
                               diffcyt.method = c("edgeR", "voom"),
                               mirkat.method = c("Opt"),
                               storeY = FALSE,
                               seedmat = NULL){
  
  #rows of Xmat denote sample units, columns denote different X variables
  #sigb is the covariance matrix of the random effects (num.cat x num.cat)
  #covadjs is the list of CODAK methods using AD, options are "SK", "alr" and "clr"
  #covpermute is a vector containing the permutation methods for the CODAK methods. It should have
  #           the same length as covadjs. Options are "shuffleY", "shuffleRes" and "FL".
  #           A value is required for the "SK" method, but it will not be used.
  
  pval.logistic <- pval.kernel <- pval.kernel.bc <- pval.euclidean <- pval.chisq <- 
    pval.logisticM <- pval.diffcyt <- pval.pm <- pval.mirkat <- pval.dcor <- NA
  
  if ("Kernel" %in% whichtests){
    pval.kernel <- matrix(NA, nrow = num.sim, ncol = length(covadjs))
    colnames(pval.kernel) <- paste(covadjs, covpermute)
  }
  
  if ("Kernel.BC" %in% whichtests){
    pval.kernel.bc <- numeric(num.sim)
  }
  
  if ("Euclidean" %in% whichtests){
    pval.euclidean <- numeric(num.sim)
  }
  
  if ("permanova" %in% whichtests){
    pval.pm <- numeric(num.sim)
  }  
  
  if ("dcor" %in% whichtests){
    pval.dcor <- numeric(num.sim)
  }
  
  if ("Chisq" %in% whichtests){
    pval.chisq <- numeric(num.sim)
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
  
  if ("mirkat" %in% whichtests){
    ncolmir <- length(mirkat.method) + 2*as.numeric("Opt" %in% mirkat.method)
    pval.mirkat <- matrix(NA, nrow = num.sim, ncol = ncolmir)
    
    if ("Opt" %in% mirkat.method){
      colnamesmir <- c("AD", "BCD", mirkat.method)
      optcols <- c(1, 2, 2 + which(mirkat.method == "Opt"))
    }else{
      colnamesmir <- mirkat.method
    }
    colnames(pval.mirkat) <- colnamesmir
  }
  
  if (storeY == T){
    ylist <- list()
  }else{
    ylist = NA
  }
  
  zsim <- F
  if (is.null(betaZ)){zsim = T}
  
  
  
  #### Simulation begins ####
  
  for (i in 1:num.sim){
    if (floor(i/printafter) == i/printafter){print(i)}
    
    if (!is.null(seedmat)){
      set.seed(seedmat[i, 1])
    }
    
    if (zsim){
      probs.addz <- add.effects(probs, maxeffect = 0.002, ind.maxeff = 1, n.nochange = 0)
      betaZ <- as.matrix(getES.mlogit(probs.addz) - betarev1)
    }
    
    y <- simProps(beta, cbind(betaX, betaZ), cbind(Xmat, Zmat), nvec, sigb)
    
    ng <- nrow(y)/2
    num.cat <- ncol(y)
    
    propmat <- y/rep(nvec, num.cat)
    
    #x is the predictor of interest. 
    #z is the covariate, works for only one covariate.
    
    #x <- as.vector(Xmat)
    z <- as.vector(Zmat)
    
    #### Scaling for continuous predictors ####
    
    #if (Xdist == "gaussian" | Xdist == "LK"){
    #  x <- (x - mean(x))/sd(x)
    #}
    
    
    #### CODAK test ####
    
    if ("Kernel" %in% whichtests){
      for (j in 1:length(covadjs)){
        
        if (!is.null(seedmat)){
          set.seed(seedmat[i, 2])
        }
        
        covadj.curr <- covadjs[j]
        
        if (covadj.curr %in% c("SK", "SK2", "SK3")){
          out <- try(dcor.comp.cov(propmat, x, Zmat, num.perm, 
                                   Xdist = Xdist, covadj = covadj.curr, 
                                   result = c("test", "cor")), silent = T)
          if (class(out) == "try-error"){
            pval.kernel[i, j] <-NA
          }else{
            pval.kernel[i, j] <- out$pval
          }
        }
        
        if (covadj.curr %in% c("alr", "clr")){
          out <- try(dcor.comp.cov(propmat, x, Zmat, num.perm, 
                                   Xdist = Xdist, covadj = covadj.curr, permutation = covpermute[j],
                                   result = c("test", "cor")), silent = T)
          if (class(out)[1] == "try-error"){
            pval.kernel[i, j] <-NA
          }else{
            pval.kernel[i, j] <- out$pval
          }
        }
        
        
      }#End covadj loop
    }#End kernel loop
    
    
    
    
    
    #PERMANOVA
    
    if ("permanova" %in% whichtests){
      
      if (!is.null(seedmat)){
        set.seed(seedmat[i, 2])
      }
      
      d <- dist(propmat, method = 'aitchison')
      
      pmfit <- try(adonis(d ~ z + x, permutations = num.perm), silent = T)
      if (class(pmfit)[1] == "try-error"){
        pval.pm[i] <- NA
      }else{
        pval.pm[i] <- pmfit$aov.tab$"Pr(>F)"[2]
      }
    }
    
    
    
    #Original dcor - skipped for covariates
    #CODAK-BC - skipped for covariates
    #CODAK-ED - skipped for covariates
    
    
    
    
    #### diffcyt tests ####
    if ("diffcyt" %in% whichtests){
      for (k in 1:length(diffcyt.method)){
        curr.method <- diffcyt.method[k]
        temp <- diffcyt.test.cov(y, x, z, method = curr.method)
        pval.diffcyt[i, k] <- min(p.adjust(temp, method = adjmethod))
      }
    }
    
    
    
    
    #mirkat
    if ("mirkat" %in% whichtests){
      typ = "C"
      if (Xdist == "categorical"){typ = "D"}
      
      for (k in 1:length(mirkat.method)){
        
        if (!is.null(seedmat)){
          set.seed(seedmat[i, 2])
        }
        
        curr.method <- mirkat.method[k]
        temp <- try(mirkat.test.cov(propmat, x, z, 
                                    method = curr.method, typ = typ, num.perm = num.perm), silent = T)
        
        if (class(temp) == "try-error"){
          temp <- data.frame(p_values = rep(NA, 1 + as.numeric("Opt" %in% mirkat.method)), 
                             omnibus_p = NA)
        }
        
        if (curr.method == "Opt"){
          pval.mirkat[i, optcols] <- c(temp$p_values, temp$omnibus_p[1])
        }else{
          pval.mirkat[i, (1:ncol(pval.mirkat))[-optcols][k]] <- temp$p_values
        }
      }
    }
    
    
    
    
    
    #### Logistic or Logistic Mixed models ####
    
    if ("LogisticM" %in% whichtests | "Logistic" %in% whichtests){
      ptemp <- rep(NA, num.cat) #num.cat-1)
      ptemp.mixed <- matrix(NA, nrow = num.cat, ncol = length(lmertest))
      
      for (k in 1:num.cat){ #(num.cat-1)){
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
                   euclidean = pval.euclidean, kernelbc = pval.kernel.bc, 
                   permanova = pval.pm, dcor = pval.dcor,
                   chisq = pval.chisq, logistic.mixed = pval.logisticM,
                   diffcyt = pval.diffcyt, mirkat = pval.mirkat,
                   ylist = ylist)
  return(retlist)
}


############################################################################






# Simulation 3 function + potential dependence in the data 

simulation_study_4 <- function(beta, betaX, Xmat, nvec, sigb, 
                               betaZ = NULL, Zmat,
                               subjind, sig.subj,
                               num.sim, num.perm = 10000, printafter = 100, 
                               Xdist = "categorical", covadjs = c("SK", "alr"),
                               covpermute = c("", "shuffleRes"),
                               whichtests = c("Kernel", "LogisticM"),
                               adjmethod = "bonferroni",
                               lmertest = c("Wald", "LRT"), lmer.nump = 100, 
                               diffcyt.method = c("edgeR", "voom"),
                               mirkat.method = c("Opt"),
                               storeY = FALSE,
                               seedmat = NULL){
  
  #rows of Xmat denote sample units, columns denote different X variables
  #sigb is the covariance matrix of the random effects (num.cat x num.cat)
  #subjind is the vector of subjects, sig.subj is the subject specific random
      # effect variance. 
  
  pval.logistic <- pval.kernel <- pval.kernel.bc <- pval.euclidean <- pval.chisq <- 
    pval.logisticM <- pval.diffcyt <- pval.pm <- pval.mirkat <- pval.dcor <- NA
  
  if ("Kernel" %in% whichtests){
    pval.kernel <- matrix(NA, nrow = num.sim, ncol = length(covadjs))
    colnames(pval.kernel) <- covadjs
  }
  
  if ("Kernel.BC" %in% whichtests){
    pval.kernel.bc <- numeric(num.sim)
  }
  
  if ("Euclidean" %in% whichtests){
    pval.euclidean <- numeric(num.sim)
  }
  
  if ("permanova" %in% whichtests){
    pval.pm <- numeric(num.sim)
  }  
  
  if ("dcor" %in% whichtests){
    pval.dcor <- numeric(num.sim)
  }
  
  if ("Chisq" %in% whichtests){
    pval.chisq <- numeric(num.sim)
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
  
  if ("mirkat" %in% whichtests){
    ncolmir <- length(mirkat.method) + 2*as.numeric("Opt" %in% mirkat.method)
    pval.mirkat <- matrix(NA, nrow = num.sim, ncol = ncolmir)
    
    if ("Opt" %in% mirkat.method){
      colnamesmir <- c("AD", "BCD", mirkat.method)
      optcols <- c(1, 2, 2 + which(mirkat.method == "Opt"))
    }else{
      colnamesmir <- mirkat.method
    }
    colnames(pval.mirkat) <- colnamesmir
  }
  
  if (storeY == T){
    ylist <- list()
  }else{
    ylist = NA
  }
  
  zsim <- F
  if (is.null(betaZ)){zsim = T}
  
  
  
  #### Simulation begins ####
  
  for (i in 1:num.sim){
    if (floor(i/printafter) == i/printafter){print(i)}
    
    if (!is.null(seedmat)){
      set.seed(seedmat[i, 1])
    }
    
    if (zsim){
      probs.addz <- add.effects(probs, maxeffect = 0.002, ind.maxeff = 1, n.nochange = 0)
      betaZ <- as.matrix(getES.mlogit(probs.addz) - betarev1)
    }
    
    y <- simProps.d(beta, cbind(betaX, betaZ), cbind(Xmat, Zmat), nvec, sigb,
                  subjind, sig.subj)
    
    ng <- nrow(y)/2
    num.cat <- ncol(y)
    
    propmat <- y/rep(nvec, num.cat)
    
    #x is the predictor of interest. 
    #z is the covariate, works for only one covariate.
    
    #x <- as.vector(Xmat)
    z <- as.vector(Zmat)
    
    #### Scaling for continuous predictors ####
    
    #if (Xdist == "gaussian" | Xdist == "LK"){
    #  x <- (x - mean(x))/sd(x)
    #}
    
    
    #### CODAK test ####
    
    if ("Kernel" %in% whichtests){
      for (j in 1:length(covadjs)){
        
        if (!is.null(seedmat)){
          set.seed(seedmat[i, 2])
        }
        
        covadj.curr <- covadjs[j]
        
        if (covadj.curr %in% c("SK", "SK2", "SK3")){
          out <- try(dcor.comp.cov(propmat, x, Zmat, num.perm, 
                                   Xdist = Xdist, covadj = covadj.curr, 
                                   result = c("test", "cor")), silent = T)
          if (class(out) == "try-error"){
            pval.kernel[i, j] <-NA
          }else{
            pval.kernel[i, j] <- out$pval
          }
        }
        
        if (covadj.curr %in% c("alr", "clr")){
          out <- try(dcor.comp.cov(propmat, x, Zmat, num.perm, 
                                   Xdist = Xdist, covadj = covadj.curr, permutation = covpermute[j],
                                   result = c("test", "cor")), silent = T)
          if (class(out)[1] == "try-error"){
            pval.kernel[i, j] <-NA
          }else{
            pval.kernel[i, j] <- out$pval
          }
        }
        
        
      }#End covadj loop
    }#End kernel loop
    
    
    
    
    
    #PERMANOVA
    
    if ("permanova" %in% whichtests){
      
      if (!is.null(seedmat)){
        set.seed(seedmat[i, 2])
      }
      
      d <- dist(propmat, method = 'aitchison')
      
      pmfit <- try(adonis(d ~ z + x, permutations = num.perm), silent = T)
      if (class(pmfit)[1] == "try-error"){
        pval.pm[i] <- NA
      }else{
        pval.pm[i] <- pmfit$aov.tab$"Pr(>F)"[2]
      }
    }
    
    
    
    #Original dcor - skipped for covariates
    #CODAK-BC - skipped for covariates
    #CODAK-ED - skipped for covariates
    
    
    
    
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
    
    
    
    
    #mirkat
    if ("mirkat" %in% whichtests){
      typ = "C"
      if (Xdist == "categorical"){typ = "D"}
      
      for (k in 1:length(mirkat.method)){
        
        if (!is.null(seedmat)){
          set.seed(seedmat[i, 2])
        }
        
        curr.method <- mirkat.method[k]
        temp <- try(mirkat.test.cov(propmat, x, z, 
                                    method = curr.method, typ = typ, num.perm = num.perm), silent = T)
        
        if (class(temp) == "try-error"){
          temp <- data.frame(p_values = rep(NA, 1 + as.numeric("Opt" %in% mirkat.method)), 
                             omnibus_p = NA)
        }
        
        if (curr.method == "Opt"){
          pval.mirkat[i, optcols] <- c(temp$p_values, temp$omnibus_p[1])
        }else{
          pval.mirkat[i, (1:ncol(pval.mirkat))[-optcols][k]] <- temp$p_values
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
          fitmixed <- try(glmer(y.bin/nvec ~ x + z + (1|id) + (1|id:z), 
                                weights = nvec, 
                                family = binomial(link = "logit")), silent = T)
          if (class(fitmixed) == "try-error"){
            ptemp.mixed[k, colind]<- NA
          }else{
            if ("Wald" %in% lmertest){
              colind <- which(lmertest == "Wald")
              ptemp.mixed[k, colind] <- summary(fitmixed)$coefficients[2, 4]
            }
            
            if ("LRT" %in% lmertest){
              colind <- which(lmertest == "LRT")
              temptry <- try(drop1(fitmixed, test="Chisq")$"Pr(Chi)"[2], silent = T)
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
                   euclidean = pval.euclidean, kernelbc = pval.kernel.bc, 
                   permanova = pval.pm, dcor = pval.dcor,
                   chisq = pval.chisq, logistic.mixed = pval.logisticM,
                   diffcyt = pval.diffcyt, mirkat = pval.mirkat,
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








##### Wrapper function to run the mirkat methods #####
mirkat.test <- function(ymat, x, method = "Opt", typ = "D", num.perm = 1e4){
  
  propmat <- ymat
  xvec <- x
  
  if (method == "BCD"){ #assuming gaussian = F
    D.BC <- as.matrix(vegdist(propmat, method="bray"))
    K.BC <- D2K(D.BC)
    out <- MiRKAT(y = xvec, 
                  #X = covar, 
                  Ks = K.BC, out_type = typ, nperm = num.perm,
                  method = "permutation", returnKRV = TRUE, returnR2 = TRUE)
  }
  
  if (method == "BCD-G"){ #assuming gaussian = T 
    D.BC <- as.matrix(vegdist(propmat, method="bray"))
    alldist <- as.numeric(D.BC)^2
    K.BC <- exp(-as.matrix(D.BC)/median(alldist))
    out <- MiRKAT(y = xvec, 
                  #X = covar, 
                  Ks = K.BC, out_type = typ, nperm = num.perm,
                  method = "permutation", returnKRV = TRUE, returnR2 = TRUE)
  }
  
  
  if (method == "AD"){ #assuming gaussian = F
    D.AD <- as.matrix(dist(propmat, method = "aitchison"))
    K.AD <- D2K(D.AD)
    out <- MiRKAT(y = xvec, 
                  #X = covar, 
                  Ks = K.AD, out_type = typ, nperm = num.perm,
                  method = "permutation", returnKRV = TRUE, returnR2 = TRUE)
  }
  
  if (method == "AD-G"){ #assuming gaussian = T
    D.AD <- as.matrix(dist(propmat, method = "aitchison"))
    alldist <- as.numeric(D.AD)
    K.AD <- exp(-as.matrix(D.AD)/median(alldist))
    out <- MiRKAT(y = xvec, 
                  #X = covar, 
                  Ks = K.AD, out_type = typ, nperm = num.perm,
                  method = "permutation", returnKRV = TRUE, returnR2 = TRUE)
  }
  
  
  if (method == "Opt"){ #only doing for gaussian = F
    D.BC <- as.matrix(vegdist(propmat, method="bray"))
    K.BC <- D2K(D.BC)
    
    D.AD <- as.matrix(dist(propmat, method = "aitchison"))
    K.AD <- D2K(D.AD)
    
    out <- MiRKAT(y = xvec, 
                  #X = covar, 
                  Ks = list(K.AD, K.BC), out_type = typ, nperm = num.perm,
                  method = "permutation", returnKRV = TRUE, returnR2 = TRUE)
  }
  
  return(out)
}




mirkat.test.cov <- function(ymat, x, z, method = "Opt", typ = "D", num.perm = 1e4){
  
  propmat <- ymat
  xvec <- x
  
  if (method == "BCD"){ #assuming gaussian = F
    D.BC <- as.matrix(vegdist(propmat, method="bray"))
    K.BC <- D2K(D.BC)
    out <- MiRKAT(y = xvec, 
                  X = as.matrix(z), 
                  Ks = K.BC, out_type = typ, nperm = num.perm,
                  method = "permutation", returnKRV = TRUE, returnR2 = TRUE)
  }
  
  if (method == "BCD-G"){ #assuming gaussian = T 
    D.BC <- as.matrix(vegdist(propmat, method="bray"))
    alldist <- as.numeric(D.BC)^2
    K.BC <- exp(-as.matrix(D.BC)/median(alldist))
    out <- MiRKAT(y = xvec, 
                  X = as.matrix(z), 
                  Ks = K.BC, out_type = typ, nperm = num.perm,
                  method = "permutation", returnKRV = TRUE, returnR2 = TRUE)
  }
  
  
  if (method == "AD"){ #assuming gaussian = F
    D.AD <- as.matrix(dist(propmat, method = "aitchison"))
    K.AD <- D2K(D.AD)
    out <- MiRKAT(y = xvec, 
                  X = as.matrix(z), 
                  Ks = K.AD, out_type = typ, nperm = num.perm,
                  method = "permutation", returnKRV = TRUE, returnR2 = TRUE)
  }
  
  if (method == "AD-G"){ #assuming gaussian = T
    D.AD <- as.matrix(dist(propmat, method = "aitchison"))
    alldist <- as.numeric(D.AD)
    K.AD <- exp(-as.matrix(D.AD)/median(alldist))
    out <- MiRKAT(y = xvec, 
                  X = as.matrix(z), 
                  Ks = K.AD, out_type = typ, nperm = num.perm,
                  method = "permutation", returnKRV = TRUE, returnR2 = TRUE)
  }
  
  
  if (method == "Opt"){ #only doing for gaussian = F
    D.BC <- as.matrix(vegdist(propmat, method="bray"))
    K.BC <- D2K(D.BC)
    
    D.AD <- as.matrix(dist(propmat, method = "aitchison"))
    K.AD <- D2K(D.AD)
    
    out <- MiRKAT(y = xvec, 
                  X = as.matrix(z), 
                  Ks = list(K.AD, K.BC), out_type = typ, nperm = num.perm,
                  method = "permutation", returnKRV = TRUE, returnR2 = TRUE)
  }
  
  return(out)
}








########## Simulations for the scenario with replacing zeros by pseudocounts ############

simulation_study_pseudo <- function(beta, betaX, Xmat, nvec, sigb, 
                                    num.sim, num.perm = 10000, printafter = 100, 
                                    Xdist = "categorical",
                                    whichtests = "Kernel.pseudo",
                                    standardize = TRUE){
  #Contains the random effect, sigb is the sd of the random effect
  #rows of Xmat denote sample units, columns denote different X variables
  #sigb is the covariance matrix of the random effects (num.cat x num.cat)
  
  
  pval.kernel <- numeric(num.sim)
  numzeros <- numeric(num.sim)
  
  #### Simulation begins ####
  
  for (i in 1:num.sim){
    if (floor(i/printafter) == i/printafter){print(i)}
    
    if (whichtests == "Kernel.pseudo"){
      N <- 0
      while(N == 0){
        y <- simProps(beta, betaX, Xmat, nvec, sigb)
        N <- sum(y == 0)
        numzeros[i] <- N
        y[which(y == 0)] <- 1
      }
    }else{
      N <- 1
      while(N > 0){
        y <- simProps(beta, betaX, Xmat, nvec, sigb)
        N <- sum(y == 0)
      }
    }
    
    
    ng <- nrow(y)/2
    num.cat <- ncol(y)
    
    propmat <- y/rep(nvec, num.cat)
    
    #x is the predictor of interest. Works for only one predictor.
    #separate function needed when adjusting for covariates
    x <- as.vector(Xmat)
    
    
    #### Scaling for continuous predictors ####
    
    if (standardize == T & Xdist %in% c("euclidean", "LK")){
      x <- (x - mean(x))/sd(x)
    }
    
    
    #### CODAK test ####
    
    out <- try(dcor.comp(propmat, x, num.perm, Xdist = Xdist, 
                         result = c("test")), silent = T) 
    if (class(out) == "try-error"){
      pval.kernel[i] <- NA
    }else{
      pval.kernel[i] <- out$pval
    }
  }  
  
  return(list(pval = pval.kernel, numzeros = numzeros))
}

##############################################################################






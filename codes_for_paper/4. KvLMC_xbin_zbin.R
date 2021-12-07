

################# Simulation in the presence of covariates ##################

##### x is binary and z is also binary #####



rm(list = ls())

plotdir <- paste0(getwd(), "/plots/")
resultdir <- paste0(getwd(), "/results/KVLMC_xbin_zbin/")


### Source files
source("./functions/Old version/CODAK_functions_old.R")
source("./functions/Old version/simulation_functions_old.R")




####### Starting simulations #######

# Let's use 20 categories
probs <- c(0.15, 0.13, 0.07, 0.05, 0.1, 0.05, 0.05,
           0.05, 0.02, 0.1, 0.02, 0.003, 0.002, 0.004,
           0.12, 0.04, 0.004, 0.007, 0.02, 0.01)


num.cat <- length(probs)

rep <- 100
num.sim <- 100

#12 observations in each group 
ng <- 6
Xmat <- matrix(rep(c(0, 1), each = 2*ng), ncol = 1)

#covariate values
Zmat <- matrix(rep(c(0, 1), 2*ng), ncol = 1)


#number of events for the individual samples
set.seed(92486)
nvec <- sample(c(30000, 40000, 50000), 4*ng, 
               replace = T)

  
#covariance matrix of random effects
set.seed(92486)
sigbvec <- rep(0.2, num.cat-1) 
sigb <- diag(num.cat-1)*sigbvec^2






#### Case 1: No differences ####

p1 <- p2 <- probs
probs.new <- probs
betarev1 <- getES.mlogit(p1)
betarev2 <- getES.mlogit(p2)

betaX <- as.matrix(betarev2 - betarev1)

#### Effect of z####
probs.addz <- add.effects(probs, maxeffect = 0.002, ind.maxeff = 1, n.nochange = 0)
betaZ <- as.matrix(getES.mlogit(probs.addz) - betarev1)

set.seed(1)
out <- simulation_study_3(betarev1, betaX, Xmat, nvec, sigb,
                          betaZ, Zmat, 
                          num.sim = 10000, num.perm = 100000, 
                          covadjs = c("SK", "alr"), #, "clr"),
                          whichtests = c("Kernel", "LogisticM", "diffcyt"), 
                          diffcyt.method = c("edgeR", "voom"),
                          printafter = 50, 
                          storeY = T)


save(out, probs, probs.new, file = paste0(resultdir, "no_diff.RObj"))








#### Case 2: Small differences in all ####

set.seed(123)
ors <- matrix(NA, nrow = rep, ncol = length(probs))
ds <- numeric(rep)
pval.kernel.sk <- pval.kernel.alr <- pval.kernel.clr <-
  pval.logistic.mixed.wald <- pval.logistic.mixed.lrt <- 
  pval.diffcyt.edger <- pval.diffcyt.voom <- matrix(NA, nrow = rep, ncol = num.sim)

for (i in 1:rep){
  print(paste("rep =", i))
  
  betarev1 <- getES.mlogit(probs)
  probs.new <- add.effects(probs, maxeffect = 0.002, ind.maxeff = 1, n.nochange = 0)
  betarev2 <- getES.mlogit(probs.new)
  
  betaX <- as.matrix(betarev2 - betarev1)
  
  out <- simulation_study_3(betarev1, betaX, Xmat, nvec, sigb, 
                            betaZ, Zmat,
                            num.sim = num.sim, num.perm = 10000,
                            covadjs = c("SK", "alr"),#, "clr"),
                            whichtests = c("Kernel", "LogisticM",
                                           "diffcyt"),
                            diffcyt.method = c("edgeR", "voom"),
                            printafter = 100)
  
  orlist <- numeric(length(probs))
  for (k in 1:length(probs)){
    orlist[k] <- (probs[k]/(1-probs[k]))/(probs.new[k]/(1-probs.new[k]))
    if (orlist[k] < 1) {orlist[k] <- 1/orlist[k]}
  }
  
  ors[i, ] <- orlist
  
  ds[i] <- dist(rbind(probs, probs.new),method = 'aitchison')
  
  pval.kernel.sk[i, ] <- out$kernel[,1] 
  pval.kernel.alr[i, ] <- out$kernel[,2] 
  pval.logistic.mixed.wald[i, ] <- out$logistic.mixed[,1]
  pval.logistic.mixed.lrt[i, ] <- out$logistic.mixed[,2]
  pval.diffcyt.edger[i, ] <- out$diffcyt[,1] 
  pval.diffcyt.voom[i, ] <- out$diffcyt[,2]
}

save(pval.kernel.sk, pval.kernel.alr, pval.logistic.mixed.wald, pval.logistic.mixed.lrt,
     pval.diffcyt.edger, pval.diffcyt.voom,
     ors, ds, probs, probs.new,
     file = paste0(resultdir, "small_diff_all.RObj"))







#### Case 3: Small differences in some ####

set.seed(123)
ors <- matrix(NA, nrow = rep, ncol = length(probs))
ds <- numeric(rep)
pval.kernel.sk <- pval.kernel.alr <- pval.kernel.clr <-
  pval.logistic.mixed.wald <- pval.logistic.mixed.lrt <- 
  pval.diffcyt.edger <- pval.diffcyt.voom <- matrix(NA, nrow = rep, ncol = num.sim)

for (i in 1:rep){
  print(paste("rep =", i))
  
  betarev1 <- getES.mlogit(probs)
  probs.new <- add.effects(probs, maxeffect = 0.002, ind.maxeff = 1, n.nochange = 10, ind.nochange = 2:11)
  betarev2 <- getES.mlogit(probs.new)
  
  betaX <- as.matrix(betarev2 - betarev1)
  
  out <- simulation_study_3(betarev1, betaX, Xmat, nvec, sigb, 
                            betaZ, Zmat,
                            num.sim = num.sim, num.perm = 10000,
                            covadjs = c("SK", "alr"),#, "clr"),
                            whichtests = c("Kernel", "LogisticM",
                                           "diffcyt"),
                            diffcyt.method = c("edgeR", "voom"),
                            printafter = 100)
  
  orlist <- numeric(length(probs))
  for (k in 1:length(probs)){
    orlist[k] <- (probs[k]/(1-probs[k]))/(probs.new[k]/(1-probs.new[k]))
    if (orlist[k] < 1) {orlist[k] <- 1/orlist[k]}
  }
  
  ors[i, ] <- orlist
  
  ds[i] <- dist(rbind(probs, probs.new),method = 'aitchison')
  
  pval.kernel.sk[i, ] <- out$kernel[,1] 
  pval.kernel.alr[i, ] <- out$kernel[,2] 
  pval.logistic.mixed.wald[i, ] <- out$logistic.mixed[,1]
  pval.logistic.mixed.lrt[i, ] <- out$logistic.mixed[,2]
  pval.diffcyt.edger[i, ] <- out$diffcyt[,1] 
  pval.diffcyt.voom[i, ] <- out$diffcyt[,2]
}

save(pval.kernel.sk, pval.kernel.alr, pval.logistic.mixed.wald, pval.logistic.mixed.lrt,
     pval.diffcyt.edger, pval.diffcyt.voom,
     ors, ds, probs, probs.new,
     file = paste0(resultdir, "small_diff_some.RObj"))






#### Case 4: large differences in few ####

set.seed(123)
ors <- matrix(NA, nrow = rep, ncol = length(probs))
ds <- numeric(rep)
pval.kernel.sk <- pval.kernel.alr <- pval.kernel.clr <-
  pval.logistic.mixed.wald <- pval.logistic.mixed.lrt <- 
  pval.diffcyt.edger <- pval.diffcyt.voom <- matrix(NA, nrow = rep, ncol = num.sim)

for (i in 1:rep){
  print(paste("rep =", i))
  
  betarev1 <- getES.mlogit(probs)
  probs.new <- add.effects(probs, maxeffect = 0.004, ind.maxeff = 1, n.nochange = 15, ind.nochange = 2:16)
  betarev2 <- getES.mlogit(probs.new)
  
  betaX <- as.matrix(betarev2 - betarev1)
  
  out <- simulation_study_3(betarev1, betaX, Xmat, nvec, sigb, 
                            betaZ, Zmat,
                            num.sim = num.sim, num.perm = 10000,
                            covadjs = c("SK", "alr"),#, "clr"),
                            whichtests = c("Kernel", "LogisticM",
                                           "diffcyt"),
                            diffcyt.method = c("edgeR", "voom"),
                            printafter = 100)
  
  orlist <- numeric(length(probs))
  for (k in 1:length(probs)){
    orlist[k] <- (probs[k]/(1-probs[k]))/(probs.new[k]/(1-probs.new[k]))
    if (orlist[k] < 1) {orlist[k] <- 1/orlist[k]}
  }
  
  ors[i, ] <- orlist
  
  ds[i] <- dist(rbind(probs, probs.new),method = 'aitchison')
  
  pval.kernel.sk[i, ] <- out$kernel[,1] 
  pval.kernel.alr[i, ] <- out$kernel[,2] 
  pval.logistic.mixed.wald[i, ] <- out$logistic.mixed[,1]
  pval.logistic.mixed.lrt[i, ] <- out$logistic.mixed[,2]
  pval.diffcyt.edger[i, ] <- out$diffcyt[,1] 
  pval.diffcyt.voom[i, ] <- out$diffcyt[,2]
}

save(pval.kernel.sk, pval.kernel.alr, pval.logistic.mixed.wald, pval.logistic.mixed.lrt,
     pval.diffcyt.edger, pval.diffcyt.voom,
     ors, ds, probs, probs.new,
     file = paste0(resultdir, "large_diff_few.RObj"))





######## Continuous predictor #########


rm(list = ls())

plotdir <- paste0(getwd(), "/plots/")
resultdir <- paste0(getwd(), "/results/LKvGK_x_continuous/")


### Source files
source("./functions/CODAK_functions.R")
source("./functions/simulation_functions.R")





####### Starting simulations #######

# Let's use 20 categories
probs <- c(0.15, 0.13, 0.07, 0.05, 0.1, 0.05, 0.05,
           0.05, 0.02, 0.1, 0.02, 0.003, 0.002, 0.004,
           0.12, 0.04, 0.004, 0.007, 0.02, 0.01)


num.cat <- length(probs)

rep <- 100 #replications
num.sim <- 100 #number of simulations per replication

ng <- 12
Xmat <- matrix(runif(2*ng, 2.5, 3.5), ncol = 1)

#number of events for the individual samples
set.seed(92486)
nvec <- sample(c(30000, 40000, 50000), 2*ng, 
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

set.seed(1)
out1 <- simulation_study_2(betarev1, betaX, Xmat, nvec, sigb,
                           num.sim = 10000, num.perm = 10000, 
                           whichtests = c("Kernel"),
                           diffcyt.method = c("edgeR", "voom"),
                           printafter = 100, 
                           Xdist = "LK",
                           storeY = T)

set.seed(1)
out2 <- simulation_study_2(betarev1, betaX, Xmat, nvec, sigb,
                           num.sim = 10000, num.perm = 10000, 
                           whichtests = c("Kernel"), 
                           diffcyt.method = c("edgeR", "voom"),
                           printafter = 100, 
                           Xdist = "gaussian",
                           storeY = T)


save(out1, out2, probs, probs.new, file = paste0(resultdir, "no_diff.RObj"))









#### Case 2: Small differences in all ####

set.seed(123)
seeds <- sample(1:10000, rep)
ors <- matrix(NA, nrow = rep, ncol = length(probs))
ds <- numeric(rep)
pval.kernel <- pval.kernel.g <-  matrix(NA, nrow = rep, ncol = num.sim)

for (i in 1:rep){
  print(paste("rep =", i))
  
  Xmat <- matrix(runif(2*ng, 2.5, 3.5), ncol = 1)
  
  betarev1 <- getES.mlogit(probs)
  probs.new <- add.effects(probs, maxeffect = 0.002, ind.maxeff = 1, n.nochange = 0)
  betarev2 <- getES.mlogit(probs.new)
  
  betaX <- sqrt(3)*as.matrix((betarev2 - betarev1)) 
  #To maintain similar variance as the binary case 
  
  set.seed(seeds[i])
  out1 <- simulation_study_2(betarev1, betaX, Xmat, nvec, sigb, 
                             num.sim = num.sim, num.perm = 10000, 
                             whichtests = c("Kernel"), Xdist = "LK",
                             diffcyt.method = c("edgeR", "voom"),
                             mirkat.method = c("Opt"),
                             printafter = 10)
  
  set.seed(seeds[i])
  out2 <- simulation_study_2(betarev1, betaX, Xmat, nvec, sigb, 
                             num.sim = num.sim, num.perm = 10000, 
                             whichtests = c("Kernel"), Xdist = "gaussian",
                             diffcyt.method = c("edgeR", "voom"),
                             mirkat.method = c("Opt"),
                             printafter = 10)
  
  orlist <- numeric(length(probs))
  for (k in 1:length(probs)){
    orlist[k] <- (probs[k]/(1-probs[k]))/(probs.new[k]/(1-probs.new[k]))
    if (orlist[k] < 1) {orlist[k] <- 1/orlist[k]}
  }
  
  ors[i, ] <- orlist
  
  ds[i] <- dist(rbind(probs, probs.new),method = 'aitchison')
  
  pval.kernel[i, ] <- out1$kernel 
  pval.kernel.g[i, ] <- out2$kernel
}

save(pval.kernel, pval.kernel.g, 
     file = paste0(resultdir, "small_diff_all.RObj"))








#### Case 3: Small differences in some ####

set.seed(123)
seeds <- sample(1:10000, rep)
ors <- matrix(NA, nrow = rep, ncol = length(probs))
ds <- numeric(rep)
pval.kernel <- pval.kernel.g <-  matrix(NA, nrow = rep, ncol = num.sim)

for (i in 1:rep){
  print(paste("rep =", i))
  
  Xmat <- matrix(runif(2*ng, 2.5, 3.5), ncol = 1)
  
  betarev1 <- getES.mlogit(probs)
  probs.new <- add.effects(probs, maxeffect = 0.002, ind.maxeff = 1, n.nochange = 10, ind.nochange = 2:11)
  betarev2 <- getES.mlogit(probs.new)
  
  betaX <- sqrt(3)*as.matrix((betarev2 - betarev1)) 
  #To maintain similar variance as the binary case 
  
  set.seed(seeds[i])
  out1 <- simulation_study_2(betarev1, betaX, Xmat, nvec, sigb, 
                             num.sim = num.sim, num.perm = 10000, 
                             whichtests = c("Kernel"), Xdist = "LK",
                             diffcyt.method = c("edgeR", "voom"),
                             mirkat.method = c("Opt"),
                             printafter = 10)
  
  set.seed(seeds[i])
  out2 <- simulation_study_2(betarev1, betaX, Xmat, nvec, sigb, 
                             num.sim = num.sim, num.perm = 10000, 
                             whichtests = c("Kernel"), Xdist = "gaussian",
                             diffcyt.method = c("edgeR", "voom"),
                             mirkat.method = c("Opt"),
                             printafter = 10)
  
  orlist <- numeric(length(probs))
  for (k in 1:length(probs)){
    orlist[k] <- (probs[k]/(1-probs[k]))/(probs.new[k]/(1-probs.new[k]))
    if (orlist[k] < 1) {orlist[k] <- 1/orlist[k]}
  }
  
  ors[i, ] <- orlist
  
  ds[i] <- dist(rbind(probs, probs.new),method = 'aitchison')
  
  pval.kernel[i, ] <- out1$kernel 
  pval.kernel.g[i, ] <- out2$kernel
}

save(pval.kernel, pval.kernel.g, 
     file = paste0(resultdir, "small_diff_some.RObj"))









#### Case 4: large differences in few ####


set.seed(123)
seeds <- sample(1:10000, rep)
ors <- matrix(NA, nrow = rep, ncol = length(probs))
ds <- numeric(rep)
pval.kernel <- pval.kernel.g <-  matrix(NA, nrow = rep, ncol = num.sim)

for (i in 1:rep){
  print(paste("rep =", i))
  
  Xmat <- matrix(runif(2*ng, 2.5, 3.5), ncol = 1)
  
  betarev1 <- getES.mlogit(probs)
  probs.new <- add.effects(probs, maxeffect = 0.004, ind.maxeff = 1, n.nochange = 15, ind.nochange = 2:16)
  betarev2 <- getES.mlogit(probs.new)
  
  betaX <- sqrt(3)*as.matrix((betarev2 - betarev1)) 
  #To maintain similar variance as the binary case 
  
  set.seed(seeds[i])
  out1 <- simulation_study_2(betarev1, betaX, Xmat, nvec, sigb, 
                             num.sim = num.sim, num.perm = 10000, 
                             whichtests = c("Kernel"), Xdist = "LK",
                             diffcyt.method = c("edgeR", "voom"),
                             mirkat.method = c("Opt"),
                             printafter = 10)
  
  set.seed(seeds[i])
  out2 <- simulation_study_2(betarev1, betaX, Xmat, nvec, sigb, 
                             num.sim = num.sim, num.perm = 10000, 
                             whichtests = c("Kernel"), Xdist = "gaussian",
                             diffcyt.method = c("edgeR", "voom"),
                             mirkat.method = c("Opt"),
                             printafter = 10)
  
  orlist <- numeric(length(probs))
  for (k in 1:length(probs)){
    orlist[k] <- (probs[k]/(1-probs[k]))/(probs.new[k]/(1-probs.new[k]))
    if (orlist[k] < 1) {orlist[k] <- 1/orlist[k]}
  }
  
  ors[i, ] <- orlist
  
  ds[i] <- dist(rbind(probs, probs.new),method = 'aitchison')
  
  pval.kernel[i, ] <- out1$kernel 
  pval.kernel.g[i, ] <- out2$kernel
}

save(pval.kernel, pval.kernel.g, 
     file = paste0(resultdir, "large_diff_few.RObj"))












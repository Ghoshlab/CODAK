

########### Simulations to understand the effect of individual components ########


rm(list = ls())


plotdir <- paste0(getwd(), "/plots/")
resultdir <- paste0(getwd(), "/results/")


library("stringr")
library("data.tree")
library("DiagrammeR")
library("Cairo")
library("RColorBrewer")
library("coda.base")


### Source files
source("./functions/CODAK_functions.R")
source("./functions/simulation_functions.R")




####### Starting simulations #######

# Let's use 20 categories
probs <- c(0.15, 0.13, 0.07, 0.05, 0.1, 0.05, 0.05,
           0.05, 0.02, 0.1, 0.02, 0.003, 0.002, 0.004,
           0.12, 0.04, 0.004, 0.007, 0.02, 0.01)


num.cat <- length(probs)

rep <- 100
num.sim <- 100

#12 observations in each group 
ng <- 12
Xmat <- matrix(rep(c(0, 1), each = ng), ncol = 1)

#number of events for the individual samples
set.seed(92486)
nvec <- sample(c(30000, 40000, 50000), 2*ng, 
               replace = T)


#covariance matrix of random effects
set.seed(92486)
sigbvec <- rep(0.2, num.cat-1) 
sigb <- diag(num.cat-1)*sigbvec^2






for (cases in 1:3){
  
  set.seed(123)
  ors <- lrs <- diffs <- matrix(NA, nrow = rep, ncol = num.cat)
  ds <- numeric(rep)
  loo.store <- wts.store <- obs.or <- obs.lr <- array(NA, c(rep,num.sim,num.cat))
  dcormat <- matrix(NA, nrow = rep, ncol = num.sim)
  
  for (i in 1:rep){
    print(paste("rep =", i))
    
    betarev1 <- getES.mlogit(probs)
    
    if (cases == 1){
      probs.new <- add.effects(probs, maxeffect = 0.002, ind.maxeff = 1, n.nochange = 0)
    }
    
    if (cases == 2){
      probs.new <- add.effects(probs, maxeffect = 0.002, ind.maxeff = 1, n.nochange = 10, ind.nochange = 2:11)
    }
    
    if (cases == 3){
      probs.new <- add.effects(probs, maxeffect = 0.004, ind.maxeff = 1, n.nochange = 15, ind.nochange = 2:16)
    }
        
    betarev2 <- getES.mlogit(probs.new)
    
    betaX <- as.matrix(betarev2 - betarev1)
    
    orlist <- (probs/(1-probs))/(probs.new/(1-probs.new))
    for (k in 1:length(probs)){
      if (orlist[k] < 1) {orlist[k] <- 1/orlist[k]}
    }
    
    ors[i, ] <- orlist
    lrs[i, ] <- abs(log(probs/probs.new))
    diffs[i, ] <- abs(probs.new - probs)
    ds[i] <- dist(rbind(probs, probs.new),method = 'aitchison')
    
    ### Data simulation ###
    for (j in 1:num.sim){
      y <- simProps(betarev1, betaX, Xmat, nvec, sigb)
      propmat <- y/rep(nvec, num.cat)
      x <- as.vector(Xmat)
      
      out <- LOO.comp(propmat, x)
      loo.store[i, j, ] <- out$loo
      dcormat[i, j] <- out$dcor
      
      out <- dcor.comp.w(propmat, x, result = "cor",
                              wt.in = NULL)
      wts.store[i, j, ] <- out$wt.out
      
      probs.obs1 <- apply(propmat[which(x == 0), ], 2, mean) 
      probs.obs2 <- apply(propmat[which(x == 1), ], 2, mean) 
      orlist <- (probs.obs1/(1-probs.obs1))/(probs.obs2/(1-probs.obs2))
      
      for (k in 1:length(orlist)){
        if (orlist[k] < 1) {orlist[k] <- 1/orlist[k]}
      }
      
      obs.or[i, j, ] <- orlist
      obs.lr[i, j, ] <- abs(log(probs.obs1/probs.obs2))
    }
  }
  
  save(loo.store, wts.store, dcormat, ors, obs.or, lrs, obs.lr,
       diffs, ds ,file = paste0(resultdir, "LOOsim_case",cases,".RObj"))
}





############ Exploration ############

loostat <- array(NA, c(3,rep,num.sim,num.cat))
ordered.loostat <- ordered.wts <- array(NA, c(3,rep,num.sim,num.cat))

rankcor.or <- rankcor.obs.or <- rankcor.diff <- 
  top.detect <- top.obs.detect <- top.detect.lr <- top.obs.detect.lr <-
  rankcor.lr <- rankcor.obs.lr <- array(NA, c(3,rep,num.sim)) 

rankcor.or.w <- rankcor.obs.or.w <- rankcor.diff.w <- 
  top.detect.w <- top.obs.detect.w <- top.detect.lr.w <- top.obs.detect.lr.w <-
  rankcor.lr.w <- rankcor.obs.lr.w <- array(NA, c(3,rep,num.sim)) 

rankcor.or_v_obsor <- top.detect.or_v_obsor <- 
  array(NA, c(3,rep,num.sim))
ordered.obsor <- array(NA, c(3,rep,num.sim,num.cat))

num.top <- 5


for (cases in 1:3){
  
  load(paste0(resultdir, "LOOsim_case",cases,".RObj"))
  
  for (i in 1:rep){
    print(i)
    for (j in 1:num.sim){
      for (k in 1:num.cat){
        loostat[cases, i, j, k] <- max(0, dcormat[i,j]-loo.store[i,j,k])
      }
      
      ordered.loostat[cases, i, j, ] <- order(loostat[cases,i, j, ], decreasing = T)
      rankcor.or[cases, i, j] <- cor(loostat[cases,i,j,], ors[i, ], 
                              method = "spearman", use = "pairwise.complete.obs")
      rankcor.obs.or[cases, i, j] <- cor(loostat[cases,i,j,], obs.or[i,j, ], 
                                  method = "spearman", use = "pairwise.complete.obs")
      rankcor.lr[cases, i, j] <- cor(loostat[cases,i,j,], lrs[i, ], 
                              method = "spearman", use = "pairwise.complete.obs")
      rankcor.obs.lr[cases, i, j] <- cor(loostat[cases,i,j,], obs.lr[i,j, ], 
                                  method = "spearman", use = "pairwise.complete.obs")
      rankcor.diff[cases, i, j] <- cor(loostat[cases,i,j,], diffs[i, ], 
                                method = "spearman", use = "pairwise.complete.obs")
      top.detect[cases, i, j] <- sum(ordered.loostat[cases,i, j, 1:num.top] %in% 
                                order(ors[i,], decreasing = T)[1:num.top])
      top.obs.detect[cases, i, j] <- sum(ordered.loostat[cases,i, j, 1:num.top] %in% 
                                    order(obs.or[i,j,], decreasing = T)[1:num.top])
      top.detect.lr[cases, i, j] <- sum(ordered.loostat[cases,i, j, 1:num.top] %in% 
                                   order(lrs[i,], decreasing = T)[1:num.top])
      top.obs.detect.lr[cases, i, j] <- sum(ordered.loostat[cases,i, j, 1:num.top] %in% 
                                       order(obs.lr[i,j,], decreasing = T)[1:num.top])
      
      ordered.wts[cases, i, j, ] <- order(wts.store[i, j, ], decreasing = T)
      rankcor.or.w[cases, i, j] <- cor(wts.store[i,j,], ors[i, ], 
                                method = "spearman", use = "pairwise.complete.obs")
      rankcor.obs.or.w[cases, i, j] <- cor(wts.store[i,j,], obs.or[i,j, ], 
                                    method = "spearman", use = "pairwise.complete.obs")
      rankcor.lr.w[cases, i, j] <- cor(wts.store[i,j,], lrs[i, ], 
                                method = "spearman", use = "pairwise.complete.obs")
      rankcor.obs.lr.w[cases, i, j] <- cor(wts.store[i,j,], obs.lr[i,j, ], 
                                    method = "spearman", use = "pairwise.complete.obs")
      rankcor.diff.w[cases, i, j] <- cor(wts.store[i,j,], diffs[i, ], 
                                  method = "spearman", use = "pairwise.complete.obs")
      top.detect.w[cases, i, j] <- sum(ordered.wts[cases, i, j, 1:num.top] %in% 
                                  order(ors[i,], decreasing = T)[1:num.top])
      top.obs.detect.w[cases, i, j] <- sum(ordered.wts[cases, i, j, 1:num.top] %in% 
                                      order(obs.or[i,j,], decreasing = T)[1:num.top])
      top.detect.lr.w[cases, i, j] <- sum(ordered.wts[cases, i, j, 1:num.top] %in% 
                                     order(lrs[i,], decreasing = T)[1:num.top])
      top.obs.detect.lr.w[cases, i, j] <- sum(ordered.wts[cases, i, j, 1:num.top] %in% 
                                         order(obs.lr[i,j,], decreasing = T)[1:num.top])
      
      
      rankcor.or_v_obsor[cases, i, j] <- cor(ors[i, ], obs.or[i,j, ], 
                                      method = "spearman", use = "pairwise.complete.obs")
      
      ordered.obsor[cases, i, j, ] <- order(obs.or[i,j, ], decreasing = T)
      top.detect.or_v_obsor[cases, i, j] <- sum(ordered.obsor[cases,i, j, 1:num.top] %in% 
                                           order(ors[i,], decreasing = T)[1:num.top])
      
    }
  }
}

save.image(paste0(resultdir, "LOOsim_final.RObj"))






hist(as.numeric(rankcor.or))
summary(as.numeric(rankcor.or))
summary(as.numeric(top.detect)) 

hist(as.numeric(rankcor.or.w))
summary(as.numeric(rankcor.or.w))
summary(as.numeric(top.detect.w)) 

summary(as.numeric(rankcor.or_v_obsor))
summary(as.numeric(top.detect.or_v_obsor))










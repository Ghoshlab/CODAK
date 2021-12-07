

################# Adding the random effects ##################

######## Exploring the effect of adding a small pseuodocount for the zeros #########


rm(list = ls())

plotdir <- paste0(getwd(), "/plots/simulations/")
resultdir <- paste0(getwd(), "/results/Pseudocounts/")


### Source files
source("./functions/CODAK_functions.R")
source("./functions/simulation_functions.R")





####### Starting simulations #######

# Let's use 20 categories
probs <- c(0.15, 0.13, 0.07, 0.05, 0.1, 0.05, 0.05,
           0.05, 0.02, 0.1, 0.02, 0.003, 0.002, 0.004,
           0.12, 0.04, 0.004, 0.007, 0.02, 0.01)

#Based on the following
#CD4+ T cells   CD8+ T cells CD27hi B cells CD27lo B cells     CD56dim NK      CD56hi NK    CD16- CD56- CD14 Monocytes CD16 Monocytes 
#0.417218430    0.245822962    0.024063537    0.113801847    0.019711307    0.002791418    0.003896324    0.121548555    0.039823659 
#pDCs           mDCs 
#0.004299898    0.007022064

num.cat <- length(probs)

rep <- 100
num.sim <- 100

#12 observations in each group 
ng <- 12
Xmat <- matrix(rep(c(0, 1), each = ng), ncol = 1)
#Xmat[13:16, ] <-  rep(0, 4)

#number of events for the individual samples
set.seed(92486)
nvec <- sample(c(30000, 40000, 50000), 2*ng, 
               replace = T)

#covariance matrix of random effects
set.seed(92486)
sigbvec <- rep(0.2, num.cat-1) 
sigb <- diag(num.cat-1)*sigbvec^2
#diag(sigb)[1:16] <- rep(0.5^2, 16)

probs[12] <- 1e-4
probs[18] <- probs[18] + 0.003 - 1e-4




#### Case 1: No differences ####

p1 <- p2 <- probs
probs.new <- probs
betarev1 <- getES.mlogit(p1)
betarev2 <- getES.mlogit(p2)

betaX <- as.matrix(betarev2 - betarev1)

set.seed(1)
out1 <- simulation_study_pseudo(betarev1, betaX, Xmat, nvec, sigb,
                          num.sim = 10000, num.perm = 10000, 
                          whichtests = "Kernel",
                          printafter = 1)

set.seed(1)
out2 <- simulation_study_pseudo(betarev1, betaX, Xmat, nvec, sigb,
                               num.sim = 10000, num.perm = 10000,  
                               whichtests = "Kernel.pseudo",
                               printafter = 1)


save(out1, out2,
    file = paste0(resultdir, "no_diff_pseudo.RObj"))



logalphavec <- seq(log10(0.00005), log10(0.05), 0.05)
alphavec <- 10^logalphavec
powervec.k <- powervec.k.p <- numeric(length(alphavec))
for(k in 1:length(alphavec)){
  powervec.k[k] <- mean(out1$pval <= alphavec[k], na.rm = T)
  powervec.k.p[k] <- mean(out2$pval <= alphavec[k], na.rm = T)
}

plot(alphavec,powervec.k, log = "x", typ = "l", pch = 15, ylim = c(0,0.1), 
     col = "red", lwd = 2, xlab = "Alpha level (log-10 scale)", ylab = "Power",
     main = "Comparison of P(Type-I error)")
lines(alphavec,powervec.k.p, typ = "l", pch = 16, col = "lightslateblue", lwd = 2)
lines(alphavec, alphavec, lwd = 2, lty = 3)










#### Case 2: Small differences in all ####

set.seed(123)
rep <- 20
num.sim <- 500
pval.kernel <- pval.kernel.pseudo <- numz.pseudo  <- matrix(NA, nrow = rep, ncol = num.sim)


for (i in 1:rep){
  print(paste("rep =", i))
  
  betarev1 <- getES.mlogit(probs)
  probs.new <- add.effects(probs, maxeffect = 0.002, ind.maxeff = 1, n.nochange = 0)
  betarev2 <- getES.mlogit(probs.new)

  betaX <- as.matrix(betarev2 - betarev1)
  
  out1 <- simulation_study_pseudo(betarev1, betaX, Xmat, nvec, sigb,
                                  num.sim = num.sim, num.perm = 10000, #1e5 and 1e6 
                                  whichtests = "Kernel",
                                  printafter = 1, 
                                  Xdist = "LK")
  
  out2 <- simulation_study_pseudo(betarev1, betaX, Xmat, nvec, sigb,
                                  num.sim = num.sim, num.perm = 10000, #1e5 and 1e6 
                                  whichtests = "Kernel.pseudo",
                                  printafter = 1, 
                                  Xdist = "LK")
  
  pval.kernel[i, ] <- out1$pval 
  pval.kernel.pseudo[i, ] <- out2$pval
  numz.pseudo[i, ] <- out2$numzeros
}

save(pval.kernel, pval.kernel.pseudo, numz.pseudo,
     file = paste0(resultdir, "small_diff_all_pseudo.RObj"))






logalphavec <- seq(log10(0.00005), log10(0.05), 0.05)
alphavec <- 10^logalphavec
powervec.k <- powervec.k.p <- numeric(length(alphavec))
for(k in 1:length(alphavec)){
  powervec.k[k] <- mean(pval.kernel <= alphavec[k], na.rm = T)
  powervec.k.p[k] <- mean(pval.kernel.pseudo <= alphavec[k], na.rm = T)
}

plot(alphavec,powervec.k, log = "x", typ = "l", pch = 15, ylim = c(0,1), 
     col = "red", lwd = 2, xlab = "Alpha level (log-10 scale)", ylab = "Power",
     main = "Comparison of P(Type-I error)")
lines(alphavec,powervec.k.p, typ = "l", pch = 16, col = "lightslateblue", lwd = 2)














#### Case 3: Small differences in some ####

set.seed(123)
rep <- 20
num.sim <- 500
pval.kernel <- pval.kernel.pseudo <- numz.pseudo  <-matrix(NA, nrow = rep, ncol = num.sim)


for (i in 1:rep){
  print(paste("rep =", i))
  
  betarev1 <- getES.mlogit(probs)
  probs.new <- add.effects(probs, maxeffect = 0.002, ind.maxeff = 1, n.nochange = 10, ind.nochange = 2:11)
  betarev2 <- getES.mlogit(probs.new)
  
  betaX <- as.matrix(betarev2 - betarev1)
  
  out1 <- simulation_study_pseudo(betarev1, betaX, Xmat, nvec, sigb,
                                  num.sim = num.sim, num.perm = 10000, #1e5 and 1e6 
                                  whichtests = "Kernel",
                                  printafter = 1, 
                                  Xdist = "LK")
  
  out2 <- simulation_study_pseudo(betarev1, betaX, Xmat, nvec, sigb,
                                  num.sim = num.sim, num.perm = 10000, #1e5 and 1e6 
                                  whichtests = "Kernel.pseudo",
                                  printafter = 1, 
                                  Xdist = "LK")
  
  pval.kernel[i, ] <- out1$pval 
  pval.kernel.pseudo[i, ] <- out2$pval
  numz.pseudo[i, ] <- out2$numzeros
}#worked until 16 repts

save(pval.kernel, pval.kernel.pseudo, numz.pseudo,
     file = paste0(resultdir, "small_diff_some_pseudo.RObj"))





logalphavec <- seq(log10(0.00005), log10(0.05), 0.05)
alphavec <- 10^logalphavec
powervec.k <- powervec.k.p <- numeric(length(alphavec))
for(k in 1:length(alphavec)){
  powervec.k[k] <- mean(pval.kernel <= alphavec[k], na.rm = T)
  powervec.k.p[k] <- mean(pval.kernel.pseudo <= alphavec[k], na.rm = T)
}

plot(alphavec,powervec.k, log = "x", typ = "l", pch = 15, ylim = c(0,1), 
     col = "red", lwd = 2, xlab = "Alpha level (log-10 scale)", ylab = "Power",
     main = "Comparison of P(Type-I error)")
lines(alphavec,powervec.k.p, typ = "l", pch = 16, col = "lightslateblue", lwd = 2)














#### Case 4: Large differences in few ####

set.seed(123)
rep <- 20
num.sim <- 500
pval.kernel <- pval.kernel.pseudo <- numz.pseudo  <-matrix(NA, nrow = rep, ncol = num.sim)


for (i in 1:rep){
  print(paste("rep =", i))
  
  betarev1 <- getES.mlogit(probs)
  probs.new <- add.effects(probs, maxeffect = 0.004, ind.maxeff = 1, n.nochange = 15, ind.nochange = 2:16)
  betarev2 <- getES.mlogit(probs.new)
  
  betaX <- as.matrix(betarev2 - betarev1)
  
  out1 <- simulation_study_pseudo(betarev1, betaX, Xmat, nvec, sigb,
                                  num.sim = num.sim, num.perm = 10000, #1e5 and 1e6 
                                  whichtests = "Kernel",
                                  printafter = 1, 
                                  Xdist = "LK")
  
  out2 <- simulation_study_pseudo(betarev1, betaX, Xmat, nvec, sigb,
                                  num.sim = num.sim, num.perm = 10000, #1e5 and 1e6 
                                  whichtests = "Kernel.pseudo",
                                  printafter = 1, 
                                  Xdist = "LK")
  
  pval.kernel[i, ] <- out1$pval 
  pval.kernel.pseudo[i, ] <- out2$pval
  numz.pseudo[i, ] <- out2$numzeros
}

save(pval.kernel, pval.kernel.pseudo, numz.pseudo,
     file = paste0(resultdir, "large_diff_few_pseudo.RObj"))


mean(pval.kernel < 0.05)
mean(pval.kernel.pseudo < 0.05)




logalphavec <- seq(log10(0.00005), log10(0.05), 0.05)
alphavec <- 10^logalphavec
powervec.k <- powervec.k.p <- numeric(length(alphavec))
for(k in 1:length(alphavec)){
  powervec.k[k] <- mean(pval.kernel <= alphavec[k], na.rm = T)
  powervec.k.p[k] <- mean(pval.kernel.pseudo <= alphavec[k], na.rm = T)
}

plot(alphavec,powervec.k, log = "x", typ = "l", pch = 15, ylim = c(0,1), 
     col = "red", lwd = 2, xlab = "Alpha level (log-10 scale)", ylab = "Power",
     main = "Comparison of P(Type-I error)")
lines(alphavec,powervec.k.p, typ = "l", pch = 16, col = "lightslateblue", lwd = 2)








Cairo(file = paste0(plotdir, "pseudocount.pdf"), typ = "pdf", dpi = 95,
      height = 1250, width = 1150, pointsize = 12)

minalpha = 0.0005
titlevec <- c("Case 1: No difference",
              "Case 2: Small difference in \n all (20) cell types",
              "Case 3: Small difference in \n some (10) cell types",
              "Case 4: Larger difference in \n a few (5) cell types")

op=par(oma=c(8.5,0,0,0))
op = par(xpd=NA)
par(mfrow = c(2,2))

for (ii in 1:4){
  if (ii == 1){
    load(paste0(resultdir, "no_diff_pseudo.RObj"))
    
    logalphavec <- c(seq(log10(minalpha), log10(0.05), 0.05), log10(0.05))
    alphavec <- 10^logalphavec
    powervec.k <- powervec.k.p <- numeric(length(alphavec))
    for(k in 1:length(alphavec)){
      powervec.k[k] <- mean(out1$pval <= alphavec[k], na.rm = T)
      powervec.k.p[k] <- mean(out2$pval <= alphavec[k], na.rm = T)
    }
    
    plot(alphavec,powervec.k, log = "x", typ = "l", pch = 15, ylim = c(0,0.1), 
         col = "red", lwd = 1.5, xlab = "Alpha level (log-10 scale)", ylab = "Power",
         main = titlevec[ii])
    lines(alphavec,powervec.k.p, typ = "l", pch = 16, col = "lightslateblue", lwd = 1.5)
    lines(alphavec, alphavec, lwd = 1, lty = 2, col = "grey25")
    lines(alphavec, 2*alphavec, lwd = 1, lty = 2, col = "darkgrey")
  }
  
  if (ii == 2){
    load(paste0(resultdir, "small_diff_all_pseudo.RObj"))
  }
  
  if (ii == 3){
    load(paste0(resultdir, "small_diff_some_pseudo.RObj"))
  }
  
  if (ii == 4){
    load(paste0(resultdir, "large_diff_few_pseudo.RObj"))
  }
  
  if (ii %in% 2:4){
    logalphavec <- c(seq(log10(minalpha), log10(0.05), 0.05), log10(0.05))
    alphavec <- 10^logalphavec
    powervec.k <- powervec.k.p <- numeric(length(alphavec))
    for(k in 1:length(alphavec)){
      powervec.k[k] <- mean(pval.kernel <= alphavec[k], na.rm = T)
      powervec.k.p[k] <- mean(pval.kernel.pseudo <= alphavec[k], na.rm = T)
    }
    
    plot(alphavec,powervec.k, log = "x", typ = "l", pch = 15, ylim = c(0,1), 
         col = "red", lwd = 1.5, xlab = "Alpha level (log-10 scale)", ylab = "Power",
         main = titlevec[ii])
    lines(alphavec,powervec.k.p, typ = "l", pch = 16, col = "lightslateblue", lwd = 1.5)
  }
}

legend("bottomright",
       legend = c("CODAK (no zeros)", "CODAK + pseudocounts"),  
       col= c("red", "lightslateblue"), box.lty = 0,
       lty= 1, bg="white", horiz=F,  
       pch= NA, inset=c(0.9,-0.8),
       cex = 1.3)

dev.off()










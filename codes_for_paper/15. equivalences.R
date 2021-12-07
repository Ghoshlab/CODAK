


######## Benchmarking - Comparing time to run for different methods #########


rm(list = ls())

plotdir <- paste0(getwd(), "/plots/simulations/")
resultdir <- paste0(getwd(), "/results/Equivalence/")


### Source files
source("./functions/CODAK_functions.R")
source("./functions/simulation_functions.R")





####### Starting simulations #######

# Let's use 20 categories
probs <- c(0.15, 0.13, 0.07, 0.05, 0.1, 0.05, 0.05,
           0.05, 0.02, 0.1, 0.02, 0.003, 0.002, 0.004,
           0.12, 0.04, 0.004, 0.007, 0.02, 0.01)


num.cat <- length(probs)

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






#### Simulations from null ####

p1 <- p2 <- probs
probs.new <- probs
betarev1 <- getES.mlogit(p1)
betarev2 <- getES.mlogit(p2)

betaX <- as.matrix(betarev2 - betarev1)
beta <- betarev1

num.perm <- 10000
num.sim <- 1000
printafter <- 10

Xdist = "categorical"
x <- as.vector(Xmat)

proplist <- ylist <- list()
pval.k <- pval.p <- pval.m <- pval.e <- pval.d <- numeric(num.sim)
r2.k <- r2.p <- r2.m <- numeric(num.sim)

set.seed(1)
for (i in 1:num.sim){
  print(i)
  
  y <- simProps(beta, betaX, Xmat, nvec, sigb)
  propmat <- y/rep(nvec, num.cat)
  proplist[[i]] <- propmat
  ylist[[i]] <- y
}




### CODAK ###
#set.seed(1)
for (i in 1:num.sim){
  if (floor(i/printafter) == i/printafter){print(i)}
  out <- try(dcor.comp(proplist[[i]], x, num.perm, Xdist = Xdist, gaussian = F,
                       result = c("test", "cor")), silent = T)
  pval.k[i] <- out$pval
  r2.k[i] <- out$dcor
}



### PERMANOVA ###
#set.seed(1)
for (i in 1:num.sim){
  if (floor(i/printafter) == i/printafter){print(i)}
  d <- dist(proplist[[i]], method = 'aitchison')
  
  pmfit <- try(adonis(d ~ x, permutations = num.perm), silent = T)
  if (class(pmfit) == "try-error"){
    pval.p[i] <- NA
  }else{
    pval.p[i] <- pmfit$aov.tab$"Pr(>F)"[1]
  }
  r2.p[i] <- pmfit$aov.tab$R2[1]
}


### MiRKAT-AD ###
#set.seed(1)
for (i in 1:num.sim){
  if (floor(i/printafter) == i/printafter){print(i)}
  
  typ = "C"
  if (Xdist == "categorical"){typ = "D"}
  
  temp <- try(mirkat.test(proplist[[i]], x, method = "AD", typ = typ, num.perm = num.perm), silent = T)

  if (class(pmfit) == "try-error"){
    pval.m[i] <- NA
  }else{
    pval.m[i] <- temp$p_values
  }
  
  r2.m[i] <- sqrt(temp$R2)
}


### CODAK - ED ###
#set.seed(1)
for (i in 1:num.sim){
  if (floor(i/printafter) == i/printafter){print(i)}
  out <- try(dcor.comp(proplist[[i]], x, num.perm, Xdist = Xdist, gaussian = F,
                       compdist = "euclidean", result = c("test", "cor")), silent = T)
  pval.e[i] <- out$pval
}



### Original dcor ###
#set.seed(1)
for (i in 1:num.sim){
  if (floor(i/printafter) == i/printafter){print(i)}
  pval.d[i] <- dcov.test(proplist[[i]], x, R = num.perm)$p.value
}



save(proplist, ylist, pval.k, pval.e, pval.d, pval.p, pval.m, 
     r2.k, r2.p, r2.m, file = paste0(resultdir, "equivalence_null.RObj"))










#### Simulations from alternative ####

set.seed(123)

betarev1 <- getES.mlogit(probs)
probs.new <- add.effects(probs, maxeffect = 0.002, ind.maxeff = 1, n.nochange = 0)
betarev2 <- getES.mlogit(probs.new)
betaX <- as.matrix(betarev2 - betarev1)

betaX <- as.matrix(betarev2 - betarev1)
beta <- betarev1

num.perm <- 10000
num.sim <- 1000
printafter <- 10

Xdist = "categorical"
x <- as.vector(Xmat)

proplist.a <- ylist.a <- list()
pval.k.a <- pval.p.a <- pval.m.a <- pval.e.a <- pval.d.a <- numeric(num.sim)
r2.k.a <- r2.p.a <- r2.m.a <- numeric(num.sim)

set.seed(1)
for (i in 1:(0.4*num.sim)){
  print(i)
  
  y <- simProps(beta, betaX, Xmat, nvec, sigb)
  propmat <- y/rep(nvec, num.cat)
  proplist.a[[i]] <- propmat
  ylist.a[[i]] <- y
}




set.seed(123)

betarev1 <- getES.mlogit(probs)
probs.new <- add.effects(probs, maxeffect = 0.001, ind.maxeff = 1, n.nochange = 0) #smaller effects
betarev2 <- getES.mlogit(probs.new)
betaX <- as.matrix(betarev2 - betarev1)

betaX <- as.matrix(betarev2 - betarev1)
beta <- betarev1

set.seed(1)
for (i in (0.4*num.sim+1):(0.7*num.sim)){
  print(i)
  
  y <- simProps(beta, betaX, Xmat, nvec, sigb)
  propmat <- y/rep(nvec, num.cat)
  proplist.a[[i]] <- propmat
  ylist.a[[i]] <- y
}


set.seed(123)

betarev1 <- getES.mlogit(probs)
probs.new <- add.effects(probs, maxeffect = 0.004, ind.maxeff = 1, n.nochange = 0) #larger effects
betarev2 <- getES.mlogit(probs.new)
betaX <- as.matrix(betarev2 - betarev1)

betaX <- as.matrix(betarev2 - betarev1)
beta <- betarev1

set.seed(1)
for (i in (0.7*num.sim+1):(num.sim)){
  print(i)
  
  y <- simProps(beta, betaX, Xmat, nvec, sigb)
  propmat <- y/rep(nvec, num.cat)
  proplist.a[[i]] <- propmat
  ylist.a[[i]] <- y
}







### CODAK ###
#set.seed(1)
for (i in 1:num.sim){
  if (floor(i/printafter) == i/printafter){print(i)}
  out <- try(dcor.comp(proplist.a[[i]], x, num.perm, Xdist = Xdist, gaussian = F,
                       result = c("test", "cor")), silent = T)
  pval.k.a[i] <- out$pval
  r2.k.a[i] <- out$dcor
}



### PERMANOVA ###
#set.seed(1)
for (i in 1:num.sim){
  if (floor(i/printafter) == i/printafter){print(i)}
  d <- dist(proplist.a[[i]], method = 'aitchison')
  
  pmfit <- try(adonis(d ~ x, permutations = num.perm), silent = T)
  if (class(pmfit) == "try-error"){
    pval.p.a[i] <- NA
  }else{
    pval.p.a[i] <- pmfit$aov.tab$"Pr(>F)"[1]
  }
  r2.p.a[i] <- pmfit$aov.tab$R2[1]
}


### MiRKAT-AD ###
#set.seed(1)
for (i in 1:num.sim){
  if (floor(i/printafter) == i/printafter){print(i)}
  
  typ = "C"
  if (Xdist == "categorical"){typ = "D"}
  
  temp <- try(mirkat.test(proplist.a[[i]], x, method = "AD", typ = typ, num.perm = num.perm), silent = T)
  
  if (class(pmfit) == "try-error"){
    pval.m.a[i] <- NA
  }else{
    pval.m.a[i] <- temp$p_values
  }
  
  r2.m.a[i] <- sqrt(temp$R2)
}


### CODAK - ED ###
#set.seed(1)
for (i in 1:num.sim){
  if (floor(i/printafter) == i/printafter){print(i)}
  out <- try(dcor.comp(proplist.a[[i]], x, num.perm, Xdist = Xdist, gaussian = F,
                       compdist = "euclidean", result = c("test", "cor")), silent = T)
  pval.e.a[i] <- out$pval
}



### Original dcor ###
#set.seed(1)
for (i in 1:num.sim){
  if (floor(i/printafter) == i/printafter){print(i)}
  pval.d.a[i] <- dcov.test(proplist.a[[i]], x, R = num.perm)$p.value
}



save(proplist.a, ylist.a, pval.k.a, pval.e.a, pval.d.a, pval.p.a, pval.m.a, 
     r2.k.a, r2.p.a, r2.m.a, file = paste0(resultdir, "equivalence_alt.RObj"))










######## Create plots ########

library(scales)

Cairo(file = paste0(plotdir, "equivalence.pdf"), typ = "pdf", dpi = 80, pointsize = 18,
      height = 800, width = 1300)
par(mfrow = c(2,3))

plot(c(pval.k, pval.k.a), c(pval.m, pval.m.a), xlab = "CODAK p-value", ylab = "MirKAT p-value",
     col = alpha(c(rep("black", num.sim),rep("darkred", num.sim)), 0.4))
text(0.25, 0.95, paste0("r = ", round(cor(c(pval.k, pval.k.a), c(pval.m, pval.m.a)), 6)), cex = 1.2)
abline(0,1, col = "grey", lwd = 1.5)

plot(c(pval.k, pval.k.a), c(pval.p, pval.p.a), xlab = "CODAK p-value", ylab = "PERMANOVA p-value",
     col = alpha(c(rep("black", num.sim),rep("darkred", num.sim)), 0.4))
text(0.25, 0.95, paste0("r = ", round(cor(c(pval.k, pval.k.a), c(pval.p, pval.p.a)), 6)), cex = 1.2)
abline(0,1, col = "grey", lwd = 1.5)

plot(c(pval.m, pval.m.a), c(pval.p, pval.p.a), xlab = "MirKAT p-value", ylab = "PERMANOVA p-value",
     col = alpha(c(rep("black", num.sim),rep("darkred", num.sim)), 0.4))
text(0.25, 0.95, paste0("r = ", round(cor(c(pval.k, pval.k.a), c(pval.m, pval.m.a)), 6)), cex = 1.2)
abline(0,1, col = "grey", lwd = 1.5)



plot(c(r2.k, r2.k.a), c(r2.m, r2.m.a), xlab = expression(CODAK~~R^2), ylab = expression(MiRKAT~~R),
     col = alpha(c(rep("black", num.sim),rep("darkred", num.sim)), 0.4),
     xlim = c(0,0.9), ylim = c(0,0.9))
text(0.25, 0.85, paste0("r = ", round(cor(c(r2.k, r2.k.a), c(r2.m, r2.m.a)), 6)), cex = 1.2)
abline(0,1, col = "grey", lwd = 1.5)

plot(c(r2.k, r2.k.a), c(r2.p, r2.p.a), xlab = expression(CODAK~~R^2), ylab = expression(PERMANOVA~~R^2),
     col = alpha(c(rep("black", num.sim),rep("darkred", num.sim)), 0.4),
     xlim = c(0,0.9), ylim = c(0,0.9))
text(0.25, 0.85, paste0("r = ", round(cor(c(r2.k, r2.k.a), c(r2.p, r2.p.a)), 6)), cex = 1.2)
abline(0,1, col = "grey", lwd = 1.5)

plot(c(r2.m, r2.m.a), c(r2.p, r2.p.a), xlab = expression(MiRKAT~~R), ylab = expression(PERMANOVA~~R^2),
     col = alpha(c(rep("black", num.sim),rep("darkred", num.sim)), 0.4),
     xlim = c(0,0.9), ylim = c(0,0.9))
text(0.25, 0.85, paste0("r = ", round(cor(c(r2.m, r2.m.a), c(r2.p, r2.p.a)), 6)), cex = 1.2)
abline(0,1, col = "grey", lwd = 1.5)

mtext(expression(Comparison~~of~~P-value), side = 3, line = -2, outer = T)
mtext(expression(Comparison~~of~~R^2), side = 3, line = -25, outer = T)
dev.off()






Cairo(file = paste0(plotdir, "equivalence_dcor_ed.pdf"), typ = "pdf", dpi = 80, pointsize = 18,
      height = 600, width = 600)

plot(c(pval.e, pval.e.a), c(pval.d, pval.d.a), xlab = "CODAK-ED p-value", ylab = "dcor p-value",
     col = alpha(c(rep("black", num.sim),rep("darkred", num.sim)), 0.4), cex = 0.7,
     main = "Original dcor vs CODAK-ED")
text(0.25, 0.95, paste0("r = ", round(cor(c(pval.e, pval.e.a), c(pval.d, pval.d.a)), 6)), cex = 1.2)
abline(0,1, col = "grey", lwd = 1.5)
dev.off()






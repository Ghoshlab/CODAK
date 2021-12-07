


######## Benchmarking - Comparing time to run for different methods #########


rm(list = ls())

plotdir <- paste0(getwd(), "/plots/")
resultdir <- paste0(getwd(), "/results/TimeToRun/")


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






#### Case 1: No differences ####

p1 <- p2 <- probs
probs.new <- probs
betarev1 <- getES.mlogit(p1)
betarev2 <- getES.mlogit(p2)

betaX <- as.matrix(betarev2 - betarev1)
beta <- betarev1
num.perm <- 1000
lmer.nump <- 1000
num.sim <- 100
Xdist = "categorical"
x <- as.vector(Xmat)
proplist <- ylist <- list()

set.seed(1)
for (i in 1:num.sim){
  print(i)
  
  y <- simProps(beta, betaX, Xmat, nvec, sigb)
  propmat <- y/rep(nvec, num.cat)
  proplist[[i]] <- propmat
  ylist[[i]] <- y
}


### CODAK ###
set.seed(1)
system.time({
for (i in 1:num.sim){
  print(i)
  out <- dcor.comp(proplist[[i]], x, num.perm, Xdist = Xdist, gaussian = F,
                       result = c("test", "cor"))
  #d <- dist(proplist[[i]], method = "aitchison")
  #out <- adonis(d ~ x, permutations = num.perm)
  
  #out <- mirkat.test(ylist[[i]], x, method = "AD", typ = "D", num.perm = num.perm)
}
})
#2.91



### Logistic LRT ###
set.seed(1)
system.time({
  for (i in 1:num.sim){
    print(i)
    
    ptemp <- rep(NA, num.cat) 
    ptemp.mixed <- matrix(NA, nrow = num.cat, ncol = 1)
    
    for (k in 1:num.cat){ 
      y.bin <- y[, k]
      id <- 1:length(x)
      fitmixed <- try(glmer(y.bin/nvec ~ x + (1|id), weights = nvec, 
                            family = binomial(link = "logit")))
      out <- try(drop1(fitmixed, test="Chisq")$"Pr(Chi)"[2])
      if (class(fitmixed) == "try-error"){
        ptemp.mixed[k, 1]<- NA
      }else{
        ptemp.mixed[k, 1] <- out
      }
    }  
    
    pvals <- apply(ptemp.mixed, 2, function(w){
      return(min(p.adjust(w, method = "bonferroni")))
    })
  }
})
#131.15



### Logistic LRT Permutation ###
set.seed(1)
system.time({
  for (i in 1:num.sim){
    print(i)
    
    ptemp <- rep(NA, num.cat) 
    ptemp.mixed <- matrix(NA, nrow = num.cat, ncol = 1)
    
    for (k in 1:num.cat){ 
      y.bin <- y[, k]
      id <- 1:length(x)
      fitmixed <- try(glmer(y.bin/nvec ~ x + (1|id), weights = nvec, 
                            family = binomial(link = "logit")))
      
      lrstat <- drop1(fitmixed, test="Chisq")$LRT[2]
      lrstat.perm <- numeric(lmer.nump)
      
      for (p in 1:lmer.nump){
        xperm <- sample(x)
        fitmixed <- glmer(y.bin/nvec ~ xperm + (1|id), weights = nvec, 
                          family = binomial(link = "logit"))
        lrstat.perm[p] <- drop1(fitmixed, test="Chisq")$LRT[2]
      }
      
      ptemp.mixed[k, 1] <- (sum(lrstat.perm >= lrstat)+1)/(lmer.nump+1)
    }  
    
    pvals <- apply(ptemp.mixed, 2, function(w){
      return(min(p.adjust(w, method = "bonferroni")))
    })
  }
})
#6522.54 for num.sim = 5




### diffcyt ###
diffcyt.method <- "edgeR"
set.seed(1)
system.time({
  for (i in 1:num.sim){
    print(i)
    
    temp <- diffcyt.test(ylist[[i]], x, method = diffcyt.method)
    pval <- min(p.adjust(temp, method = "bonferroni"))
  }
})
#3.76






#save(out, probs, probs.new, file = paste0(resultdir, "no_diff.RObj"))






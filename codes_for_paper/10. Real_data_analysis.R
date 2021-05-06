
########### Real data analysis #############


rm(list = ls())


plotdir <- paste0(getwd(), "/plots/")
resultdir <- paste0(getwd(), "/results/")


library("Cairo")
library("RColorBrewer")

### Source files
source("./functions/CODAK_functions.R")



########### Loading data ###########

load("./data/propready.RObj")




############# Kernel dcov method for testing and plot ###############

d <- dist(propmat, method = 'aitchison')
alldist <- as.numeric(d)
K <- as.matrix(exp(-d/median(alldist)))


##### Distance correlation test for groups ######

d2 <- outer(as.numeric(casecontrolvec == "SLE"), as.numeric(casecontrolvec == "SLE"), 
            FUN = function(x,y){abs(x-y)})
K2 <- exp(-d2)

#Do permutation test
set.seed(123)

a <- Sys.time()
dcor.comp(propmat, casecontrolvec, 1e6, result = c("cor", "test")) #0.6967, 1e-6
Sys.time() - a


##### Distance correlation test for conditions, don't report this ######
d2 <- outer(as.numeric(stimcond == "T0"), as.numeric(stimcond == "T0"),  #should be T0 and T0
            FUN = function(x,y){abs(x-y)})
K2 <- exp(-d2)


#Do permutation
set.seed(123)
dcor.comp(propmat, as.numeric(stimcond=="T0"), 1e6, result = c("cor", "test")) 
#0.0839, 0.9260


#### Covariate adjustment ####
set.seed(123)
dcor.comp.cov(propmat, casecontrolvec, 
              zmat =  as.matrix(as.numeric(stimcond=="T0"), ncol = 1), 
              1e6, result = c("cor", "test"), covadj = "SK")







####### dcor test for T0 #######
hold <- expression(which(infotable$Group %in% c("SLE", "Healthy Control") &
                                infotable$Condition == "T0"))
propmat.curr <- propmat[eval(hold), ]
infotable.curr <- infotable[eval(hold), ]
xvec <- infotable.curr$Group

set.seed(123)
dcor.comp(propmat.curr, xvec, num.perm = 1e6, result = c("cor", "test"))
#0.66, 0.0006

out <- LOO.comp(propmat.curr, xvec, Xdist = "categorical")
loo.t0 <- out$dcor - out$loo
loo.t0[which(loo.t0 < 0)] <- 0




####### dcor test for T6 #######
hold <- expression(which(infotable$Group %in% c("SLE", "Healthy Control") &
                           infotable$Condition == "T6"))
propmat.curr <- propmat[eval(hold), ]
infotable.curr <- infotable[eval(hold), ]
xvec <- infotable.curr$Group

set.seed(123)
dcor.comp(propmat.curr, xvec, num.perm = 1e6, result = c("cor", "test"))
#0.71, 0.0006

out <- LOO.comp(propmat.curr, xvec, Xdist = "categorical")
loo.t6 <- out$dcor - out$loo
loo.t6[which(loo.t6 < 0)] <- 0




######### Loo plot ##########

Cairo(paste0(plotdir, "realdata/looplot.pdf"),
      typ = "pdf", dpi = 85, height = 800, width = 900)
par(mfrow = c(2,1), mar = c(6,3,4,1))
barplot(loo.t0, names.arg = colnames(propmat), las = 2, col = "darkorange", 
        cex.names = 0.85)
legend("topright", legend = "T0     ", inset = c(0.1, 0), cex = 1.2)
barplot(loo.t6, names.arg = colnames(propmat), las = 2, col = "darkorange",
        cex.names = 0.85)
legend("topright", legend = "T6     ", inset = c(0.1, 0), cex = 1.2)
dev.off()






#### Kernel Plot for motivation ####

Cairo(file = paste0(plotdir,"realdata/kernelplot.pdf"), typ = "pdf", dpi = 85, height = 400, width = 900)

par(mfrow = c(1,2))

xx <- c(K[which(casecontrolvec == "SLE"), which(casecontrolvec == "SLE")], 
        K[which(casecontrolvec == "Healthy Control"), which(stimcond == "Healthy Control")])
yy <- K[which(casecontrolvec == "SLE"), which(casecontrolvec == "Healthy Control")]
xx <- xx[-which(xx == 0)]
plot(density(xx),
     col = "red", xlim = c(0,0.9), ylim = c(0, 5), xlab = "Similarity", 
     main = bquote("Comparison across disease groups"), # (p-value <"~10^-6~")"), 
     font.main = 1, cex.main = 0.9)
lines(density(yy), 
      col = "blue")
legend("topright", legend = c("Same group", "Different group"), col = c("red", "blue"), 
       lty = c(1, 1), cex = 0.7)



xx <- c(K[which(stimcond == "T0"), which(stimcond == "T0")], K[which(stimcond == "T6"), which(stimcond == "T6")])
yy <- K[which(stimcond == "T0"), which(stimcond == "T6")]
xx <- xx[-which(xx == 0)]

plot(density(xx), col = "red", xlim = c(0,0.9), ylim = c(0, 4), xlab = "Similarity", 
     main = "Comparison across conditions", # (p-value = 0.9266)", 
     font.main = 1, cex.main = 0.9)
lines(density(yy), col = "blue")
legend("topright", legend = c("Same condition", "Different condition"), col = c("red", "blue"), 
       lty = c(1, 1), cex = 0.7)

dev.off()









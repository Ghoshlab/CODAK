


################# Code for creating effect size plots ####################


rm(list = ls())

library("Cairo")
library("plotfunctions")




plotdir <- paste0(getwd(), "/plots/")
resultdir <- paste0(getwd(), "/results/")


restype <- "KvLM_binary"



effsize_powerplot <- function(compare.pvals = c(expression(pval.kernel),
                                                expression(pval.logistic.mixed.lrt),
                                                expression(pval.diffcyt.voom)), 
                              alpha = 0.001, 
                              restype = "KvLM_binary", 
                              cases = "all", gradient = F){
  
  #compare.pvals: expression vector with name of the pvalue objects to compare
                  # Should always be of length 3
  #alpha: level of test
  #restype: simulation study type
  #cases: "all" if combined result only, "allplus" if combined and individual
  
  pathvec <- c(paste0(resultdir, restype, "/no_diff.RObj"),
               paste0(resultdir, restype, "/small_diff_all.RObj"),
               paste0(resultdir, restype, "/small_diff_some.RObj"),
               paste0(resultdir, restype, "/large_diff_few.RObj"))
  
  power1 <- power2 <- power3 <- 
    max.logor <- dsmat <- matrix(NA, nrow = 3, ncol = 100)
  
  for (j in 1:3){
    ii <- j + 1
    load(pathvec[ii])
    power1[j,] <- apply(eval(compare.pvals[1]), 1, function(x){mean(x <= alpha, na.rm = T)})
    power2[j,] <- apply(eval(compare.pvals[2]), 1, function(x){mean(x <= alpha, na.rm = T)})
    power3[j,] <- apply(eval(compare.pvals[3]), 1, function(x){mean(x <= alpha, na.rm = T)})
    
    max.logor[j,] <- apply(abs(log10(ors)), 1, max)
    dsmat[j,] <- ds
  }
  
  
  grad <- c("nograd", "grad")[1+gradient]
  filename <- paste0(plotdir, "simulations/", 
                     restype, "_", alpha, "_", cases, "_", grad, ".pdf")
  
  
  if (cases == "allplus"){
    Cairo(height = 1600, width = 900, dpi = 85, typ = "pdf", file = filename)
    par(mfrow = c(4,2))
    
    for (l in 1:3){
      
      power1.all <- as.numeric(power1[l,])
      power2.all <- as.numeric(power2[l,])
      power3.all <- as.numeric(power3[l,])
      
      #Logistic
      better <- numeric(length(ds))
      better[which(power1.all > power2.all)] <- 1
      better[which(power1.all < power2.all)] <- 2
      better[which(power1.all == power2.all)] <- 3
      colvec <- c("red", "blue", "black")[better]
      
      if (gradient == T){
        r <- 0.33
        powerdiff <- abs(power1.all - power3.all)
        k <- 1/max(powerdiff^r)
        colvec[which(better == 1)] <- 
          rgb(k*powerdiff[which(better == 1)]^r, 0, 0)
        colvec[which(better == 2)] <- 
          rgb(0, 0, k*powerdiff[which(better == 2)]^r)
      }
      
      colvec <- alpha(colvec, 0.5)
      plot(as.numeric(dsmat[l,]), as.numeric(max.logor[l,]), 
           col = colvec, pch = 16, xlim = c(0.15,1), ylim = c(0,0.5),
           xlab = "Aitchison Distance", ylab = "Maximum log(OR)",
           main = "CODAK vs GLMM", cex.main = 0.85)
      
      if (gradient == F){
        legend("bottomright", legend = c("CODAK", "GLMM", "tie"), 
               fill = alpha(c("red", "blue", "black"), 0.5), box.lty = 0, cex = 0.7)
      }
      
      if (gradient == T){
        gradcol <- alpha(c(rgb(0, 0,rev(seq(0, max(k*powerdiff[which(better == 2)]^r),0.01))),
                           rgb(seq(0, max(k*powerdiff[which(better == 2)]^r),0.01), 0, 0)), 0.5)
        gradientLegend(valRange=c(-1, 1), side=1, inside=TRUE, cex = 0.5, 
                       border.col = NA, col = gradcol,
                       pos=c(.75,-0.005,0.95,.02), coords=TRUE, pos.num = c(2,3,4))
        legend(0.57, 0.07, "Power difference: CODAK - GLMM", 
               box.col = "white", bg = "white", cex = 0.5)
      }
      
      
      #Diffcyt
      better <- numeric(length(ds))
      better[which(power1.all > power3.all)] <- 1
      better[which(power1.all < power3.all)] <- 2
      better[which(power1.all == power3.all)] <- 3
      colvec <- c("red", "darkgreen", "black")[better]
      
      if (gradient == T){
        r <- 0.33
        powerdiff <- abs(power1.all - power3.all)
        k <- 1/max(powerdiff^r)
        colvec[which(better == 1)] <- 
          rgb(k*powerdiff[which(better == 1)]^r, 0, 0)
        colvec[which(better == 2)] <- 
          rgb(0, k*powerdiff[which(better == 2)]^r, 0)
      }  
      
      colvec <- alpha(colvec, 0.5)
      plot(as.numeric(dsmat[l,]), as.numeric(max.logor[l,]), 
           col = colvec, pch = 16, xlim = c(0.15,1), ylim = c(0,0.5),
           xlab = "Aitchison Distance", ylab = "Maximum log(OR)",
           main = "CODAK vs diffcyt", cex.main = 0.85)
      
      if (gradient == F){
        legend("bottomright", legend = c("CODAK", "diffcyt", "tie"), 
               fill = alpha(c("red", "darkgreen", "black"), 0.5), 
               box.lty = 0, cex = 0.7)
      }
      
      if (gradient == T){
        gradcol <- alpha(c(rgb(0,rev(seq(0, max(k*powerdiff[which(better == 2)]^r),0.01)), 0),
                           rgb(seq(0, max(k*powerdiff[which(better == 2)]^r),0.01), 0, 0)), 0.5)
        gradientLegend(valRange=c(-1, 1), side=1, inside=TRUE, cex = 0.5, 
                       border.col = NA, col = gradcol,
                       pos=c(.75,-0.005,0.95,.02), coords=TRUE, pos.num = c(2,3,4))
        legend(0.57, 0.07, "Power difference: CODAK - diffcyt", 
               box.col = "white", bg = "white", cex = 0.5)
      }
    }
  }
  
  if (cases == "all"){
    Cairo(height = 400, width = 900, dpi = 85, typ = "pdf", file = filename)
    par(mfrow = c(1,2))
  }  
  
  
  power1.all <- as.numeric(power1)
  power2.all <- as.numeric(power2)
  power3.all <- as.numeric(power3)
  
  
  #Logistic
  better <- numeric(length(ds))
  better[which(power1.all > power2.all)] <- 1
  better[which(power1.all < power2.all)] <- 2
  better[which(power1.all == power2.all)] <- 3
  colvec <- c("red", "blue", "black")[better]
  
  if (gradient == T){
    r <- 0.33
    powerdiff <- abs(power1.all - power3.all)
    k <- 1/max(powerdiff^r)
    colvec[which(better == 1)] <- 
      rgb(k*powerdiff[which(better == 1)]^r, 0, 0)
    colvec[which(better == 2)] <- 
      rgb(0, 0, k*powerdiff[which(better == 2)]^r)
  }
  
  colvec <- alpha(colvec, 0.5)
  plot(as.numeric(dsmat), as.numeric(max.logor), 
       col = colvec, pch = 16, xlim = c(0.15,1), ylim = c(0,0.5),
       xlab = "Aitchison Distance", ylab = "Maximum log(OR)",
       main = "CODAK vs GLMM", cex.main = 0.85)
  
  if (gradient == F){
    legend("bottomright", legend = c("CODAK", "GLMM", "tie"), 
           fill = alpha(c("red", "blue", "black"), 0.5), box.lty = 0, cex = 0.7)
  }
  
  if (gradient == T){
    gradcol <- alpha(c(rgb(0, 0,rev(seq(0, max(k*powerdiff[which(better == 2)]^r),0.01))),
                       rgb(seq(0, max(k*powerdiff[which(better == 2)]^r),0.01), 0, 0)), 0.5)
    gradientLegend(valRange=c(-1, 1), side=1, inside=TRUE, cex = 0.5, 
                   border.col = NA, col = gradcol,
                   pos=c(.75,-0.005,0.95,.02), coords=TRUE, pos.num = c(2,3,4))
    legend(0.57, 0.07, "Power difference: CODAK - GLMM", 
           box.col = "white", bg = "white", cex = 0.5)
  }
  
    
  #Diffcyt
  better <- numeric(length(ds))
  better[which(power1.all > power3.all)] <- 1
  better[which(power1.all < power3.all)] <- 2
  better[which(power1.all == power3.all)] <- 3
  colvec <- c("red", "darkgreen", "black")[better]
  
  if (gradient == T){
    r <- 0.33
    powerdiff <- abs(power1.all - power3.all)
    k <- 1/max(powerdiff^r)
    colvec[which(better == 1)] <- 
      rgb(k*powerdiff[which(better == 1)]^r, 0, 0)
    colvec[which(better == 2)] <- 
      rgb(0, k*powerdiff[which(better == 2)]^r, 0)
  }  
  
  colvec <- alpha(colvec, 0.5)
  plot(as.numeric(dsmat), as.numeric(max.logor), 
       col = colvec, pch = 16, xlim = c(0.15,1), ylim = c(0,0.5),
       xlab = "Aitchison Distance", ylab = "Maximum log(OR)",
       main = "CODAK vs diffcyt", cex.main = 0.85)
  
  if (gradient == F){
    legend("bottomright", legend = c("CODAK", "diffcyt", "tie"), 
           fill = alpha(c("red", "darkgreen", "black"), 0.5), 
           box.lty = 0, cex = 0.7)
  }
  
  if (gradient == T){
    gradcol <- alpha(c(rgb(0,rev(seq(0, max(k*powerdiff[which(better == 2)]^r),0.01)), 0),
                 rgb(seq(0, max(k*powerdiff[which(better == 2)]^r),0.01), 0, 0)), 0.5)
    gradientLegend(valRange=c(-1, 1), side=1, inside=TRUE, cex = 0.5, 
                   border.col = NA, col = gradcol,
                   pos=c(.75,-0.005,0.95,.02), coords=TRUE, pos.num = c(2,3,4))
    legend(0.57, 0.07, "Power difference: CODAK - diffcyt", 
           box.col = "white", bg = "white", cex = 0.5)
  }
  
  dev.off()  
}


effsize_powerplot()

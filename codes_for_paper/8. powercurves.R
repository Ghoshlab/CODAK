

############ Combining power curve simulations ##############



rm(list = ls())

library("Cairo")

plotdir <- paste0(getwd(), "/plots/")
resultdir <- paste0(getwd(), "/results/")




powercurves <- function(resultdir, plotdir, restype,
                        methodlist = c("CODAK", "GLMM-Wald", "GLMM-LRT",
                                    "diffcyt-edger", "diffcyt-voom"),
                        method.power = c(expression(pval.kernel),
                                         expression(pval.logistic.mixed.wald),
                                         expression(pval.logistic.mixed.lrt),
                                         expression(pval.diffcyt.edger),
                                         expression(pval.diffcyt.voom)),
                        method.size = c(expression(out$kernel),
                                        expression(out$logistic.mixed[,1]),
                                        expression(out$logistic.mixed[,2]),
                                        expression(out$diffcyt[,1]),
                                        expression(out$diffcyt[,2])),
                        methodcols = c("red", "magenta", "blue", "green", "darkgreen"),
                        lmerperm = F,
                        minalpha = 0.0005,
                        ymax.size = 0.17, ylim.power =c(0,1),
                        drop.after.size = c(2,4),
                        legend2 = T){
  
  pathvec <- c(paste0(resultdir, restype, "/no_diff.RObj"),
               paste0(resultdir, restype, "/small_diff_all.RObj"),
               paste0(resultdir, restype, "/small_diff_some.RObj"),
               paste0(resultdir, restype, "/large_diff_few.RObj"))
  
  pathvec.lmer <- c(paste0(resultdir, restype, "/no_diff_lmerperm.RObj"),
               paste0(resultdir, restype, "/small_diff_all_lmerperm.RObj"),
               paste0(resultdir, restype, "/small_diff_some_lmerperm.RObj"),
               paste0(resultdir, restype, "/large_diff_few_lmerperm.RObj"))
  
  titlevec <- c("Case 1: No difference",
                "Case 2: Small difference in \n all (20) cell types",
                "Case 3: Small difference in \n some (10) cell types",
                "Case 4: Larger difference in \n a few (5) cell types")
  
  
  Cairo(file = paste0(plotdir, "simulations/powercurves_",restype,".pdf"), 
        height = 550, width = 1800, 
        dpi = 95, typ = "pdf", pointsize = 15)
  
  op=par(oma=c(6.5,0,0,0))
  par(mfrow=c(1,4))
  
  
  ######### Size ##########
  
  ii <- 1
  load(pathvec[ii])
  
  logalphavec <- c(seq(log10(minalpha), log10(0.05), 0.05), log10(0.05))
  alphavec <- 10^logalphavec
  
  for (j in 1:length(methodlist)){
    powervec <- numeric(length(alphavec))
    
    for(k in 1:length(alphavec)){
      powervec[k] <- mean(eval(method.size[j]) <= alphavec[k], na.rm = T)
    }
    
    if (j == 1){
      plot(alphavec,powervec, typ = "l", log = "x", pch = 15, ylim = c(0,ymax.size), 
           col = methodcols[j], lwd = 2, xlab = "Alpha level  (log-10 scale)", 
           ylab = "Power", main = titlevec[ii], 
           cex.main = 1.2, cex.lab = 1.5)
    }else{
      lines(alphavec,powervec, pch = 16, 
            col = methodcols[j], lwd = 2)
    }
  }
  
  lines(alphavec, alphavec, lwd = 1, lty = 3)
  lines(alphavec, 2*alphavec, lwd = 1, lty = 3, col = "darkgrey")
  
  if (lmerperm == T){
    load(pathvec.lmer[ii])
    power.lm.perm.05 <- mean(out$logistic.mixed <= 0.05, na.rm = T)
    points(0.05, power.lm.perm.05, col = "black", pch = "*", cex = 2)
  }
  
  
  
  ######## Power ########
  
  if (!is.null(drop.after.size)){
    method.power <- method.power[-drop.after.size]
    methodcols.power <- methodcols[-drop.after.size]
  }
  
  for (ii in 2:4){
    load(pathvec[ii])
    
    for (j in 1:length(method.power)){
      powervec <- numeric(length(alphavec))
      
      for(k in 1:length(alphavec)){
        powervec[k] <- mean(eval(method.power[j]) <= alphavec[k], na.rm = T)
      }
      
      if (j == 1){
        plot(alphavec,powervec, typ = "l", log = "x", pch = 15, ylim = ylim.power, 
             col = methodcols.power[j], lwd = 2, xlab = "Alpha level  (log-10 scale)", 
             ylab = "Power", main = titlevec[ii], 
             cex.main = 1.2, cex.lab = 1.5)
      }else{
        lines(alphavec,powervec, pch = 16, 
              col = methodcols.power[j], lwd = 2)
      }
    }
    
    if (lmerperm == T){
      load(pathvec.lmer[ii])
      power.lm.perm.05 <- mean(pval.logistic.mixed.perm <= 0.05, na.rm = T)
      points(0.05, power.lm.perm.05, col = "black", pch = "*", cex = 2)
    }
  }
  
  op = par(xpd=NA)
  
  ltyvec <- rep(1, length(methodlist))
  pchvec <- rep(NA, length(methodlist))
  
  if (lmerperm == T){
    ltyvec <- c(ltyvec, NA)
    pchvec <- c(pchvec, "*")
    methodlist <- c(methodlist,"LRT Permutation")
    methodcols <- c(methodcols, "black")
  }
  
  
  if (legend2 == T){
    legend("bottomright",
           legend = methodlist[1:3],  
           col= methodcols[1:3], box.lty = 0,
           lty= ltyvec[1:3], bg="white", horiz=F, 
           pch= pchvec[1:3], inset=c(2.8,-1.05),
           cex = 1.3)
    
    legend("bottomright",
           legend = methodlist[-(1:3)],  
           col= methodcols[-(1:3)], box.lty = 0,
           lty= ltyvec[-(1:3)], bg="white", horiz=F,  
           pch= pchvec[-(1:3)], inset=c(1.8,-1.05),
           cex = 1.3)
  }else{
    legend("bottomright",
           legend = methodlist,  
           col= methodcols, box.lty = 0,
           lty= ltyvec, bg="white", horiz=F,  
           pch= pchvec, inset=c(2.3,-1.05),
           cex = 1.3)
  }
  
  
  dev.off()
}






################ Create the power curves ################


methods_KvL_dc <- list(methods = c("CODAK", "GLMM-Wald", "GLMM-LRT",
                                   "diffcyt-edger", "diffcyt-voom"),
                       method.power = c(expression(pval.kernel),
                                        expression(pval.logistic.mixed.wald),
                                        expression(pval.logistic.mixed.lrt),
                                        expression(pval.diffcyt.edger),
                                        expression(pval.diffcyt.voom)),
                       method.size = c(expression(out$kernel),
                                       expression(out$logistic.mixed[,1]),
                                       expression(out$logistic.mixed[,2]),
                                       expression(out$diffcyt[,1]),
                                       expression(out$diffcyt[,2])),
                       methodcols = c("red", "magenta", "blue", "green", "darkgreen"))

#1
powercurves(resultdir, plotdir, "KvLM_binary",
              methods_KvL_dc[[1]], methods_KvL_dc[[2]], 
              methods_KvL_dc[[3]], methods_KvL_dc[[4]],
              lmerperm = T,
              minalpha = 0.0005,
              ymax.size = 0.17, ylim.power =c(0.3,1))

#2
powercurves(resultdir, plotdir, "KvLM_binary_bigsig",
            methods_KvL_dc[[1]], methods_KvL_dc[[2]], 
            methods_KvL_dc[[3]], methods_KvL_dc[[4]],
            lmerperm = F,
            minalpha = 0.0005,
            ymax.size = 0.17, ylim.power =c(0, 0.6))

#3
powercurves(resultdir, plotdir, "KvLM_continuous",
            methods_KvL_dc[[1]][-c(4,5)], methods_KvL_dc[[2]][-c(4,5)], 
            methods_KvL_dc[[3]][-c(4,5)], methods_KvL_dc[[4]][-c(4,5)],
            lmerperm = F,
            minalpha = 0.0005, drop.after.size = c(2),
            ymax.size = 0.17, ylim.power =c(0, 0.9), legend2 = F)



methods_KvL_dc <- list(methods = c("CODAK-sk", "CODAK-alr", "GLMM-Wald", "GLMM-LRT",
                                   "diffcyt-edger", "diffcyt-voom"),
                       method.power = c(expression(pval.kernel.sk),
                                        expression(pval.kernel.alr),
                                        expression(pval.logistic.mixed.wald),
                                        expression(pval.logistic.mixed.lrt),
                                        expression(pval.diffcyt.edger),
                                        expression(pval.diffcyt.voom)),
                       method.size = c(expression(out$kernel[,1]),
                                       expression(out$kernel[,2]),
                                       expression(out$logistic.mixed[,1]),
                                       expression(out$logistic.mixed[,2]),
                                       expression(out$diffcyt[,1]),
                                       expression(out$diffcyt[,2])),
                       methodcols = c("red", "orange", "magenta", 
                                      "blue", "green", "darkgreen"))


#4
powercurves(resultdir, plotdir, "KvLMC_xbin_zbin",
            methods_KvL_dc[[1]], methods_KvL_dc[[2]], 
            methods_KvL_dc[[3]], methods_KvL_dc[[4]],
            lmerperm = F,
            minalpha = 0.0005,
            ymax.size = 0.2, ylim.power =c(0.3, 1),
            drop.after.size = c(3,5))

#5
powercurves(resultdir, plotdir, "KvLMC_dependent",
            methods_KvL_dc[[1]], methods_KvL_dc[[2]], 
            methods_KvL_dc[[3]], methods_KvL_dc[[4]],
            lmerperm = F,
            minalpha = 0.0005,
            ymax.size = 0.25, ylim.power =c(0.2, 1),
            drop.after.size = c(2,3,5))


#6
powercurves(resultdir, plotdir, "KvLMC_xcont_zbin",
            methods_KvL_dc[[1]][-c(5,6)], methods_KvL_dc[[2]][-c(5,6)], 
            methods_KvL_dc[[3]][-c(5,6)], methods_KvL_dc[[4]][-c(5,6)],
            lmerperm = F,
            minalpha = 0.0005,
            ymax.size = 0.2, ylim.power =c(0, 0.9),
            drop.after.size = c(3), legend2 = F)



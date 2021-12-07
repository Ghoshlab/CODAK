

############ Combining power curve simulations ##############



rm(list = ls())

library("Cairo")
library("scales")

plotdir <- paste0(getwd(), "/plots/")
resultdir <- paste0(getwd(), "/results/")




powercurves <- function(resultdir, plotdir, restype, resrev,
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
                        legend2 = T, l2num = 4,
                        mylwd = 1.5,
                        methorder = NULL, plotname = NULL, ltyvec = NULL){
  
  if (is.null(methorder)){
    methorder <- 1:length(methodlist)
  }
  
  if (is.null(plotname)){
    plotname <- restype
  }
  
  if (is.null(ltyvec)){
    ltyvec <- rep(1, length(methodlist))
  }
  
  pathvec <- c(paste0(resultdir, restype, "/no_diff.RObj"),
               paste0(resultdir, restype, "/small_diff_all.RObj"),
               paste0(resultdir, restype, "/small_diff_some.RObj"),
               paste0(resultdir, restype, "/large_diff_few.RObj"))
  
  pathvec.lmer <- c(paste0(resultdir, restype, "/no_diff_lmerperm.RObj"),
               paste0(resultdir, restype, "/small_diff_all_lmerperm.RObj"),
               paste0(resultdir, restype, "/small_diff_some_lmerperm.RObj"),
               paste0(resultdir, restype, "/large_diff_few_lmerperm.RObj"))
  
  pathvec.rev <- c(paste0(resultdir, resrev, "/no_diff.RObj"),
               paste0(resultdir, resrev, "/small_diff_all.RObj"),
               paste0(resultdir, resrev, "/small_diff_some.RObj"),
               paste0(resultdir, resrev, "/large_diff_few.RObj"))
  
  titlevec <- c("Case 1: No difference",
                "Case 2: Small difference in \n all (20) cell types",
                "Case 3: Small difference in \n some (10) cell types",
                "Case 4: Larger difference in \n a few (5) cell types")
  
  
  Cairo(file = paste0(plotdir, "simulations/powercurves_",plotname,".pdf"), 
        height = 1250, width = 1150, 
        dpi = 95, typ = "pdf", pointsize = 12)
  
  op=par(oma=c(8.5,0,0,0))
  par(mfrow=c(2,2))
  
  
  ######### Size ##########
  
  ii <- 1
  load(pathvec[ii])
  if (exists("out")){out1 <- out}
  load(pathvec.rev[ii])
  
  logalphavec <- c(seq(log10(minalpha), log10(0.05), 0.05), log10(0.05))
  alphavec <- 10^logalphavec
  
  for (j in 1:length(methodlist)){
    powervec <- numeric(length(alphavec))
    
    for(k in 1:length(alphavec)){
      powervec[k] <- mean(eval(method.size[methorder[j]]) <= alphavec[k], na.rm = T)
    }
    
    if (j == 1){
      plot(alphavec,powervec, typ = "l", log = "x", pch = 15, ylim = c(0,ymax.size), 
           col = methodcols[methorder[j]], lwd = mylwd, xlab = "Alpha level  (log-10 scale)", 
           ylab = "Power", main = titlevec[ii], 
           cex.main = 1.2, cex.lab = 1.5, lty = ltyvec[methorder[j]])
      #powervec1 <- powervec
    }else{
      lines(alphavec,powervec, pch = 16, 
            col = methodcols[methorder[j]], lwd = mylwd, lty = ltyvec[methorder[j]])
    }
    #lines(alphavec, powervec1, pch = 16, 
    #      col = methodcols[1], lwd = mylwd)
  }
  
  lines(alphavec, alphavec, lwd = 1, lty = 2, col = "grey25")
  lines(alphavec, 2*alphavec, lwd = 1, lty = 2, col = "darkgrey")
  
  if (lmerperm == T){
    load(pathvec.lmer[ii])
    power.lm.perm.05 <- mean(out$logistic.mixed <= 0.05, na.rm = T)
    points(0.05, power.lm.perm.05, col = "black", pch = "*", cex = 2)
  }
  
  
  
  ######## Power ########
  
  if (!is.null(drop.after.size)){
    method.power <- method.power[-drop.after.size]
    methodcols.power <- methodcols[-drop.after.size]
    methorder <- rank(methorder[-which(methorder %in% drop.after.size)])
  }else{
    method.power <- method.power
    methodcols.power <- methodcols
  }
  
  for (ii in 2:4){
    load(pathvec[ii])
    load(pathvec.rev[ii])
    
    for (j in 1:length(method.power)){
      powervec <- numeric(length(alphavec))
      
      for(k in 1:length(alphavec)){
        powervec[k] <- mean(eval(method.power[methorder[j]]) <= alphavec[k], na.rm = T)
      }
      
      if (j == 1){
        plot(alphavec,powervec, typ = "l", log = "x", pch = 15, ylim = ylim.power, 
             col = methodcols.power[methorder[j]], lwd = mylwd, xlab = "Alpha level  (log-10 scale)", 
             ylab = "Power", main = titlevec[ii], 
             cex.main = 1.2, cex.lab = 1.5, lty = ltyvec[methorder[j]])
        #powervec1 <- powervec
      }else{
        lines(alphavec,powervec, pch = 16, 
              col = methodcols.power[methorder[j]], lwd = mylwd, lty = ltyvec[methorder[j]])
      }
      #lines(alphavec, powervec1, pch = 16, 
      #      col = methodcols[1], lwd = mylwd)
    }
    
    if (lmerperm == T){
      load(pathvec.lmer[ii])
      power.lm.perm.05 <- mean(pval.logistic.mixed.perm <= 0.05, na.rm = T)
      points(0.05, power.lm.perm.05, col = "black", pch = "*", cex = 2)
    }
  }
  
  op = par(xpd=NA)
  
  pchvec <- rep(NA, length(methodlist))
  
  if (lmerperm == T){
    ltyvec <- c(ltyvec, NA)
    pchvec <- c(pchvec, "*")
    methodlist <- c(methodlist,"LRT Permutation")
    methodcols <- c(methodcols, "black")
  }
  
  
  if (legend2 == T){
    legend("bottomright",
           legend = methodlist[1:l2num],  
           col= methodcols[1:l2num], box.lty = 0,
           lty= ltyvec[1:l2num], bg="white", horiz=F, 
           pch= pchvec[1:l2num], inset=c(1.4,-0.8),
           cex = 1.3)
    
    legend("bottomright",
           legend = methodlist[-(1:l2num)],  
           col= methodcols[-(1:l2num)], box.lty = 0,
           lty= ltyvec[-(1:l2num)], bg="white", horiz=F,  
           pch= pchvec[-(1:l2num)], inset=c(0.4,-0.8),
           cex = 1.3)
  }else{
    legend("bottomright",
           legend = methodlist,  
           col= methodcols, box.lty = 0,
           lty= ltyvec, bg="white", horiz=F,  
           pch= pchvec, inset=c(0.9,-0.8),
           cex = 1.3)
  }
  
  
  dev.off()
}






################ Create the power curves ################


####### No covariates #######

methods_KvL_dc <- list(methods = c("CODAK", "GLMM-Wald", "GLMM-LRT",
                                   "diffcyt-edger", "diffcyt-voom", 
                                   "dcor", "CODAK-ED", "CODAK-BC", "MiRKAT-Opt"),
                       method.power = c(expression(pval.kernel),
                                        expression(pval.logistic.mixed.wald),
                                        expression(pval.logistic.mixed.lrt),
                                        expression(pval.diffcyt.edger),
                                        expression(pval.diffcyt.voom),
                                        expression(pval.dcor),
                                        expression(pval.euclidean),
                                        expression(pval.kernel.bc),
                                        expression(pval.mirkat.opt)),
                       method.size = c(expression(out1$kernel),
                                       expression(out1$logistic.mixed[,1]),
                                       expression(out1$logistic.mixed[,2]),
                                       expression(out1$diffcyt[,1]),
                                       expression(out1$diffcyt[,2]),
                                       expression(out$dcor),
                                       expression(out$euclidean),
                                       expression(out$kernelbc),
                                       expression(out$mirkat[,3])),
                       methodcols = c("red", "magenta", "blue", "green", "darkgreen",
                                      "purple", "gold", "cyan", "darkseagreen"))

#1
powercurves(resultdir, plotdir, restype = "KvLM_binary", resrev = "KvLM_binary - rev",
              methodlist =  methods_KvL_dc[[1]][-c(6, 9)], method.power =  methods_KvL_dc[[2]][-c(6, 9)], 
              method.size = methods_KvL_dc[[3]][-c(6, 9)], methodcols = methods_KvL_dc[[4]][-c(6, 9)],
              lmerperm = T,
              minalpha = 0.0005,
              ymax.size = 0.17, ylim.power =c(0,1),
              methorder = 7:1)

#2
powercurves(resultdir, plotdir, "KvLM_binary_bigsig", "KvLM_binary_bigsig - rev",
            methods_KvL_dc[[1]][-c(6, 9)], methods_KvL_dc[[2]][-c(6, 9)], 
            methods_KvL_dc[[3]][-c(6, 9)], methods_KvL_dc[[4]][-c(6, 9)],
            lmerperm = F,
            minalpha = 0.0005,
            ymax.size = 0.17, ylim.power =c(0, 0.6),
            methorder = 7:1)

#3
powercurves(resultdir, plotdir, "KvLM_continuous", "KvLM_continuous - rev",
            methods_KvL_dc[[1]][-c(4,5,6,9)], methods_KvL_dc[[2]][-c(4,5,6,9)], 
            methods_KvL_dc[[3]][-c(4,5,6,9)], methods_KvL_dc[[4]][-c(4,5,6,9)],
            lmerperm = F,
            minalpha = 0.0005, drop.after.size = c(2),
            ymax.size = 0.17, ylim.power =c(0, 0.95), legend2 = F,
            methorder = 5:1)

#4
powercurves(resultdir, plotdir, restype = "KvLM_continuous", resrev = "KvLM_continuous - rev",
            methodlist =  methods_KvL_dc[[1]][c(6, 7)], method.power =  methods_KvL_dc[[2]][c(6, 7)], 
            method.size = methods_KvL_dc[[3]][c(6, 7)], methodcols = methods_KvL_dc[[4]][c(6, 7)],
            lmerperm = F,
            minalpha = 0.0005,
            ymax.size = 0.12, ylim.power =c(0,0.09),
            methorder = 2:1,
            plotname = "KvLM_dcor_ed", drop.after.size = NULL,
            legend2 = F, ltyvec = c(2, 1))





###### With covariates ######

methods_KvL_dc <- list(methods = c("CODAK-sk", "CODAK-alr (KC)", "CODAK-alr (FL)",
                                   "GLMM-Wald", "GLMM-LRT",
                                   "diffcyt-edger", "diffcyt-voom", "PERMANOVA",
                                   "dcor", "euclidean", "CODAK-BC", "MiRKAT-AD"),
                       method.power = c(expression(pval.kernel.sk),
                                        expression(pval.kernel.alr1),
                                        expression(pval.kernel.alr2),
                                        expression(pval.logistic.mixed.wald),
                                        expression(pval.logistic.mixed.lrt),
                                        expression(pval.diffcyt.edger),
                                        expression(pval.diffcyt.voom),
                                        expression(pval.permanova),
                                        expression(pval.dcor),
                                        expression(pval.euclidean),
                                        expression(pval.kernel.bc),
                                        expression(pval.mirkat.opt)),
                       method.size = c(expression(out$kernel[,1]),
                                       expression(out$kernel[,2]),
                                       expression(out$kernel[,3]),
                                       expression(out1$logistic.mixed[,1]),
                                       expression(out1$logistic.mixed[,2]),
                                       expression(out1$diffcyt[,1]),
                                       expression(out1$diffcyt[,2]),
                                       expression(out$permanova),
                                       expression(out$dcor),
                                       expression(out$euclidean),
                                       expression(out$kernelbc),
                                       expression(out$mirkat[,1])),
                       methodcols = c("red", "orange", "yellowgreen",
                                      "magenta", 
                                      "blue", "green", "darkgreen", "cornflowerblue",
                                      "gold", "cyan", "darkred", "palevioletred1"))


#1
powercurves(resultdir, plotdir, restype = "KvLMC_xbin_zbin", resrev = "KvLMC_xbin_zbin - rev",
            methodlist = methods_KvL_dc[[1]][-c(9,10,11)], 
            method.power = methods_KvL_dc[[2]][-c(9,10,11)], 
            method.size = methods_KvL_dc[[3]][-c(9,10,11)], 
            methodcols = methods_KvL_dc[[4]][-c(9,10,11)],
            lmerperm = F,
            minalpha = 0.0005,
            ymax.size = 0.2, ylim.power =c(0.3, 1),
            drop.after.size = c(4,6),
            methorder = c(9, 8, 7, 6, 5, 4, 3, 2, 1))

#2
powercurves(resultdir, plotdir, restype = "KvLMC_dependent", resrev = "KvLMC_dependent - rev",
            methodlist = methods_KvL_dc[[1]][-c(9,10,11)], 
            method.power = methods_KvL_dc[[2]][-c(9,10,11)], 
            method.size = methods_KvL_dc[[3]][-c(9,10,11)], 
            methodcols = methods_KvL_dc[[4]][-c(9,10,11)],
            lmerperm = F,
            minalpha = 0.0005,
            ymax.size = 0.25, ylim.power =c(0.25, 1),
            drop.after.size = c(2,3,4,6,7,8),
            methorder = 9:1)


#3
powercurves(resultdir, plotdir, restype = "KvLMC_xcont_zbin", resrev = "KvLMC_xcont_zbin - rev",
            methodlist = methods_KvL_dc[[1]][-c(6,7,9,10,11)], 
            method.power = methods_KvL_dc[[2]][-c(6,7,9,10,11)], 
            method.size = methods_KvL_dc[[3]][-c(6,7,9,10,11)], 
            methodcols = methods_KvL_dc[[4]][-c(6,7,9,10,11)],
            lmerperm = F,
            minalpha = 0.0005,
            ymax.size = 0.2, ylim.power =c(0, 0.9),
            drop.after.size = c(4),
            l2num = 3,
            methorder = 7:1)

#4
powercurves(resultdir, plotdir, restype = "KvLMC_dependent_bigsig", resrev = "KvLMC_dependent_bigsig - rev",
            methodlist = methods_KvL_dc[[1]][-c(9,10,11)], 
            method.power = methods_KvL_dc[[2]][-c(9,10,11)], 
            method.size = methods_KvL_dc[[3]][-c(9,10,11)], 
            methodcols = methods_KvL_dc[[4]][-c(9,10,11)],
            lmerperm = F,
            minalpha = 0.0005,
            ymax.size = 0.55, ylim.power =c(0, 0.85),
            drop.after.size = c(2,3,4,6,7,8,9),
            methorder = 9:1)


#5
# powercurves(resultdir, plotdir, restype = "KvLMC_xbin_zbin_confounder", 
#             resrev = "KvLMC_xbin_zbin_confounder - rev",
#             methodlist = methods_KvL_dc[[1]][c(1,2,7,11)], 
#             method.power = methods_KvL_dc[[2]][c(1,2,7,11)], 
#             method.size = methods_KvL_dc[[3]][c(1,2,7,11)], 
#             methodcols = methods_KvL_dc[[4]][c(1,2,7,11)],
#             lmerperm = F,
#             minalpha = 0.0005,
#             ymax.size = 0.2, ylim.power =c(0.2, 1),
#             drop.after.size = NULL, l2num = 2)


#6
# powercurves(resultdir, plotdir, restype = "KvLMC_xcont_zbin_confounder", 
#             resrev = "KvLMC_xcont_zbin_confounder - rev",
#             methodlist = methods_KvL_dc[[1]][c(1,2,7,11)], 
#             method.power = methods_KvL_dc[[2]][c(1,2,7,11)], 
#             method.size = methods_KvL_dc[[3]][c(1,2,7,11)], 
#             methodcols = methods_KvL_dc[[4]][c(1,2,7,11)],
#             lmerperm = F,
#             minalpha = 0.0005,
#             ymax.size = 0.2, ylim.power =c(0, 0.9),
#             drop.after.size = NULL, l2num = 2)




######## Comparison of different kernels #######

methods_KvL_diffK <- list(methods = c("CODAK", "CODAK-G", "CODAK-xL",
                                   "CODAK-xG"),
                       method.power = c(expression(pval.kernel),
                                        expression(pval.kernel.g),
                                        expression(pval.kernel),
                                        expression(pval.kernel.g)),
                       method.size = c(expression(out1$kernel),
                                       expression(out2$kernel),
                                       expression(out1$kernel),
                                       expression(out2$kernel)),
                       methodcols = c("red", "darksalmon", "red", "darksalmon"))

#1
powercurves(resultdir, plotdir, "KvKG_binary", "KvKG_binary",
            methods_KvL_diffK[[1]][1:2], methods_KvL_diffK[[2]][1:2], 
            methods_KvL_diffK[[3]][1:2], methods_KvL_diffK[[4]][1:2],
            lmerperm = F,
            minalpha = 0.0005,
            ymax.size = 0.1, ylim.power =c(0.5,1), 
            drop.after.size = NULL, legend2 = F,
            methorder = 2:1)


#2
powercurves(resultdir, plotdir, "LKvGK_x_continuous", "LKvGK_x_continuous",
            methods_KvL_diffK[[1]][3:4], methods_KvL_diffK[[2]][3:4], 
            methods_KvL_diffK[[3]][3:4], methods_KvL_diffK[[4]][3:4],
            lmerperm = F,
            minalpha = 0.0005,
            ymax.size = 0.1, ylim.power =c(0,1), 
            drop.after.size = NULL, legend2 = F,
            methorder = 2:1)


#### preprocessing for the observed count matrix Y ####
#### drop off the feature has less than 2 observations for each sample group
Y.filter = function(Ycount, zvec, min.number = 2){
  n = dim(Ycount)[1];p = dim(Ycount)[2]
  K = 2
  filter.vec = rep(1, p)
  cut.point = table(zvec) - min.number
  for(j in 1:p){
    for(k in 1:K){
      group_size = 0; zero_count = 0;
      for(i in 1:n){
        if(zvec[i] == k){
          group_size = group_size + 1
          if(Ycount[i,j] == 0) {zero_count = zero_count + 1}
        }
      }
      if(zero_count > cut.point[k]){
        filter.vec[j] = 0
      }
    }
  }
  cat(sum(filter.vec), "out of", p, "features are kept. \n")
  cat("Feature", which(filter.vec == 0), "are dropped. \n")
  return(list(feature.keep = filter.vec, 
              Y.filter = Ycount[, filter.vec == 1]))
}




#### estimate size factor s(vector) ####
sizefactor.estimator = function(Ycount){
  n = dim(Ycount)[1]; p = dim(Ycount)[2]
  # get sample 50% quantiles aftering removing all 0's
  Y_rm0 = Ycount
  Y_rm0[Y_rm0 == 0] = NA
  sample_Q50 = apply(Y_rm0, 1, quantile, probs = 0.5, na.rm = T)
  si.vec = rep(0, n)
  for(i in 1:n){
    for(j in 1:p){
      if(Ycount[i,j] <= sample_Q50[i]){
        si.vec[i] = si.vec[i] + Ycount[i,j]
      }
    }
  }
  # normalize to have sum of log(si) equals 0
  log.si = log(si.vec) - mean(log(si.vec))
  si.est = round(exp(log.si), 4)
  return(si.est)
}




#### calculate the threshold for gamma / delta for FDR control ####
BayFDR <- function(PPI, alpha){
  PPI_sorted = sort(PPI,decreasing = TRUE)
  k = 1
  fdr = 0
  while(fdr < alpha){
    fdr = mean(1 - PPI_sorted[1:k])
    k = k+1
    if(k > length(PPI_sorted)){
      k = length(PPI_sorted);
      break;
    }
  }
  return.value = PPI_sorted[k]
  return.value = ifelse(is.na(return.value), 0, return.value)
  return(return.value)
}




#### Visualizition functions ####
gamma_VS = function(gammaPPI.raw, gamma.true, sig.level = 0.05){
  th.gamma = BayFDR(gammaPPI.raw, sig.level)
  #### plot ####
  layout(rbind(1,2), heights=c(7,1))  # put legend on bottom 1/8th of the chart
  plot(round(gammaPPI.raw,2), type='h', ylim = c(0,1), 
       ylab = "Posterior Probability of Inclusion",xlab = "Feature Index", 
       main = expression(paste("Evaluation of the Feature Selection for ", gamma)))
  abline(h = th.gamma,lty = 2,col = 'darkred')
  # ppi for disc features:
  disc.ppi = gammaPPI.raw[gamma.true==1]
  # plots:
  points(gammaPPI.raw,col = ifelse(gamma.true==1,'red','grey'), 
         pch = ifelse(gamma.true==1,20,16), cex = ifelse(gamma.true==1,1, 0.5))
  # setup for no margins on the legend
  par(mar=c(0, 0, 0, 0))
  #c(bottom, left, top, right)
  plot.new()
  legend(x = 'center',c("True Discriminating Feature","Non-discriminating Feature"),
        pch = c(20, 16),col=c('red','grey'),ncol=2,bty ="n" ,cex = 0.75)
}
####
delta_ROC = function(deltaPPI.raw, delta.true){
  library(pROC)
  roc.delta = roc(response = delta.true, predictor = deltaPPI.raw, 
                  auc = T)
  plot.roc(roc.delta,print.auc=TRUE,main = expression(paste("ROC plot for estimating ", delta)))
}




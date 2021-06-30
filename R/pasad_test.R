pasad_test = function(obj, 
                      test_idx,
                      newdata,
                      r = 3, 
                      calib = 1,
                      thres = NULL,
                      movn = 10, 
                      plot = TRUE){
  
  #### class
  output = list()
  class(output) = 'training: svd_based '
  
  newdata = newdata[test_idx]
  
  if(length(newdata) < obj$ws){
    warning("data length must be longer than ws")
    
  }
  
  #### input
  obj = obj
  L = obj$L
  N = obj$N
  ws = obj$ws
  C = apply(obj$X, 1, mean)
  P = t(obj$singulars) %*% (obj$singulars)
  C_hat = P %*% C
  train_idx = obj$train_idx
  
  ### make extracted signal 
  s = svd(x = obj$X)
  U = s$u
  Sigma = s$d
  V = s$v
  
  X_element = NULL
  for(i in c(1:(r+1))){
    re = outer(U[,i], V[, i])
    X_element[[i]] = Sigma[i] * re
  }
  X_train_extracted = Reduce("+", lapply(X_element, data.frame)) # list dataframe sum
  X_train_extracted_data = c(X_train_extracted[, 1], 
                             X_train_extracted[, dim(X_train_extracted)[2]]) 
  
  
  
  
  ### make anomaly score
  score = NULL
  for(i in 1:c(length(newdata)-obj$L)){
    val = newdata[i:(i+L-1)]
    temp_m = C_hat - P%*%val
    score[i] = norm(temp_m, type ="2")
  }
  
  threshold = max(score[1:c(obj$N-obj$L+1)])
  threshold = threshold*calib
  
  if(!is.null(thres)){
    threshold = thres
  }
  
  
  #### plot
  old_pars = par(mfrow = c(2,1), mar = c(2,2,2,2))
  on.exit(par(old_pars), add = TRUE)
  
  if(plot == TRUE){
    x_idx = 1:length(obj$x)
    
    ### signal 
    plot(x = test_idx, y = newdata,
         col = "black", type = "l",
         main = "Test Signal")
    abline(v = N, col  = "green", lty = 2)
    
    lines(x = train_idx,
          y = X_train_extracted_data, col = "green", lwd = 1)
    legend("topleft", 
           legend = c("signal", "extracted signal", "training range"), 
           col = c("black", "green", "green"), lty = c(1,1,2))
    
    ### anomaly score
    plot(x = test_idx[-c(1:obj$ws)], y = score, 
         type = "l", col = "blue", pch = 20, cex = .5,
         xlab = "index", ylab = "value", 
         xlim = c(test_idx[1], tail(test_idx,1)),
         ylim = c(0, max(score ,threshold)),
         main = "Anomaly score of Test")
    lines(x = test_idx[-c(1:obj$ws)],
          y = movavg(score, n = movn), col = "red")
    abline(v = N, col  = "green", lty = 2)
    abline(h = threshold, col = "dark grey", lty = 2, lwd = 3)
    legend("topleft", 
           legend = c("score", "mvavg_score"), 
           col = c("blue", "red"), lty = 1)
  }
  
  #### returns
  output$total_scores = score 
  output$tr_score = score[1:c(obj$N-obj$L+1)]
  output$te_score = score[-(1:c(obj$N-obj$L+1))]
  output$extraced = X_train_extracted_data
  output$threshold = threshold
  output$outidx = which(score > threshold)
  return(output) 
  
}
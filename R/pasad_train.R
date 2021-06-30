pasad_train = function(x, 
                       train_idx,
                       r = 3,
                       ws,
                       scree_plot = FALSE
){
  
  
  #### class
  output = list()
  class(output) = 'class object of svd_based'
  
  
  
  ### data input
  x_train = x[train_idx] 
  
  N = length(x_train)
  ws = ifelse(is.null(ws), 
              floor(N/2),
              ws)
  L = ws #  floor(N/2) or ws
  
  
  
  ### Creating trajectory matrix: Trajectory matrix
  X = matrix(NA, nrow = ws, ncol = (N-L+1))
  for(i in 1:(N-L+1)){
    X[, i] = x_train[i:(i+L-1)]
  }

  
  
  ### SVD of log covariance matrix: Solve the SVD problem
  z = svd(X)
  U = z$u
  s = z$d
  
  

  ### scree plot: The statistical dimension
  if(scree_plot == TRUE){
    old_pars = par(mfrow = c(r,1), mar = c(2,2,2,2))
    on.exit(par(old_pars), add = TRUE)
    
    plot(s, type = "o", col = scales::alpha("dark blue", 0.7),
         pch = 20, cex = 1,
         xlim = c(0, 10),
         xlab = "cardinal # of eigen value_", ylab = "eigen value",
         main = paste0("scree plot_r:", 1))
    abline(v = r, col = "Red")
    if(r > 1){
      for(i in 2:r){
        plot(s[-c(1:i)], type = "o", col = scales::alpha("dark blue", 0.7),
             pch = 20, cex = 1,
             xlim = c(0, 10),
             xlab = "cardinal # of eigen value_", ylab = "eigen value",
             main = paste0("scree plot_r:", i))
        abline(v = r, col = "Red")
      }
    }
    
  }
  
  singulars = NULL 
  for(i in 1:r){
    singulars = cbind(singulars, U[,i])
  }
  
  
  
  ### returns
  output$N = N
  output$L = L
  output$U = U
  output$X = X
  output$x = x
  output$ws = ws
  output$train_idx = train_idx
  output$x_train = x_train
  output$singulars = t(singulars)
  
  return(output)
}
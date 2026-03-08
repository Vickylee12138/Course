#Likelihood function
L <- function(theta,y,X) {
  p <- 1/(1 + exp(-X %*% theta))  # logistic function
  likelihood <- prod(dbinom(y, size = 1, prob = p, log = FALSE)) # likelihood 
  return(likelihood)
}

#Log-likelihood function
l <- function(theta,y,X){
  p <- 1/(1+exp(-X %*% theta))
  log_likelihood <- sum((dbinom(y,size=1,prob=p,log= TRUE)))
  return(log_likelihood)
}

#Score function
S <- function(theta,y,X){
  p <- 1/(1+exp(-X %*% theta))
  score <- t(X) %*% (y-p)
  return(score)
}

#Fisher Information I function
I <- function(theta,y,X){
  p <- 1/(1+exp(-X %*% theta))
  v <- p * (1-p)
  D = diag(as.vector(v))
  fisher_info = t(X) %*% D %*% X
  return(fisher_info)
}

#Newton-Raphson
NR <- function(theta0,niter,y,X){
  theta <- theta0
  for (i in 1:niter){
    score <- S(theta0, y, X)
    fisher_info <- I(theta0,y,X)
    theta <- theta + solve(fisher_info) %*% score
  }
  return(theta) # Return the final estimate after niter iterations
}

NR <- function(theta0, niter, y, X, tol = 0.01) {
  theta <- theta0  # 初始化参数
  theta_history <- matrix(NA, niter + 1, length(theta0))  
  theta_history[1, ] <- theta0 
  
  for (i in 1:niter) {
    score <- S(theta, y, X)  
    fisher_info <- I(theta, y, X)  
    theta_new <- theta + solve(fisher_info) %*% score 
    
    
    theta_history[i + 1, ] <- theta_new
    
    
    if (all(abs(theta_new - theta) < tol)) {
      cat("Converged after", i, "iterations\n")
      return(list(final_theta = theta_new, theta_history = theta_history[1:(i + 1), ]))
    }
    
    
    theta <- theta_new
  }
  
  cat("Did not fully converge after", niter, "iterations\n")
  return(list(final_theta = theta, theta_history = theta_history))
}


theta0 <- c(0, 0, 0, 0)  # 初始参数值
niter <- 50  # 最大迭代次数


result <- NR(theta0, niter, y, X)

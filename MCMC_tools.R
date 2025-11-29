library(mvtnorm)

# banana_func: 输入原始高斯的中心，协方差，弯曲强度，返回banana-shape的pdf
banana_func <- function(center, Sigma, b = 0.03){
  d <- length(center) 
  Sigma_inv <- solve(Sigma)
  const <- -0.5 * ( d * log(2*pi) + log(100) )
  
  # 返回的密度函数：输入 y（长度 d），输出香蕉分布密度
  banana_density <- function(y){
    f <- y
    f[2] <- y[2] - b * (y[1]^2 - 10)
    log_density <- -0.5 * t(f - center) %*% Sigma_inv %*% (f - center) + const
    return( exp(log_density) )
  }
  banana_density
}

# final_target: 根据指定，拼接最后的加权多峰
final_target <- function(Centers, Sigmas, weights, ifbanana, b = 0.03){
  d <- length(Centers[[1]]) # Dimension
  n <- length(Centers) # how many density in total
  final_density <- function(x){
    density <- c(rep(0,n)) 
    for (i in 1:n){
      if (ifbanana[i]){
        density[i] <- banana_func(Centers[[i]], Sigmas[[i]])(x)
      }else{
        density[i] <- dmvnorm(x, mean = Centers[[i]],sigma = Sigmas[[i]])
      }
    }
    final_den <- weights %*% density
    return(final_den[[1]])
  }
  return(final_density)
}

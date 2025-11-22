
# 初始输入

# # 维度与模式结构
# N                # 模式个数
# mean_modes[i]    # 每个 mode 的中心 μ_i
# Sigma_modes[i]   # 每个 mode 的初始 Σ_i
# 
# # MCMC 基础设定
# epsilon          # jump move 频率 P(jump)
# steps            # 总迭代次数
# 
# # target
# target(x)        # 目标分布 π(x) 的密度（或 log-density）
# 
# # 初始状态
# state = { x0 , i0 }
# 
# # gamma 初值
# gamma = {
#     w[i],        # 初始权重
#     Sigma[i],    # mode 协方差（可与 Σ_modes 相同）
#     Sigma_til[i] # 等同于Sigma[i],用于早期经验判断
# }
# 
# # 用于记录每一步输出的日志
# log = empty table or list
# mode_count[i] = 0


# 辅助函数

library(mvtnorm)

# local_move: return newproposal point & acceptance rate
local_move <- function(state, mean_modes, target, gamma){
  x          <- state$x          # x
  i          <- state$i          # i
  mean_i     <- mean_modes[[i]]
  Sigma_i    <- gamma$Sigma[[i]]
  w          <- gamma$w
  N          <- length(w)        # number of modes
  
  # --- propose y from N(x, Sigma_i)
  y <- as.numeric(rmvnorm(1, mean = x, sigma = Sigma_i))
  
  # Q_i(x) and Q_i(y)
  Qi_x <- dmvnorm(x,  mean = mean_i, sigma = Sigma_i)
  Qi_y <- dmvnorm(y,  mean = mean_i, sigma = Sigma_i)
  
  # S(x) = sum_j w_j Q_j(x)
  S_x <- 0
  S_y <- 0
  for(j in 1:N){
    S_x <- S_x + w[[j]] * dmvnorm(x, mean = mean_modes[[j]], sigma = gamma$Sigma[[j]])
    S_y <- S_y + w[[j]] * dmvnorm(y, mean = mean_modes[[j]], sigma = gamma$Sigma[[j]])
  }
  
  # numerator: π(y) * Q_i(y) * S(x)
  num <- target(y) * Qi_y * S_x
  
  # denominator: π(x) * Q_i(x) * S(y)
  den <- target(x) * Qi_x * S_y
  
  acc_rate <- min(1, num / den)
  
  return(list(x = y, acc = acc_rate, mode = i))
}



# jump_move: return newproposal point & acceptance rate
jump_move <- function(state, mean_modes, target, gamma){
  
  x          <- state$x
  i          <- state$i     # current mode
  w          <- gamma$w
  Sigma_list <- gamma$Sigma
  N          <- length(w)
  
  # --- choose new mode k
  k <- sample(setdiff(1:N, i), 1)
  
  # --- propose y from mode k
  y <- as.numeric(rmvnorm(1, mean = mean_modes[[k]], sigma = Sigma_list[[k]]))
  
  # -----------------------------
  # define tilde-pi function
  # -----------------------------
  p_tilde <- function(x, mode){
    Qi_x <- dmvnorm(x, mean = mean_modes[[mode]], sigma = Sigma_list[[mode]])
    
    S_x <- 0
    for(j in 1:N){
      S_x <- S_x + w[j] * dmvnorm(x, mean = mean_modes[[j]], sigma = Sigma_list[[j]])
    }
    return( target(x) * w[mode] * Qi_x / S_x )
  }
  
  # R_i(x) and R_k(y)
  R_i_x <- dmvnorm(x, mean = mean_modes[[i]], sigma = Sigma_list[[i]])
  R_k_y <- dmvnorm(y, mean = mean_modes[[k]], sigma = Sigma_list[[k]])
  
  # a_ik = a_ki = 1/(N-1)
  a_ik <- 1/(N-1)
  a_ki <- 1/(N-1)
  
  # acceptance ratio
  num <- p_tilde(y, k) * a_ki * R_i_x
  den <- p_tilde(x, i) * a_ik * R_k_y
  
  acc_rate <- min(1, num / den)
  
  return(list(x = y, acc = acc_rate, mode = k))
}





# 内层循环

JAMS_step <- function(state, gamma, mean_modes, target, epsilon){
  
  x <- state$x
  i <- state$i
  
  # --- 选择 local 或 jump ---
  u <- runif(1)
  
  if (u > epsilon) {
    ## ---- local move ----
    proposal <- local_move(state, mean_modes, target, gamma)
    
    y   <- proposal$x
    acc <- proposal$acc
    k   <- i                     
    
    # 标记提议类型
    proposal_type <- "local"
    
    # 接受或拒绝
    u2 <- runif(1)
    if (u2 < acc) {
      x_new <- y
      i_new <- i
      accepted <- TRUE
    } else {
      x_new <- x
      i_new <- i
      accepted <- FALSE
    }
    
  } else {
    ## ---- jump move ----
    proposal <- jump_move(state, mean_modes, target, gamma)
    
    y   <- proposal$x
    acc <- proposal$acc
    k   <- proposal$mode        
    
    proposal_type <- "jump"
    
    u3 <- runif(1)
    if (u3 < acc) {
      x_new <- y
      i_new <- k
      accepted <- TRUE
    } else {
      x_new <- x
      i_new <- i
      accepted <- FALSE
    }
  }
  
  # --- 更新 state ---
  state_new <- list(x = x_new, i = i_new)
  
  
  # --- 返回所有信息（供外层记录日志） ---
  return(list(
    state          = state_new,
    proposal_type  = proposal_type,
    y              = y,
    k              = k,
    acc            = acc,
    accepted       = accepted
  ))
}





# update_empirical (更新每个mode的经验位置/协方差) 

update_empirical <- function(state, mode_count, mean_obs, Sigma_obs){
  i_new <- state$i
  x_new <- state$x
  
  mode_count[i_new] <- mode_count[i_new] + 1
  n <- mode_count[i_new]
  
  # ---- 更新经验均值 ----
  old_mean <- mean_obs[[i_new]]
  new_mean <- old_mean + (x_new - old_mean)/n
  mean_obs[[i_new]] <- new_mean
  
  # ---- 更新经验协方差 S_i ----
  # Welford 递推：基于旧均值和新均值
  Sigma_obs[[i_new]] <-
    Sigma_obs[[i_new]] +
    (x_new - old_mean) %*% t(x_new - new_mean)
  
  return(list(
    mode_count = mode_count,
    mean_obs   = mean_obs,
    Sigma_obs  = Sigma_obs
  ))
} 


# update_gamma

update_gamma <- function(out, gamma, mode_count, Sigma_obs,
                         AC1, AC2, alpha, alpha_opt, beta, epsilon_w, d,
                         update_Sigma,
                         update_w)
{
  i_new <- out$state$i
  count <- mode_count[i_new]
  acc   <- out$acc
  
  N     <- length(mode_count)
  n_tot <- sum(mode_count)
  
  ## ------ 1. Sigma_tilde / Sigma 的缩放阶段 ------
  if (count < AC1) {
    if (out$proposal_type == "local" && update_Sigma) {
      
      gamma$Sigma_tilde[[i_new]] <-
        exp(count^(-alpha) * (acc - alpha_opt)) * gamma$Sigma_tilde[[i_new]]
      
      gamma$Sigma[[i_new]] <-
        gamma$Sigma_tilde[[i_new]] + beta * diag(d)
    }
    
    ## ------ 2. 协方差=经验协方差的阶段 ------
  } else {
    if (count %% AC2 == 0) {
      
      ## --- 更新 Sigma ---
      if (update_Sigma) {
        if (mode_count[i_new] > 1) {
          S_i <- Sigma_obs[[i_new]] / (mode_count[i_new] - 1)
          gamma$Sigma[[i_new]] <- S_i + beta * diag(d)
        }
      }
      
      ## --- 更新 w ---
      if (update_w) {
        w_add <- n_tot / (1/epsilon_w - N)
        gamma$w <- (mode_count + w_add) / (n_tot + N * w_add)
      }
    }
  }
  
  return(gamma)
}




# 外层循环

JAMS_MCMC <- function(steps, start_state, start_gamma,
                      mean_modes, target, epsilon, d, AC1, AC2,alpha, alpha_opt, beta, epsilon_w,
                      update_Sigma = TRUE,
                      update_w     = TRUE){
  
  state <- start_state
  gamma <- start_gamma
  N     <- length(mean_modes)
  
  mode_count <- rep(0, N)
  
  mean_obs  <- vector("list", N)     # 每个 mode 一个 d 维向量
  Sigma_obs <- vector("list", N)     # 每个 mode 一个 d×d 协方差矩阵
  
  for(i in 1:N){
    mean_obs[[i]]  <- rep(0, d)
    Sigma_obs[[i]] <- matrix(0, d, d)
  }
  
  log <- vector("list", steps)   
  
  for (s in 1:steps){
    
    # ---- 内层：执行一步 JAMS ----
    out <- JAMS_step(state, gamma,
                     mean_modes, target,
                     epsilon)
    
    # ---- 更新状态 ----
    state <- out$state         # { x_new , i_new }
    
    # ---- 更新 mode_count, mean_obs, Sigma_obs ----
    res <- update_empirical(state, mode_count, mean_obs, Sigma_obs)
    
    mode_count <- res$mode_count
    mean_obs   <- res$mean_obs
    Sigma_obs  <- res$Sigma_obs
    
    # ---- 更新 gamma（基于新状态 + mode_count） ----
    gamma <- update_gamma(out, gamma, mode_count, Sigma_obs,
                          AC1, AC2, alpha, alpha_opt, beta, epsilon_w, d, 
                          update_Sigma = update_Sigma,
                          update_w = update_w)
    
    # ---- 记录日志 ----
    log[[s]] <- list(
      step = s,
      x = state$x,
      i = state$i,
      proposal_type = out$proposal_type,
      proposed_point = out$y,
      proposed_mode  = out$k,
      acc_rate = out$acc,
      accepted = out$accepted,
      gamma = gamma
    )
  }
  
  return(list(
    mode_count = mode_count,
    log = log
  ))
}

























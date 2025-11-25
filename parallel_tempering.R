# ---------------Iteration--------------- #
parallel_tempering<- function(beta, temper_num, target, start_where,Sigma,step,d){
  record          <- array(NA_real_, dim = c(temper_num,(step+1),d))
  proposal_points <- array(NA_real_, dim = c(temper_num,step    ,d))
  prob.all        <- matrix(0      , nrow = temper_num, ncol = (2*step))
  
  for (i in 1:temper_num){
    record[i,1,] <- start_where[[i]]
  }
  
  for (s in 1:step){
    for (temper in 1:temper_num){
      # ---------------Local jump--------------- #
      original_point              <- record[temper, s,]
      # 先把proposal_points的对应位置改成local提议点
      proposal_points[temper,s,]  <- rmvnorm(1, original_point, Sigma[[temper]])
      # prob.all 2s-1列记录第s次更新，初始点        的target()
      # prob.all 2s  列记录第s次更新，local_proposal的target()
      prob.all[temper,2*s-1]     <- target(original_point)
      prob.all[temper,2*s]       <- target(proposal_points[temper,s,])
      
      # 计算acceptance rate
      acc_rate_local <- min(1, exp(beta[temper]*(log(prob.all[temper,2*s])-log(prob.all[temper,2*s-1]))))
      
      u_local <- runif(1)
      if (u_local<acc_rate_local){
        proposal_points[temper,s,] <- proposal_points[temper,s,]
      }else{
        # 如果拒绝更换位置
        # proposal_points的对应位置改回第s次更新初始点
        # prob.all 2s列对应修改，保证始终反应proposal_points的s列的target()
        proposal_points[temper,s,] <- original_point
        prob.all[temper,2*s]       <- prob.all[temper,2*s-1]
      }
    }
    
    # ---------------Swith decision--------------- #
    if (s%%2 == 1){
      # 更换1<->2; 3<->4.
      pairs <- cbind(
        seq(1,temper_num-1,2),
        seq(2,temper_num  ,2)
      )}else{
        # 更换2<->3; 4<->5...
      pairs <- cbind(
        seq(2,temper_num-1,2),
        seq(3,temper_num  ,2))
      }
      
    for (j in 1:nrow(pairs)){
      original_point2 <- proposal_points[pairs[j,1],s,]
      propose_switch2 <- proposal_points[pairs[j,2],s,]
      
      acc_rate_switch <- min(1, exp(log(prob.all[pairs[j,2],2*s])*(beta[pairs[j,1]]-beta[pairs[j,2]])+
                                      log(prob.all[pairs[j,1],2*s])*(beta[pairs[j,2]]-beta[pairs[j,1]])))
      
      u_switch1 <- runif(1)
      if (u_switch1<acc_rate_switch){
        record[pairs[j,1],s+1,] <- proposal_points[pairs[j,2],s,]
        record[pairs[j,2],s+1,] <- proposal_points[pairs[j,1],s,]
      }
    }
    
    # 将没有写入的用proposal-point填补 包含：1.未参与交换 2. 参与了但被拒绝
    idx <- apply(record[, s+1, ], 1, function(v) all(is.na(v)))
    record[idx, s+1, ] <- proposal_points[idx, s, ]
  }
  return(list(record = record))
}
      
      
    
      
      
      
      
      
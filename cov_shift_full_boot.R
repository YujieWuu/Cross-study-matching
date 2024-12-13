rm(list=ls())
# Simulation setup 1:
# X have different distribution because of study selection probs depend on X
# Same treatment effect
library(mvtnorm)
library(ggplot2)
library(MatchIt)
library(cluster)
library(nnls)
library(marginaleffects)
library(glmnet)
library(alabama)
d     <- 4
K     <- 4
size  <- 1200
set.seed(2122)
alpha_fix    <- matrix(nrow = K-1, ncol = d, rnorm(d*(K-1), 0, sqrt(0.1)))
beta_X_fix    <- rnorm(d, 0.1, 0.1)
beta_Tx_fix   <- 1
match_method = "nearest"
replacement = TRUE; alpha = alpha_fix
beta_X = beta_X_fix ; beta_Tx = beta_Tx_fix 
true_cluster <- c(1,1,1,1)

logsumexp     <- function(x){
  c = max(x)
  temp = c + log(sum(exp(x-c)))
  return(temp)
}
fn <- function(x, Sigma, study_num){
  weights <- matrix(x[1:study_num], ncol=1)
  obj     <- t(weights)%*%Sigma%*%weights/2
  return(obj)
}
gr <- function(x, Sigma, study_num){
  weights <- matrix(x[1:study_num], ncol=1)
  gr <- Sigma%*%weights 
  return(gr)
}
hin <- function(x, Sigma,study_num){
  weights <- x[1:study_num]
  return(weights)
}
hin.jac <- function(x, Sigma,study_num){
  h.jac <- diag(rep(1,study_num))
  return(h.jac)
}
heq    <- function(x, Sigma,study_num){
  weights <- x[1:study_num]
  h       <- sum(weights) -1
  return(h)
}
heq.jac <- function(x, Sigma,study_num){
  j <- matrix(rep(1,study_num), nrow=1)
  return(j)
}


simfun <- function(d = 4, K = 4,size = 4000, boot_num=1000,
                   match_method = "nearest", 
                   replacement = TRUE, alpha = alpha_fix ,
                   beta_X = beta_X_fix , beta_Tx = beta_Tx_fix ){
  X        <- rmvnorm(size, mean = rep(1, d), sigma = diag(2, d))
  #alpha    <- matrix(nrow = K-1, ncol = d, rnorm(d*(K-1), 0.1, sqrt(0.1)))
  exp_temp <- matrix(NA, nrow = size, ncol = K-1)
  for(i in 1:(K-1)){
    exp_temp[,i] <- exp(X%*%alpha[i,])
  }
  denom    <- rowSums(exp_temp) + 1
  exp_temp <- cbind(exp_temp, denom)
  prob_all <- matrix(NA, nrow = size, ncol = K-1) # prob for k=2,...,K.
  for(i in 1:(K-1)){
    prob_all[,i] <- exp_temp[,i]/exp_temp[,K]
  }
  prob1    <- 1 - rowSums(prob_all) # prob for k=1
  prob_all <- cbind(prob1, prob_all)
  stud     <- sapply(1:size, function(x)which(rmultinom(1,1,prob_all[x,])==1))
  sample_sizes <- table(stud)
  Tx       <- rbinom(size, 1, 0.5)
  
  Y        <- 1 + X%*%beta_X + Tx*beta_Tx + rnorm(size,0,1)
  
  df_all   <- data.frame(X, Tx, stud, Y)
  
  
  # Summary statistics
  mean_all <- c()
  sd_all   <- c()
  control.mean <- c()
  treated.mean <- c()
  for(i in 1:K){
    df_temp  <- df_all[df_all$stud==i,]
    mean_X   <- sapply(1:d, function(x)mean(df_temp[,x]))
    sd_X     <- sapply(1:d, function(x)sd(df_temp[,x]))
    mean_all <- rbind(mean_all, mean_X)
    sd_all   <- rbind(sd_all, sd_X)
    control.mean[i] <- mean(df_temp$Y[df_temp$Tx==0])
    treated.mean[i] <- mean(df_temp$Y[df_temp$Tx==1])
  }
  rownames(mean_all) <- paste0("Study",1:4)
  colnames(mean_all) <- paste0("X",1:4)
  rownames(sd_all) <- paste0("Study",1:4)
  colnames(sd_all) <- paste0("X",1:4)
  
  ######################################################
  # ATE/ATT in each study without cross-study matching
  ######################################################
  ATT_within_reg_avg <- c()
  ATT_within_var_reg_avg <- c()
  
  for(i in 1:K){
    df_temp       <- df_all[df_all$stud==i,]
    
    fit1 <- lm(Y ~ Tx + X1 + X2 + X3 + X4, df_temp)
    # ATT_within_reg[i] <- fit1$coefficients[2]
    # ATT_within_var_reg[i] <- vcov(fit1)[2,2]
    re_avg <- avg_comparisons(fit1, variables = "Tx")
    ATT_within_reg_avg[i] <- re_avg$estimate
    ATT_within_var_reg_avg[i] <- (re_avg$std.error)^2
  }
  
  
  
  # ATT_within
  ######################################################
  # Estimate ATT through Cross-study matching
  ######################################################
  ctl_mat_ps       <- matrix(NA, nrow = K, ncol = K)
  rownames(ctl_mat_ps) <- paste0("Study ", 1:4)
  colnames(ctl_mat_ps) <- paste0("Study ", 1:4)
  
  ctl_mat_ps_reg_avg       <- matrix(NA, nrow = K, ncol = K)
  rownames(ctl_mat_ps_reg_avg) <- paste0("Study ", 1:4)
  colnames(ctl_mat_ps_reg_avg) <- paste0("Study ", 1:4)
  
  ctl_mat_ps_reg_var_avg   <- matrix(NA, nrow = K, ncol = K)
  rownames(ctl_mat_ps_reg_var_avg) <- paste0("Study ", 1:4)
  colnames(ctl_mat_ps_reg_var_avg) <- paste0("Study ", 1:4)
  
  ATT_mat_ps_reg_avg     <- matrix(NA, nrow = K, ncol = K)
  diag(ATT_mat_ps_reg_avg) <- ATT_within_reg_avg
  ATT_mat_ps_reg_var_avg     <- matrix(NA, nrow = K, ncol = K)
  diag(ATT_mat_ps_reg_var_avg) <- ATT_within_var_reg_avg
  
  ATT_mat_ps_reg_avg_source_only           <- matrix(NA, nrow = K, ncol = K)
  diag(ATT_mat_ps_reg_avg_source_only)     <- ATT_within_reg_avg
  ATT_mat_ps_reg_var_avg_source_only       <- matrix(NA, nrow = K, ncol = K)
  diag(ATT_mat_ps_reg_var_avg_source_only) <- ATT_within_var_reg_avg
  
  ATT_mat_ps_reg_avg_match_ctrl_only         <- matrix(NA, nrow = K, ncol = K)
  diag(ATT_mat_ps_reg_avg_match_ctrl_only)  <- ATT_within_reg_avg
  ATT_mat_ps_reg_var_avg_match_ctrl_only     <- matrix(NA, nrow = K, ncol = K)
  diag(ATT_mat_ps_reg_var_avg_match_ctrl_only)   <- ATT_within_var_reg_avg
  
  for(index_target in 1:K){
    for(index_source in 1:K){
      if(index_target!=index_source){
        # Match the two controls
        target_ctl      <- df_all[df_all$stud==index_target & df_all$Tx==0,]
        target_ctl$stud <- 1 # stud here is used to indicate which is the target study (reference for matching)
        source_ctl      <- df_all[df_all$stud==index_source & df_all$Tx==0,]
        source_ctl$stud <- 0
        df_for_mat      <- rbind(target_ctl, source_ctl)
        
        match_ctls      <- matchit(stud ~ X1+X2+X3+X4,
                                   data = df_for_mat, method = match_method,
                                   replace = replacement,
                                   distance = "glm")
        source_ctl_matched <- source_ctl[match_ctls$match.matrix,]
        ctrls_data         <- rbind(target_ctl, source_ctl_matched)
        ctl_mat_ps[index_target,index_source] <- mean(ctrls_data$Y[ctrls_data$stud==1]) - mean(ctrls_data$Y[ctrls_data$stud==0])
        
        matched_ctrls <- match.data(match_ctls)
        fit1 <- lm(Y ~ stud + X1 + X2 + X3 + X4,
                   data = matched_ctrls, weights = weights)
        # ATT_mat_ps_reg[index_target, index_source] <- fit1$coefficients[2]
        # ATT_mat_ps_reg_var[index_target, index_source] <- vcov(fit1)[2,2]
        re_avg <- avg_comparisons(fit1, variables = "stud",
                                  vcov = "HC3",
                                  wts = "weights")
        ctl_mat_ps_reg_avg[index_target,index_source] <- re_avg$estimate
        ctl_mat_ps_reg_var_avg[index_target, index_source] <- (re_avg$std.error)^2
        
        
        # Match treated in study 1 to controls in study 2
        target_tx      <- df_all[df_all$stud==index_target & df_all$Tx==1,]
        target_tx$stud <- 1 # stud here is used to indicate which is the target study (reference for matching)
        source_ct      <- df_all[df_all$stud==index_source & df_all$Tx==0,]
        source_ct$stud <- 0
        df_for_mat     <- rbind(target_tx, source_ct)
        match_tx       <- matchit(stud ~ X1+X2+X3+X4,
                                  data = df_for_mat, method = match_method,
                                  replace = replacement,
                                  distance = "glm", estimand = "ATT")
        # source_ct_matched <- source_ct[match_tx$match.matrix,]
        # tx.data           <- rbind(target_tx, source_ct_matched)
        
        # Match controls in study 1 to treated in study 2
        target_ct       <- df_all[df_all$stud==index_target & df_all$Tx==0,]
        target_ct$stud  <- 1
        source_tx       <- df_all[df_all$stud==index_source & df_all$Tx==1,]
        source_tx$stud  <- 0
        df_for_mat2     <- rbind(target_ct, source_tx)
        match_ct        <- matchit(stud ~ X1+X2+X3+X4,
                                   data = df_for_mat2, method = match_method,
                                   replace = replacement,
                                   distance = "glm", estimand = "ATT")
        
        # run regression
        match_data_treat  <- match.data(match_tx)
        match_data_ct     <- match.data(match_ct)
        match_comb        <- rbind(match_data_treat, match_data_ct)
        fit1 <- lm(Y ~ Tx + X1 + X2 + X3 + X4,
                   data = match_comb, weights = weights)
        re_avg <- avg_comparisons(fit1, variables = "Tx",
                                  vcov = "HC3",
                                  wts = "weights")
        ATT_mat_ps_reg_avg[index_target, index_source] <- re_avg$estimate
        ATT_mat_ps_reg_var_avg[index_target, index_source] <- (re_avg$std.error)^2
        
        # Get the matched source data
        source_ctrl    <- match_data_treat[match_data_treat$stud==0,]
        source_tx      <- match_data_ct[match_data_ct$stud==0,]
        source_matched <- rbind(source_ctrl, source_tx)
        fit2 <- lm(Y ~ Tx + X1 + X2 + X3 + X4,
                   data = source_matched, weights = weights)
        re_avg <- avg_comparisons(fit2, variables = "Tx",
                                  vcov = "HC3",
                                  wts = "weights")
        ATT_mat_ps_reg_avg_source_only[index_target, index_source] <- re_avg$estimate
        ATT_mat_ps_reg_var_avg_source_only[index_target, index_source] <- (re_avg$std.error)^2
        
        # Use only the controls
        match_data_treat  <- match.data(match_tx)
        fit3 <- lm(Y ~ Tx + X1 + X2 + X3 + X4,
                   data = match_data_treat, weights = weights)
        re_avg <- avg_comparisons(fit3, variables = "Tx",
                                  vcov = "HC3",
                                  wts = "weights")
        ATT_mat_ps_reg_avg_match_ctrl_only[index_target, index_source] <- re_avg$estimate
        ATT_mat_ps_reg_var_avg_match_ctrl_only[index_target, index_source] <- (re_avg$std.error)^2
        
      }
    }
  }
  
  
  rownames(ATT_mat_ps_reg_avg) <- paste0("Study ", 1:4)
  colnames(ATT_mat_ps_reg_avg) <- paste0("Study ", 1:4)
  rownames(ATT_mat_ps_reg_var_avg) <- paste0("Study ", 1:4)
  colnames(ATT_mat_ps_reg_var_avg) <- paste0("Study ", 1:4)
  
  rownames(ATT_mat_ps_reg_avg_source_only) <- paste0("Study ", 1:4)
  colnames(ATT_mat_ps_reg_avg_source_only) <- paste0("Study ", 1:4)
  rownames(ATT_mat_ps_reg_var_avg_source_only) <- paste0("Study ", 1:4)
  colnames(ATT_mat_ps_reg_var_avg_source_only) <- paste0("Study ", 1:4)
  
  rownames(ATT_mat_ps_reg_avg_match_ctrl_only) <- paste0("Study ", 1:4)
  colnames(ATT_mat_ps_reg_avg_match_ctrl_only) <- paste0("Study ", 1:4)
  rownames(ATT_mat_ps_reg_var_avg_match_ctrl_only) <- paste0("Study ", 1:4)
  colnames(ATT_mat_ps_reg_var_avg_match_ctrl_only) <- paste0("Study ", 1:4)
  
  df <- list(mean_X = mean_all, sd_X = sd_all, sample_sizes = sample_sizes, 
             ctl_mat_ps_reg_avg = ctl_mat_ps_reg_avg,
             ctl_mat_ps_reg_var_avg = ctl_mat_ps_reg_var_avg, data_gen = df_all, 
             ATT_mat_ps_reg_avg = ATT_mat_ps_reg_avg, ATT_mat_ps_reg_var_avg = ATT_mat_ps_reg_var_avg,
             ATT_mat_ps_reg_avg_source_only = ATT_mat_ps_reg_avg_source_only, 
             ATT_mat_ps_reg_var_avg_source_only = ATT_mat_ps_reg_var_avg_source_only,
             ATT_mat_ps_reg_avg_match_ctrl_only = ATT_mat_ps_reg_avg_match_ctrl_only, 
             ATT_mat_ps_reg_var_avg_match_ctrl_only = ATT_mat_ps_reg_var_avg_match_ctrl_only,
             ATT_within_reg_avg = ATT_within_reg_avg, ATT_within_var_reg_avg = ATT_within_var_reg_avg)
  
  
  
  ATT_mat_ps_reg_avg_boot <- matrix(NA, nrow = boot_num, ncol = 16)
  ATT_mat_ps_reg_avg_source_only_boot <- matrix(NA, nrow = boot_num, ncol = 16)
  ATT_mat_ps_reg_avg_match_ctrl_only_boot <- matrix(NA, nrow = boot_num, ncol = 16)
  for(boot in 1:boot_num){
    print(boot)
    df_all_boot <- c()
    for(k in 1:4){
      index_s <- sample(1:nrow(df_all[df_all$stud==k,]), 
                        nrow(df_all[df_all$stud==k,]), replace = TRUE)
      df_all_boot <- rbind(df_all_boot, df_all[df_all$stud==k,][index_s,])
    }
    
    ATT_within_reg_avg <- c()
    for(i in 1:K){
      df_temp       <- df_all_boot[df_all_boot$stud==i,]
      
      fit1 <- lm(Y ~ Tx + X1 + X2 + X3 + X4, df_temp)
      # ATT_within_reg[i] <- fit1$coefficients[2]
      # ATT_within_var_reg[i] <- vcov(fit1)[2,2]
      re_avg <- avg_comparisons(fit1, variables = "Tx")
      ATT_within_reg_avg[i] <- re_avg$estimate
    }
    
    ATT_mat_ps_reg_avg     <- matrix(NA, nrow = K, ncol = K)
    diag(ATT_mat_ps_reg_avg) <- ATT_within_reg_avg
    
    ATT_mat_ps_reg_avg_source_only           <- matrix(NA, nrow = K, ncol = K)
    diag(ATT_mat_ps_reg_avg_source_only)     <- ATT_within_reg_avg
    
    ATT_mat_ps_reg_avg_match_ctrl_only         <- matrix(NA, nrow = K, ncol = K)
    diag(ATT_mat_ps_reg_avg_match_ctrl_only)  <- ATT_within_reg_avg
    
    for(index_target in 1:K){
      for(index_source in 1:K){
        if(index_target!=index_source){
          # Match treated in study 1 to controls in study 2
          target_tx      <- df_all_boot[df_all_boot$stud==index_target & df_all_boot$Tx==1,]
          target_tx$stud <- 1 # stud here is used to indicate which is the target study (reference for matching)
          source_ct      <- df_all_boot[df_all_boot$stud==index_source & df_all_boot$Tx==0,]
          source_ct$stud <- 0
          df_for_mat     <- rbind(target_tx, source_ct)
          match_tx       <- matchit(stud ~ X1+X2+X3+X4,
                                    data = df_for_mat, method = match_method,
                                    replace = replacement,
                                    distance = "glm", estimand = "ATT")
          
          # Match controls in study 1 to treated in study 2
          target_ct       <- df_all_boot[df_all_boot$stud==index_target & df_all_boot$Tx==0,]
          target_ct$stud  <- 1
          source_tx       <- df_all_boot[df_all_boot$stud==index_source & df_all_boot$Tx==1,]
          source_tx$stud  <- 0
          df_for_mat2     <- rbind(target_ct, source_tx)
          match_ct        <- matchit(stud ~ X1+X2+X3+X4,
                                     data = df_for_mat2, method = match_method,
                                     replace = replacement,
                                     distance = "glm", estimand = "ATT")
          
          # run regression
          match_data_treat  <- match.data(match_tx)
          match_data_ct     <- match.data(match_ct)
          match_comb        <- rbind(match_data_treat, match_data_ct)
          fit1 <- lm(Y ~ Tx + X1 + X2 + X3 + X4,
                     data = match_comb, weights = weights)
          re_avg <- avg_comparisons(fit1, variables = "Tx",
                                    vcov = "HC3",
                                    wts = "weights")
          ATT_mat_ps_reg_avg[index_target, index_source] <- re_avg$estimate
          
          # Get the matched source data
          source_ctrl    <- match_data_treat[match_data_treat$stud==0,]
          source_tx      <- match_data_ct[match_data_ct$stud==0,]
          source_matched <- rbind(source_ctrl, source_tx)
          fit2 <- lm(Y ~ Tx + X1 + X2 + X3 + X4,
                     data = source_matched, weights = weights)
          re_avg <- avg_comparisons(fit2, variables = "Tx",
                                    vcov = "HC3",
                                    wts = "weights")
          ATT_mat_ps_reg_avg_source_only[index_target, index_source] <- re_avg$estimate
          
          # Use only the controls
          match_data_treat  <- match.data(match_tx)
          fit3 <- lm(Y ~ Tx + X1 + X2 + X3 + X4,
                     data = match_data_treat, weights = weights)
          re_avg <- avg_comparisons(fit3, variables = "Tx",
                                    vcov = "HC3",
                                    wts = "weights")
          ATT_mat_ps_reg_avg_match_ctrl_only[index_target, index_source] <- re_avg$estimate
        }
      }
    }
    ATT_mat_ps_reg_avg_boot[boot,] <- c(t(ATT_mat_ps_reg_avg))
    ATT_mat_ps_reg_avg_source_only_boot[boot,] <- c(t(ATT_mat_ps_reg_avg_source_only))
    ATT_mat_ps_reg_avg_match_ctrl_only_boot[boot,] <- c(t(ATT_mat_ps_reg_avg_match_ctrl_only))
  }
  
  
  return(list(df_point = df,
              ATT_mat_ps_reg_avg_boot = ATT_mat_ps_reg_avg_boot,
              ATT_mat_ps_reg_avg_source_only_boot = ATT_mat_ps_reg_avg_source_only_boot,
              ATT_mat_ps_reg_avg_match_ctrl_only_boot = ATT_mat_ps_reg_avg_match_ctrl_only_boot))
}


DPM_general   <- function(theta, Sigma_theta_star, Sigma_theta, theta_var,lam1=1,
                          a = 0.5, b = 0.5){ # shape, scale for gamma
  #### DPM on theta
  # get the theta* matrix
  theta_mat <- theta
  theta_star_mat <- theta_mat - matrix(rep(diag(theta_mat), each = 4),
                                       nrow = 4, ncol = 4, byrow=TRUE )
  # plug in values for priors
  alpha <- 0.1
  eta   <- mean((theta_star_mat[row(theta_star_mat)!=col(theta_star_mat)]))
  gamma <- lam1*sd((theta_star_mat[row(theta_star_mat)!=col(theta_star_mat)]))
  
  # get all partitions
  load("partitions4.RData")
  # all number of partitions
  total_part <- length(out)
  poste_prob <- c()
  A_mat_all  <- list()
  mu_ij_all  <- list()
  log_like   <- c()
  for(which_part in 1:total_part){
    # iterate to get PMF for each partition
    part_re    <- out[[which_part]]
    temp       <- prod(gamma(sapply(1:length(part_re), function(x)length(part_re[[x]]))))
    P_Pi       <- ((alpha^(length(part_re))*gamma(alpha))/gamma(alpha+4))*temp
    after_inte <- log(P_Pi)
    # mu_ij is used to record all unique mu's 
    if(length(part_re)==1){
      mu_ij <- matrix(c(1,1), ncol = 1)
    }else{
      # get all mu_ij's
      mu_ij      <- combn(length(part_re),2)
      mu_ij2     <- rbind(mu_ij[2,], mu_ij[1,])
      mu_ij      <- cbind(mu_ij, mu_ij2)
      mu_ij      <- cbind(matrix(rep(1:length(part_re), each =2),
                                 nrow = 2), mu_ij)
    }
    to_remove <- c() # for example, (1)(2)(34), then mu_{A,A} should not exist
    for(ij_pair in 1:ncol(mu_ij)){
      # iterate over all mu_ij 
      select_clust <- mu_ij[,ij_pair]
      target_stud  <- part_re[[select_clust[1]]]
      source_stud  <- part_re[[select_clust[2]]]
      
      if(select_clust[1]==select_clust[2]){
        if(length(target_stud)==1){
          to_remove <- append(to_remove, ij_pair)
        }
      }
    }
    if(length(to_remove)!=0){
      mu_ij <- mu_ij[,-to_remove]
    }
    mu_ij_all[[which_part]] <- mu_ij
    # get the A matrix, organized by row
    A_mat <- matrix(0, nrow = 4*(4-1), ncol = ncol(mu_ij))
    name_A <- c()
    for(i in 1:4){
      for(j in 1:4){
        if(i!=j){
          name_A <- append(name_A, paste0(i,j))
        }
      }
    }
    rownames(A_mat) <- name_A
    colnames(A_mat) <- paste0(mu_ij[1,], mu_ij[2,])
    for(ij_pair in 1:ncol(mu_ij)){
      # iterate over all mu_ij 
      select_clust <- mu_ij[,ij_pair]
      target_stud  <- part_re[[select_clust[1]]]
      source_stud  <- part_re[[select_clust[2]]]
      
      if(select_clust[1]!=select_clust[2]){
        for(i in target_stud){
          for(j in source_stud){
            A_mat[paste0(i,j), ij_pair] <- 1
          }
        }
      }else if(length(target_stud)!=1){
        for(i in target_stud){
          for(j in source_stud){
            if(i != j){
              A_mat[paste0(i,j), ij_pair] <- 1
            }
          }
        }
      }
    }
    A_mat_all[[which_part]] <- A_mat
    # calculate the posterior
    D_Pi                 <- ncol(mu_ij)
    theta_star_vec       <- c(t(theta_star_mat)[col(t(theta_star_mat))!=row(t(theta_star_mat))])
    eta_vec              <- matrix(eta, ncol = 1, nrow = D_Pi)
    Sigma_theta_star     <- Sigma_theta_star
    Sigma_mu             <- diag(gamma^2, D_Pi)
    Sigma_mu_inv         <- solve(Sigma_mu)
    Sigma_theta_star_inv <- solve(Sigma_theta_star)
    log_C_1_simple       <- (-D_Pi/2)*log(2*pi) - (1/2)*log(det(Sigma_mu)) - (1/2)*t(eta_vec)%*%Sigma_mu_inv%*%eta_vec
    Phi                  <- solve(Sigma_mu_inv + t(A_mat)%*%Sigma_theta_star_inv%*%A_mat)
    lambda               <- Phi%*%(t(A_mat)%*%Sigma_theta_star_inv%*%theta_star_vec + Sigma_mu_inv%*%eta_vec)
    log_post_temp        <- after_inte + log_C_1_simple + 0.5*t(lambda)%*%solve(Phi)%*%lambda + (D_Pi/2)*log(2*pi) + 0.5*log(det(Phi))
    log_like[which_part] <- log_C_1_simple + 0.5*t(lambda)%*%solve(Phi)%*%lambda + (D_Pi/2)*log(2*pi) + 0.5*log(det(Phi))
    poste_prob[which_part] <- log_post_temp
    
  }
  
  poste_prob <- exp(poste_prob - logsumexp(poste_prob))
  
  ### begin iterate
  
  alpha_ite <- alpha
  Pi_ite    <- sample(1:length(out), prob = poste_prob, size=1)
  for(iterr in 2:15000){
    Pi_size   <- length(out[[Pi_ite[iterr-1]]])
    
    # sample x
    x <- rbeta(1,alpha_ite[iterr-1]+1, 4)
    # sample alpha
    pi_x <- (a + Pi_size - 1)/(a + Pi_size - 1 + 4*(b-log(x)) )
    which_gamma <- rbinom(1,1,pi_x) 
    if(which_gamma == 1){
      alpha_ite[iterr] <- rgamma(1,shape = a+Pi_size, rate = b-log(x))
    }else{
      alpha_ite[iterr] <- rgamma(1,shape = a+Pi_size-1, rate = b-log(x))
    }
    # sample Pi
    poste_prob <- c()
    for(which_part in 1:total_part){
      # iterate to get PMF for each partition
      part_re    <- out[[which_part]]
      temp       <- prod(gamma(sapply(1:length(part_re), function(x)length(part_re[[x]]))))
      P_Pi       <- ((alpha_ite[iterr]^(length(part_re))*gamma(alpha_ite[iterr]))/gamma(alpha_ite[iterr]+4))*temp
      after_inte <- log(P_Pi)
      
      # calculate the posterior
      poste_prob[which_part] <- log_like[which_part] + after_inte
      
    }
    poste_prob <- exp(poste_prob - logsumexp(poste_prob))
    Pi_ite[iterr] <- sample(1:length(out), prob = poste_prob, size=1)
  }
  
  poste_prob <- rep(0,15)
  names(poste_prob) <- 1:15
  cal_post <- table(Pi_ite[5001:15000])/10000
  poste_prob[names(cal_post)] <- cal_post
  
  # use the inverse of estimated variance
  DPM_partition_est_inv_var    <- matrix(NA, nrow = 15, ncol = 4)
  # use the inverse of bootstrap variance
  DPM_partition_est_inv_var_boot    <- matrix(NA, nrow = 15, ncol = 4)
  # use simple average
  DPM_partition_est_simple_avg <- matrix(NA, nrow = 15, ncol = 4)
  # minimize w^T*Sigma*w, st. w>0 and sum_i w=1
  DPM_partition_est_minimize_var <- matrix(NA, nrow = 15, ncol = 4)
  
  ATE     <- theta
  ATE_var <- theta_var
  ATE_var_boot <- matrix(diag(Sigma_theta), nrow=4,ncol=4,
                         byrow = TRUE)
  for(partition in 1:length(out)){
    temp <- out[[partition]]
    cluste_re <- rep(NULL, 4)
    for(i in 1:length(temp)){
      cluste_re[temp[[i]]] <- i
    }
    for(s in 1:4){
      index.same <- which(cluste_re == cluste_re[s])
      point_est <- ATE[s,index.same]
      
      # inverse variance
      var_est   <- ATE_var[s,index.same]
      weights   <- 1/var_est
      DPM_partition_est_inv_var[partition, s] <- sum((weights*((sum(weights))^(-1)))*point_est)
      # inverse bootstrap variance
      var_est_boot   <- ATE_var_boot[s,index.same]
      weights_boot   <- 1/var_est_boot
      DPM_partition_est_inv_var_boot[partition, s] <- sum((weights_boot*((sum(weights_boot))^(-1)))*point_est)
      # simple average
      DPM_partition_est_simple_avg[partition, s] <- mean(point_est)
      # minimizing variance
      study_num <- length(index.same)
      name_find <- paste0(s, index.same)
      Sigma_boot <- Sigma_theta[name_find, name_find]
      stack  <- auglag(par=rep(1/study_num,study_num), 
                       fn = fn, gr = gr,
                       hin = hin, hin.jac = hin.jac,
                       heq = heq, heq.jac = heq.jac,
                       Sigma = Sigma_boot, study_num = study_num)
      weights_minimize <- stack$par
      DPM_partition_est_minimize_var[partition, s] <- sum((weights_minimize*((sum(weights_minimize))^(-1)))*point_est)
      
    }
  }
  # use the inverse of estimated variance
  DPM_est_inv_var      <- t(DPM_partition_est_inv_var)%*%matrix(poste_prob, ncol = 1)
  # use the inverse of bootstrap variance
  DPM_est_inv_var_boot <- t(DPM_partition_est_inv_var_boot)%*%matrix(poste_prob, ncol = 1)
  # use simple average
  DPM_est_simple_avg   <- t(DPM_partition_est_simple_avg)%*%matrix(poste_prob, ncol = 1)
  # minimize w^T*Sigma*w, st. w>0 and sum_i w=1
  DPM_est_minimize_var <- t(DPM_partition_est_minimize_var)%*%matrix(poste_prob, ncol = 1)
  
  names <- c()
  for(i in 1:length(out)){
    temp <- NULL
    for(j in 1:length(out[[i]])){
      temp <- append(temp, paste0("(",paste(out[[i]][[j]],collapse=""),")"))
    }
    names <- append(names, paste(temp, collapse=""))
  }
  names(poste_prob) <- names
  return(list(DPM_est_inv_var = DPM_est_inv_var, 
              DPM_est_inv_var_boot = DPM_est_inv_var_boot,
              DPM_est_simple_avg = DPM_est_simple_avg,
              DPM_est_minimize_var = DPM_est_minimize_var,
              poste_prob = round(poste_prob, digits=3),
              theta_star_mat = theta_star_mat))
}
DPM_general_fix_alpha   <- function(theta, Sigma_theta_star, Sigma_theta, theta_var, 
                                    alpha=.1, lam1=1){
  #### DPM on theta
  # get the theta* matrix
  theta_mat <- theta
  theta_star_mat <- theta_mat - matrix(rep(diag(theta_mat), each = 4),
                                       nrow = 4, ncol = 4, byrow=TRUE )
  # plug in values for priors
  alpha <- alpha
  eta   <- mean((theta_star_mat[row(theta_star_mat)!=col(theta_star_mat)]))
  gamma <- lam1*sd((theta_star_mat[row(theta_star_mat)!=col(theta_star_mat)]))
  
  # get all partitions
  load("partitions4.RData")
  # all number of partitions
  total_part <- length(out)
  poste_prob <- c()
  for(which_part in 1:total_part){
    # iterate to get PMF for each partition
    part_re    <- out[[which_part]]
    temp       <- prod(gamma(sapply(1:length(part_re), function(x)length(part_re[[x]]))))
    P_Pi       <- ((alpha^(length(part_re))*gamma(alpha))/gamma(alpha+4))*temp
    after_inte <- log(P_Pi)
    # mu_ij is used to record all unique mu's 
    if(length(part_re)==1){
      mu_ij <- matrix(c(1,1), ncol = 1)
    }else{
      # get all mu_ij's
      mu_ij      <- combn(length(part_re),2)
      mu_ij2     <- rbind(mu_ij[2,], mu_ij[1,])
      mu_ij      <- cbind(mu_ij, mu_ij2)
      mu_ij      <- cbind(matrix(rep(1:length(part_re), each =2),
                                 nrow = 2), mu_ij)
    }
    to_remove <- c() # for example, (1)(2)(34), then mu_{A,A} should not exist
    for(ij_pair in 1:ncol(mu_ij)){
      # iterate over all mu_ij 
      select_clust <- mu_ij[,ij_pair]
      target_stud  <- part_re[[select_clust[1]]]
      source_stud  <- part_re[[select_clust[2]]]
      
      if(select_clust[1]==select_clust[2]){
        if(length(target_stud)==1){
          to_remove <- append(to_remove, ij_pair)
        }
      }
    }
    if(length(to_remove)!=0){
      mu_ij <- mu_ij[,-to_remove]
    }
    # get the A matrix, organized by row
    A_mat <- matrix(0, nrow = 4*(4-1), ncol = ncol(mu_ij))
    name_A <- c()
    for(i in 1:4){
      for(j in 1:4){
        if(i!=j){
          name_A <- append(name_A, paste0(i,j))
        }
      }
    }
    rownames(A_mat) <- name_A
    colnames(A_mat) <- paste0(mu_ij[1,], mu_ij[2,])
    for(ij_pair in 1:ncol(mu_ij)){
      # iterate over all mu_ij 
      select_clust <- mu_ij[,ij_pair]
      target_stud  <- part_re[[select_clust[1]]]
      source_stud  <- part_re[[select_clust[2]]]
      
      if(select_clust[1]!=select_clust[2]){
        for(i in target_stud){
          for(j in source_stud){
            A_mat[paste0(i,j), ij_pair] <- 1
          }
        }
      }else if(length(target_stud)!=1){
        for(i in target_stud){
          for(j in source_stud){
            if(i != j){
              A_mat[paste0(i,j), ij_pair] <- 1
            }
          }
        }
      }
    }
    
    # calculate the posterior
    D_Pi                 <- ncol(mu_ij)
    theta_star_vec       <- c(t(theta_star_mat)[col(t(theta_star_mat))!=row(t(theta_star_mat))])
    eta_vec              <- matrix(eta, ncol = 1, nrow = D_Pi)
    Sigma_theta_star     <- Sigma_theta_star
    Sigma_mu             <- diag(gamma^2, D_Pi)
    Sigma_mu_inv         <- solve(Sigma_mu)
    Sigma_theta_star_inv <- solve(Sigma_theta_star)
    log_C_1_simple       <- (-D_Pi/2)*log(2*pi) - (1/2)*log(det(Sigma_mu)) - (1/2)*t(eta_vec)%*%Sigma_mu_inv%*%eta_vec
    Phi                  <- solve(Sigma_mu_inv + t(A_mat)%*%Sigma_theta_star_inv%*%A_mat)
    lambda               <- Phi%*%(t(A_mat)%*%Sigma_theta_star_inv%*%theta_star_vec + Sigma_mu_inv%*%eta_vec)
    log_post_temp        <- after_inte + log_C_1_simple + 0.5*t(lambda)%*%solve(Phi)%*%lambda + (D_Pi/2)*log(2*pi) + 0.5*log(det(Phi))
    
    poste_prob[which_part] <- log_post_temp
    
  }
  
  poste_prob <- exp(poste_prob - logsumexp(poste_prob))
  
  # use the inverse of estimated variance
  DPM_partition_est_inv_var    <- matrix(NA, nrow = 15, ncol = 4)
  # use the inverse of bootstrap variance
  DPM_partition_est_inv_var_boot    <- matrix(NA, nrow = 15, ncol = 4)
  # use simple average
  DPM_partition_est_simple_avg <- matrix(NA, nrow = 15, ncol = 4)
  # minimize w^T*Sigma*w, st. w>0 and sum_i w=1
  DPM_partition_est_minimize_var <- matrix(NA, nrow = 15, ncol = 4)
  
  ATE     <- theta
  ATE_var <- theta_var
  ATE_var_boot <- matrix(diag(Sigma_theta), nrow=4,ncol=4,
                         byrow = TRUE)
  for(partition in 1:length(out)){
    temp <- out[[partition]]
    cluste_re <- rep(NULL, 4)
    for(i in 1:length(temp)){
      cluste_re[temp[[i]]] <- i
    }
    for(s in 1:4){
      index.same <- which(cluste_re == cluste_re[s])
      point_est <- ATE[s,index.same]
      
      # inverse variance
      var_est   <- ATE_var[s,index.same]
      weights   <- 1/var_est
      DPM_partition_est_inv_var[partition, s] <- sum((weights*((sum(weights))^(-1)))*point_est)
      # inverse bootstrap variance
      var_est_boot   <- ATE_var_boot[s,index.same]
      weights_boot   <- 1/var_est_boot
      DPM_partition_est_inv_var_boot[partition, s] <- sum((weights_boot*((sum(weights_boot))^(-1)))*point_est)
      # simple average
      DPM_partition_est_simple_avg[partition, s] <- mean(point_est)
      # minimizing variance
      study_num <- length(index.same)
      name_find <- paste0(s, index.same)
      Sigma_boot <- Sigma_theta[name_find, name_find]
      stack  <- auglag(par=rep(1/study_num,study_num), 
                       fn = fn, gr = gr,
                       hin = hin, hin.jac = hin.jac,
                       heq = heq, heq.jac = heq.jac,
                       Sigma = Sigma_boot, study_num = study_num)
      weights_minimize <- stack$par
      DPM_partition_est_minimize_var[partition, s] <- sum((weights_minimize*((sum(weights_minimize))^(-1)))*point_est)
      
    }
  }
  # use the inverse of estimated variance
  DPM_est_inv_var      <- t(DPM_partition_est_inv_var)%*%matrix(poste_prob, ncol = 1)
  # use the inverse of bootstrap variance
  DPM_est_inv_var_boot <- t(DPM_partition_est_inv_var_boot)%*%matrix(poste_prob, ncol = 1)
  # use simple average
  DPM_est_simple_avg   <- t(DPM_partition_est_simple_avg)%*%matrix(poste_prob, ncol = 1)
  # minimize w^T*Sigma*w, st. w>0 and sum_i w=1
  DPM_est_minimize_var <- t(DPM_partition_est_minimize_var)%*%matrix(poste_prob, ncol = 1)
  
  names <- c()
  for(i in 1:length(out)){
    temp <- NULL
    for(j in 1:length(out[[i]])){
      temp <- append(temp, paste0("(",paste(out[[i]][[j]],collapse=""),")"))
    }
    names <- append(names, paste(temp, collapse=""))
  }
  names(poste_prob) <- names
  return(list(DPM_est_inv_var = DPM_est_inv_var, 
              DPM_est_inv_var_boot = DPM_est_inv_var_boot,
              DPM_est_simple_avg = DPM_est_simple_avg,
              DPM_est_minimize_var = DPM_est_minimize_var,
              poste_prob = round(poste_prob, digits=3),
              theta_star_mat = theta_star_mat))
}

within_study_est_avg      <- matrix(NA, nrow = 300, ncol = 4)
within_study_est_avg_var  <- matrix(NA, nrow = 300, ncol = 4)

weighted_ave_reg_avg_ps               <- matrix(NA, nrow = 300, ncol = 4)
weighted_ave_reg_avg_ps_source_only   <- matrix(NA, nrow = 300, ncol = 4)
weighted_ave_reg_avg_ps_ctrl_only     <- matrix(NA, nrow = 300, ncol = 4)

weighted_ave_reg_avg_ps_minvar               <- matrix(NA, nrow = 300, ncol = 4)
weighted_ave_reg_avg_ps_source_only_minvar   <- matrix(NA, nrow = 300, ncol = 4)
weighted_ave_reg_avg_ps_ctrl_only_minvar     <- matrix(NA, nrow = 300, ncol = 4)

simple_ave_reg_avg_ps              <- matrix(NA, nrow = 300, ncol = 4)
simple_ave_reg_avg_ps_source_only  <- matrix(NA, nrow = 300, ncol = 4)
simple_ave_reg_avg_ps_ctrl_only    <- matrix(NA, nrow = 300, ncol = 4)

cross_s1_estimates_avg_ps  <- matrix(NA, nrow = 300, ncol = 3)
cross_s2_estimates_avg_ps  <- matrix(NA, nrow = 300, ncol = 3)
cross_s3_estimates_avg_ps  <- matrix(NA, nrow = 300, ncol = 3)
cross_s4_estimates_avg_ps  <- matrix(NA, nrow = 300, ncol = 3)

cross_s1_estimates_avg_ps_source_only  <- matrix(NA, nrow = 300, ncol = 3)
cross_s2_estimates_avg_ps_source_only  <- matrix(NA, nrow = 300, ncol = 3)
cross_s3_estimates_avg_ps_source_only  <- matrix(NA, nrow = 300, ncol = 3)
cross_s4_estimates_avg_ps_source_only  <- matrix(NA, nrow = 300, ncol = 3)

cross_s1_estimates_avg_ps_ctrl_only  <- matrix(NA, nrow = 300, ncol = 3)
cross_s2_estimates_avg_ps_ctrl_only  <- matrix(NA, nrow = 300, ncol = 3)
cross_s3_estimates_avg_ps_ctrl_only  <- matrix(NA, nrow = 300, ncol = 3)
cross_s4_estimates_avg_ps_ctrl_only  <- matrix(NA, nrow = 300, ncol = 3)

cross_s1_estimates_avg_ps_var  <- matrix(NA, nrow = 300, ncol = 3)
cross_s2_estimates_avg_ps_var  <- matrix(NA, nrow = 300, ncol = 3)
cross_s3_estimates_avg_ps_var  <- matrix(NA, nrow = 300, ncol = 3)
cross_s4_estimates_avg_ps_var  <- matrix(NA, nrow = 300, ncol = 3)

cross_s1_estimates_avg_ps_var_source_only  <- matrix(NA, nrow = 300, ncol = 3)
cross_s2_estimates_avg_ps_var_source_only  <- matrix(NA, nrow = 300, ncol = 3)
cross_s3_estimates_avg_ps_var_source_only  <- matrix(NA, nrow = 300, ncol = 3)
cross_s4_estimates_avg_ps_var_source_only  <- matrix(NA, nrow = 300, ncol = 3)

cross_s1_estimates_avg_ps_var_ctrl_only  <- matrix(NA, nrow = 300, ncol = 3)
cross_s2_estimates_avg_ps_var_ctrl_only  <- matrix(NA, nrow = 300, ncol = 3)
cross_s3_estimates_avg_ps_var_ctrl_only  <- matrix(NA, nrow = 300, ncol = 3)
cross_s4_estimates_avg_ps_var_ctrl_only  <- matrix(NA, nrow = 300, ncol = 3)

DPM_est_inv_var_current      <- matrix(NA, nrow = 300, ncol = 4)
DPM_est_inv_var_boot_current <- matrix(NA, nrow = 300, ncol = 4)
DPM_est_simple_avg_current   <- matrix(NA, nrow = 300, ncol = 4)
DPM_est_minimize_var_current <- matrix(NA, nrow = 300, ncol = 4)
DPM_theta_post <- matrix(NA, nrow = 300, ncol = 15)

DPM_est_inv_var_source      <- matrix(NA, nrow = 300, ncol = 4)
DPM_est_inv_var_boot_source <- matrix(NA, nrow = 300, ncol = 4)
DPM_est_simple_avg_source   <- matrix(NA, nrow = 300, ncol = 4)
DPM_est_minimize_var_source <- matrix(NA, nrow = 300, ncol = 4)
DPM_theta_post_source_only <- matrix(NA, nrow = 300, ncol = 15)

DPM_est_inv_var_ctrl      <- matrix(NA, nrow = 300, ncol = 4)
DPM_est_inv_var_boot_ctrl <- matrix(NA, nrow = 300, ncol = 4)
DPM_est_simple_avg_ctrl   <- matrix(NA, nrow = 300, ncol = 4)
DPM_est_minimize_var_ctrl <- matrix(NA, nrow = 300, ncol = 4)
DPM_theta_post_ctrl_only  <- matrix(NA, nrow = 300, ncol = 15)


re_all <- list()

cov_all <- list()

params       <- commandArgs(trailingOnly=TRUE)
size_sim     <- as.numeric(params[[1]]) 
sim          <- as.numeric(params[[2]]) 
set.seed(sim+sim^2)

re1 <- simfun(size = size_sim, boot_num = 1000)
paste0("This is the ", sim,"-th iteration")

temp2 <- re1$ATT_mat_ps_reg_avg_boot
colnames(temp2) <- c("11", "12", "13", "14",
                     "21", "22", "23", "24",
                     "31", "32", "33", "34",
                     "41", "42", "43", "44")
Sigma_theta_current <- cov(temp2)
temp2[,c("12","13","14")] <- temp2[,c("12","13","14")]-temp2[,c("11")]
temp2[,c("21","23","24")] <- temp2[,c("21","23","24")]-temp2[,c("22")]
temp2[,c("31","32","34")] <- temp2[,c("31","32","34")]-temp2[,c("33")]
temp2[,c("41","42","43")] <- temp2[,c("41","42","43")]-temp2[,c("44")]
temp2 <- temp2[,!colnames(temp2)%in%c("11", "22",
                                      "33", "44")]
Sigma_theta_star_current <- cov(temp2)

temp2 <- re1$ATT_mat_ps_reg_avg_source_only_boot
colnames(temp2) <- c("11", "12", "13", "14",
                     "21", "22", "23", "24",
                     "31", "32", "33", "34",
                     "41", "42", "43", "44")
Sigma_theta_source <- cov(temp2)
temp2[,c("12","13","14")] <- temp2[,c("12","13","14")]-temp2[,c("11")]
temp2[,c("21","23","24")] <- temp2[,c("21","23","24")]-temp2[,c("22")]
temp2[,c("31","32","34")] <- temp2[,c("31","32","34")]-temp2[,c("33")]
temp2[,c("41","42","43")] <- temp2[,c("41","42","43")]-temp2[,c("44")]
temp2 <- temp2[,!colnames(temp2)%in%c("11", "22",
                                      "33", "44")]
Sigma_theta_star_source <- cov(temp2)

temp2 <- re1$ATT_mat_ps_reg_avg_match_ctrl_only_boot
colnames(temp2) <- c("11", "12", "13", "14",
                     "21", "22", "23", "24",
                     "31", "32", "33", "34",
                     "41", "42", "43", "44")
Sigma_theta_ctrl <- cov(temp2)
temp2[,c("12","13","14")] <- temp2[,c("12","13","14")]-temp2[,c("11")]
temp2[,c("21","23","24")] <- temp2[,c("21","23","24")]-temp2[,c("22")]
temp2[,c("31","32","34")] <- temp2[,c("31","32","34")]-temp2[,c("33")]
temp2[,c("41","42","43")] <- temp2[,c("41","42","43")]-temp2[,c("44")]
temp2 <- temp2[,!colnames(temp2)%in%c("11", "22",
                                      "33", "44")]
Sigma_theta_star_ctrl <- cov(temp2)

cov_all[[sim]] <- list(Sigma_theta_current = Sigma_theta_current,
                       Sigma_theta_source = Sigma_theta_source,
                       Sigma_theta_ctrl = Sigma_theta_ctrl)


# PS: weighted average, after avg_comparisons()
within_study_est_avg[sim,] <- re1$df_point$ATT_within_reg_avg
within_study_est_avg_var[sim,] <- re1$df_point$ATT_within_var_reg_avg
cross_est <- re1$df_point$ATT_mat_ps_reg_avg
cross_s1_estimates_avg_ps[sim,] <- cross_est[1,c(2,3,4)]
cross_s2_estimates_avg_ps[sim,] <- cross_est[2,c(1,3,4)]
cross_s3_estimates_avg_ps[sim,] <- cross_est[3,c(1,2,4)]
cross_s4_estimates_avg_ps[sim,] <- cross_est[4,c(1,2,3)]
cross_var <- re1$df_point$ATT_mat_ps_reg_var_avg
cross_s1_estimates_avg_ps_var[sim,] <- cross_var[1,c(2,3,4)]
cross_s2_estimates_avg_ps_var[sim,] <- cross_var[2,c(1,3,4)]
cross_s3_estimates_avg_ps_var[sim,] <- cross_var[3,c(1,2,4)]
cross_s4_estimates_avg_ps_var[sim,] <- cross_var[4,c(1,2,3)]
for(k in 1:K){
  index.same <- which(true_cluster == true_cluster[k])
  point_est <- cross_est[k,index.same]
  var_est   <- cross_var[k,index.same]
  weights   <- 1/var_est
  weighted_ave_reg_avg_ps[sim,k] <- sum((weights*((sum(weights))^(-1)))*point_est)
  
  # minimizing variance
  study_num <- length(index.same)
  name_find <- paste0(k, index.same)
  Sigma_boot <- Sigma_theta_current[name_find, name_find]
  stack  <- auglag(par=rep(1/study_num,study_num), 
                   fn = fn, gr = gr,
                   hin = hin, hin.jac = hin.jac,
                   heq = heq, heq.jac = heq.jac,
                   Sigma = Sigma_boot, study_num = study_num)
  weights_minimize <- stack$par
  weighted_ave_reg_avg_ps_minvar[sim,k] <- sum((weights_minimize*((sum(weights_minimize))^(-1)))*point_est)
  simple_ave_reg_avg_ps[sim,k] <- mean(point_est)
}

cross_est2 <- re1$df_point$ATT_mat_ps_reg_avg_source_only
cross_s1_estimates_avg_ps_source_only[sim,] <- cross_est2[1,c(2,3,4)]
cross_s2_estimates_avg_ps_source_only[sim,] <- cross_est2[2,c(1,3,4)]
cross_s3_estimates_avg_ps_source_only[sim,] <- cross_est2[3,c(1,2,4)]
cross_s4_estimates_avg_ps_source_only[sim,] <- cross_est2[4,c(1,2,3)]
cross_var2 <- re1$df_point$ATT_mat_ps_reg_var_avg_source_only
cross_s1_estimates_avg_ps_var_source_only[sim,] <- cross_var2[1,c(2,3,4)]
cross_s2_estimates_avg_ps_var_source_only[sim,] <- cross_var2[2,c(1,3,4)]
cross_s3_estimates_avg_ps_var_source_only[sim,] <- cross_var2[3,c(1,2,4)]
cross_s4_estimates_avg_ps_var_source_only[sim,] <- cross_var2[4,c(1,2,3)]
for(k in 1:K){
  index.same <- which(true_cluster == true_cluster[k])
  point_est <- cross_est2[k,index.same]
  var_est   <- cross_var2[k,index.same]
  weights   <- 1/var_est
  weighted_ave_reg_avg_ps_source_only[sim,k] <- sum((weights*((sum(weights))^(-1)))*point_est)
  
  # minimizing variance
  study_num <- length(index.same)
  name_find <- paste0(k, index.same)
  Sigma_boot <- Sigma_theta_source[name_find, name_find]
  stack  <- auglag(par=rep(1/study_num,study_num), 
                   fn = fn, gr = gr,
                   hin = hin, hin.jac = hin.jac,
                   heq = heq, heq.jac = heq.jac,
                   Sigma = Sigma_boot, study_num = study_num)
  weights_minimize <- stack$par
  weighted_ave_reg_avg_ps_source_only_minvar[sim,k] <- sum((weights_minimize*((sum(weights_minimize))^(-1)))*point_est)
  simple_ave_reg_avg_ps_source_only[sim,k] <- mean(point_est)
}

cross_est3 <- re1$df_point$ATT_mat_ps_reg_avg_match_ctrl_only
cross_s1_estimates_avg_ps_ctrl_only[sim,] <- cross_est3[1,c(2,3,4)]
cross_s2_estimates_avg_ps_ctrl_only[sim,] <- cross_est3[2,c(1,3,4)]
cross_s3_estimates_avg_ps_ctrl_only[sim,] <- cross_est3[3,c(1,2,4)]
cross_s4_estimates_avg_ps_ctrl_only[sim,] <- cross_est3[4,c(1,2,3)]
cross_var3 <- re1$df_point$ATT_mat_ps_reg_var_avg_match_ctrl_only
cross_s1_estimates_avg_ps_var_ctrl_only[sim,] <- cross_var3[1,c(2,3,4)]
cross_s2_estimates_avg_ps_var_ctrl_only[sim,] <- cross_var3[2,c(1,3,4)]
cross_s3_estimates_avg_ps_var_ctrl_only[sim,] <- cross_var3[3,c(1,2,4)]
cross_s4_estimates_avg_ps_var_ctrl_only[sim,] <- cross_var3[4,c(1,2,3)]
for(k in 1:K){
  index.same <- which(true_cluster == true_cluster[k])
  point_est <- cross_est3[k,index.same]
  var_est   <- cross_var3[k,index.same]
  weights   <- 1/var_est
  weighted_ave_reg_avg_ps_ctrl_only[sim,k] <- sum((weights*((sum(weights))^(-1)))*point_est)
  
  # minimizing variance
  study_num <- length(index.same)
  name_find <- paste0(k, index.same)
  Sigma_boot <- Sigma_theta_ctrl[name_find, name_find]
  stack  <- auglag(par=rep(1/study_num,study_num), 
                   fn = fn, gr = gr,
                   hin = hin, hin.jac = hin.jac,
                   heq = heq, heq.jac = heq.jac,
                   Sigma = Sigma_boot, study_num = study_num)
  weights_minimize <- stack$par
  weighted_ave_reg_avg_ps_ctrl_only_minvar[sim, k] <- sum((weights_minimize*((sum(weights_minimize))^(-1)))*point_est)
  
  simple_ave_reg_avg_ps_ctrl_only[sim,k] <- mean(point_est)
  
}



# DPM
# DPM
re_dpm <- DPM_general(theta = re1$df_point$ATT_mat_ps_reg_avg, 
                      Sigma_theta_star = Sigma_theta_star_current, 
                      Sigma_theta = Sigma_theta_current,
                      theta_var = re1$df_point$ATT_mat_ps_reg_var_avg)
re_dpm_11 <- DPM_general(theta = re1$df_point$ATT_mat_ps_reg_avg, 
                         Sigma_theta_star = Sigma_theta_star_current, 
                         Sigma_theta = Sigma_theta_current,
                         theta_var = re1$df_point$ATT_mat_ps_reg_var_avg,
                         a=1,b=1)
re_dpm_fixalpha01 <- DPM_general_fix_alpha(theta = re1$df_point$ATT_mat_ps_reg_avg, 
                                           Sigma_theta_star = Sigma_theta_star_current, 
                                           Sigma_theta = Sigma_theta_current,
                                           theta_var = re1$df_point$ATT_mat_ps_reg_var_avg,
                                           alpha=0.1)
re_dpm_fixalpha1 <- DPM_general_fix_alpha(theta = re1$df_point$ATT_mat_ps_reg_avg, 
                                          Sigma_theta_star = Sigma_theta_star_current, 
                                          Sigma_theta = Sigma_theta_current,
                                          theta_var = re1$df_point$ATT_mat_ps_reg_var_avg,
                                          alpha=1)
DPM_est_inv_var_current[sim,]      <- re_dpm$DPM_est_inv_var
DPM_est_inv_var_boot_current[sim,] <- re_dpm$DPM_est_inv_var_boot
DPM_est_simple_avg_current[sim,]   <- re_dpm$DPM_est_simple_avg
DPM_est_minimize_var_current[sim,] <- re_dpm$DPM_est_minimize_var
DPM_theta_post[sim,] <- re_dpm$poste_prob

re_dpm2 <- DPM_general(theta = re1$df_point$ATT_mat_ps_reg_avg_source_only, 
                       Sigma_theta_star = Sigma_theta_star_source, 
                       Sigma_theta = Sigma_theta_source,
                       theta_var = re1$df_point$ATT_mat_ps_reg_var_avg_source_only)
re_dpm2_11 <- DPM_general(theta = re1$df_point$ATT_mat_ps_reg_avg_source_only, 
                       Sigma_theta_star = Sigma_theta_star_source, 
                       Sigma_theta = Sigma_theta_source,
                       theta_var = re1$df_point$ATT_mat_ps_reg_var_avg_source_only,
                       a=1,b=1)
re_dpm2_fixalpha01 <- DPM_general_fix_alpha(theta = re1$df_point$ATT_mat_ps_reg_avg_source_only, 
                                            Sigma_theta_star = Sigma_theta_star_source, 
                                            Sigma_theta = Sigma_theta_source,
                                            theta_var = re1$df_point$ATT_mat_ps_reg_var_avg_source_only,
                                            alpha=0.1)
re_dpm2_fixalpha1 <- DPM_general_fix_alpha(theta = re1$df_point$ATT_mat_ps_reg_avg_source_only, 
                                           Sigma_theta_star = Sigma_theta_star_source, 
                                           Sigma_theta = Sigma_theta_source,
                                           theta_var = re1$df_point$ATT_mat_ps_reg_var_avg_source_only,
                                           alpha=1)
DPM_est_inv_var_source[sim,]      <- re_dpm2$DPM_est_inv_var
DPM_est_inv_var_boot_source[sim,] <- re_dpm2$DPM_est_inv_var_boot
DPM_est_simple_avg_source[sim,]   <- re_dpm2$DPM_est_simple_avg
DPM_est_minimize_var_source[sim,] <- re_dpm2$DPM_est_minimize_var
DPM_theta_post_source_only[sim,]  <- re_dpm2$poste_prob

re_dpm3 <- DPM_general(theta = re1$df_point$ATT_mat_ps_reg_avg_match_ctrl_only, 
                       Sigma_theta_star = Sigma_theta_star_ctrl, 
                       Sigma_theta = Sigma_theta_ctrl,
                       theta_var = re1$df_point$ATT_mat_ps_reg_var_avg_match_ctrl_only)
re_dpm3_11 <- DPM_general(theta = re1$df_point$ATT_mat_ps_reg_avg_match_ctrl_only, 
                          Sigma_theta_star = Sigma_theta_star_ctrl, 
                          Sigma_theta = Sigma_theta_ctrl,
                          theta_var = re1$df_point$ATT_mat_ps_reg_var_avg_match_ctrl_only,
                          a=1,b=1)
re_dpm3_fixalpha01 <- DPM_general_fix_alpha(theta = re1$df_point$ATT_mat_ps_reg_avg_match_ctrl_only, 
                                            Sigma_theta_star = Sigma_theta_star_ctrl, 
                                            Sigma_theta = Sigma_theta_ctrl,
                                            theta_var = re1$df_point$ATT_mat_ps_reg_var_avg_match_ctrl_only,
                                            alpha=0.1)
re_dpm3_fixalpha1 <- DPM_general_fix_alpha(theta = re1$df_point$ATT_mat_ps_reg_avg_match_ctrl_only, 
                                           Sigma_theta_star = Sigma_theta_star_ctrl, 
                                           Sigma_theta = Sigma_theta_ctrl,
                                           theta_var = re1$df_point$ATT_mat_ps_reg_var_avg_match_ctrl_only,
                                           alpha=1)
DPM_est_inv_var_ctrl[sim,]      <- re_dpm3$DPM_est_inv_var
DPM_est_inv_var_boot_ctrl[sim,] <- re_dpm3$DPM_est_inv_var_boot
DPM_est_simple_avg_ctrl[sim,]   <- re_dpm3$DPM_est_simple_avg
DPM_est_minimize_var_ctrl[sim,] <- re_dpm3$DPM_est_minimize_var
DPM_theta_post_ctrl_only[sim,]  <- re_dpm3$poste_prob

save.image(file=paste0("size=",size_sim,"_sim=", sim,"_sim_re.RData"))



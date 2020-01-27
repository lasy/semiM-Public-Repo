source("Scripts/00_functions_viz.R")

lu = function(x)length(unique(x))


my_acf = function(x, lag.min = 1, lag.max = 10){cc = acf(x, lag.max = lag.max, plot = FALSE); return(max(cc$acf[lag.min:lag.max]))}

# x1 = c(rep(c(rep(1, 5), rep(0, 23)), 10), rep(0, 100))
# macf = frollapply(x1, n = 100, FUN = my_acf, lag.min = 23, lag.max = 45,
#                   align = "center")
# macf[is.nan(macf)] = 0
# plot(x1, type = "l")
# points(macf, type = "l", col = "red", lwd = 2)


create_observation_scores = function(d, features = c("bleeding")){
  
  user_id = unique(d$user_id)
  if(length(user_id)>1){stop("several users in the input\n")}
  
  # we expand the feature matrix d so that is has one row per day
  
  obs = data.frame(user_id = user_id,rel_date = min(d$rel_date):max(d$rel_date))
  m = match(obs$rel_date, d$rel_date)
  
  # First days
  
  obs$first_day = d$first_day[m]
  obs$d1_score = obs$first_day * 1
  
  
  # bleeding score
  if("bleeding" %in% features){
    obs$bleeding = d$bleeding[m]
    #obs$bleeding_score = (obs$bleeding/3)^0.4
    obs$bleeding_score = (pmin(2,obs$bleeding)/3)^0.4
    obs$bleeding_score[is.na(obs$bleeding_score)] = 0
    #cens$bleeding = 1
    
    # We also convert a "first day" into bleeding
    obs$bleeding_score[which(obs$first_day)] = max(obs$bleeding_score)
    
    obs$bleeding_density_5d_score = frollmean(obs$bleeding_score, n = 5, align = "center")
    obs$bleeding_density_90d_score = frollmean(obs$bleeding_score, n = 90, align = "left")
    
    obs$acf_score = frollapply(replace_na(obs$bleeding_density_5d_score,0), n = 100, FUN = my_acf, lag.min = 23, lag.max = 45,
                               align = "center")
    obs$acf_score[is.nan(obs$acf_score)] = 0
  }
  
  # LH score
  if("LH" %in% features){
    obs$LH = d$LH[m]
    obs$LH_score = obs$LH
    obs$LH_score[obs$LH == 0] = NA # any days with nor a positive or a negative LH test is missing
    obs$LH_score = (obs$LH_score+1)/2 # = 1 if positive; = 0 if negative
  }
  
  # Mucus score
  if("mucus" %in% features){
    obs$mucus = d$mucus_type[m]
    obs$mucus_score = mucus.dict$score[match(obs$mucus, mucus.dict$names)]
    obs$mucus_score[obs$bleeding >= 1] = NA
  }
  
  # Temperature
  if("temp" %in% features){
    obs$temp = d$temperature[m]
    obs$quest_temp = d$questionable_temp[m]
    j = which(obs$quest_temp); if(length(j)>0) obs$temp[j] = NA  # questionable temperatures are set as missing data points
    
    # detect odd values (e.g. way too low, way too high, or oddly repeated values)
    if(sum(!is.na(obs$temp))>1){ # if there is at least one temperature
      ttemp = table(obs$temp) # histogram of temperatures
      if((length(unique(obs$temp))==1)|(max(ttemp)> 5*(sort(ttemp,decreasing = TRUE)[2]))){weird_temp = as.numeric(names(ttemp)[which.max(ttemp)])}else{weird_temp = -9999}
      odd =  which((obs$temp < 95) | (obs$temp > 101) | (obs$temp == weird_temp))
      if(length(odd)>0){obs$temp[odd] = NA}
    }
    
    # Scale temperature
    median_temp = median(obs$temp, na.rm = TRUE)
    obs$temp_score = obs$temp-median_temp
    # removing extreme values
    obs$temp_score = pmin(2,pmax(-2,obs$temp_score))
    # removing temperatures when bleeding >= heavy
    obs$temp_score[obs$bleeding == 3] = NA
  }
  
  # Pregnancy tests
  if("preg" %in% features){
    #cat("preg\n")
    obs$preg = d$preg_test[m]
    obs$preg_score = obs$preg
    obs$preg_score[obs$preg == 0] = NA # any days with nor a positive or a negative test is missing
    obs$preg_score = (obs$preg_score+1)/2 # = 1 if positive; = 0 if negative
  }
  

  
  # Pregnancy hints
  if(all(c("preg_hint","bleeding","preg","temp") %in% features)){
    ns = 25*7
    # Any long stretches of high temperature?
    temp_score_mod = obs$temp_score; temp_score_mod[obs$bleeding_score>0] = NA
    tmp = (frollapply(temp_score_mod, n = 45, FUN = function(x){sum(x>0.2, na.rm = TRUE)}) >= 30) %>% replace_na(.,0)
    j = which(tmp & c(0, diff(tmp)));
    if(any(diff(j)<20)){j = j[-(which(diff(j)<20)+1)]}
    obs$preg_from_temp = 0*tmp; if(length(j)>0){k = rep(j,each = ns) + (1:ns); k = k[k %in% 1:nrow(obs)]; obs$preg_from_temp[k] = 1}
    # Any pregnancy tests ?
    j = which(obs$preg_score==1)
    obs$preg_from_preg_test = 0*tmp; if(length(j)>0){k = rep(j,each = ns) + (1:ns); k = k[k %in% 1:nrow(obs)]; obs$preg_from_preg_test[k] = 1}
    # combined with low bleeding density over 90 days
    obs$preg_hint_score = pmax(obs$preg_from_temp, obs$preg_from_preg_test ) * (obs$acf_score < 0.2)
    
  }
  
  #obs_mat = as.matrix(obs[,match(paste0(features,"_score"),colnames(obs))])
  #colnames(obs_mat) = features

  return(obs)
}




create_observation_scores_obs_cens = function(d, features = c("bleeding")){   # DEPRECATED
  
  user_id = unique(d$user_id)
  if(length(user_id)>1){stop("several users in the input\n")}
  
  # we expand the feature matrix d so that is has one row per day
  
  obs = data.frame(user_id = user_id,rel_date = min(d$rel_date):max(d$rel_date))
  cens = data.frame(user_id = user_id,rel_date = min(d$rel_date):max(d$rel_date))
  m = match(obs$rel_date, d$rel_date)
  
  # bleeding score
  if("bleeding" %in% features){
    obs$bleeding = d$bleeding[m]
    obs$bleeding_score = (obs$bleeding/3)^0.4
    obs$bleeding_score[is.na(obs$bleeding_score)] = 0
    cens$bleeding = 1
  }
  
  # LH score
  if("LH" %in% features){
    obs$LH = d$LH[m]
    cens$LH = 1*(obs$LH !=0) # = 1 if reported; = 0 if not reported
    cens$LH[is.na(cens$LH)] = 0 
    obs$LH_score = (obs$LH+1)/2 # = 1 if positive; = 0 if negative
    obs$LH_score[is.na(obs$LH_score)] = 0.5 # value does not matter when unreported
  }
  
  # Mucus score
  if("mucus" %in% features){
    obs$mucus = d$mucus_type[m]
    cens$mucus = 1*(!is.na(obs$mucus))
    obs$mucus_score = mucus.dict$score[match(obs$mucus, mucus.dict$names)]
    obs$mucus_score[is.na(obs$mucus_score)] = 0 # value does not matter when unreported
  }
  
  # Temperature
  if("temp" %in% features){
    obs$temp = d$temperature[m]
    obs$quest_temp = d$questionable_temp[m]
    if(sum(obs$quest_temp, na.rm = TRUE)>1) obs$temp[which(obs$quest_temp)] = NA
    
    # censoring
    cens$temp = 1*(!is.na(obs$temp))
    cens$temp[which(obs$quest_temp)] = 0
    
    # detect odd values (e.g. way too low, way too high, or oddly repeated values)
    if(any(!is.na(obs$temp))){
      ttemp = table(obs$temp)
      if(max(ttemp)> 5*(sort(ttemp,decreasing = TRUE)[2])){weird_temp = as.numeric(names(ttemp)[which.max(ttemp)])}else{weird_temp = -9999}
      odd =  which((obs$temp < 95) | (obs$temp > 101) | (obs$temp == weird_temp))
      if(length(odd)>1){
        obs$temp[odd] = NA
        cens$temp[odd] = 0
      }
    }
    
    median_temp = median(obs$temp, na.rm = TRUE)
    obs$temp_score = obs$temp-median_temp
    obs$temp_score[which(cens$temp == 0)] = 0 # value does not matter when unreported
    # removing extreme values
    obs$temp_score = pmin(2,pmax(-2,obs$temp_score))
  }
  
  # Pregnancy tests
  if("preg" %in% features){
    #cat("preg\n")
    obs$preg = d$preg_test[m]
    cens$preg = 1*(obs$preg !=0) # = 1 if reported; = 0 if not reported
    cens$preg[is.na(cens$preg)] = 0 
    
    obs$preg_score = (obs$preg+1)/2 # = 1 if positive; = 0 if negative
    obs$preg_score[is.na(obs$preg_score)] = 0.5 # value does not matter when unreported
  }
  
  # Pregnancy tests
  if("d1" %in% features){
    #cat("d1\n")
    obs$first_day = d$first_day[m]
    cens$d1 = 1*(obs$first_day) # = 1 if reported; = 0 if not reported
    cens$d1[is.na(cens$d1)] = 0 
    
    obs$d1_score = cens$d1 # = 1 if first day ; = 0 otherwise
  }
  
  obs_mat = as.matrix(obs[,match(paste0(features,"_score"),colnames(obs))])
  colnames(obs_mat) = features
  cens_mat = as.matrix(cens[,match(features,colnames(cens))])
  
  obs_and_cens = list(obs = obs, cens = cens, obs_mat = obs_mat, cens_mat = cens_mat)
  return(obs_and_cens)
}




compute_log_likelihood = function(viterbi = viterbi_states_fitted, model = hsmm_fitted$model, obsdata = obsdata, N = 100){
  
  decoded_seq = viterbi$state_num
  M = sum(obsdata$N)
  
  seq_rle = rle(str_c(decoded_seq,"_",viterbi$user_id))
  state_seq = seq_rle$values %>% str_replace(., "_.*","") %>%  as.numeric()
  state_sojourn = seq_rle$lengths
  
  loglik_base = rep(0, M)
  # loglik_init
  loglik_init = loglik_base
  ii = c(0, cumsum(obsdata$N) %>% head(.,-1))+1
  loglik_init[ii] = hsmm$init[decoded_seq[ii]] %>% log()
  # loglik_trans
  it = (c(0, cumsum(state_sojourn))+1)  %>% head(.,-1) %>% tail(.,-1)
  state_trans = data.frame(from = state_seq, to = dplyr::lead(state_seq)) %>% dplyr::filter(!is.na(to)) %>% as.matrix()
  loglik_trans = loglik_base
  loglik_trans[it] = hsmm$trans[state_trans] %>% log()
  loglik_trans[ii] = 0
  
  # log_lik_sojourn
  is = cumsum(state_sojourn) %>% head(.,-1)
  sojourn_i = cbind(state_seq, state_sojourn) %>% head(.,-1)
  log_lik_sojourn = loglik_base
  log_lik_sojourn[is] = hsmm$sojourn[sojourn_i] %>%  log() 
  log_lik_sojourn[ii[-1]-1] = 0
  
  ideal_sojourn_per_state = apply(hsmm$sojourn, 1, which.max)
  sojourn_i_b = cbind(state_seq, ideal_sojourn_per_state[state_seq]) %>% head(.,-1)
  log_lik_sojourn_b = loglik_base
  log_lik_sojourn_b[is] = hsmm$sojourn[sojourn_i_b] %>%  log() 
  log_lik_sojourn_b[ii[-1]-1] = 0
  
  
  
  
  # loglik_obs
  x = obsdata$x
  cens_o = 1*(!is.na(x))
  p_o = sapply(1:hsmm$n_states,function(state) model$dens.emission(obs = x,state = state,parem = model$parms.emission, w = model$weights[state,]))
  loglik_obs = p_o[cbind(1:M,decoded_seq)] %>% log() 
  
  # most likely observations given the decoded sequence
  mus = sapply(1:hsmm$n_obs, function(obs) {
    sapply(1:hsmm$n_states, function(state) most_probable_value(par = model$parms.emission[[obs]], state = state))})
  
  x_b = mus[decoded_seq,]
  p_b = sapply(1:hsmm$n_states,function(state) model$dens.emission(obs = x_b,state = state,parem = model$parms.emission, w = model$weights[state,]))
  loglik_obs_b = p_b[cbind(1:M,decoded_seq)] %>% log() 
  
  
  # distribution of likelihood
  cat("Bootstrap starts.....")
  dims = c(M, hsmm$n_obs,N)
  # creating the "random" observations
  mus_array = array(mus, dim = dims)
  sigmas_diag = sapply(1:hsmm$n_states, function(state) model$parms.emission$sigma[[state]] %>% diag()) %>% t()
  sigmas = sigmas_diag[decoded_seq,] %>% array(.,dim = dims)
  d = seq(-10,10,by = 0.05); dd = dnorm(d); dd = dd/sum(dd) * 10000; dd = round(dd); ddd = rep(d,dd)
  noise = array(sample(ddd), dim = dims)
  x_r = mus_array + noise * sigmas
  # creating the "random" missingness
  f_missing = runif(N); r_missing = sapply(1:N, function(i) rbinom(1000,size = 1, prob = 1-f_missing[i]))
  cens_r = array(1, dim = dims); j = which(!(hsmm$obs_names %in% c("bleeding","d1","bleeding_density_5d","bleeding_density_90d")))
  for(k in 1:N){cens_r[,j,k] = sample(r_missing[,k], M, replace = TRUE)}
  cat("....and....")
  # computing the likelihood of these sequences
  loglik_obs_r = sapply(1:N, function(k) {
    p_obs_r = rep(1,M)
    ok = foreach(state = 1:hsmm$n_states) %do%{ 
      i = which(decoded_seq == state);
      if(length(i)>0){p_obs_r[i] = model$dens.emission(x = x_r[i,,k],c = cens_r[i,,k],state,model = model)}
    }
    return(p_obs_r)
  }) %>%  log()
  cat(".... is done\n")
  
  
  
  loglik = loglik_init + loglik_trans + log_lik_sojourn/2 + loglik_obs
  loglik_b = loglik_init + loglik_trans + log_lik_sojourn_b/2 + loglik_obs_b
  
  LR = exp(loglik_obs - loglik_obs_b)
  
  cs = rowSums(loglik_obs>loglik_obs_r)/N
  
  return(list(loglik = loglik, loglik_b = loglik_b, LR = LR, cs = cs))
  
}





bluid_hsmm = function(hsmm, M){
  
  SOJOURN = t(hsmm$sojourn)
  if(M>nrow(SOJOURN)){
    SOJOURN = rbind(SOJOURN, matrix(0, nrow = M-nrow(SOJOURN), ncol = ncol(SOJOURN)))
  }else{
    SOJOURN = SOJOURN[1:M,]
  }
  SOJOURN = t(t(SOJOURN)/apply(SOJOURN, 2, sum))
  
  this_user_hsmm =  hsmmspec(
    init = hsmm$init, 
    trans = hsmm$trans_no_names, 
    parms.emission = hsmm$emission_par, 
    sojourn = list(d = SOJOURN, type = "nonparametric"), 
    dens.emission = compute_prob,
    rand.emission = generate_random_obs, 
    mstep = compute_new_em_parms,
    weights = hsmm$weights)
  
  return(this_user_hsmm)
}


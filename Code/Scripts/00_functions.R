source("Scripts/00_functions_viz.R")

lu = function(x)length(unique(x))


sigmoid = function(x, x0, slope){
  y = 1/(1+exp(-slope*(x-x0)))
  return(y)
}



event_accuracy = function(decoding, manual_labels, events){
  
  event_accuracy = foreach(event = events, .combine = bind_rows) %do% {
    if(event == "Ovu"){
      j = which(manual_labels$state_name == "Ovu")
    }else if(event %in% c("Birth","Loss")){
      jj = which(manual_labels$state_name == event);
      j = jj[which(c(0, diff(manual_labels$rel_date[jj]))!=1)]
    }else{stop = "unknown event\n"}
    
    ml = manual_labels[j,]
    
    this_event_accuracy = foreach(i = 1:nrow(ml), .combine = bind_rows) %do% {
      D = decoding %>% dplyr::filter(user_id == ml$user_id[i])
      if(event %in% c("Birth","Loss")){
        jj = which(D$state_name == event)
        j = jj[which(c(0, diff(D$rel_date[jj]))!=1)]
        D = D[j,]
      }else{
        D = D %>% dplyr::filter(state_name == event)
      }
      d = D %>% dplyr::filter(rel_date %in% (-15:15+ml$rel_date[i]))
      
      if(nrow(d) == 0){
        decoded_rel_date = NA
        time_diff = - 17
      }else if(nrow(d)>1){
        decoded_rel_date = d$rel_date[which.min(abs(d$rel_date - ml$rel_date[i]))]
        time_diff = decoded_rel_date - ml$rel_date[i] 
      }else{
        decoded_rel_date = d$rel_date
        time_diff =  decoded_rel_date - ml$rel_date[i]
      }
      acc = data.frame(event = event, 
                       state_name = ml$state_name[i],
                       user_id = ml$user_id[i], 
                       labelled_rel_date = ml$rel_date[i], 
                       decoded_rel_date = decoded_rel_date,
                       time_diff = time_diff, 
                       stringsAsFactors = FALSE)
      return(acc)
    }
    return(this_event_accuracy)
  }
  return(event_accuracy)
}






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
    obs$bleeding_score = pmin(2.5,obs$bleeding)/2.5
    obs$bleeding_score[is.na(obs$bleeding_score)] = 0
    #cens$bleeding = 1
    
    # We also convert a "first day" into bleeding
    j = which(obs$first_day & (obs$bleeding_score == 0))
    obs$bleeding_score[j] = pmax(obs$bleeding_score[j],rep(1/2.5, length(j)))
    
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
      # looking for oddly repeated values
      ttemp = table(obs$temp) # histogram of temperatures
      if((length(unique(obs$temp))==1)|(max(ttemp)> 5*(sort(ttemp,decreasing = TRUE)[2]))){weird_temp = as.numeric(names(ttemp)[which.max(ttemp)])}else{weird_temp = -9999}
      # and removing too high or too low
      odd =  which((obs$temp < 95) | (obs$temp > 101) | (obs$temp == weird_temp))
      if(length(odd)>0){obs$temp[odd] = NA}
    }
    
    # Scale temperature
    median_temp = median(obs$temp, na.rm = TRUE)
    obs$temp_score = obs$temp-median_temp
    # removing extreme values
    obs$temp_score = pmin(1.5,pmax(-1.5,obs$temp_score))
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
    # duration of pregnancy hint after a pregnancy hint triggering event
    ns = 50*7
    
    # pregnancy hint triggering events: either a long stretch of high temperatures OR a positive pregnancy test
     # Any long stretches of high temperature?
    temp_score_mod = obs$temp_score; temp_score_mod[obs$bleeding_score>0] = NA # we remove temperature measurements when bleeding is logged
    tmp = (frollapply(temp_score_mod, n = 35, FUN = function(x){x2 = x; x2 = x2[!is.na(x)];sum(x2>0.2, na.rm = TRUE)}, align = "left") >= 20) & # a long stretch of high temperature is at least 20/35 days of relative temperature above 0.2
      (frollapply(temp_score_mod, n = 15, FUN = function(x){median(x, na.rm = TRUE)}, align = "left") >= 0.2) # and a median temperature over 0.2 for 15 days
    tmp = tmp %>% replace_na(.,0)
    j_temp = which(tmp & c(0, diff(tmp))); j_temp = j_temp+14
    if(any(diff(j_temp)<20)){j_temp = j_temp[-(which(diff(j_temp)<20)+1)]}
    obs$preg_from_temp = 0*tmp; if(length(j_temp)>0){k = rep(j_temp,each = ns) + (1:ns); k = k[k %in% 1:nrow(obs)]; obs$preg_from_temp[k] = 1}
     # Any pregnancy tests ?
    j_pregt = which(obs$preg_score==1)
    obs$preg_from_preg_test = 0*tmp; if(length(j_pregt)>0){k = rep(j_pregt,each = ns) + (1:ns); k = k[k %in% 1:nrow(obs)]; obs$preg_from_preg_test[k] = 1}
    if(any(diff(j_pregt)<20)){j_pregt = j_pregt[-(which(diff(j_pregt)<20)+1)]}
     # combining triggering events
    j = sort(c(j_temp, j_pregt))
    if(any(diff(j)<20)){j = j[-(which(diff(j)<20)+1)]}
    
    obs$preg_hint_score = 0
    for(jj in j){
      ix = jj:min((jj+ns-1),nrow(obs))
      rx = 1:length(ix)
      preg_hint = 1+0*rx
      # the pregnancy hint variable is modulated by different factors
      # 1. the bleeding patterns autocorrelation: 
      #    if the autocorrelation is high, the user is unlikely to be pregnant; except at the beginning of the pregnancy
      acf = obs$acf_score[ix] * sigmoid(rx,x0 = 30,slope = 0.3)
      acf[is.na(acf)] = 0
      preg_hint = 1-pmax(0,acf)
      # 2. any bleeding turns this score off
      preg_hint = preg_hint * (obs$bleeding_score[ix] == 0)
      
      # 3. adverse events such as period-like events or negative pregnancy tests decrease the intensity of the pregnancy hint score
      neg_preg_test = (obs$preg_score[ix] == -1)*1
      period_like_bleeding = (obs$bleeding_density_5d_score[ix]>0.5)*1
      period_like_bleeding = ((period_like_bleeding == 1) & c(TRUE, diff(period_like_bleeding) == 1) )* 1
      adverse_events = cumsum(pmax(0,neg_preg_test, period_like_bleeding, na.rm = TRUE))
      preg_hint = preg_hint / (2^adverse_events)
      
      # the pregnancy hint is triggers by any of the events:
      obs$preg_hint_score[ix] = pmax(obs$preg_hint_score[ix], preg_hint) 
    }
  }
  obs$preg_hint_score = round(obs$preg_hint_score/0.2)*0.2
  
  
  # tracking density
  if("obs_dens" %in% features){
    obs$tracking = !is.na(m)
    obs$obs_dens_score = frollmean(obs$tracking, n = 5, align = "left")
  }
  
  
  #obs_mat = as.matrix(obs[,match(paste0(features,"_score"),colnames(obs))])
  #colnames(obs_mat) = features
  
  return(obs)
}



compute_log_likelihood = function(state_decoded_seq = viterbi_states_fitted, model = hsmm_fitted$model, obsdata = obsdata, N = 200){
  
  decoded_seq = state_decoded_seq$state_num
  M = sum(obsdata$N)
  
  seq_rle = rle(str_c(decoded_seq,"_",state_decoded_seq$user_id))
  state_seq = seq_rle$values %>% str_replace(., "_.*","") %>%  as.numeric()
  state_sojourn = seq_rle$lengths
  
  loglik_base = rep(0, M)
  # loglik_init
  loglik_init = loglik_base
  ii = c(0, cumsum(obsdata$N) %>% head(.,-1))+1
  loglik_init[ii] = model$init[decoded_seq[ii]] %>% log()
  # loglik_trans
  it = (c(0, cumsum(state_sojourn))+1)  %>% head(.,-1) %>% tail(.,-1)
  state_trans = data.frame(from = state_seq, to = dplyr::lead(state_seq)) %>% dplyr::filter(!is.na(to)) %>% as.matrix()
  loglik_trans = loglik_base
  loglik_trans[it] = hsmm$trans[state_trans] %>% log()
  loglik_trans[ii] = 0
  
  # log_lik_sojourn
  is = cumsum(state_sojourn) %>% head(.,-1)
  sojourn_i = cbind(state_sojourn,state_seq) %>% head(.,-1)
  log_lik_sojourn = loglik_base
  log_lik_sojourn[is] = model$sojourn$d[sojourn_i] %>%  log() 
  log_lik_sojourn[ii[-1]-1] = 0
  
  ideal_sojourn_per_state = apply(model$sojourn$d, 2, which.max)
  sojourn_i_b = cbind(ideal_sojourn_per_state[state_seq],state_seq) %>% head(.,-1)
  log_lik_sojourn_b = loglik_base
  log_lik_sojourn_b[is] = model$sojourn$d[sojourn_i_b] %>%  log() 
  log_lik_sojourn_b[ii[-1]-1] = 0
  
  
  
  
  # loglik_obs
  x = obsdata$x
  cens_o = 1*(!is.na(x))
  p_o = sapply(1:hsmm$n_states,function(state) model$dens.emission(obs = x,state = state,parem = model$parms.emission, w = model$weights[state,]))
  loglik_obs = p_o[cbind(1:M,decoded_seq)] %>% log() 
  
  # most likely observations given the decoded sequence
  mus = sapply(1:hsmm$n_obs, function(obs) {
    sapply(1:hsmm$n_states, function(state) most_probable_value(par = model$parms.emission[[obs]], state = state))})
  p_mus_state = sapply(1:hsmm$n_states,function(state) model$dens.emission(obs = mus[state,],state = state,parem = model$parms.emission, w = model$weights[state,]))
  loglik_obs_b = p_mus_state[decoded_seq] %>% log()
  
  # distribution of likelihood
  cat("Bootstrap starts.....")
  n_obs = model$parms.emission %>%  length()
  n_states = model$J
  dims = c(M, n_obs,N)
  reduced_dims = c(n_states, n_obs, N)
  # creating the "random" observations
  x_r = array(NA, dim = reduced_dims)
  for(state in 1:n_states){
    x_r[state,,] = generate_random_obs(n = N, parem = model$parms.emission, state = state) %>% t()
  }
  # creating the "random" missingness
  f_missing = runif(N); r_missing = sapply(1:N, function(i) rbinom(1000,size = 1, prob = 1-f_missing[i]))
  cens_r = array(1, dim = reduced_dims);
  for(k in 1:N){
    cens_r[,,k] = sample(r_missing[,k], n_states * n_obs , replace = TRUE);
    j = which(rowSums(cens_r[,,k]) == 0); 
    if(length(j)>0){cens_r[cbind(j,sample(1:n_obs, length(j), replace = TRUE),rep(k, length(j)))] =1 } # making sure that at least one variable is not missing at each time-point
  }
  x_r[cens_r == 0] = NA
  # computing the likelihood of these sequences
  loglik_obs_r1 = sapply(1:N, function(n) {
    p_obs_r = sapply(1:n_states, function(state) model$dens.emission(obs = x_r[state,,n],
                                                                          state = state,
                                                                          parem = model$parms.emission, 
                                                                          w = model$weights[state,]))
    return(p_obs_r)
  }) %>%  log()
  cat("....and....")
  loglik_obs_r = loglik_obs_r1[decoded_seq, ]
  cat(".... is done\n")
  
  loglik = loglik_init + loglik_trans + log_lik_sojourn/2 + loglik_obs
  loglik_b = loglik_init + loglik_trans + log_lik_sojourn_b/2 + loglik_obs_b
  
  LR = exp(loglik_obs - loglik_obs_b)
  
  cs = rowSums(loglik_obs>loglik_obs_r)/N
  
  return(list(loglik = loglik, loglik_b = loglik_b, LR = LR, cs = cs, loglik_obs_r1 = loglik_obs_r1))
  
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



source("Scripts/00_functions_viz.R")
source("Scripts/00_functions_feature_engineering.R")
source("Scripts/00_function_prepare_obs.R")



FAM_data_augmentation = function(X, features = c("dtemp","f.bleeding.last.week"), get_var_names_only = FALSE, get_var_types_only = FALSE){

  if(get_var_names_only) return(features)

  if(get_var_types_only) return(
    list(
      dtemp = list(
        type = "cat",
        values = c("missing","increase","no_change","decrease")
      ),
      f.bleeding.last.week = list(
        type = "continuous"
      )#,
      # bleeding_density_5d = list(
      #   type = "continuous"
      # )
    )
  )

  # the functions for each of these features are in "Scripts/00_functions_feature_engineering.R"
  ordered_features = features

  if("seq_id" %in% colnames(X)){X = X %>% rename(user_id = seq_id, rel_date = t)}

  if(any(diff(order(X$user_id, X$rel_date)) != 1)) stop("X must be arranged by seq_id (user_id) and t (rel_date)")

  E = X
  for(f in ordered_features){
    fun = eval(parse(text = str_c("compute_",f)))
    E = fun(E = E)
  }

  if(any(E$user_id != X$user_id) | any(E$rel_date != X$rel_date)) stop("Order of E is not the same as order of X")

  E = E %>%  select(all_of(features)) %>%  ungroup()
  E
}



prepare_obs_with_missing_level = function(d, observations = c("bleeding","mucus","temp","LH","preg","extra.log")){
  
  user_id = unique(d$user_id)
  if(length(user_id)>1){stop("several users in the input\n")}
  
  # we expand the feature matrix d so that is has one row per day
  X = data.frame(user_id = user_id,rel_date = min(d$rel_date):max(d$rel_date))
  m = match(X$rel_date, d$rel_date)
  
  if("bleeding" %in% observations){
    # we keep the 0, 0.5, 1, 2, 3 code
    X = X %>% mutate(bleeding = d$bleeding[m],
                     first_day = d$first_day[m] %>% replace_na(FALSE))
    # and replace all missing values by 0 (no bleeding)
    X = X %>% mutate(bleeding = bleeding %>% replace_na(0))
    # and we also convert a "first day" into bleeding
    X  = X %>%
      mutate(first_day_bleeding = ifelse(first_day & (bleeding < 1),2,0),
             bleeding = pmax(bleeding, first_day_bleeding)) %>%
      select(-first_day, -first_day_bleeding)
    X = X %>% rename(bleeding_num = bleeding)
    X = X %>% mutate(bleeding = c("none","spotting","light","medium","heavy")[match(bleeding_num,c(0,0.5,1:3))]) %>%
      select(-bleeding_num)
  }
  
  if("mucus" %in% observations)
    X = X %>%
    mutate(mucus = mucus.dict$category[match(d$mucus_type[m], mucus.dict$names)]%>%
             fct_expand("missing") %>%
             replace_na("missing") %>%
             factor(., levels = c("missing",levels(mucus.dict$category))))
  
  
  if("temp" %in% observations){
    # questionable temperatures are transformed into missing data
    d = d %>% mutate(
      quest_temp = ifelse(is.na(questionable_temp), FALSE, questionable_temp),
      temp = ifelse(quest_temp, NA, temperature) )
    
    # we remove temperature values that are oddly repeated
    if(any(!is.na(d$temp))){ # if there is at least one temperature
      # looking for oddly repeated values
      ttemp = sort(table(d$temp), decreasing = TRUE) # histogram of temperatures
      if((length(ttemp) == 1) | (ttemp[1]> (5*ttemp[2]))) {
        weird_temp = names(ttemp)[1] %>% as.numeric()} else{weird_temp = 999}
      d = d %>% mutate(temp = ifelse(temp == weird_temp, NA, temp))
    }
    
    # we scale the temperature
    median_temp = median(d$temp, na.rm = TRUE)
    d = d %>%  mutate(temp =  (temp - median_temp) %>% pmin(.,1.5) %>% pmax(.,-1.5))
    
    #
    X = X %>% mutate(temp = d$temp[m])
  }
  
  
  test_values = c("neg","missing","pos")
  test_levels = c("pos","neg","missing")
  
  if("LH" %in% observations)
    X$LH = test_values[d$LH[m] %>% replace_na(0)+2] %>% factor(.,levels = test_levels)
  
  if("preg" %in% observations)
    X$preg = test_values[d$preg_test[m] %>% replace_na(0)+2] %>% factor(.,levels = test_levels)
  
  if("sex" %in% observations) X = X %>%  mutate(sex = d$sex_short[m])
  
  if("extra.log" %in% observations)  X = X %>% mutate(extra.log = (!is.na(d$sex[m]))*1)
  
  X
}




###------


add_missing_data = function(X, m){
  if(m == 0) return(X)
  if((m < 0) | (m > 1)) stop("m must be between 0 and 1\n")

  # extra.log
  j = which(X$extra.log == 1)
  n_to_remove = round(m*length(j))
  X$extra.log[sample(j,n_to_remove)] = 0

  # preg
  j = which(!is.na(X$preg)) # != "missing")
  n_to_remove = round(m*length(j))
  X$preg[sample(j,n_to_remove)] = NA #"missing"

  # LH
  j = which(!is.na(X$LH)) # j = which(X$LH != "missing")
  n_to_remove = round(m*length(j))
  X$LH[sample(j,n_to_remove)] = NA # "missing"

  # mucus
  j = which(!is.na(X$mucus)) # j = which(X$mucus != "missing")
  n_to_remove = round(m*length(j))
  X$mucus[sample(j,n_to_remove)] = NA #"missing"

  # temp
  j = which(!is.na(X$temp))
  n_to_remove = round(m*length(j))
  X$temp[sample(j,n_to_remove)] = NA

  X
}


###-------


run_fitting_and_decoding = function(Xa = Xa, ground_truth = ground_truth, model_init = FAM_init, experiment_name = NULL, verbose = FALSE){

  X_ad = Xa %>% left_join(., ground_truth %>%  rename(state_GT = state), by = c("seq_id","t"))

  viterbi_init = predict_states_hsmm(model = model_init, X = Xa, method = "viterbi", verbose = FALSE)
  X_ad$state_VI = viterbi_init$state_seq$state

  smoothed_init = predict_states_hsmm(model = model_init, X = Xa, method = "smoothed", verbose = FALSE)
  X_ad$state_SI = smoothed_init$state_seq$state

  model_fitted = fit_hsmm(model = model_init, X = Xa, ground_truth = ground_truth, n_iter = 10,
                          lock.transition = TRUE, lock.sojourn = TRUE, lock.emission = FALSE,
                          verbose = verbose)

  viterbi_fitted = predict_states_hsmm(model = model_fitted$model, X = Xa, method = "viterbi")
  X_ad$state_VF = viterbi_fitted$state_seq$state

  smoothed_fitted = predict_states_hsmm(model = model_fitted$model, X = Xa, method = "smoothed")
  X_ad = X_ad %>% full_join(.,
                            smoothed_fitted$state_seq %>% select(seq_id, t, state) %>% rename(state_SF = state),
                            by=  c("seq_id","t"))

  fitting_par = data.frame(
    seq_id = "all",
    iter = 1:length(model_fitted$fit_param$ll),
    ll = model_fitted$fit_param$ll,
    message = model_fitted$fit_param$message,
    stringsAsFactors = FALSE)


  if(is.null(experiment_name)) experiment_name = str_c("Exp ",Sys.time())
  output = list(name = experiment_name,
                model = model_init,
                X_ad = X_ad, fitting_par = fitting_par ,
                smoothed_prob = list(init = smoothed_init$state_probs, fitted = smoothed_fitted$state_probs))
}



run_fitting_and_decoding_seq_by_seq = function(Xa = Xa, ground_truth = ground_truth, model_init = FAM_init, experiment_name = NULL, verbose = FALSE){

  X_ad = Xa %>% left_join(., ground_truth %>%  rename(state_GT = state), by = c("seq_id","t"))

  viterbi_init = predict_states_hsmm(model = model_init, X = Xa, method = "viterbi")
  X_ad$state_VI = viterbi_init$state_seq$state

  smoothed_init = predict_states_hsmm(model = model_init, X = Xa, method = "smoothed")
  X_ad$state_SI = smoothed_init$state_seq$state


  X_adf = data.frame()
  smoothed_prob_fitted = data.frame()
  fitting_par = data.frame()

  for(sid in unique(Xa$seq)){
    cat(sid, "\n")
    this_Xa = Xa %>% filter(seq_id == sid)
    this_ground_truth = ground_truth %>% filter(seq_id == sid)

    this_X_adf = this_Xa %>% select(seq_id, t)

    #fitting
    model_fitted = fit_hsmm(model = model_init,
                            X = this_Xa, ground_truth = this_ground_truth,
                            n_iter = 10,
                            lock.transition = TRUE, lock.sojourn = TRUE, lock.emission = FALSE,
                            verbose = verbose)


    fitting_par = rbind(fitting_par,
                        data.frame(
                          seq_id = sid,
                          iter = 1:length(model_fitted$fit_param$ll),
                          ll = model_fitted$fit_param$ll,
                          message = model_fitted$fit_param$message,
                          stringsAsFactors = FALSE))

    # decoding with fitted model
    viterbi_fitted = predict_states_hsmm(model = model_fitted$model, X = this_Xa, method = "viterbi")
    this_X_adf$state_VF = viterbi_fitted$state_seq$state

    smoothed_fitted = predict_states_hsmm(model = model_fitted$model, X = this_Xa, method = "smoothed")
    this_X_adf = this_X_adf %>% full_join(.,
                                          smoothed_fitted$state_seq %>% select(seq_id, t, state) %>% rename(state_SF = state),
                                          by=  c("seq_id","t"))

    X_adf = rbind(X_adf, this_X_adf)
    smoothed_prob_fitted = rbind(smoothed_prob_fitted, smoothed_fitted$state_probs)
  }

  tmp = full_join(X_ad, X_adf, by =  c("seq_id","t"))
  X_ad = tmp


  if(is.null(experiment_name)) experiment_name = str_c("Exp ",Sys.time())
  output = list(name = experiment_name,
                model = model_init,
                X_ad = X_ad, fitting_par = fitting_par ,
                smoothed_prob = list(init = smoothed_init$state_probs, fitted = smoothed_fitted$state_probs))

}




subset_model = function(hsmm = hsmm, model = 1){
  i = which(hsmm$states$model <= model)
  if(length(i)==0){stop("this model does not exist or models have not been defined")}
  sub_hsmm = hsmm
  sub_hsmm$states = hsmm$states[i,]
  sub_hsmm$n_states = length(i)
  sub_hsmm$init =  hsmm$init[i]
  sub_hsmm$trans =  hsmm$trans[i,i];  sub_hsmm$trans =  sub_hsmm$trans/ rowSums(sub_hsmm$trans)
  sub_hsmm$trans_no_names =  sub_hsmm$trans; colnames(sub_hsmm$trans_no_names) = NULL;  rownames(sub_hsmm$trans_no_names) = NULL;
  sub_hsmm$sojourn =  hsmm$sojourn[i,]
  for(o in 1:hsmm$n_obs){
    par = hsmm$emission_par[[o]]
    if(par$type == "norm"){par$param$mean = par$param$mean[i];par$param$sd = par$param$sd[i]}
    if(par$type == "binom"){par$param$size = par$param$size[i];par$param$prob = par$param$prob[i]}
    if(par$type == "non-par"){par$param$probs = par$param$probs[,i]}
    sub_hsmm$emission_par[[o]] = par
  }

  sub_hsmm$glm_models =  hsmm$glm_models[i]
  sub_hsmm$weights =  hsmm$weights[i,]
  return(sub_hsmm)
}




refresh_mhsmm_LSY = function(){
  mhsmm_LSY_dir = "Scripts/mhsmm_LSY/R/"
  mhsmm_LSY_functions = list.files(mhsmm_LSY_dir)
  for(f in mhsmm_LSY_functions){
    cat(f,"\n")
    source(str_c(mhsmm_LSY_dir, f))
  }
}



lu = function(x)length(unique(x))


sigmoid = function(x, x0, slope){
  y = 1/(1+exp(-slope*(x-x0)))
  return(y)
}



compute_accuracy_by_group = function(decoding = rule_based_decoding, users = users, decoding_name = "rule-based"){
  tmp = full_join(decoding, users %>% dplyr::select(user_id, ends_with("_group")), by = "user_id")
  group_names = colnames(users %>% dplyr::select(ends_with("_group")))

  accuracy_by_group = foreach(group_name = group_names, .combine = bind_rows) %do% {
    tmp$group = tmp %>% select(matches(group_name)) %>% unlist()
    df = data.frame(decoding = decoding_name, group_type = group_name,
                    tmp %>% group_by(group) %>%
                      dplyr::summarize(accuracy = mean(correct_decoding, na.rm = TRUE)),
                    stringsAsFactors = FALSE)
    return(df)
  }

  accuracy_by_group = accuracy_by_group %>%
    dplyr::mutate(group = factor(group, levels = group))

  return(accuracy_by_group)

}




# event_accuracy = function(decoding, manual_labels, events){
#
#   event_accuracy = foreach(event = events, .combine = bind_rows) %do% {
#
#     jj = which(manual_labels$state_name == event);
#     j = jj[which(c(0, diff(manual_labels$rel_date[jj]))!=1)]
#
#     if(length(j)>0){
#       ml = manual_labels[j,]
#
#       this_event_accuracy = foreach(i = 1:nrow(ml), .combine = bind_rows) %do% {
#         D = decoding %>% dplyr::filter(user_id == ml$user_id[i])
#         if(event %in% c("Birth","Loss")){
#           jj = which(D$state_name == event)
#           j = jj[which(c(0, diff(D$rel_date[jj]))!=1)]
#           D = D[j,]
#         }else{
#           D = D %>% dplyr::filter(state_name == event)
#         }
#         d = D %>% dplyr::filter(rel_date %in% (-15:15+ml$rel_date[i]))
#
#         if(nrow(d) == 0){
#           decoded_rel_date = NA
#           time_diff = - 17
#         }else if(nrow(d)>1){
#           decoded_rel_date = d$rel_date[which.min(abs(d$rel_date - ml$rel_date[i]))]
#           time_diff = decoded_rel_date - ml$rel_date[i]
#         }else{
#           decoded_rel_date = d$rel_date
#           time_diff =  decoded_rel_date - ml$rel_date[i]
#         }
#         acc = data.frame(event = event,
#                          state_name = ml$state_name[i],
#                          user_id = ml$user_id[i],
#                          labelled_rel_date = ml$rel_date[i],
#                          decoded_rel_date = decoded_rel_date,
#                          time_diff = time_diff,
#                          stringsAsFactors = FALSE)
#
#       }else{
#         acc = data.frame()
#       }
#       return(acc)
#     }
#     return(this_event_accuracy)
#   }
#   return(event_accuracy)
# }



prepare_obs_deprecated = function(d, observations = c("bleeding","mucus","temperature","LH","preg","extra.log")){

  user_id = unique(d$user_id)
  if(length(user_id)>1){stop("several users in the input\n")}

  # we expand the feature matrix d so that is has one row per day
  #
  X = data.frame(user_id = user_id,rel_date = min(d$rel_date):max(d$rel_date))
  m = match(X$rel_date, d$rel_date)

  if("bleeding" %in% observations){
    # we keep the 0, 0.5, 1, 2, 3 code
    X = X %>% mutate(bleeding = d$bleeding[m],
                     first_day = d$first_day[m] %>% replace_na(FALSE))
    # and replace all missing values by 0 (no bleeding)
    X = X %>% mutate(bleeding = bleeding %>% replace_na(0))
    # and we also convert a "first day" into bleeding
    X  = X %>%
      mutate(first_day_bleeding = ifelse(first_day & (bleeding < 1),2,0),
             bleeding = pmax(bleeding, first_day_bleeding)) %>%
      select(-first_day)


  }

  if("mucus" %in% observations){
    X = X %>%
      mutate(mucus = mucus.dict$category[match(d$mucus_type[m], mucus.dict$names)] %>%  factor(., levels = levels(mucus.dict$category)))
  }

  if("temperature" %in% observations){

    # questionable temperatures are transformed into missing data
    d = d %>% mutate(
      quest_temp = ifelse(is.na(questionable_temp), FALSE, questionable_temp),
      temp = ifelse(quest_temp, NA, temperature) )

    # we remove temperature values that are oddly repeated
    if(any(!is.na(d$temp))){ # if there is at least one temperature
      # looking for oddly repeated values
      ttemp = sort(table(d$temp), decreasing = TRUE) # histogram of temperatures
      if((length(ttemp) == 1) | (ttemp[1]> (5*ttemp[2]))) {weird_temp = names(ttemp)[1] %>% as.numeric()} else{weird_temp = 999}
      d = d %>% mutate(temp = ifelse(temp == weird_temp, NA, temp))
    }

    # we scale the temperature
    median_temp = median(d$temp, na.rm = TRUE)
    d = d %>%  mutate(temp =  (temp - median_temp) %>% pmin(.,1.5) %>% pmax(.,-1.5))

    #
    X = X %>% mutate(temp = d$temp[m])
  }

  if("LH" %in% observations){
    X = X %>%
      mutate(
        LH_test = d$LH[m],
        LH = case_when(
          (LH_test == 0) ~ NA_real_,
          (LH_test == -1) ~ 0,
          (LH_test == 1) ~ 1,
          TRUE ~ NA_real_)
      ) %>% select(-LH_test)
  }

  if("preg" %in% observations){
    X = X %>%
      mutate(
        preg_test = d$preg_test[m],
        preg = case_when(
          (preg_test == 0) ~ NA_real_,
          (preg_test == -1) ~ 0,
          (preg_test == 1) ~ 1,
          TRUE ~ NA_real_)
      ) %>% select(-preg_test)
  }

  if("sex" %in% observations){
    X = X %>%  mutate(sex = d$sex_short[m])
  }

  if("extra.log" %in% observations){
    X = X %>% mutate(extra.log = (!is.na(d$sex_short[m]))*1)
  }

  X
}




source("Scripts/00_functions_feature_engineering.R")

create_observation_scores = function(d, features = c("bleeding"), return_tmp = FALSE){

  user_id = unique(d$user_id)
  if(length(user_id)>1){stop("several users in the input\n")}

  # we expand the feature matrix d so that is has one row per day
  # tmp is a data.frame which will store temporary variables needed to compute features
  tmp = data.frame(user_id = user_id,rel_date = min(d$rel_date):max(d$rel_date))

  # we get the list of the d columns that we need to compute the features
  input_cols = foreach(f = features)%do%{
    fun = eval(parse(text = str_c("compute_",f)))
    ic = fun(get_col_names_only = TRUE)$input_cols
  } %>% unlist() %>% unique()

  # we create the input data.frame, i.e. the d data.frame expanded to the size of the output.
  input = tmp
  d$tracking = TRUE
  input = full_join(input, d %>%  dplyr::select(rel_date,all_of(input_cols)), by = c("rel_date")) %>%  arrange(rel_date)
  input = input %>% mutate(tracking = tracking %>% replace_na(FALSE))

  # we re-order the features and the tmp_cols to account for the dependencies
  # first we need the tmp_cols
  tmp_cols = foreach(f = features)%do%{
    fun = eval(parse(text = str_c("compute_",f)))
    ic = fun(get_col_names_only = TRUE)$tmp_cols
  } %>% unlist() %>% unique()
  # then join them with the features
  all_features = unique(c(tmp_cols, features))
  # and order them
  ordered_features = c("bleeding","mucus","temp","LH","preg",
                       "bleeding_in_5d","susp_ano",
                       "long_high_temp",
                       "M_hint","bleeding_density_5d","acf",
                       "m_return","preg_hint",
                       "fmucus","gap_Xd")
  ordered_features = intersect(ordered_features, all_features)


  # Now we can compute the different features
  for(f in ordered_features){
    #cat(f,"\n")
    fun = eval(parse(text = str_c("compute_",f)))
    output = fun(input = input, tmp = tmp)
    output_cols = fun(input = input, tmp = tmp, get_col_names_only = TRUE)$output_cols
    # add to tmp
    cols_to_remove = intersect(output_cols,colnames(tmp))
    tmp = tmp %>% dplyr::select(-all_of(cols_to_remove))
    tmp = full_join(tmp, output, by = "rel_date")
  }
  dim(tmp)
  colnames(tmp)

  if(return_tmp){
    return(tmp)
  }else{
    obs = tmp %>% select(user_id, rel_date, all_of(features))
    return(obs)
  }
}


compute_log_likelihood = function(seq = rep(1,10), model = hsmm_fitted$model, obsdata, init_trans_sojourn_at_transitions = TRUE){

  L = length(seq)
  N = length(obsdata$N)

  seq_rle = rle(str_c(rep(1:N, obsdata$N) ,"_",seq))
  state_seq = seq_rle$values %>% str_replace(.,".*_","") %>% as.numeric()
  state_sojourn = seq_rle$lengths
  seq_id = seq_rle$values %>% str_replace(.,"_.*","") %>% as.numeric()

  loglik_base = rep(0, L)
  # Likelihood of the initial state : loglik_init
  loglik_init = loglik_base
  ii = c(0, cumsum(obsdata$N) %>% head(.,-1))+1; names(ii) = NULL
  if(init_trans_sojourn_at_transitions){
    loglik_init[ii] = model$init[seq[ii]] %>% log()
  }else{
    isj = c(1,which(diff(seq_id) != 0)+1); sojourn_first_state_of_sequences = state_sojourn[isj];
    iii = rep(ii, sojourn_first_state_of_sequences)
    i_first_state_of_sequences = iii+unlist(sapply(sojourn_first_state_of_sequences,function(x)0:(x-1)))
    nn = rep(sojourn_first_state_of_sequences, sojourn_first_state_of_sequences)
    loglik_init[i_first_state_of_sequences] = log(model$init[seq[iii]])/nn
  }


  # Likelihood from the transitions : loglik_trans
  loglik_trans = loglik_base
  if(length(state_seq)>1){
    it = (cumsum(state_sojourn)+1) %>% head(.,-1) # index of the sequence when the state changes (at the start of the new state)
    state_trans = data.frame(from = state_seq, to = dplyr::lead(state_seq)) %>% dplyr::filter(!is.na(to)) %>% as.matrix()
    if(init_trans_sojourn_at_transitions){
      loglik_trans[it] = hsmm$trans[state_trans] %>% log()
      loglik_trans[ii] = 0 # we remove the values at the start of a user sequence
    }else{
      it_seq = rep(it, state_sojourn[-1])
      state_trans_seq = state_trans[rep(1:nrow(state_trans),state_sojourn[-1]),]
      nn = rep(state_sojourn[-1], state_sojourn[-1])
      loglik_trans[it_seq] = log(hsmm$trans[state_trans_seq])/nn
      loglik_trans[ii_seq] = 0
    }
  }

  # Likelihood from the states sojourn: loglik_sojourn
  loglik_sojourn = loglik_base
  is = cumsum(state_sojourn) # index of the sequences just before a state change
  sojourn_i = cbind(state_sojourn,state_seq) # table with the length of each state along the sequence
  if(init_trans_sojourn_at_transitions){
    loglik_sojourn[is] = model$sojourn$d[sojourn_i] %>%  log()
    loglik_sojourn[ii[-1] - 1] = 0 # we ignore the duration of the last state of a sequence
    loglik_sojourn[L] = 0
    loglik_sojourn[ii] = 0 # and we also ignore the duration of the first state of a sequence
  }else{
    is = c(0, cumsum(state_sojourn) %>% head(.,-1))+1
    is_seq = rep(1:length(state_sojourn), state_sojourn)
    nn = rep(state_sojourn, state_sojourn)
    loglik_sojourn[1:L] = log(model$sojourn$d[sojourn_i[is_seq,]])/nn

    sojourn_of_last_state_of_sequences = state_sojourn[c(which(diff(seq_id) != 0)+1, length(seq_id))]
    i_last_state_of_sequences = rep(cumsum(obsdata$N), sojourn_of_last_state_of_sequences) +
      unlist(sapply(sojourn_of_last_state_of_sequences,function(x) -(x-1):0)); names(i_last_state_of_sequences) = NULL
    loglik_sojourn[i_last_state_of_sequences] = 0
    loglik_sojourn[i_first_state_of_sequences] = 0
  }

  # Likelihood of the observations: loglik_obs
  x = obsdata$x
  p_o = sapply(1:model$J,function(state) model$dens.emission(obs = x,state = state, glm_model = model$glm_models[[state]], parem = model$parms.emission))
  loglik_obs = (p_o[cbind(1:L,seq)]) %>% log()

  # putting it together
  loglik = loglik_init + loglik_trans + loglik_sojourn + loglik_obs

  # K = 300
  # plot(loglik_init[1:K], type = "l", ylim = range(loglik[1:K], na.rm = TRUE))
  # points(loglik_trans[1:K], type = "l", col = "red")
  # points(loglik_sojourn[1:K], type = "l", col = "blue")
  # points(loglik_obs[1:K], type = "l", col = "green")

  LLs = list(loglik = loglik, loglik_init = loglik_init, loglik_trans = loglik_trans, loglik_sojourn = loglik_sojourn,  loglik_obs = loglik_obs)

  return(LLs)
}



compute_confidence_score = function(seq = rep(1,10), model = hsmm_fitted$model, obsdata = obsdata){

  LL = compute_log_likelihood(seq = seq, model = model, obsdata = obsdata)



}



compute_confidence_score_2 = function(seq = rep(1,10), model = hsmm_fitted$model, obsdata = obsdata, N_bootstrap = 200){

  LL = compute_log_likelihood(seq = seq, model = model, obsdata = obsdata)


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
  loglik_obs_r1 = sapply(1:N, function(n) { # for each bootstrap iteration
    p_obs_r = sapply(1:n_states, function(state) model$dens.emission(obs = x_r[state,,n],
                                                                     state = state,
                                                                     parem = model$parms.emission,
                                                                     w = model$weights[state,]))
    w = model$weights
    ww = rowSums(w * (!is.na(x_r[,,n])))
    return(ww * p_obs_r)
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





my_hsmm_bluid = function(hsmm){

  SOJOURN = hsmm$sojourn/rowSums(hsmm$sojourn)
  SOJOURN = t(SOJOURN)

  this_user_hsmm =  specify_hsmm(
    J = hsmm$n_states,
    init = hsmm$init,
    trans = hsmm$trans_no_names,
    marg_em_probs = hsmm$emission_par,
    censoring_probs = hsmm$censoring_probs,
    sojourn = list(d = SOJOURN, type = "nonparametric"),
    augment_data_fun = hsmm$augment_data_fun,
    state_names = hsmm$states$abbr,
    state_colors = hsmm$states$colors)

  return(this_user_hsmm)
}



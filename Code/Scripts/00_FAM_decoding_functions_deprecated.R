
get_most_likely_sequence_with_prob = 
  function(X, 
           R_model, 
           M_model, T_model,
           is_hmm = FALSE,
           fit_models = FALSE, 
           verbose = FALSE){
  decoding = decode_user_timeseries(X, R_model = R_model, M_model = M_model, T_model = T_model, is_hmm = is_hmm, fit_models = fit_models, verbose = verbose)
  decoding %>% filter(is_on_most_likely_path) %>% select(-is_on_most_likely_path, -common_overlap)
}


decode_user_timeseries = function(X, R_model, M_model, T_model, is_hmm = FALSE, fit_models = FALSE, verbose = FALSE){
  res_M = decode_with_M_model(X, M_model = M_model, verbose = verbose)
  stretches = identify_tracking_categories(X, res_M = res_M, T_model = T_model, verbose = verbose)
  decodings = purrr::map_dfr(.x = 1:nrow(stretches), 
                             .f = decode_stretch, 
                             stretches = stretches, 
                             X = X, 
                             res_M = res_M,
                             fit_model = fit_models,
                             is_hmm = is_hmm,
                             R_model = R_model,
                             verbose = verbose)
  decoding = stitch_stretches(decodings, verbose = verbose)
  if(any(!decoding$common_overlap, na.rm = TRUE)) decoding = fix_stitches(decoding = decoding, X = X, res_M = res_M, R_model = R_model, T_model = T_model,is_hmm = is_hmm, fit_models = fit_models, verbose = verbose)
  decoding
}


decode_with_M_model = function(X, M_model, verbose = FALSE){
  if(verbose) cat("identifying menses\n")
  
  # 1. prepare observations for menses-identifying model
  XM = X %>% select(seq_id, t, bleeding)
  
  # 2. decode with menses-identifying model
  vit_M = predict_states_hsmm(model = M_model, X = XM, method = "Viterbi")
  smoo_M = predict_states_hsmm(model = M_model, X = XM, method = "FwBw")
  if(nrow(smoo_M$state_seq) == 0) smoo_M = predict_states_hsmm(model = M_model, X = XM, method = "FwBw", ground_truth = vit_M$state_seq)
  
  # 3. res_M
  res_M =  left_join(vit_M$state_seq %>% select(seq_id, t, state),
                     smoo_M$probabilities %>% select(seq_id, t, state, state_prob) %>% rename(prob = state_prob),
                     by = c("seq_id","t","state"))
   
  res_M
}


identify_tracking_categories = function(X, res_M, T_model, verbose = FALSE){
  if(verbose) cat("Identifying stretches of specific tracking behaviors\n")

  # 1. prepare observations for tracking categories model
  XT = X %>% 
    left_join(., res_M, by = c("seq_id", "t")) %>% 
    mutate(any.menses = ((state == 1) & (prob > 0.75)) * 1,
           any.preg = (!is.na(preg))*1,
           any.LH = (!is.na(LH))*1,
           any.mucus = (!is.na(mucus))*1,
           any.temp = (!is.na(temp))*1) %>% 
    select(seq_id, t, starts_with("any"))
  
  # 2. decode with T_model
  smoo =  predict_states_hsmm(model = T_model, X = XT, method = "FwBw")
  
  # 3. identify the stretches
  tn = which(T_model$state_names == "transition")
  stretches = smoo$state_seq %>%  
    # define a stretch_id
    mutate(change = (state != lag(state)) %>% replace_na(TRUE),
           stretch_id = cumsum(change)) %>% 
    group_by(stretch_id) %>% 
    # identify starts and ends
    mutate(start = min(t),
           end = max(t)) %>%
    filter(t == start) %>% 
    select(seq_id, state, stretch_id, start, end) %>% 
    ungroup() %>% 
    # make sure there is a transition state between each other states
    mutate(has_transition_state = (((state != tn) & lead(state) == tn)) | ((state == tn) & (lag(state) != tn)),
           has_transition_state = has_transition_state %>% replace_na(TRUE))
  
  stretches = stretches %>% 
    bind_rows(., stretches %>% filter(!has_transition_state) %>% mutate(state = tn, start = end-1, end = end+1)) %>% 
    select(-has_transition_state) %>% 
    arrange(seq_id, start) %>% 
    mutate(stretch_id = row_number()) %>% 
    # remove transition state
    mutate(start2 = lag(start),
           end2 = lead(end)) %>% 
    filter(state != tn) %>% 
    mutate(start = ifelse(is.na(start2),start, start2),
           end = ifelse(is.na(end2), end, end2)) %>% 
    select(seq_id,state, start, end) %>% 
    # remove consecutive duplicates
    mutate(stretch_id = (state != lag(state)) %>% replace_na(TRUE) %>% cumsum()) %>% 
    group_by(seq_id, stretch_id, state) %>% 
    summarize(start = min(start),
              end = max(end),
              length = end - start,
              .groups = "drop") %>% 
    # add the state/model name
    mutate(tracking_behavior = T_model$state_names[state]) %>% 
    select(-state)
  
  #4. returns the stretches
  stretches
}



decode_stretch = function(i, stretches, X, res_M, R_model, is_hmm = FALSE, fit_model = FALSE, verbose = FALSE){
  stretch = stretches[i,]
  if(verbose) cat("Decoding stretch",i," -  model:",stretch$tracking_behavior," -  [",stretch$start,":",stretch$end,"] length =",stretch$length," ,\n")
  
  
  # 1. prepare data
  X = X %>% dplyr::filter(t %in% stretch$start:stretch$end)
  
  # 2. customize model
  model = adjust_model(tracking_behavior = stretch$tracking_behavior, X = X, R_model = R_model, is_hmm = is_hmm, verbose = verbose)
  
  # 3. remove any observation that is not part of the model
  X = X %>% select(seq_id, t, all_of(names(model$marg_em_probs)))
  
  # 4. use res_M as ground_truth
  GT = res_M %>% filter(state %in% which(R_model$state_names %in% c("M","B")), prob > 0.75, t %in% stretch$start:stretch$end)
  
  # 4. fit model
  if(fit_model){
    if(verbose) cat("Fitting the model to this stretch\n")
    model_fit = fit_hsmm(model = model, X = X, ground_truth = GT,  lock_transition = TRUE, lock_sojourn = FALSE, verbose = FALSE, N0 = 100)
    model = model_fit$model
  }
  
  # 5. decode observations
  if(verbose) cat("Decoding the stretch\n")
  vit = predict_states_hsmm(X = X, model = model, ground_truth = GT, method = "Viterbi")
  smoo = predict_states_hsmm(X = X, model = model, ground_truth = GT, method = "FwBw")
  
  # return results
  res = smoo$probabilities %>%
    select(seq_id, t, state, state_prob) %>% rename(prob = state_prob) %>% 
    left_join(.,vit$state_seq %>% select(seq_id, t, state) %>% mutate(is_on_most_likely_path = TRUE), by = c("seq_id", "t", "state")) %>% 
    mutate(is_on_most_likely_path = is_on_most_likely_path %>% replace_na(FALSE),
           tracking_behavior = stretch$tracking_behavior,
           stretch_id = stretch$stretch_id)
  res
}



adjust_model = function(tracking_behavior, X, R_model, is_hmm = FALSE, verbose = FALSE){
  if(verbose) cat("Adjusting model to observations\n")
  
  # we start from the general FAM model
  model = R_model

  # we remove variables that are never observed
  var_names = names(model$marg_em_probs)
  Xmissing = X %>% select(all_of(var_names)) %>% mutate(across(everything(), is.na))
  missing_frequencies = Xmissing %>% summarize(across(everything(), mean))
  tracking_frequencies = 1-missing_frequencies
  if(any(tracking_frequencies == 0)){
    vars = var_names[which(tracking_frequencies == 0)]
    for(var in vars){
      model$marg_em_probs[[var]] = NULL
    }
  }
  
  # Do we need to fix the sojourn of Lut ?
  if(!is_hmm & (tracking_behavior %in% c("b","bp"))) model$sojourn$Lut$d = c(rep(0,10),1,rep(0, length(model$sojourn$Lut$d)-11))
  
  
  # Do we need to fix the sojourn of hE ?
  if(!is_hmm & (tracking_behavior %in% c("b","bp","b_tests","no_mucus"))) model$sojourn$hE$d = c(rep(0,2),1,rep(0, length(model$sojourn$hE$d)-3))
  
  # Do we need to modify the transition from hE to lE?
  if(tracking_behavior %in% c("b")) model$transition["hE","lE"] = 0
  
  # Do we need to cancel the transition to Ano?
  if(!(tracking_behavior %in% c("full","no_mucus"))) model$transition[,"Ano"] = 0
  
  model$transition = model$transition/rowSums(model$transition)
  
  # censoring probs
  model$censoring_probs = list(
    p = rep(0.5, model$J),
    q = matrix(0.5, ncol = model$J, nrow = length(model$marg_em_probs))
  )
  model$censoring_probs$q[1,] = 0
  model$censoring_probs$p[model$state_names %in% c("M","Ano","L")] = 0.45
  model$censoring_probs$p[model$state_names %in% c("AB")] = 0.15
  model$censoring_probs$p[model$state_names %in% c("PP","BF")] = 0.55
  
  # specify the model
  model = specify_hsmm(J = model$J, init = model$init, transition = model$transition, 
                       sojourn = model$sojourn, marg_em_probs = model$marg_em_probs, 
                       # censoring_probs = NULL,
                       # censoring_probs = list(
                       #   p = rep(0,model$J),
                       #   q = matrix(0.5, ncol = model$J, nrow = length(model$marg_em_probs))
                       # ),
                       censoring_probs = model$censoring_probs,
                       state_names = model$state_names, state_colors = model$state_colors)
  
 
  # return adjusted model
  model
}


stitch_stretches = function(decodings, verbose = FALSE){
  if(verbose) cat("Stitching stretches\n")
  
  decoding = decodings %>% 
    arrange(seq_id, t, state) %>% 
    group_by(seq_id, t, state) %>% 
    mutate(is_transition = n() > 1)
  
  if(!any(decoding$is_transition)) return(decodings %>%  mutate(common_overlap = NA))
  
  decoding_outside_transitions = decoding %>%  filter(!is_transition) %>% select(-is_transition)
  
  decoding_in_transitions = 
    decoding %>% filter(is_transition) %>% 
    group_by(seq_id, t) %>% 
    mutate(transition_id = mean(stretch_id),
           tracking_behavior = str_c(sort(unique(tracking_behavior)), collapse = " - ")) %>% 
    group_by(seq_id, transition_id) %>% 
    mutate(transition_start = min(t),
           transition_end = max(t))
  
  overlaps_in_transitions = 
    decoding_in_transitions %>% 
    filter(is_on_most_likely_path) %>% 
    group_by(seq_id, transition_id, t) %>% 
    summarize(same_state = (length(unique(state)) == 1),
              .groups = "drop") %>% 
    filter(same_state) %>% 
    group_by(seq_id, transition_id) %>% 
    summarize(overlap_start = suppressWarnings(min(t)),
              overlap_end = suppressWarnings(max(t)), 
              .groups = "drop")
  
  most_likely_path_in_transitions = 
    decoding_in_transitions %>% 
    filter(is_on_most_likely_path) %>%  
    select(seq_id, t, state, stretch_id, transition_id)  %>% 
    left_join(., overlaps_in_transitions, by = c("seq_id", "transition_id")) %>% 
    group_by(seq_id, t) %>% 
    mutate(selected_stretch = ifelse(t <= overlap_end, min(stretch_id), max(stretch_id)),
           state = ifelse(is.na(selected_stretch), NA, state),
           selected_stretch = selected_stretch %>% replace_na(min(stretch_id))) %>% 
    filter(stretch_id == selected_stretch) %>% 
    select(seq_id, t, state) %>% 
    mutate(is_on_most_likely_path = TRUE,
           common_overlap = ifelse(any(is.na(state)), FALSE, TRUE))
  
  decoding_in_transitions = decoding_in_transitions %>% 
    select(seq_id, t, state, prob, tracking_behavior, transition_id) %>% 
    rename(stretch_id = transition_id) %>% 
    group_by(seq_id, t, state, tracking_behavior, stretch_id) %>% 
    summarize(prob = mean(prob), .groups = "drop") %>% 
    left_join(. , most_likely_path_in_transitions, by = c("seq_id","t","state")) %>% 
    group_by(seq_id, stretch_id) %>% 
    mutate(is_on_most_likely_path = is_on_most_likely_path %>% replace_na(FALSE),
           common_overlap = any(common_overlap, na.rm = TRUE))
  
  decoding = bind_rows(decoding_outside_transitions %>% mutate(common_overlap = NA),
                       decoding_in_transitions) %>% 
    ungroup()
  
  decoding
}


fix_stitches = function(decoding, X, res_M = res_M, R_model, T_model, is_hmm = FALSE, fit_models = FALSE, verbose = FALSE){
  
  if(verbose) cat("Fixing stitches \n")
  
  # 1. Identify if there are any overlap problems
  any_overlap_problem = any(!decoding$common_overlap, na.rm = TRUE)
  
  # 2. If no overlap problem, return decoding as is.
  if(!any_overlap_problem) return(decoding)
  
  # 3. If there are overlap problems, identify the new stretches that need to be decoded.
  
  cycles = decoding %>% 
    filter(is_on_most_likely_path) %>% 
    arrange(seq_id, t) %>% 
    mutate(first_cycleday = (prob > 0.75) & (state == 1) & ((lag(state) != 1) | (t == 1)),
           cycle_nb = cumsum(first_cycleday)) %>%  
    filter(first_cycleday, prob > 0.75) %>% 
    select(seq_id, t, stretch_id, cycle_nb)
  
  overlaps_to_fix = decoding %>% 
    filter(!common_overlap) %>% 
    select(seq_id, t, stretch_id, common_overlap, tracking_behavior) %>% 
    distinct() %>% 
    group_by(seq_id, stretch_id, tracking_behavior) %>% 
    summarize(overlap_start = min(t),
              overlap_end = max(t),
              .groups = "drop") %>% 
    mutate(o_id = row_number())
  
  stretches_to_fix = purrr::map_dfr(
    .x = overlaps_to_fix$o_id,
    .f = function(o){
      cbind(
        overlaps_to_fix %>% filter(o_id == o) %>%  select(seq_id, stretch_id, o_id, tracking_behavior) %>% 
          mutate(tracking_behavior = tracking_behavior %>% 
                   str_split_fixed(.," - ", 2) %>% 
                   factor(., levels = T_model$state_names) %>%  
                   sort(., decreasing = TRUE) %>% 
                   head(1) %>% 
                   as.character()
                   ),
        cycles %>% filter(t < overlaps_to_fix$overlap_start[o]) %>% 
          summarize(start = suppressWarnings(max(t))) %>% mutate(start = ifelse(is.infinite(start), 1, start)),
        cycles %>% filter(t > overlaps_to_fix$overlap_end[o]) %>% 
          summarize(end = suppressWarnings(min(t))) %>%  mutate(end = ifelse(is.infinite(end), max(X$t), end))
      )
    }
  )
  
  # 3b. We make sure each stretch is unique
  if(nrow(stretches_to_fix) > nrow(stretches_to_fix %>% select(seq_id, start, end) %>% distinct())){
    stretches_to_fix = stretches_to_fix %>% 
      mutate(tracking_behavior_num = tracking_behavior %>% factor(., levels = T_model$state_names) %>%  as.numeric()) %>% 
      arrange(seq_id, start, end, tracking_behavior_num) %>% 
      group_by(seq_id, start, end) %>% 
      summarize(tracking_behavior_num = max(tracking_behavior_num),
                stretch_id = max(stretch_id),
                o_id = max(o_id),
                .groups = "drop") %>% 
      mutate(tracking_behavior = T_model$state_names[tracking_behavior_num]) %>% 
      select(-tracking_behavior_num)
  }
  
  # 4. For each stretch, decode with decode_stretch() function
  decodings_in_stretches_to_fix =
    purrr::map_dfr(.x = 1:nrow(stretches_to_fix), 
                   .f = decode_stretch, 
                   stretches = stretches_to_fix, X = X, 
                   res_M = res_M,
                   R_model = R_model,
                   is_hmm = is_hmm,
                   fit_model = fit_models,
                   verbose = verbose)

  # 5. Remove the time-points that are overlapping with the new stretches, append the new decodings (from step 4)
  decoding = decoding %>% filter(!(t %in% decodings_in_stretches_to_fix$t)) %>% 
    bind_rows(., decodings_in_stretches_to_fix)
  
  # 6. Sort by seq_id, t, state
  decoding = decoding %>% arrange(seq_id, t, state)
    
  # 7. Return results
  decoding
}


decode_with_model = function(m, verbose = FALSE) {
  cat(m, "\n")
  this_model = models[m,]
  res_file = paste0(output_dir, this_model$file_name)
  
  tic()
  if (!file.exists(res_file)){
    eval(parse(text = paste0("model = ",this_model$model)))
    if (this_model$approach == "adaptative"){
      res = 
        decode_with_adaptation(
          X_all_users, model, is_hmm = ifelse(this_model$model == "R_hmm", TRUE, FALSE), 
          this_model$fit_model, verbose
          )
    } else {
      res = decode_without_adaptation(X_all_users, model, this_model$fit_model, verbose)
    }
    write_feather(res, path = res_file)
  }
  a = toc(quiet = TRUE)
  this_model$exec_time = (a$toc - a$tic) %>% as.numeric()
  this_model
}



decode_with_adaptation = function(X, model, is_hmm, fit_model, verbose = FALSE){
  purrr::map_dfr(
    .x = unique(X$seq_id),
    .f = function(sid){
      if (verbose) cat(sid, "\n")
      get_most_likely_sequence_with_prob(
        X = X %>% filter(seq_id == sid), 
        model = model,
        is_hmm = is_hmm,
        fit_model = fit_model,
        verbose = FALSE)
    }
  )
}



get_most_likely_sequence_with_prob = 
  function(X, 
           model, 
           is_hmm = FALSE,
           fit_model = FALSE, 
           verbose = FALSE){
    decoding = decode_user_timeseries(X, model = model, is_hmm = is_hmm, fit_model = fit_model, verbose = verbose)
    decoding %>% filter(is_on_most_likely_path) %>% select(-is_on_most_likely_path)
  }


decode_user_timeseries = function(X, model, is_hmm = FALSE, fit_model = FALSE, verbose = FALSE){
  res_M =  identify_cycles_from_bleeding(X, model, is_hmm, verbose)
  GT = ground_truth_from_res_M(res_M, model)
  stretches = identify_tracking_categories(X, res_M = res_M, verbose = verbose)
  decodings = decode_stretches(X, model, is_hmm, GT = GT, stretches, fit_model = fit_model, verbose = verbose)
  # ggplot(decodings %>%  filter(is_on_most_likely_path) %>% mutate(stretch_number = stretch_number %>%  factor()),
  #        aes(x = t, y = state, col = stretch_number)) + 
  #   geom_line()
  decoding = stitch_stretches(decodings, X, model, is_hmm, fit_model,  verbose = verbose)
  # ggplot(decoding %>%  filter(is_on_most_likely_path) %>% mutate(stretch_number = stretch_number %>%  factor()),
  #        aes(x = t, y = state)) + 
  #   geom_line() + geom_point(size = 0.5, aes(col = stretch_number))
  decoding
}

ground_truth_from_res_M = function(res_M, model){
  res_M %>% 
    filter(
      state %in% which(model$state_names %in% c("M","B")), 
      prob > 0.75, 
    ) %>% 
    mutate(
      trust = 0.75
    ) %>% 
    select(seq_id, t, state, trust)
}

decode_stretches = function(X, model, is_hmm = FALSE, GT = data.frame(), stretches, fit_model = FALSE, verbose = FALSE) {
  purrr::map_dfr(.x = 1:nrow(stretches), 
                 .f = decode_stretch, 
                 stretches = stretches, 
                 X = X, 
                 GT = GT,
                 fit_model = fit_model,
                 is_hmm = is_hmm,
                 model = model,
                 verbose = verbose)
}


identify_cycles_from_bleeding = function(X, model, is_hmm, verbose = FALSE){
  if(verbose) cat("identifying menses\n")
  
  # 1. prepare observations and menses-identifying model
  XM = X %>% mutate(LH = NA, preg = NA, mucus = NA, temp = NA)
  model = adjust_model(tracking_behavior = "B", X = XM, model = model, is_hmm = is_hmm, verbose = verbose)
  
  # 2. decode with menses-identifying model
  vit_M = predict_states_hsmm(model = model, X = XM, method = "Viterbi")
  smoo_M = predict_states_hsmm(model = model, X = XM, method = "FwBw")
  if(nrow(smoo_M$state_seq) == 0) smoo_M = predict_states_hsmm(model = model, X = XM, method = "FwBw", ground_truth = vit_M$state_seq)
  
  # 3. res_M
  res_M =  
    left_join(vit_M$state_seq %>% select(seq_id, t, state),
              smoo_M$probabilities %>% 
                select(seq_id, t, state, state_prob) %>% 
                rename(prob = state_prob),
              by = c("seq_id","t","state"))
  
  # 4. identify cycle starts
  res_M = 
    res_M %>% 
    group_by(seq_id) %>% 
    mutate(
      menses = ((state == 1) & (prob > 0.75)) * 1,
      is_day_1 = menses & (!lag(menses) | is.na(lag(menses))),
      cycle_number = cumsum(is_day_1)
    ) %>% 
    ungroup() 
  
  
  res_M
}


identify_tracking_categories = function(X, res_M, verbose = FALSE){
  if(verbose) cat("Identifying stretches of specific tracking behaviors\n")
  
  # 1. identify cycle starts
  XT = X %>% 
    left_join(., res_M, by = c("seq_id", "t")) %>% 
    group_by(seq_id, cycle_number) %>% 
    mutate(
      cycle_start = min(t),
      cycle_end = max(t),
      cycle_length = n(),
      period_length = sum(menses),
      group = str_c(seq_id,"_",cycle_number)
    )
  
  # 2. for each cycle, determine the tracking category
  tracking_behavior_per_cycle =
    compute_tracking_behavior(XT)
  
  # 3. group cycles with the same category together
  stretches = 
    XT %>% 
    group_by(group) %>%
    slice_head(n = 1) %>% 
    arrange(t) %>% 
    left_join(
      .,
      tracking_behavior_per_cycle %>%  select(group, tracking_behavior), 
      by = "group") %>%     
    group_by(seq_id) %>% 
    mutate(
      stretch_start = 
        (tracking_behavior != lag(tracking_behavior)) | 
        (lag(tracking_behavior) %>% is.na()),
      stretch_number = cumsum(stretch_start)) %>% 
    group_by(seq_id, stretch_number, tracking_behavior) %>% 
    summarize(
      start = min(cycle_start),
      end = max(cycle_end),
      first_period_duration = period_length[1],
      .groups = "drop"
    ) %>% 
    mutate(
      end = end + ifelse(is.na(lead(first_period_duration)), 0, lead(first_period_duration)),
      length = end - start
    ) %>% 
    select(seq_id, stretch_number, start, end, length, tracking_behavior) %>% 
    filter(length > 0)
  
  #4. returns the stretches
  stretches
}



compute_tracking_behavior <- function(df){
  df %>% 
    group_by(group) %>% 
    summarize(
      start = min(t),
      end = max(t), 
      length = n(),
      across(.cols = c("mucus", "temp","LH", "preg"),
             .fns = function(x) sum(!is.na(x)),
             .names = "n_{.col}"),
      n_cycles = length(unique(str_c(seq_id,"_",cycle_number))),
      .groups = "drop") %>% 
    mutate(
      tracking_behavior = 
        case_when(
          (n_mucus >= 5*n_cycles) & (n_temp >= 8*n_cycles) ~ "BTM",
          (n_mucus >= 5*n_cycles) ~ "BM",
          (n_temp >= 8*n_cycles) ~ "BT",
          (n_LH >= n_cycles) ~ "BLH",
          (n_preg >= n_cycles) ~ "BP",
          TRUE ~ "B"
        )
    )
}


decode_stretch = function(i, stretches, X, GT = data.frame(), model, is_hmm = FALSE, fit_model = FALSE, start_end_on_menses = FALSE, verbose = FALSE){
  stretch = stretches[i,]
  if(verbose) cat("Decoding stretch",stretch$stretch_number," -  model:",stretch$tracking_behavior," -  [",stretch$start,":",stretch$end,"] length =",stretch$length," ,\n")
  
  
  # 1. prepare data
  X = X %>% dplyr::filter(t %in% stretch$start:stretch$end)
  GT = GT %>% dplyr::filter(t %in% stretch$start:stretch$end)
  trust_in_ground_truth = min(GT$trust)
  
  # 2. customize model
  model = adjust_model(tracking_behavior = stretch$tracking_behavior, X = X, model = model, is_hmm = is_hmm, verbose = verbose)
  
  # 3. remove any observation that is not part of the model
  X = X %>% select(seq_id, t, all_of(names(model$marg_em_probs)))
  
  # 4. fit model
  if(fit_model){
    if(verbose) cat("Fitting the model to this stretch\n")
    model_fit = fit_hsmm(model = model, X = X, ground_truth = GT, trust_in_ground_truth = trust_in_ground_truth,  lock_transition = TRUE, lock_sojourn = FALSE, verbose = FALSE, N0_emission = 100)
    model = model_fit$model
  }
  
  # 5. decode observations
  if(verbose) cat("Decoding the stretch\n")
  vit = predict_states_hsmm(X = X, model = model, ground_truth = GT, trust_in_ground_truth = trust_in_ground_truth, method = "Viterbi")
  smoo = predict_states_hsmm(X = X, model = model, ground_truth = GT, trust_in_ground_truth = trust_in_ground_truth, method = "FwBw")
  
  # return results
  res = smoo$probabilities %>%
    select(seq_id, t, state, state_prob) %>% rename(prob = state_prob) %>% 
    left_join(.,vit$state_seq %>% select(seq_id, t, state) %>% mutate(is_on_most_likely_path = TRUE), by = c("seq_id", "t", "state")) %>% 
    mutate(is_on_most_likely_path = is_on_most_likely_path %>% replace_na(FALSE),
           tracking_behavior = stretch$tracking_behavior,
           stretch_number = stretch$stretch_number)
  res
}



adjust_model = function(tracking_behavior, X, model, is_hmm = FALSE, verbose = FALSE){
  if(verbose) cat("Adjusting model to observations\n")
  
  # we remove variables that are never observed
  var_names = names(model$marg_em_probs)
  Xmissing = X %>% select(all_of(var_names)) %>% mutate(across(everything(), is.na))
  missing_frequencies = Xmissing %>% summarize(across(everything(), mean))
  all_tracking_frequencies = 1-missing_frequencies
  tracking_frequencies = all_tracking_frequencies
  if(any(tracking_frequencies == 0)){
    j = which(tracking_frequencies == 0)
    vars = var_names[j]
    tracking_frequencies = tracking_frequencies[-j]
    for(var in vars){
      model$marg_em_probs[[var]] = NULL
    }
  }
  tracking_frequency = 1 - (apply(Xmissing, 1, prod) %>% mean())
  
  # Do we need to fix the sojourn of Lut ?
  if(!is_hmm & (tracking_behavior %in% c("B","BP"))) 
    model$sojourn$Lut$d = c(rep(0,10),1,rep(0, length(model$sojourn$Lut$d)-11))
  
  # Do we need to fix the sojourn of hE ?
  if(!is_hmm & (tracking_behavior %in% c("B","BP","BLH","BT"))) 
    model$sojourn$hE$d = c(rep(0,2),1,rep(0, length(model$sojourn$hE$d)-3))
  
  # Do we need to modify the transition from hE to lE?
  if(all_tracking_frequencies$mucus > 0)
    model$transition["hE","lE"] = model$transition["hE","lE"] * all_tracking_frequencies$mucus
  
  if(tracking_behavior %in% c("BP"))
    model$transition["hE","lE"] = 0.001

  if(tracking_behavior %in% c("B"))
    model$transition["hE","lE"] = 0
  

  # Do we need to cancel the transition to Ano?
  if(!(tracking_behavior %in% c("BTM","BT"))) 
    model$transition[,"Ano"] = 0
  
  model$transition = model$transition/rowSums(model$transition)
  
  # censoring probs
  p = 1-tracking_frequency
  q = 1-unlist(tracking_frequencies)
  model$censoring_probs = list(
    p = rep(p, model$J) %>% set_names(model$state_names),
    q = matrix(q, 
               ncol = model$J, nrow = length(model$marg_em_probs), 
               dimnames = list(names(model$marg_em_probs), model$state_names))
  )
  # probability of all variables missing
  model$censoring_probs$p["M"] = min(0.2,p*0.4)
  model$censoring_probs$p[c("Ano","L")] = min(0.3, p*0.7)
  model$censoring_probs$p[c("AB")] = 0.2
  #model$censoring_probs$p[c("PP","BF")] = p+(1-p)*0.1
  
  # probability of specific variables missing
  model$censoring_probs$q["bleeding",] = 0
  if ("preg" %in% names(model$marg_em_probs)) {
    model$censoring_probs$q["preg",c("P")] = model$censoring_probs$q["preg","P"]*0.5
    model$censoring_probs$q["preg",c("Lut","PB1","PL")] = model$censoring_probs$q["preg","Lut"]*0.8
  }
  if ("temp" %in% names(model$marg_em_probs)){
    model$censoring_probs$q["temp","P"] = model$censoring_probs$q["temp","P"]*0.8
    model$censoring_probs$q["temp","Ano"] = model$censoring_probs$q["temp","Ano"]*0.5
  }
  
  
  # specify the model
  model = specify_hsmm(J = model$J, init = model$init, transition = model$transition, 
                       sojourn = model$sojourn, 
                       marg_em_probs = model$marg_em_probs, 
                       censoring_probs = model$censoring_probs,
                       state_names = model$state_names, state_colors = model$state_colors)
  
  
  # return adjusted model
  model
}


stitch_stretches = function(decodings, X, model, is_hmm = FALSE, fit_model = FALSE, verbose = FALSE){
  if(verbose) cat("Stitching stretches\n")
  
  decodings = 
    decodings %>% 
    arrange(seq_id, t, state) %>% 
    group_by(seq_id, t, state) %>% 
    mutate(is_transition = n() > 1,
           transition_number = mean(stretch_number)) %>% 
    ungroup()
  
  if(!any(decodings$is_transition)) return(decodings)
  
  
  # there might be discrepancies within the transitions 
  # if one of them does not have "menses" in their most likely path
  # so, for each transition, we check if "menses" is on the most likely path 
  # for both stretch on each side of the transition
  
  transitions = 
    decodings %>% 
    filter(is_transition, is_on_most_likely_path) %>% 
    group_by(seq_id, transition_number, stretch_number) %>% 
    summarize(any_menses = any(state == 1), .groups = "drop") %>% 
    group_by(seq_id, transition_number) %>% 
    summarize(has_menses_in_both_stretches = all(any_menses), .groups = "drop") %>% 
    mutate(needs_fixing = !has_menses_in_both_stretches)
  
  decoding_outside_transitions = 
    decodings %>% 
    filter(!is_transition) %>% 
    mutate(needs_fixing = FALSE)
  
  decoding_in_transitions = 
    decodings %>% 
    filter(is_transition) %>% 
    left_join(., 
              transitions %>% 
                select(seq_id, transition_number, needs_fixing),
              by = c("seq_id", "transition_number")) %>% 
    group_by(seq_id, t, transition_number, state) %>% 
    summarize(prob = mean(prob),
              tracking_behavior = str_c(unique(tracking_behavior), collapse = "-"),
              needs_fixing = any(needs_fixing),
              .groups = "drop")  %>% 
    arrange(seq_id, t, -prob) %>%
    group_by(seq_id, t) %>% 
    mutate(is_on_most_likely_path = c(TRUE, rep(FALSE, n()-1))) %>% 
    ungroup()
  
  decoding =
    bind_rows(
      decoding_outside_transitions,
      decoding_in_transitions
    ) %>% 
    mutate(stretch_number = transition_number) %>% 
    arrange(seq_id, t, state) %>% 
    select(-is_transition, -transition_number)
  
  if(any(decoding$needs_fixing, na.rm = TRUE)) 
    decoding = fix_stitches(decoding, X, model, is_hmm, fit_model, verbose)
  
  decoding  = decoding %>%  select(-needs_fixing)
  decoding
}


fix_stitches = function(decoding, X, model, is_hmm = FALSE, fit_model = FALSE, verbose = FALSE){
  
  if(verbose) cat("Fixing stitches \n")
  
  # 1. Identify the new stretches that need to be re-decoded.
  to_fix = 
    decoding %>% 
    filter(needs_fixing) %>% 
    select(seq_id, t, stretch_number, tracking_behavior) %>% 
    distinct() %>% 
    group_by(seq_id, stretch_number, tracking_behavior) %>% 
    summarize(stretch_start = min(t),
              stretch_end = max(t),
              .groups = "drop") 
  
  # 2. extend the stretch 
  # from the last menses of the previous stretch
  # to the first menses of the next stretch
  
  cycles = 
    decoding %>% 
    filter(is_on_most_likely_path) %>% 
    arrange(seq_id, t) %>% 
    mutate(day_1 = (prob > 0.75) & (state == 1) & ((lag(state) != 1) | (t == 1)),
           cycle_number = cumsum(day_1),
           cycle_start = t) %>%  
    filter(day_1) %>% 
    select(seq_id, cycle_start, stretch_number, cycle_number)
  
  
  stretches_to_fix = purrr::map_dfr(
    .x = 1:nrow(to_fix),
    .f = function(i){
      prev_cycles = 
        cycles %>% 
        filter(cycle_start < to_fix$stretch_start[i]) %>% 
        slice_tail(n = 1)
      min_t = max(prev_cycles$cycle_start, 1, na.rm = TRUE)
      
      next_cycles =
        cycles %>% 
        filter(cycle_start > to_fix$stretch_end[i]) %>% 
        slice_head(n = 1)
      max_t = min(next_cycles$cycle_start + 5, max(X$t), na.rm = TRUE)
      
      data.frame(
        seq_id = to_fix$seq_id[i],
        stretch_number = to_fix$stretch_number[i], 
        start = min_t, end = max_t, length = max_t - min_t + 1)
    }
  )
  
  # 3. We make sure each stretch is unique & that there isn't any overlap
  stretches_to_fix = 
    stretches_to_fix %>% 
    mutate(
      condition = !is.na(lag(end)) & (lag(end) >= start),
      start = ifelse(condition, lag(start), start),
    ) %>% 
    group_by(seq_id, start) %>% 
    summarize(
      end = max(end),
      stretch_number = min(stretch_number),
      .groups = "drop") %>% 
    mutate(
      length = end - start + 1,
      group = str_c(seq_id, "_",stretch_number)
      )
  
  # 4. We determine the tracking behavior in these stretches
  X_in_stretches = 
    X %>% 
    inner_join(
      ., 
      stretches_to_fix[rep(1:nrow(stretches_to_fix), stretches_to_fix$length),] %>% 
        group_by(stretch_number) %>% 
        mutate(t = start + row_number() - 1) %>% 
        select(seq_id, t, stretch_number, group),
      by = c("seq_id", "t")
    ) %>%
    inner_join(
      .,
      cycles %>%  select(seq_id, cycle_start, cycle_number) %>%  dplyr::rename(t = cycle_start),
      by = c("seq_id","t")
    )
  
  tracking_in_stretches_to_fix = 
    compute_tracking_behavior(X_in_stretches)
  stretches_to_fix = 
    stretches_to_fix %>% 
    left_join(
      .,
      tracking_in_stretches_to_fix %>% select(group, tracking_behavior), 
      by = "group") %>% 
    select(-group)
  
  GT = 
    bind_rows(
      stretches_to_fix %>%  mutate(state = 1, t = start),
      stretches_to_fix %>%  mutate(state = 1, t = end),
    ) %>% 
    select(seq_id, t, state) %>% 
    mutate(trust = 1)
  
  
  # 4. For each stretch, decode with decode_stretch() function
  decodings_in_stretches_to_fix =
    decode_stretches(X = X, GT = GT, 
                     model = model, is_hmm = is_hmm, fit_model = fit_model,
                     stretches = stretches_to_fix, 
                     verbose = verbose)
  
  
  # 5. Remove the time-points that are overlapping with the new stretches, append the new decodings (from step 4)
  decoding = decoding %>% filter(!(t %in% decodings_in_stretches_to_fix$t)) %>% 
    bind_rows(., decodings_in_stretches_to_fix)
  
  # 6. Sort by seq_id, t, state
  decoding = decoding %>% arrange(seq_id, t, state)
  
  # 7. Return results
  decoding
}




decode_without_adaptation = function(X, model, fit_model, verbose = FALSE){
  
  if(fit_model){
    if(verbose) cat("fitting model to observations\n")
    fitted_res = fit_hsmm(model = model, X = X, lock_transition = TRUE, N0_emission = 100)
    model = fitted_res$model
  }
  
  if(verbose) cat("predicting most likely sequence\n")
  vit = predict_states_hsmm(model = R_hsmm, X = X, method = "Viterbi")
  if(verbose) cat("predicting state probabilities at each time step\n")
  fwbw = predict_states_hsmm(model = R_hsmm, X = X, method = "FwBw")
  
  if(verbose) cat("preparing output\n")
  vit$state_seq %>% 
    select(seq_id, t, state) %>% 
    left_join(
      fwbw$probabilities %>%  
        select(seq_id, t, state, state_prob) %>% 
        rename(prob = state_prob),
      by = c("seq_id","t","state")
    ) %>% 
    mutate(tracking_behavior = "unspecified",
           stretch_id = 1)
  
}

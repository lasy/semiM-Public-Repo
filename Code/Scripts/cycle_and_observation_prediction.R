
identify_cycles = function(RES, X, ovu_state = 5, ovu_states = 1:7){
  RES = 
    RES %>% 
    arrange(seq_id, t) %>% 
    group_by(seq_id) %>% 
    mutate(cycle_start = (state == 1) & ((lag(state) != 1) | (t == 1)),
           cycle_end = (state != 1) & (lead(state) == 1),
           cycle_number = cumsum(cycle_start)) %>% 
    group_by(seq_id, cycle_number) %>% 
    mutate(cycle_length = n()) %>% 
    ungroup() %>% 
    select(-tracking_behavior, -stretch_number) %>% 
    left_join(X, by = c("seq_id", "t")) %>% 
    mutate(group = str_c(seq_id, "_", cycle_number))
  
  tracking_categories = compute_tracking_behavior(RES)
  
  RES = 
    RES %>% 
    left_join(tracking_categories, by = "group") %>% 
    select(-group, -all_of(colnames(X %>% select(-seq_id,-t))) )
  
  if("is_transition" %in% colnames(RES)){
    RES = RES %>% select(-is_transition, -transition_number)
  }
  
  cycles = 
    RES %>% 
    group_by(seq_id, cycle_number, cycle_length, tracking_behavior) %>% 
    summarize(
      cycle_start = min(t),
      cycle_end = max(t),
      min_prob = min(prob),
      max_prob = max(prob),
      ovu_prob = max(prob[state == ovu_state]),
      n_ovu = sum(state == ovu_state),
      ovulatory_cycle = all(state %in% ovu_states),
      .groups = "drop")
  
  list(RES = RES, cycles = cycles)
}


cycle_characteristics_per_user = function(cycles){
  cycles %>%
    group_by(seq_id) %>% 
    summarize(n_cycles = n(),
              median_cycle_length_all = median(cycle_length, na.rm = TRUE),
              median_cycle_length = median(cycle_length[cycle_length < 50], 
                                           na.rm = TRUE),
              sd_cycle_length_all = sd(cycle_length, na.rm = TRUE),
              sd_cycle_length = sd(cycle_length[cycle_length < 50], na.rm = TRUE),
              .groups = "drop")
}




predict_obs = function(
  input, fifth_cycle, model, file, reset = FALSE
){
  if(reset | !file.exists(file)){
    if(model != "baseline")
      preds = predict_obs_hsmm(input = input, fifth_cycle = fifth_cycle, model = model)
    else
      preds = predict_obs_baseline(input = input, fifth_cycle = fifth_cycle)
    write_feather(preds, path = file )
  }else{
    preds = read_feather(path = file)
  }
  preds
}



predict_obs_hsmm = function(input, fifth_cycle, model){
  map_dfr(
    .x = unique(fifth_cycle$seq_id),
    .f = function(sid){
      cat(sid,"\t")
      if(model == "fitted"){
        load(file = str_c(fitted_model_dir, sid,".Rdata"))
        hsmm_model = trained_model$model
      }else{
        hsmm_model = P_hsmm
      }
      
      pred = predict_next_time_point(input = input %>% filter(seq_id == sid),
                                     model = hsmm_model,
                                     fifth_cycle = fifth_cycle %>%  filter(seq_id == sid))
    }
  )
}




predict_next_time_point = function(input, model, fifth_cycle){
  pred_table = build_pred_table_from_past_cycles(input, model)
  pred = predict_next_time_point_for_5th_cycle(fifth_cycle, pred_table)
  pred
}


build_pred_table_from_past_cycles = function(input, model){
  
  res = predict_states_hsmm(model = model, X = input, method = "Viterbi")
  
  res = 
    res$state_seq %>% 
    mutate(
      is_first_day = 
        (state == 1) & 
        ((lag(state) != 1) | is.na(lag(state))),
      cycle_nb = cumsum(is_first_day),
      is_ovu = (state == 5),
      cycle_nb_o = cumsum(is_ovu)
    ) %>% 
    group_by(cycle_nb) %>% 
    mutate(
      cycle_length = n(),
      cycleday_fw = row_number(),
      cycleday_bw = row_number() - cycle_length - 1
    ) %>% 
    group_by(cycle_nb_o) %>% 
    mutate(
      cycle_length_o = n(),
      cycleday_fw_o = row_number() - 1,
      cycleday_bw_o = row_number() - cycle_length_o -1
    )  %>% 
    ungroup() %>% 
    mutate(
      phase = ifelse((cycle_nb == cycle_nb_o) & (cycleday_fw_o > 0), 
                     "luteal","follicular"),
      cycleday_m = ifelse(phase == "luteal", cycleday_bw, cycleday_fw),
      cycleday_o = ifelse(phase == "luteal" | is_ovu , cycleday_fw_o, cycleday_bw_o)
    ) %>% 
    filter(cycle_nb %in% 1:4)
  
  
  phase_length = 
    res %>% 
    group_by(cycle_nb, phase) %>% 
    summarize(n = n(), .groups = "drop") %>% 
    group_by(phase) %>% 
    summarize(mean_duration = mean(n) %>% round(), .groups = "drop")
  
  pred_table = 
    res %>% 
    left_join(., phase_length, by = "phase") %>% 
    mutate(
      cycleday_phase = 
        case_when(
          (phase == "follicular") & (cycleday_o >= -3) ~ mean_duration + cycleday_o,
          (phase == "follicular") & (cycleday_o < -3) & (cycleday_m <= (mean_duration - 4)) ~ cycleday_m,
          (phase == "luteal") & (cycleday_o <= 5) ~  cycleday_o,
          (phase == "luteal") & (cycleday_o > 5) & (cycleday_m >= -(mean_duration - 5)) ~ mean_duration + cycleday_m + 1,
          TRUE ~ NA_real_
        )) %>% 
    filter(!is.na(cycleday_phase))
  
  pred_table = 
    pred_table %>% 
    select(seq_id, t, phase, cycleday_phase) %>% 
    left_join(., input, by = c("seq_id", "t")) %>% 
    arrange(phase, cycleday_phase)
  
  pred_table
}




predict_next_time_point_for_5th_cycle = function(fifth_cycle, pred_table){
  
  fifth_cycle = fifth_cycle %>% mutate(cycleday = row_number())
  this_cycle_length = nrow(fifth_cycle)
  
  
  pred = 
    map_dfr(
      .x = 2:this_cycle_length,
      .f = function(cd){
        this_X = fifth_cycle %>% filter(t < min(fifth_cycle$t) + cd - 1)
        dec = predict_states_hsmm(
          model = trained_model$model, X = this_X, method = "Viterbi",
          ground_truth = data.frame(seq_id = sid, t = min(this_X$t), state = 1))$state_seq
        if(any(dec$state == 5)){
          cd_since_ovu = sum(dec$state > 5)
          this_pred = 
            pred_table %>% 
            filter(phase == "luteal") %>% 
            filter(cycleday_phase == min(cd_since_ovu+1, max(cycleday_phase)))
        }else{
          this_pred = 
            pred_table %>% 
            filter(phase == "follicular") %>% 
            filter(cycleday_phase == min(cd, max(cycleday_phase)))
        }
        this_pred %>% 
          select(seq_id, all_of(names(trained_model$model$marg_em_probs)), phase, cycleday_phase) %>% 
          mutate(t = max(this_X$t)+1,cycleday = cd)
      }
    )
  pred = 
    bind_rows(
      pred, 
      pred_table %>% 
        filter(phase == "follicular", cycleday_phase == 1) %>% 
        select(seq_id, all_of(names(trained_model$model$marg_em_probs)), phase, cycleday_phase) %>% 
        mutate(t = min(fifth_cycle$t),cycleday = 1)
    ) %>% 
    arrange(t)
  
}


predict_obs_baseline = function(input, fifth_cycle){
  
  input_with_cycles = 
    input %>% 
    left_join(., 
              cycles %>% rename(user_id = seq_id) %>% 
                select(user_id, cycle_start, cycle_number, cycle_length) %>% 
                mutate(t = cycle_start), 
              by = c("user_id","t")) %>% 
    mutate(cycle_start = !is.na(cycle_number),
           cycle_number = cycle_number %>% replace_NAs_with_latest_value(),
           cycle_length = cycle_length %>% replace_NAs_with_latest_value()
    ) %>% 
    group_by(user_id, cycle_number) %>% 
    mutate(cycleday_fw = row_number(),
           cycleday_bw = row_number() - cycle_length - 1,
           cycleday = ifelse(cycleday_bw >= -18, cycleday_bw , cycleday_fw)) %>% 
    ungroup()
  
  preds_baseline = 
    map_dfr(
      .x = unique(fifth_cycle$seq_id),
      .f = function(sid){
        cat(sid,"\t")
        
        # output
        this_fifth_cycle = fifth_cycle %>% filter(seq_id == sid)
        this_cycle_length = nrow(this_fifth_cycle)
        
        #input
        this_input = input_with_cycles %>% filter(seq_id == sid, t < min(this_fifth_cycle$t))
        expected_cycle_length = 
          this_input %>% filter(cycle_start) %>% 
          arrange(cycle_number) %>% 
          select(cycle_length) %>% unlist() %>% mean() %>% round()
        cycledays = c(-18:-1, 1:(expected_cycle_length-18))
        
        preds = 
          this_input %>% 
          filter(cycleday %in% cycledays) %>% 
          mutate(pred_cycleday = ifelse(cycleday < 0, cycleday + expected_cycle_length + 1, cycleday),
                 t = min(this_fifth_cycle$t) - 1 + pred_cycleday) %>% 
          select(seq_id, user_id, t, bleeding, mucus, temp, LH, preg)
        
        days_without_preds =  # that's if the cycle length of the fifth cycle is longer than the expected cycle length
          preds %>%  filter(is.na(seq_id))
        
        if(nrow(days_without_preds) > 0){
          j = which(preds$pred_cycleday == expected_cycle_length)
          preds = 
            bind_rows(
              preds %>% filter(!is.na(seq_id)), 
              preds[rep(j, nrow(days_without_preds)),] %>% 
                mutate(t = rep(days_without_preds$t, length(j)),
                       pred_cycleday = rep(days_without_preds$pred_cycleday, length(j)))
            )
        }
        
        preds
      }
    )
  preds_baseline
}




# simulate_next_time_points = function(X, model, n_time_points, ground_truth = data.frame(), n_samples = 1){
#   # first we decode X with viterbi
#   vit = predict_states_hsmm(model = model, X = X, method = "Viterbi", ground_truth = ground_truth)
#   
#   # then we identify the last state transition
#   last_state = .identify_last_state_transition(state_seq = vit$state_seq)
#   
#   # we modify the initial conditions so that they start on the last state
#   model$init = rep(0, model$J)
#   model$init[last_state$state] = 1
#   
#   pred = 
#     map_dfr(.x = 1:n_samples,
#             .f = function(i){
#               simulate_hsmm(model = model, seq_id = unique(X$seq_id), 
#                             n_timepoints = n_time_points + last_state$n_time_points_in_last_state) %>% 
#                 mutate(sim_nb = i,
#                        t = t + last_state$t - 1)
#             }
#     ) %>% 
#     filter(t > max(X$t))
#   
#   pred 
# }

# predict_next_time_points = function(X, model, n_time_points, ground_truth = data.frame()){
#   # first we decode X with viterbi
#   vit = predict_states_hsmm(model = model, X = X, method = "Viterbi", ground_truth = ground_truth)
#   
#   # then we identify the last state transition
#   last_state = .identify_last_state_transition(state_seq = vit$state_seq)
#   
#   # then we get the most likely state sequence
#   next_time_points_state_seq = 
#     .get_most_likely_state_seq(current_state = last_state$state, 
#                                model = model, 
#                                n_time_points = n_time_points + last_state$n_time_points_in_last_state)
#   
#   # then we get the most likely observation for each state
#   pred = .get_most_likely_obs(state_seq = next_time_points_state_seq, model = model)
#   
#   
#   # we return the dataframe
#   pred %>% 
#     mutate(seq_id = unique(X$seq_id),
#            t = last_state$t + (1:nrow(pred)) - 1) %>% 
#     filter(t > max(X$t))
# }
# 
# 
# .identify_last_state_transition = function(state_seq){
#   state_seq %>% 
#     mutate(state_change = (state != lag(state)) %>% replace_na(TRUE)) %>% 
#     filter(state_change) %>% 
#     slice_tail(n = 1) %>% 
#     mutate(n_time_points_in_last_state = max(state_seq$t) - t + 1)
# }
# 
# .get_most_likely_state_seq = function(current_state, model, n_time_points){
#   most_likely_state_transitions = .get_most_likely_state_transitions(current_state, model, n_time_points)
#   most_likely_sojourns = apply(.build_d_from_sojourn_dist(model = model, M = 10000), 2, which.max)[most_likely_state_transitions]
#   state_seq = rep(most_likely_state_transitions, times = most_likely_sojourns)[1:n_time_points]
#   state_seq
# }
# 
# .get_most_likely_state_transitions = function(current_state, model, n_time_points){
#   state_seq = c(current_state, rep(NA, n_time_points))
#   for(i in 1:n_time_points){
#     state_seq[i+1] = which.max(model$transition[state_seq[i],])
#   }
#   state_seq
# }

.get_most_likely_obs = function(state_seq, model){
  obs_probs = 
    model$b %>%  
    pivot_longer(., cols = starts_with("p_"),
                 names_to = "state", values_to = "p",
                 names_prefix = "p_") %>% 
    mutate(state = state %>%  as.integer())
  obs_probs = 
    obs_probs %>% 
    group_by(state) %>% 
    filter(p == max(p))
  
  data.frame(state = state_seq) %>% 
    left_join(., obs_probs, by = "state")
}

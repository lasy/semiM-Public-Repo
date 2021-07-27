




compute_accuracy = function(data_file, GT_file, decoding_file){
  
  X = read_feather(path = data_file)
  ML = read_feather(path = GT_file) %>% rename(state_GT = state) 
  RES = read_feather(path = decoding_file)
  
  
  XP = full_join(X, RES, by = c("seq_id", "t")) %>% 
    left_join(., ML, by = c("seq_id", "t") ) 
  
  
  Accuracy = mean(XP$state == XP$state_GT, na.rm = TRUE)
  
  
  Weighted_mean_of_sample_accuracy = weighted.mean(x = XP$state == XP$state_GT, w = XP$prob, na.rm = TRUE)
  
  Accuracy_per_state = 
    XP %>% 
    filter(!is.na(state_GT)) %>% 
    group_by(state_GT) %>% 
    summarize(accuracy = mean(state == state_GT, na.rm = TRUE), 
              weighted_mean_of_sample_accuracy =  weighted.mean(x = (state == state_GT), w = prob, na.rm = TRUE),
              .groups = "drop")
  
  confusion_matrix_df = 
    XP %>% 
    filter(!is.na(state), !is.na(state_GT)) %>% 
    group_by(state_GT, state) %>% 
    summarize(n = n(),
              wn = sum(prob),
              .groups = "drop") %>% 
    group_by(state_GT) %>% 
    mutate(tot = sum(n),
           wtot = sum(wn)) %>% 
    ungroup() %>% 
    mutate(perc = n/tot,
           wperc = wn/wtot) %>% 
    select(-n, -wn, -tot, -wtot) %>% 
    pivot_longer(cols = c("perc","wperc"), names_to = "proportion_type", values_to = "proportion") %>% 
    mutate(proportion_type = ifelse(proportion_type == "perc", "Proportion","Weighted proportion"),
           GT_state_name = R_hsmm$state_names[state_GT] %>%  factor(., levels = R_hsmm$state_names),
           decoded_state_name = R_hsmm$state_names[state] %>%  factor(., levels = R_hsmm$state_names))
  
  XP = 
    XP %>% 
    filter(!is.na(state_GT)) %>% 
    mutate(correct_prediction = (state_GT == state))
  
  
  list(
    Accuracy = Accuracy, 
    Weighted_mean_of_sample_accuracy = Weighted_mean_of_sample_accuracy, 
    Accuracy_per_state = Accuracy_per_state, 
    Conf_Matrix = confusion_matrix_df,
    XP = XP)
}





compute_pregnancy_duration = function(decoding_file){
  RES = read_feather(decoding_file)
  
  RES %>% 
    arrange(seq_id, t) %>% 
    filter(state %in% c(1, 12, 17), # we keep the menses, losses and births
           prob >= 0.7) %>% 
    group_by(seq_id) %>% 
    mutate(new_event = ((t-1) != lag(t)) %>% replace_na(TRUE)) %>% 
    filter(new_event) %>% 
    select(-new_event) %>% 
    mutate(duration = t - lag(t),
           type = 
             case_when(
               (state == 1) & (lag(state) == 1) ~ "cycle",
               (state == 1) & (lag(state) == 17) ~ "post-partum",
               (state == 1) & (lag(state) == 12) ~ "post-loss",
               (state == 17) ~ "pregnancy with birth",
               (state == 12) ~ "pregnancy with loss",
               TRUE ~ "undefined")
    ) %>% 
    select(seq_id, t, type, duration, tracking_behavior) %>% 
    ungroup()
  
}




create_X_for_decoding_viz =  function(data_file, GT_file, decoding_files){
  
  X = read_feather(path = data_file)
  ML = read_feather(path = GT_file) %>% rename(state_GT = state) 
  
  X = full_join(X, ML, by = c("seq_id", "t"))
  
  for(i in 1:nrow(decoding_files)){
    RES = read_feather(path = str_c(IO$output_data,"decodings/", decoding_files$file_name[i]))
    
    RES = RES %>% select(seq_id, t, state, prob) %>% 
      set_colnames(
        c("seq_id","t",
          str_c("state_",decoding_files$name[i]),
          str_c("state_prob_",decoding_files$name[i]))
      )
    
    X = X %>% full_join(., RES, by=  c("seq_id", "t"))
  }
  X
}

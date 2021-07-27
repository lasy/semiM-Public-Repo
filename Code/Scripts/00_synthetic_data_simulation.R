
realistic_censoring_probabilities <- function(model){
  
  p = model$censoring_probs$p
  p_ = p
  p_[c("M")] = 0.01
  p_[c("lE","Lut","Ano","AB","P","L","lEpL")] = 0.05
  p_[c("hE","preO","O","postO")] = 0.025
  p_[c("PL","PB1","PP")] = 0.9
  p_[c("PB2","PB3","BF")] = 0.98
  p_[c("B")] = 0.2
  
  q = model$censoring_probs$q
  q_ = q
  q_["LH", c("hE","preO","O","postO")] = c(0.3,0.1,0.1,0.1)
  q_["mucus", c("hE","preO","O","postO")] = 0.1
  q_["temp", c("hE","preO","O","postO","Ano")] = 0.1
  q_["temp", c("Lut","P")] = 0.3
  q_["temp", c("PL","PB1")] = 0.75
  
  list(p = p_, q = q_)
  
}




modulate_censoring_probs <-  function(rcp, alpha_ref, alphas, alpha_levels){
  
  missing_probs = 
    map_dfr(
      .x = 1:R_hsmm$J,
      .f = function(state){
        modulated_p = spline(x = alpha_ref, y = c(1,rcp$p[state],0), method = "hyman", xout = alphas)
        modulated_qs = 
          map_dfr(
            .x = R_hsmm$marg_em_probs %>% names(),
            .f = function(var){
              modulated_q = spline(x = alpha_ref, y = c(1,rcp$q[var,state],0), method = "hyman", xout = alphas)
              data.frame(variable = var,
                         alpha = alphas,
                         q = modulated_q$y)
            }
          )
        missing_prob_this_state = 
          left_join(
            data.frame(state = state,
                       alpha_level = alpha_levels,
                       alpha = alphas,
                       p =  modulated_p$y
            ),
            modulated_qs,
            by = c("alpha")
          )
      }
    )
  
  missing_probs = 
    missing_probs %>% 
    mutate(
      q = ifelse(variable == "bleeding", 0, q), # bleeding is never missing if the app is open
      missing_prob = p + (1-p)*q,
      alpha_level = alpha_level %>% factor(., levels = alpha_levels),
      state_name = R_hsmm$state_names[state] %>% factor(., levels = R_hsmm$state_names),
      state_color = R_hsmm$state_colors[state])
  
  
  missing_probs
}





simulate_individual_characteristics <- function(N_per_alpha, missing_probs){
  
  missing_probs = 
    missing_probs %>% 
    arrange(alpha)
  
  alphas = unique(missing_probs$alpha)
  alpha_levels = unique(missing_probs$alpha_level)
  
  N_users = length(alphas) * N_per_alpha
  
  s_users = data.frame(
    seq_id = 1:N_users, 
    alpha = rep(alphas, N_per_alpha),
    alpha_level = rep(alpha_levels,N_per_alpha),
    temp_sd = runif(N_users, min = 0.05, max = 0.3),
    mean_lE_sojourn = rnorm(N_users, mean = 10, sd = 3),
    sd_lE_sojourn = 0.2 + rpois(N_users, lambda = 2),
    mean_Lut_sojourn = rnorm(N_users, mean = 11, sd = 1),
    sd_Lut_sojourn = 0.2 + rpois(N_users, lambda = 0.75)
  )
  
  s_users
}



generate_time_series <- function(s_users, missing_probs, model){
  Xsim = purrr::map_dfr(
    .x = s_users$seq_id,
    .f = function(sid){
      #cat("\n",sid, "\n")
      modified_hsmm = model
      
      # modifying censoring probabilities
      this_user_alpha = s_users$alpha[s_users$seq_id == sid]
      modified_hsmm$censoring_probs$p = 
        missing_probs %>% filter(alpha == this_user_alpha) %>% 
        select(state, p) %>% arrange(state) %>% 
        distinct() %>% 
        select(p) %>% unlist() 
      modified_hsmm$censoring_probs$q = 
        missing_probs %>% filter(alpha == this_user_alpha) %>% 
        select(state, variable, q) %>% 
        pivot_wider(id_cols = variable, names_from = state, values_from = q) %>% 
        select(-variable) %>% as.matrix()
      
      # modifying emission distributions
      modified_hsmm$marg_em_probs$temp$params$sd = rep(s_users$temp_sd,R_hsmm$J)
      # modifying sojourns
      M = length(modified_hsmm$sojourn$M$d)
      modified_hsmm$sojourn$lE$d = dnorm(1:M, 
                                         mean = s_users$mean_lE_sojourn, 
                                         sd = s_users$sd_lE_sojourn)
      modified_hsmm$sojourn$lE$d = modified_hsmm$sojourn$lE$d/sum(modified_hsmm$sojourn$lE$d)
      modified_hsmm$sojourn$Lut$d = dnorm(1:M, 
                                          mean = s_users$mean_Lut_sojourn, 
                                          sd = s_users$sd_Lut_sojourn)
      modified_hsmm$sojourn$Lut$d = modified_hsmm$sojourn$Lut$d/sum(modified_hsmm$sojourn$Lut$d)
      
      # specifying modified model
      modified_hsmm = specify_hsmm(J = modified_hsmm$J,
                                   init = modified_hsmm$init, trans = modified_hsmm$transition,
                                   sojourn = modified_hsmm$sojourn, 
                                   marg_em_probs = modified_hsmm$marg_em_probs,
                                   censoring_probs = modified_hsmm$censoring_probs)
      
      # simulating sequence 
      Xsim_i = simulate_hsmm(model = modified_hsmm, seq_id = as.character(sid), n_state_transitions = 100)
    })
  Xsim %>% 
    left_join(s_users %>% mutate(seq_id = as.character(seq_id)), by = "seq_id")
}


modify_set_of_tracked_variables <-  function(Xsim){
  Xsim  = 
    Xsim %>% 
    group_by(seq_id) %>% 
    mutate(transitions = rbinom(n(), size = 1, prob = 0.5*(state == 1)*replace_na(lag(state)!= 1, 0)),
           tracking_sequence = 1+cumsum(transitions)) %>% 
    group_by(seq_id, tracking_sequence) %>% 
    mutate(tracking_category = sample(c("b","bp","btm","full"),1)) %>% 
    mutate(temp = ifelse(tracking_category %in% c("b","bp"), NA, temp),
           mucus = ifelse(tracking_category %in% c("b","bp"), NA, mucus),
           LH = ifelse(tracking_category %in% c("b","bp","btm"), NA, LH),
           preg = ifelse(tracking_category %in% c("b","btm"), NA, preg)) %>%
    ungroup() %>% 
    select(-transitions, -tracking_sequence)
  Xsim
}


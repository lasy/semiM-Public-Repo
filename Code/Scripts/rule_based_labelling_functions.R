



rule_based_labelling = function(days){
  
  users = unique(days$user_id)
  
  all_users_labels = foreach(user = users, .combine = bind_rows) %do%{
    this_user_day = days %>% dplyr::filter(user_id == user) %>% arrange(date)
    cycle_starts = find_cycle_starts(this_user_day = this_user_day, debug = FALSE)
    
    this_user_day = this_user_day %>% 
      dplyr::mutate(cycle_start = ifelse(rel_date %in% cycle_starts, TRUE, FALSE))
    
    periods = identify_period(this_user_day = this_user_day)
    if(any(this_user_day$preg_test == 1, na.rm = TRUE)){
      pregnancies = identify_pregnancies(this_user_day = this_user_day)
      pregnancy = pregnancies$pregnancy_starts
      post_partum = c()
      breast_feeding = c()
      pregnancies_with_losses = c()
      pregnancies_PB1 = c()
      pregnancies_PB2 = c()
      pregnancies_PB3 = c()
      births = c()
      
      ### End of pregnancies
      ## pregnancies with clear live births
      i_births = which(pregnancies$pregnancy_lengths %in% (30*7):(39*7))
      if(length(i_births)>0){
        # 1. we need to replace the period with a birth
        births = pregnancies$pregnancy_ends[i_births]
        j = which(periods %in% births)
        if(length(j)>0){
          for(jj in 1:length(j)){
            k = j[jj]; i = 0; d = periods[j[jj]+i]
            while((d+1) %in% periods){i = i+1; k = c(k,j[jj]+i); d = d+1}
          }
          births = periods[k]
          periods = periods[-k]
        }
        # 2. we need to classify the next cycle as PP or BF
        for(i in i_births){
          ics = which(cycle_starts == pregnancies$pregnancy_ends[i])
          if(ics == length(cycle_starts)){clength = max(this_user_day$rel_date) - pregnancies$pregnancy_ends[i]}else{clength = cycle_starts[ics+1]-cycle_starts[ics]-1}
          if(clength>0){
            pp = (pregnancies$pregnancy_ends[i]+1):(pregnancies$pregnancy_ends[i]+clength)
            if(clength > 8*7){breast_feeding = c(breast_feeding,pp)}else{post_partum = c(post_partum,pp)}
          }
        }
        # 3. we need to label the trimesters
        pregnancies_PB1 = c( pregnancies_PB1, sapply(i_births,function(i){(pregnancies$pregnancy_starts[i]+1):(pregnancies$pregnancy_starts[i]+12*7)}) %>% unlist())
        pregnancies_PB2 = c( pregnancies_PB2, sapply(i_births,function(i){(pregnancies$pregnancy_starts[i]+12*7+1):(pregnancies$pregnancy_starts[i]+24*7)}) %>% unlist())
        pregnancies_PB3 = c( pregnancies_PB3, sapply(i_births,function(i){(pregnancies$pregnancy_starts[i]+24*7+1):(pregnancies$pregnancy_end[i]-1)}) %>% unlist())
      }
      
      ## pregnancies with loss
      # 1. we need to replace the period with a loss
      # 2. we need to re-label the pregnancies as a "pregnancy with loss"
      i_losses = which(pregnancies$pregnancy_lengths %in% 0:(30*7-1))
      if(length(i_losses)>0){
        losses = pregnancies$pregnancy_ends[i_losses]
        j = which(periods %in% losses)
        if(length(j)>0){
          k = c()
          for(jj in 1:length(j)){
            k = c(k,j[jj]); i = 0; d = periods[j[jj]+i]
            while((d+1) %in% periods){i = i+1; k = c(k,j[jj]+i); d = d+1}
          }
          losses = periods[k]
          periods = periods[-k]
        }
        jl = which(losses > max(this_user_day$rel_date))
        if(length(jl)>0){losses = losses[-jl]}
        pregnancies_with_losses = sapply(i_losses,function(i){(pregnancies$pregnancy_starts[i]+1):(pregnancies$pregnancy_ends[i]-1)}) %>% unlist() %>% as.numeric()
        jl = which(pregnancies_with_losses > max(this_user_day$rel_date))
        if(length(jl)>0){pregnancies_with_losses = losses[-jl]}
      }
      
      
      ## pregnancies with unclear ends
      # 1. we need to create a birth at gestation week = 36
      # 2. and we need to define if pp or bf happened after
      i_u = which(pregnancies$pregnancy_lengths > (39*7))
      if(length(i_u)>0){
        undefined = pregnancies$pregnancy_starts[i_u]
        births = c(births, undefined + 36*7)
        for(u in undefined){
          j = pregnancies$pregnancy_starts == u
          plength = pregnancies$pregnancy_lengths[j]
          if(plength > 46*7){
            breast_feeding = c(breast_feeding, (u+36*7+1) :  (pregnancies$pregnancy_ends[j]-1))
          }else{
            post_partum = c(post_partum, (u+36*7+1) :  (pregnancies$pregnancy_ends[j]-1))
          }
        }
        # 3. we need to label the pregnancy trimesters
        pregnancies_PB1 = c( pregnancies_PB1, sapply(i_u,function(i){(pregnancies$pregnancy_starts[i]+1):(pregnancies$pregnancy_starts[i]+12*7)}) %>% unlist())
        pregnancies_PB2 = c( pregnancies_PB2, sapply(i_u,function(i){(pregnancies$pregnancy_starts[i]+12*7+1):(pregnancies$pregnancy_starts[i]+24*7)}) %>% unlist())
        pregnancies_PB3 = c( pregnancies_PB3, sapply(i_u,function(i){(pregnancies$pregnancy_starts[i]+24*7+1):(pregnancies$pregnancy_starts[i]+36*7-1)}) %>% unlist())
      }
      
    }else{
      pregnancy = c(); births = c(); losses = c(); pregnancies_with_losses = c();
      post_partum = c(); breast_feeding = c(); 
      pregnancies_PB1 = c();pregnancies_PB2 = c();pregnancies_PB3 = c()
    }
    
    
    
    labels = data.frame(user_id = user,
                        rel_date = 1:max(this_user_day$rel_date), 
                        label_id = as.character(now()),
                        state_name = NA,
                        stringsAsFactors = FALSE)
    
    ovulatory_cycle_pattern = c(rep("Early Foll",1000), rep("high E",2), "Ovu - 1","Ovu","Ovu + 2","Ovu + 2",rep("Luteal",11))
    for(i in 1:length(cycle_starts)){
      cs = cycle_starts[i]
      if(cs > 1){
        if(i == 1){ k = 1:(cs-1)}else{k = cycle_starts[i-1]:(cs-1)}
        labels$state_name[rev(k)] = rev(ovulatory_cycle_pattern)[1:length(k)]
      }
      if(i == length(cycle_starts)){
        k = cs:nrow(labels)
        labels$state_name[k] = "Early Foll"
      }
    }
    
    labels$state_name[periods] = "Menses"
    labels$state_name[post_partum] = "Post-Partum"
    labels$state_name[breast_feeding] = "Breast-Feeding"
    labels$state_name[births] = "Birth"
    labels$state_name[losses] = "Loss"
    labels$state_name[pregnancy] = "Pregnancy"
    labels$state_name[pregnancies_PB1] = "Pregnancy with Birth - 1st trimester"
    labels$state_name[pregnancies_PB2] = "Pregnancy with Birth - 2nd trimester"
    labels$state_name[pregnancies_PB3] = "Pregnancy with Birth - 3rd trimester"
    labels$state_name[pregnancies_with_losses] = "Pregnancy with Loss"
    
    
    m = match(labels$state_name, hsmm$states$names)
    labels$color = hsmm$states$colors[m]
    labels$type = "rule-based"
    labels$state_num = m
    
    #plot(labels$state_num, type = "l")
    
    return(labels)
  }
  return(all_users_labels)
}



find_cycle_starts = function(this_user_day = this_user_day, debug = FALSE){
  
  this_user_day$day = this_user_day$rel_date
  this_user_day$preg_test[is.na(this_user_day$preg_test)] = 0
  this_user_day$menstruation[is.na(this_user_day$menstruation)] = 0
  sum_bleeding_last_7_days = my_cumsum(t = this_user_day$day, x = this_user_day$menstruation, lags = -7:-1)
  sum_bleeding_2_days = my_cumsum(t = this_user_day$day, x = this_user_day$menstruation, lags = 0:1)
  
  
  day_last_bleeding =  -10
  day_last_cycle_start = -29
  median_cycle_length = 29
  is_pregnant = FALSE
  potential_pregnancy_end_type = "none"
  cycle_starts = c()
  pregnancies = c()
  for(i in 1:nrow(this_user_day)){ # nrow(this_user_day)
    today = this_user_day$day[i]
    if(debug){cat(today,"\n")}
    current_cycle_length = today  - day_last_cycle_start
    
    period_start = (this_user_day$menstruation[i] > 0) & 
      (sum_bleeding_2_days[i] > 1) & 
      (sum_bleeding_last_7_days[i] <= 1) &
      (current_cycle_length >= suppressWarnings(min(0.6*median_cycle_length,12)))
    
    user_defined_cycle_start = this_user_day$first_day[i] & (current_cycle_length >= 3*7) & (sum_bleeding_last_7_days[i] <= 2)
    
    # users get pregnant
    if(!is_pregnant & (this_user_day$preg_test[i] == 1)){
      if(debug){cat("\t got pregnant\n")}
      is_pregnant = TRUE
      day_first_pos_preg_test = today
      if(current_cycle_length < suppressWarnings(min(23,(median_cycle_length-6)))){ 
        # positive pregnancy test early in the cycle > probably got pregnant at previous cycle
        # we cancel the last cycle start
        if(length(cycle_starts)>0){
          cycle_starts = cycle_starts[-length(cycle_starts)]
          day_last_cycle_start = max(-29,cycle_starts[length(cycle_starts)])
        }
      }
    }
    
    # users not pregnant & period starts or user-defined new cycle
    if(!is_pregnant & (period_start | user_defined_cycle_start)){
      if(debug){cat("\t new cycle\n")}
      day_last_cycle_start = today
      cycle_starts = c(cycle_starts,day_last_cycle_start)
      pregnancies = c(pregnancies, 0)
    }
    
    # user is pregnant and all went well
    if(is_pregnant & (potential_pregnancy_end_type == "none")){
      
      time_since_first_positive_preg_test = today - day_first_pos_preg_test
      # pregnancy may end
      # either with a negative pregnancy test
      # or with a period
      # or with a user-defined cycle start
      potential_pregnancy_end_type = case_when(
        (this_user_day$first_day[i]) ~ "first day",
        (this_user_day$preg_test[i] == -1) ~ "neg preg test",
        period_start ~ "period",
        TRUE ~ potential_pregnancy_end_type
      ) 
      
      if(potential_pregnancy_end_type != "none"){
        day_at_potential_pregnancy_end = today
        if(debug){cat("\t potential pregnancy end: ",potential_pregnancy_end_type,"\n")}
      }
      
      # if current pregnancy has been going on for more than 2 cycle; pregnancy ends
      if((potential_pregnancy_end_type != "none")& 
         ((time_since_first_positive_preg_test >= (2*median_cycle_length)) |
          (current_cycle_length >= (3*median_cycle_length))
         )){
        day_last_cycle_start = today
        cycle_starts = c(cycle_starts,day_last_cycle_start)
        pregnancies = c(pregnancies, 1)
        is_pregnant = FALSE
        potential_pregnancy_end_type = "none"
        if(debug){cat("\t pregnancy end\n")}
      }
    }
    
    # user is pregnant and something happened
    if(is_pregnant && (potential_pregnancy_end_type != "none") && (day_at_potential_pregnancy_end != today)){
      # something happens again! => need to decide if to end pregnancy now or at previous event
      
      # new event is a period start and previous event happened within the last 3 weeks > end now
      # new event is a period start and previous event happened more than 3 weeks ago but less than 1.75*median cycle length ago > end then + end new cycle now
      if(period_start){
        if((today - day_at_potential_pregnancy_end ) %in% (3*7):round(1.75*median_cycle_length)){
          cycle_starts = c(cycle_starts,day_at_potential_pregnancy_end)
          pregnancies = c(pregnancies, 1,0)
        }else{pregnancies = c(pregnancies, 1)}
        day_last_cycle_start = today
        cycle_starts = c(cycle_starts,day_last_cycle_start)
        is_pregnant = FALSE
        if(debug){cat("\t period started; pregnancy end\n")}
        
      }
      
      # new event is not a period start 
      #   - but previous event was a period start > end then + continue new cycle
      #   - previous event was not a period start > continue "pregnancy"
      if((potential_pregnancy_end_type == "period")&
         ( (this_user_day$first_day[i]) | (this_user_day$preg_test[i] == -1) )){
        cycle_starts = c(cycle_starts,day_at_potential_pregnancy_end)
        pregnancies = c(pregnancies, 1)
        day_last_cycle_start = day_at_potential_pregnancy_end
        is_pregnant = FALSE
        potential_pregnancy_end_type = "none"
        if(debug){cat("\t previous period was actually a pregnancy end\n")}
      }
      
      # if a positive pregnancy test is logged, previous event is cancelled
      if(this_user_day$preg_test[i] == 1){ 
        potential_pregnancy_end_type = "none"
        if(debug){cat("\t pregnancy continues\n")}
      }
      
    }
    
    if(this_user_day$menstruation[i] > 0){
      day_last_bleeding =  today
      if(debug){cat("\t\t bleeding\n")}
    }
    
    if(sum(pregnancies == 0)>4){
      prev_median_cycle_length = median_cycle_length
      cycle_starts_no_preg = cycle_starts[pregnancies == 0]
      median_cycle_length = median(diff(cycle_starts_no_preg[(length(cycle_starts_no_preg)-3):length(cycle_starts_no_preg)]), na.rm = TRUE)
      if(is.na(median_cycle_length)){median_cycle_length = prev_median_cycle_length}
      if(debug & (prev_median_cycle_length != median_cycle_length)){cat("\t median cycle length : ",median_cycle_length,"\n")}
    }
    
    
    
  }
  #if(length(cycle_starts) == 0){ cycle_starts = 0}
  return(cycle_starts)
}



my_cumsum = function(t = d$days, x = d$menstruation, lags = -7:-1){
  if(any(is.na(t))){stop("NAs in t are not allowed\n")}
  tt = seq(suppressWarnings(min(min(t),min(t)+min(lags))), suppressWarnings(max(max(t)+ max(lags), max(t))), by = 1)
  x[is.na(x)] = 0
  xx = 0*tt
  m = match(t, tt)
  xx[m] = x
  yy = frollapply(xx, n = length(lags), FUN = sum, align = "left")
  y = yy[m+min(lags)]
  return(y)
}



identify_period = function(this_user_day = this_user_day){
  cycle_starts = this_user_day$rel_date[this_user_day$cycle_start]
  this_user_day$bleeding_mod = this_user_day$bleeding %>% replace_na(0)
  periods = c()
  for(cs in cycle_starts){
    #cat(cs,"\n")
    j = which(this_user_day$rel_date == cs)
    if(this_user_day$bleeding_mod[j] == 0){
      periods = c(periods, this_user_day$rel_date[j]+0:4) # if user did not report bleeding, period is 5 days
    }else{ # if user reported bleeding, we extend the period for as long as there is light bleeding
      periods = c(periods, this_user_day$rel_date[j])
      i = j+1; 
      while( 
        (i <= nrow(this_user_day)) &
        (this_user_day$bleeding_mod[i]>0.5) & 
        ((this_user_day$rel_date[i] - this_user_day$rel_date[i-1])<=2)){
        if((this_user_day$rel_date[i] - this_user_day$rel_date[i-1])==2){extra_day = c(this_user_day$rel_date[i]-1)}else{extra_day = c()}
        periods = c(periods, extra_day, this_user_day$rel_date[i]); i = i+1; 
      }
    }
  }
  k = which(periods > max(this_user_day$rel_date))
  if(length(k)>0){periods = periods[-k]}
  return(periods)
}


identify_pregnancies = function(this_user_day = this_user_day){
  pos_preg_test_days = this_user_day$rel_date[which(this_user_day$preg_test == 1)]
  cycle_starts = this_user_day$rel_date[this_user_day$cycle_start]
  pregnancies = c()
  pregnancy_starts = c()
  pregnancy_ends = c()
  last_pregnancy_start = pos_preg_test_days[1]
  for(i in 1:length(pos_preg_test_days)){
    preg_day = pos_preg_test_days[i]
    next_cycle_start = suppressWarnings(min(cycle_starts[cycle_starts>preg_day]))
    next_pos_preg_test = suppressWarnings(min(pos_preg_test_days[-(1:i)]))
    last_tracked_day = suppressWarnings(max(this_user_day$rel_date)+1)
    days_possible_next_events = c(next_cycle_start,next_pos_preg_test,last_tracked_day )
    j = which.min(days_possible_next_events)
    day_next_event = days_possible_next_events[j]
    pregnancies = c(pregnancies, preg_day:(day_next_event-1))
    if(j != 2){pregnancy_starts = c(pregnancy_starts,last_pregnancy_start); pregnancy_ends = c(pregnancy_ends, day_next_event)}
    last_pregnancy_start = pos_preg_test_days[i+1]
  }
  pregnancy_lengths = pregnancy_ends - pregnancy_starts
  return(list(pregnancies = pregnancies, pregnancy_starts = pregnancy_starts, pregnancy_ends = pregnancy_ends, pregnancy_lengths = pregnancy_lengths))
}


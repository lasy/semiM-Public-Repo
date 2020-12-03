

compute_dtemp = function(E){
  colnames_to_keep = colnames(E)

  E = E %>%
    group_by(user_id) %>%
    mutate(mean_temp_7d = frollmean(temp, n = 7, align = "right", na.rm = TRUE, hasNA = TRUE),
           mean_temp_7d_lagged = lag(mean_temp_7d, 13),
           dtemp_num = mean_temp_7d - mean_temp_7d_lagged,
           dtemp_char = case_when(
             is.na(dtemp_num) ~ "missing",
             dtemp_num >= 0.5 ~ "increase",
             dtemp_num <= -0.5 ~ "decrease",
             TRUE ~ "no_change"
           ),
           dtemp = dtemp_char %>% factor(.,levels = c("missing","increase","no_change","decrease"))
    ) %>%
    select(all_of(colnames_to_keep), dtemp) %>%
    ungroup()

  E
}



compute_bleeding_score = function(E){

  bleeding_dict = data.frame(code =  c("none","spotting" ,"light"  ,"medium"   ,"heavy"),
                             score = c(0,0.25,0.5,0.75,1))
  E = E  %>%
    mutate(bleeding_score = bleeding_dict$score[match(bleeding, bleeding_dict$code)] %>%  replace_na(0))
  E
}



compute_f.bleeding.last.week= function(E){

  colnames_to_keep = colnames(E)

  E = compute_bleeding_score(E)

  E = E  %>%
    group_by(user_id) %>%
    mutate(f.bleeding.last.week =
             frollmean(bleeding_score>0, n = 7, align = "right", fill = 0, na.rm = TRUE)) %>%
    ungroup() %>%
    select(all_of(colnames_to_keep), f.bleeding.last.week)

  E
}




compute_bleeding_density_5d = function(E){

  colnames_to_keep = colnames(E)

  E = compute_bleeding_score(E)

  E = E  %>%
    group_by(user_id) %>%
    mutate(bleeding_density_5d =
             frollmean(bleeding_score, n = 5, align = "center", fill = 0, na.rm = TRUE)) %>%
    ungroup() %>%
    select(all_of(colnames_to_keep), bleeding_density_5d)

  E
}




#####-----------



compute_dtemp_continuous = function(E){
  colnames_to_keep = colnames(E)

  E = E %>%
    group_by(user_id) %>%
    mutate(mean_temp_7d = frollmean(temp, n = 7, align = "right", na.rm = TRUE, hasNA = TRUE),
           sd_temp_7d = frollapply(x = temp, n = 7, fill = NA, align = "right", FUN = sd, na.rm = TRUE),
           mean_temp_7d_lagged = lag(mean_temp_7d, 6) %>% replace_na(-0.25),
           sd_temp_7d_lagged = lag(sd_temp_7d, 6) %>% replace_na(0.2) %>% pmax(.,0.1),
           dtemp = ((mean_temp_7d - mean_temp_7d_lagged)/sd_temp_7d_lagged) %>% asinh(),
           dtemp = dtemp %>% replace_na(0)
    ) %>%
    select(all_of(colnames_to_keep), dtemp) %>%
    ungroup()

  E
}






compute_interbleeding_id = function(E){
  colnames_to_keep = colnames(E)
  E = E %>%
    group_by(user_id) %>%
    # we define "interbleeding intervals" and give them an id so that we can aggregate by them.
    mutate(q_bleeding = quantile(bleeding_score[bleeding_score>0], p = 0.8),
           bv = ifelse((!is.na(q_bleeding)) && (q_bleeding > 0.5), 0.5, 0.25),
           bleeding_episode = (bleeding_score >= bv),
           start_of_interbleeding = (!bleeding_episode & lead(bleeding_episode)) %>% replace_na(0),
           interbleeding_id = cumsum(start_of_interbleeding)) %>%
    select(all_of(colnames_to_keep), interbleeding_id)

  E
}



compute_susp_ano = function(E){

  colnames_to_keep = colnames(E)

  E = E %>%
    group_by(user_id, interbleeding_id) %>%
    mutate(
      cd = row_number(),
      n_days = max(cd),
      cd_back = cd - n_days -1,
      interbleeding_id2 = ifelse( (bleeding_score > 0.25) | (cd <= 5) , NA, interbleeding_id)
    ) %>%
    group_by(user_id, interbleeding_id2) %>%
    # then, we look at the temperature and the bleeding patterns to detect hints of an anovulatory cycle.
    mutate(
      # an anovulatory cycle is a cycle in which the user is constantly bleeding lightly
      n_light_bleeding = sum((bleeding_score > 0) & (bleeding_score <= 0.5)),
      frac_light_bleeding = n_light_bleeding/n_days,
      fbleed = rescale(frac_light_bleeding, from = c(0.25,1), to = c(0,1)) %>% pmax(0), # we rescale to ignore "normal" amounts of light bleeding (normal spotting at the end of a cycle for example).
      #fbleed = pmax(fbleed, bleeding_density_40days), # why???
      # or a cycle in which the temperature is low, especially in the last 10 days.
      n_low_temp = sum((temp < 0) & (cd_back %in% -10:-1), na.rm = TRUE),
      n_high_temp = sum((temp > 0.3) & (cd_back %in% -10:-1), na.rm = TRUE),
      n_temp = pmax(0, n_low_temp -  n_high_temp), # count the number of low temperatures in the last 10 days (- the number of high temp)
      temperature = sigmoid(n_temp / 10,x0 = 0.4, slope = 8), # rescale with a threshold # the higher the better
      n_high_temp_over_whole_interval = sum(temp > 0.3, na.rm = TRUE),
      suspected_ano = (!is.na(interbleeding_id2)) * (n_days >= 3) * pmax(temperature, fbleed) *(n_high_temp_over_whole_interval <= 10)
    ) %>% ungroup() %>%
    mutate(susp_ano = suspected_ano %>% replace_na(0)) %>%
    ungroup() %>%
    select(all_of(colnames_to_keep), susp_ano)

  # plot(1:nrow(E),E$bleeding_score, type = "l", col = "red", xlim = c(3000,3220))
  # points(E$frac_light_bleeding, type = "l", col = "orange")
  # points(E$fbleed, type = "l", col = "red4")
  # points(E$temperature, type = "l", col = "blue")
  # points(E$susp_ano, type = "l", col = "black", lwd = 2)

  E
}



is_long_high_temp = function(user_switch, temp, user_temp_range, user_temp_threshold, gap_size ){

  # no switch of users
  same_user = sum(user_switch) == 0
  if(!same_user) return(0)

  # gaps no longer than 2 days
  gap_ok = max(gap_size, na.rm = TRUE) <= 2
  if(!gap_ok) return(0)

  # there are at least 4 distinct temperatures (it's suspect otherwise).
  distinct_temp = length(unique(temp[!is.na(temp)]))>=4
  if(!distinct_temp) return(0)

  # the 25 percentile of the temperature over the Ndays is higher than the user_temp_threshold = user_baseline_temp + 0.7*user_temp_range
  high_temp = quantile(temp, p = 0.25, na.rm = TRUE) >= user_temp_threshold[1]
  if(!high_temp) return(0)

  # the temperature remains high throughout the Ndays
  remains_high = abs(diff(quantile(temp,p = c(0.15,0.85), na.rm = TRUE))) <= (0.5 * user_temp_range[1])
  if(!remains_high) return(0)

  if(gap_ok & distinct_temp & high_temp & remains_high) return(1)

  0
}



compute_long_high_temp = function(E){
  colnames_to_keep = colnames(E)

  # This function returns 1 after Ndays of elevated temperatures, 0 otherwise.
  Ndays = 22


  ###### long_high_temp = 1 if
  # the gap between two temperature measurements is not longer than 2 days
  # there are at least 4 distinct temperatures (it's suspect otherwise).
  # the 25 percentile of the temperature over the Ndays is higher than the user_temp_threshold = user_baseline_temp + 0.7*user_temp_range
  # the temperature remains high throughout the Ndays (i.e. the  median temperature in the last days of the Ndays is as high as the median temperature at the start of these Ndays)

  E = E %>%
    group_by(user_id) %>%
    # defining user_baseline_temp and user_temp_range
    mutate(user_baseline_temp = quantile(temp, p = 0.2, na.rm = TRUE),
           user_temp_range = quantile(temp, p = 0.8, na.rm = TRUE) - user_baseline_temp,
           user_temp_threshold = user_baseline_temp + 0.7*user_temp_range) %>%
    # gaps
    mutate(is_gap = is.na(temp),
           gap_id = cumsum(c(0,diff(is_gap) == 1)),
           gap_id = ifelse(is_gap, gap_id, -1)) %>%
    group_by(user_id, gap_id) %>%
    mutate(gap_size = n() * is_gap) %>%
    ungroup() %>%
    # user_switch
    mutate(user_switch = (user_id != lead(user_id)) %>% replace_na(0))

  is_long_high_temp_rolling = rollify(is_long_high_temp, window = Ndays)
  E = E %>%
    mutate(long_high_temp = is_long_high_temp_rolling(user_switch, temp, user_temp_range, user_temp_threshold, gap_size)) %>%
    select(all_of(colnames_to_keep), long_high_temp)

  E
}



compute_long_high_temp_deprecated = function(E){
  # This function returns 1 after Ndays of elevated temperatures, 0 otherwise.

  colnames_to_keep = colnames(E)

  Ndays = 22

  ###### long_high_temp = 1 if
  # the gap between two temperature measurements is not longer than 2 days
  # there are at least 4 distinct temperatures (it's suspect otherwise).
  # the 25 percentile of the temperature over the Ndays is higher than the user_temp_threshold = user_baseline_temp + 0.7*user_temp_range
  # the temperature remains high throughout the Ndays (i.e. the  median temperature in the last days of the Ndays is as high as the median temperature at the start of these Ndays)

  E = E %>%
    group_by(user_id) %>%
    # defining user_baseline_temp and user_temp_range
    mutate(user_baseline_temp = quantile(temp, p = 0.2, na.rm = TRUE),
           user_temp_range = quantile(temp, p = 0.8, na.rm = TRUE) - user_baseline_temp,
           user_temp_threshold = user_baseline_temp + 0.7*user_temp_range) %>%
    # gaps
    mutate(is_gap = is.na(temp),
           gap_id = cumsum(c(0,diff(is_gap) == 1)),
           gap_id = ifelse(is_gap, gap_id, -1)) %>%
    group_by(user_id, gap_id) %>%
    mutate(gap_size = n()) %>%
    ungroup() %>% group_by(user_id) %>%
    mutate(x = gap_size*is_gap,
           no_large_gaps = frollapply(x = x, n = Ndays, align = "right", FUN = function(x) all(x<=2)),
           ok = no_large_gaps) %>%
    # at least 4 distinct temps
    mutate(at_least_4_distinct_temp = frollapply(x = temp, n = Ndays, align = "right", FUN = function(x) length(unique(x[!is.na(x)]))>=4),
           ok = ok & at_least_4_distinct_temp) %>%
    # the 25 percentile of the temperature over the Ndays is higher than the user_temp_threshold = user_baseline_temp + 0.7*user_temp_range
    mutate(temp_25pc = frollapply(x = temp, n = Ndays, align = "right", FUN = function(x) quantile(x, p = 0.25, na.rm = TRUE)),
           is_high_temp = temp_25pc >= user_temp_threshold,
           ok = ok & is_high_temp) %>%
    # the temperature stays high
    mutate(
      temp_range = ifelse(ok,
                          frollapply(x = temp, n = Ndays, align = "right",
                                     FUN = function(x) abs(diff(quantile(x,p = c(0.15,0.85), na.rm = TRUE)))),
                          NA),
      temp_remains_high = (temp_range <= user_temp_range*0.5) ) %>%
    # combining the conditions
    mutate(long_high_temp = no_large_gaps & at_least_4_distinct_temp & is_high_temp & temp_remains_high) %>%
    ungroup()


  E = E %>% select(all_of(colnames_to_keep), long_high_temp)

  E
}








compute_m_return = function(E){

  E_cycle_agg = E %>%
    group_by(user_id, interbleeding_id) %>%
    summarize(n = n())

  E_agg = E_cycle_agg %>%
    group_by(user_id) %>%
    summarise(median_cycle_length = median(n[n %in% 15:100]))

  E_cycle_agg = E_cycle_agg %>%
    mutate(median_cycle_length = E_agg$median_cycle_length[match(user_id, E_agg$user_id)],
           perc_diff_with_mcl = (n - median_cycle_length)/median_cycle_length,
           m_return = frollapply(x = perc_diff_with_mcl, n = 2, align = "left", FUN = function(x) sum(abs(x)) <= 0.5))

  E = E %>% full_join(.,E_cycle_agg %>% select(user_id, interbleeding_id, m_return), by = c("user_id","interbleeding_id"))

  E
}


my_filter = function(x, filter){
  x = c(rep(0,length(filter)-1), x)
  y = stats::filter(x, filter = filter, sides = 1)
  y = y[!is.na(y)]
  y
}




compute_preg_hint = function(E){
  columns_to_keep = colnames(E)

  # duration of pregnancy hint after a pregnancy hint triggering event
  ns = 80*7

  # pregnancy hint triggering events: either a long stretch of high temperatures OR a positive pregnancy test
  E = E %>%
    mutate(trigger = (preg | long_high_temp) %>%  replace_na(0)) %>%
    group_by(user_id) %>%
    mutate(trigger_clean = frollapply(trigger, n = 21, align = "right", FUN = function(x) last(x) * (sum(x)==1) ) %>% replace_na(0),
           preg_id = cumsum(trigger_clean)) %>%
    group_by(user_id, preg_id) %>%
    mutate(preg_hint = my_filter(x = trigger_clean, filter = rep(1, ns)) %>%  pmin(.,1))

  ##### adverse events
  # return of menstruations
  E = E %>%
    mutate(
      m_return_cs = cumsum(m_return %>% replace_na(0)) %>% pmin(.,1),
      preg_hint = pmax(0, preg_hint - m_return_cs))
  # a positive LH test after birth
  E = E %>% mutate(
    pd = row_number(),
    LH_after_birth = ((pd > 39*7) & (LH == 1)) %>% replace_na(0),
    LH_after_birth_cs = cumsum(LH_after_birth %>% replace_na(0)) %>% pmin(.,1),
    preg_hint = pmax(0, preg_hint - LH_after_birth_cs)
  )
  # any bleeding turns this score off momentarily except in the first 2 weeks (this helps detecting births or losses)
  E = E %>% mutate(
    preg_hint = ifelse((bleeding_score >= 0.5) & (pd <= 2*7), 0, preg_hint)
  )
  # adverse events such as period-like events or negative pregnancy tests decrease the intensity of the pregnancy hint score
  E = E %>% mutate(
    neg_preg_test = (preg == 0) %>% replace_na(0),
    period_like_bleeding = (bleeding_density_5d >= 0.5) %>% replace_na(0),
    period_like_bleeding = period_like_bleeding & !lag(period_like_bleeding),
    adverse_event = cumsum(pmax(0,neg_preg_test, period_like_bleeding, na.rm = TRUE)),
    preg_hint = preg_hint/(2^adverse_event)
  )

  # at the trigger, preg_hint has a higher value
  E = E %>% mutate(
    preg_hint = ifelse(trigger_clean, 1.2, preg_hint)
  )

  #
  E = E %>%  ungroup() %>% select(all_of(columns_to_keep), preg_hint)

  E
}



compute_tracking = function(E){

  E = E %>% mutate(
    tracking = case_when(bleeding != 0 ~ TRUE,
                         !is.na(mucus) ~ TRUE,
                         !is.na(temp) ~ TRUE,
                         !is.na(LH) ~ TRUE,
                         !is.na(preg) ~ TRUE,
                         TRUE ~ FALSE)
  )
  E
}


compute_gap_Xd = function(E, gap_size = 70){
  colnames_to_keep = colnames(E)

  Xd = gap_size # 70 days

  E = compute_tracking(E = E %>% ungroup())

  E = E %>%  group_by(user_id) %>%
    mutate(gap_start = !tracking & lag(tracking),
           gap_start = gap_start %>% replace_na(0),
           gap_id = cumsum(gap_start),
           gap_id = ifelse(tracking, NA, gap_id)) %>%
    group_by(user_id, gap_id) %>%
    mutate(gap_size = n(),
           gap_size = ifelse(tracking, 0, gap_size),
           gap_Xd = (gap_size >= Xd) %>% replace_na(0)) %>%
    ungroup() %>%
    select(all_of(colnames_to_keep), gap_Xd)

  E
}




#####------------------


my_acf = function(x, lag.min = 1, lag.max = 10){cc = acf(x, lag.max = lag.max, plot = FALSE); return(max(cc$acf[lag.min:lag.max]))}

compute_acf = function(E){

  E = E %>%
    group_by(user_id) %>%
    mutate(bleeding_range = quantile(bleeding_score[bleeding_score>0], p = 0.8),
           bleeding_lower_lim = ifelse(bleeding_range > 0.5, 0.5, 0.25),
           any_bleeding = bleeding_score >= bleeding_lower_lim,
           y = 1*any_bleeding,
           acf = frollapply(y, n = 100, FUN = my_acf, lag.min = 23, lag.max = 45,align = "center") %>% replace_na(0)) %>%
    ungroup()
  E
}



compute_preg_hint_v1 = function(input, tmp, get_col_names_only = FALSE){
  input_cols = c("preg_test","bleeding")
  tmp_cols = c("temp", "bleeding_density_5d", "acf")
  output_cols = c("preg_hint")

  if(get_col_names_only){return(list(input_cols = input_cols, tmp_cols = tmp_cols, output_cols = output_cols))}

  output = input %>%  dplyr::select(rel_date)
  output$preg_hint = 0

  # duration of pregnancy hint after a pregnancy hint triggering event
  ns = 50*7

  # pregnancy hint triggering events: either a long stretch of high temperatures OR a positive pregnancy test

  # Any long stretches of high temperature?
  temp_mod = tmp$temp; temp_mod[tmp$bleeding_score>0] = NA # we remove temperature measurements when bleeding is logged
  high_temp = (frollapply(temp_mod, n = 35, FUN = function(x){x2 = x; x2 = x2[!is.na(x)];sum(x2>0.2, na.rm = TRUE)}, align = "left") >= 20) & # a long stretch of high temperature is at least 20/35 days of relative temperature above 0.2
    (frollapply(temp_mod, n = 15, FUN = function(x){median(x, na.rm = TRUE)}, align = "left") >= 0.2) # and a median temperature over 0.2 for 15 days
  high_temp = high_temp %>% replace_na(.,0)
  j_temp = which(high_temp & c(0, diff(high_temp))); j_temp = j_temp+14
  if(any(diff(j_temp)<20)){j_temp = j_temp[-(which(diff(j_temp)<20)+1)]}
  preg_from_temp = 0*high_temp; if(length(j_temp)>0){k = rep(j_temp,each = ns) + (1:ns); k = k[k %in% 1:nrow(tmp)]; preg_from_temp[k] = 1}

  # Any pregnancy tests ?
  j_pregt = which(tmp$preg==1)
  preg_from_preg_test = 0*high_temp; if(length(j_pregt)>0){k = rep(j_pregt,each = ns) + (1:ns); k = k[k %in% 1:nrow(tmp)]; preg_from_preg_test[k] = 1}
  if(any(diff(j_pregt)<20)){j_pregt = j_pregt[-(which(diff(j_pregt)<20)+1)]}
  # combining triggering events
  j = sort(c(j_temp, j_pregt))
  if(any(diff(j)<20)){j = j[-(which(diff(j)<20)+1)]}

  # for each of the pregnancy trigger:
  for(jj in j){
    ix = jj:min((jj+ns-1),nrow(tmp))
    rx = 1:length(ix)
    preg_hint = 1+0*rx
    # the pregnancy hint variable is modulated by different factors
    # 1. the bleeding patterns autocorrelation:
    #    if the autocorrelation is high, the user is unlikely to be pregnant; except at the beginning of the pregnancy
    acf = tmp$acf[ix] * sigmoid(rx,x0 = 30,slope = 0.3)
    acf[is.na(acf)] = 0
    preg_hint = 1-pmax(0,acf)
    # 2. any bleeding turns this score off (this helps detecting births or losses)
    preg_hint = preg_hint * (tmp$bleeding_score[ix] == 0)

    # 3. adverse events such as period-like events or negative pregnancy tests decrease the intensity of the pregnancy hint score
    neg_preg_test = ((tmp$preg[ix] == -1)*1) %>%  replace_na(0)
    period_like_bleeding = (tmp$bleeding_density_5d[ix]>0.5)*1
    period_like_bleeding = ((period_like_bleeding == 1) & c(TRUE, diff(period_like_bleeding) == 1) )* 1 # we detect only the starts of period-like bleeding
    adverse_events = cumsum(pmax(0,neg_preg_test, period_like_bleeding, na.rm = TRUE))
    preg_hint = preg_hint / (2^adverse_events)

    # the pregnancy hint is triggers by any of the events:
    output$preg_hint[ix] = pmax(output$preg_hint[ix], preg_hint)
  }
  output$preg_hint = round(output$preg_hint/0.2)*0.2 # we make sure the values of the preg_hint matches the parameter emission definition
  return(output)
}

compute_M_hint = function(input, tmp, get_col_names_only = FALSE){
  input_cols = c("first_day")
  tmp_cols = c("mucus", "bleeding")
  output_cols = c("M_hint","mucus","LH","temp")

  if(get_col_names_only){return(list(input_cols = input_cols, tmp_cols = tmp_cols, output_cols = output_cols))}

  output = tmp %>% dplyr::select(rel_date, mucus, LH, temp)
  output$M_hint = 0


  # computing the range of cycle length we suspect new cycles starts based on the first_day defined by the app
  starts = which(input$first_day)
  lengths = diff(starts)
  Mcl = ceiling(median(lengths))
  if(diff(range(lengths))>10){R = Mcl+(-5:5)}else{R = min(lengths):max(lengths)}


  for(st in starts){
    #cat(st, "\n")
    for(dir in c("forward","backward")){
      if(dir == "forward"){r = st + R }else{r = st - R }
      r = intersect(r,1:nrow(tmp))
      j_any_b = which(tmp$bleeding_score[r]>0.25); j_any_b = r[j_any_b]
      if(length(j_any_b)>0){
        j_b_no_spot = r[which(tmp$bleeding_score[r]> 0.25)];
        if(length(j_b_no_spot) == 0){j = j_any_b}else{j = j_b_no_spot}
        if (length(j)==1){j = intersect(c(j, j+1), 1:nrow(tmp))}
        # we check there is no fertile mucus just before
        j_mucus_check = (min(j)-4):(min(j)-1); j_mucus_check = intersect(j_mucus_check,1:nrow(tmp))
        mucus_cond = (length(j_mucus_check)==0)|all(tmp$mucus[j_mucus_check]<0.5, na.rm = TRUE)
        # we also check that there is no start or as much bleeding in between this start and the next bleeding
        bleeding_sum = sum(tmp$bleeding_score[j], na.rm = TRUE)
        if(dir == "forward"){j_in_between = intersect((st+5):(min(j)-3),1:nrow(tmp))}else{j_in_between = intersect((max(j)+3):(st-5),1:nrow(tmp))}
        if(length(j_in_between)>0){
          bleeding_sum_in_between = sum(tmp$bleeding_score[j_in_between], na.rm =TRUE)
          inter_bleeding_cond = bleeding_sum_in_between < bleeding_sum
        }else{inter_bleeding_cond = TRUE}

        if(mucus_cond & inter_bleeding_cond){
          output$M_hint[c(j,st)] = 1
          # We will further impute mucus, LH and temp around M_hint
          before = intersect((min(j)-5):(min(j)-1), 1:nrow(output))
          after = intersect((max(j)+1):(max(j)+5), 1:nrow(output))
          around = c(before, j, after)
          output$mucus[(output$rel_date %in% around) & (is.na(output$mucus))] = 0
          output$LH[(output$rel_date %in% around) & (is.na(output$LH))] = 0
          output$temp[(output$rel_date %in% j) & (is.na(output$temp))] = 0
          output$temp[(output$rel_date %in% before) & (is.na(output$temp))] = 0.5
          output$temp[(output$rel_date %in% after) & (is.na(output$temp))] = -0.5
        }
      }
    }
  }
  return(output)
}

compute_fmucus = function(input, tmp, get_col_names_only = FALSE){
  input_cols = c()
  tmp_cols = c("mucus")
  output_cols = c("fmucus")

  if(get_col_names_only){return(list(input_cols = input_cols, tmp_cols = tmp_cols, output_cols = output_cols))}

  output = input %>% dplyr::select(rel_date);
  output$fmucus = frollapply(tmp$mucus, n = 3, align = "center", FUN = function(x) any(x>=0.5, na.rm = TRUE)) %>% replace_na(0)
  return(output)
}









####------


compute_LH = function(input, tmp, get_col_names_only = FALSE){
  input_cols = c("LH")
  tmp_cols = c()
  output_cols = c("LH")

  if(get_col_names_only){return(list(input_cols = input_cols, tmp_cols = tmp_cols, output_cols = output_cols))}

  output = input %>% dplyr::select(rel_date, LH) %>%
    rename(LH_test = LH) %>%
    dplyr::mutate(LH =
                    case_when(
                      is.na(LH_test) ~ NA_real_,
                      (LH_test == 0) ~ NA_real_,
                      (LH_test == -1) ~ 0,
                      (LH_test == 1) ~ 1,
                      TRUE ~ NA_real_)
    ) %>%
    dplyr::select(rel_date, LH)

  return(output)
}

compute_mucus = function(input, tmp, get_col_names_only = FALSE){
  input_cols = c("mucus_type", "bleeding")
  tmp_cols = c()
  output_cols = c("mucus")

  if(get_col_names_only){return(list(input_cols = input_cols, tmp_cols = tmp_cols, output_cols = output_cols))}

  output = input %>% dplyr::select(rel_date)
  output$mucus =  mucus.dict$score[match(input$mucus_type, mucus.dict$names)]
  output$mucus[input$bleeding_score >= 1] = 0
  return(output)
}


compute_temp = function(input, tmp, get_col_names_only = FALSE){
  input_cols = c("temperature","questionable_temp", "bleeding")
  tmp_cols = c()
  output_cols = c("temp")

  if(get_col_names_only){return(list(input_cols = input_cols, tmp_cols = tmp_cols, output_cols = output_cols))}

  output = input %>% dplyr::select(rel_date)
  output$temp = input$temperature
  j = which(input$questionable_temp)
  if(length(j)>0) output$temp[j] = NA  # questionable temperatures are set as missing data points

  # detect odd values (e.g. way too low, way too high, or oddly repeated values)
  if(any(!is.na(output$temp))){ # if there is at least one temperature
    # looking for oddly repeated values
    ttemp = table(output$temp) # histogram of temperatures
    if((length(unique(output$temp[!is.na(output$temp)]))==1)|(max(ttemp)> 5*(sort(ttemp,decreasing = TRUE)[2]))){weird_temp = as.numeric(names(ttemp)[which.max(ttemp)])}else{weird_temp = -9999}
    # and removing too high or too low
    odd =  which((output$temp < 95) | (output$temp > 101) | (output$temp == weird_temp))
    if(length(odd)>0) output$temp[odd] = NA
  }

  # Scale temperature
  median_temp = median(output$temp, na.rm = TRUE)
  output$temp = output$temp - median_temp
  # removing extreme values
  output$temp = output$temp %>% pmin(.,1.5) %>% pmax(., -1.5)
  # removing temperatures when bleeding >= heavy
  output$temp[input$bleeding == 3] = NA

  return(output)
}

compute_preg = function(input, tmp, get_col_names_only = FALSE){
  input_cols = c("preg_test")
  tmp_cols = c()
  output_cols = c("preg")

  if(get_col_names_only){return(list(input_cols = input_cols, tmp_cols = tmp_cols, output_cols = output_cols))}

  output = input %>% dplyr::select(rel_date, preg_test) %>%
    dplyr::mutate(preg =
                    case_when(
                      is.na(preg_test) ~ NA_real_,
                      (preg_test == 0) ~ NA_real_,
                      (preg_test == -1) ~ 0,
                      (preg_test == 1) ~ 1,
                      TRUE ~ NA_real_)
    ) %>%
    dplyr::select(rel_date, preg)

  return(output)
}




## ---- 00_function_prepare_obs.R

prepare_obs = function(d, observations = c("bleeding","mucus","temp","LH","preg")){ 
  
  user_id = unique(d$user_id)
  if(length(user_id)>1){stop("several users in the input\n")}
  
  # we expand the feature matrix d so that is has one row per day
  X = data.frame(user_id = user_id,rel_date = min(d$rel_date):max(d$rel_date))
  m = match(X$rel_date, d$rel_date)
  d$bleeding = d$bleeding %>% replace_na(0)
  
  if("bleeding" %in% observations){
    # we keep the 0, 0.5, 1, 2, 3 code
    X = X %>% mutate(bleeding = d$bleeding[m],
                     first_day = d$first_day[m] %>% replace_na(FALSE))
    # and replace all missing values by 0 (no bleeding)
    # X = X %>% mutate(bleeding = bleeding %>% replace_na(0))
    # and we also convert a "first day" into bleeding
    X  = X %>%
      mutate(first_day_bleeding = ifelse(first_day & ((bleeding < 1)|is.na(bleeding)),2,0),
             bleeding = pmax(bleeding, first_day_bleeding)) %>%
      select(-first_day, -first_day_bleeding)
    X = X %>% rename(bleeding_num = bleeding)
    X = X %>% mutate(bleeding = 
                       c("none","spotting","light","medium","heavy")[match(bleeding_num,c(0,0.5,1:3))]) %>%
      select(-bleeding_num)
  }
  
  if("mucus" %in% observations)
    X = X %>%
    mutate(mucus = mucus.dict$category[match(d$mucus_type[m], mucus.dict$names)] %>% 
             factor(., levels = levels(mucus.dict$category)))
  
  
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
  
  
  test_levels = c("pos","neg")
  
  if("LH" %in% observations)
    X = X %>% 
    mutate(LH = 
             case_when(d$LH[m] == -1 ~ "neg", d$LH[m] == 1 ~ "pos", TRUE ~ NA_character_) %>% 
             factor(.,levels = test_levels))
  
  if("preg" %in% observations)
    X = X %>% 
    mutate(preg = 
             case_when( d$preg_test[m] == -1 ~ "neg", d$preg_test[m] == 1 ~ "pos", TRUE ~ NA_character_) %>% 
             factor(.,levels = test_levels))
  
  X
}



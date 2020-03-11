


plot_user_history= function(d = data.frame(),
                            obs = data.frame(),
                            states_labels_list = list(), 
                            current_selection = data.frame(), 
                            show_all_variables = FALSE){
  
  if((lu(d$user_id)>1)|(lu(obs$user_id)>1)){stop("there is more than one user in the data frame\n")}
  user_id = unique(c(d$user_id, obs$user_id))
  if((nrow(d) + nrow(obs)) == 0){stop("both 'd' and 'obs' are empty \n")}
  
  # AXES
  start = min(d$rel_date, obs$rel_date)-0.55
  end = max(d$rel_date, obs$rel_date)+0.55
  n_days =  end - start
  scales = c("Y","Q","M")
  if(n_days > 2*365) i = 1 else if(n_days>365) i = 2 else i = 3
  scale = scales[i]
  by_break = c(365,365/4,365/12)[i]
  by_break_min = c(365/4, 365/12,365/12/4)[i]
  x_axis_breaks = seq(start,end, by = by_break)
  x_axis_minor_breaks = seq(start,end, by = by_break_min)
  x_axis_labels = paste0(scale, x_axis_breaks %/% by_break)
  
  df_base = data.frame(rel_date = c(start, end))
  
  # X axis
  g_axis = ggplot(df_base)+
    scale_x_continuous(breaks = x_axis_breaks, minor_breaks = x_axis_minor_breaks,
                       labels = x_axis_labels, limits = c(start, end))+
    xlab("Relative date")
  
  # base plot
  g = g_axis+
    scale_fill_identity()+
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y =  element_text(angle = 0))+
    guides(col = FALSE, fill = FALSE, shape = FALSE, size = FALSE, alpha = FALSE)
  
  
  # current selection (for interactive shiny app)
  
  if(nrow(current_selection)>0){
    geom_select = geom_rect(data = current_selection, aes(xmin = start, xmax = end,ymin = -Inf, ymax = Inf, fill = color), col = NA, alpha = 0.5)
  }else{
    geom_select = geom_blank()
  }
  
  # initialising the plotlist
  plotlist = list()
  rel_heights = c()
  
  
  # title = user_id
  g_title = 
    g+geom_text(label = user_id, x = start, y = 0, hjust = 0, vjust = 0)+
    theme(panel.grid = element_blank())
  plotlist$title = g_title
  rel_heights = c(rel_heights, 0.15)
  
  
  if(nrow(d)>0){
    
    # data transformation
    d$temperature[which(d$questionable_temp)] = NA
    d$med_temp = median(d$temperature, na.rm = TRUE)
    d$rel_temp = d$temperature - d$med_temp
    
    d$LH_cols = c("red","transparent","steelblue")[d$opk+2]
    d$PT_cols = c("red","transparent","steelblue")[d$preg_test+2]
    d$mucus_y = as.numeric(d$mucus_category) + (pmax(1,d$mucus_amount)-1)/5
    d$sex_y = as.numeric(d$sex_num)
    d$sex_y[is.na(d$sex_y)] = 0
    
    
    # Temperature
    if(show_all_variables | any(!is.na(d$temperature))){
      g_temp = g + 
        geom_select+
        geom_line(data = d, aes(x = rel_date, y = temperature), col = "gray90")+
        geom_point(data = d, aes(x = rel_date, y = temperature, col = rel_temp))+
        scale_color_gradient2(low = "cornflowerblue",mid = "gray90",high = "red", limits = c(-1.5,1.5))+
        ylim(c(96,100))+
        ylab("Temperature")
      plotlist$temp = g_temp
      rel_heights = c(rel_heights,0.4)
    }
    
    
    # Mucus
    if(show_all_variables | any(!is.na(d$mucus_y))){
      g_mucus = g + 
        geom_select+
        geom_line(data = d, aes(x = rel_date, y = mucus_y), col = "gray90")+
        geom_point(data = d, aes(x = rel_date, y = mucus_y, col = mucus_category))+ 
        scale_y_continuous(breaks = c(1:5), minor_breaks = NULL, labels = levels(d$mucus_category))+
        scale_color_manual(values= c("gray","lightgoldenrod2","lightyellow3","lightskyblue","skyblue3"))+
        ylab("Mucus")
      plotlist$mucus =  g_mucus
      rel_heights = c(rel_heights,0.3)
    }
    
    
    # Bleeding
    if(show_all_variables | any(!is.na(d$bleeding))){
      g_bleeding = g+
        geom_select+
        geom_vline(data = d[which(d$first_day),], aes(xintercept = rel_date), linetype = 1, col = "gray")+
        geom_segment(data = d,aes(x = rel_date, xend = rel_date, y = 0, yend = bleeding, color = bleeding))+
        scale_color_gradient(low = "orange1", high = "orangered")+
        scale_y_continuous(labels = NULL, minor_breaks = NULL)+
        ylab("Bleeding")
      plotlist$bleeding =  g_bleeding
      rel_heights = c(rel_heights,0.25)
    }
    
    # LH tests
    if(show_all_variables | any((d$opk !=0)&(!is.na(d$opk)))){
      g_LH = g+
        geom_select+
        geom_point(data = d,aes(x = rel_date, y = opk, col = LH_cols))+
        scale_color_identity()+
        scale_y_continuous(breaks = c(-1,1), minor_breaks = NULL, labels = c("neg","pos"), limits = c(-1.5,1.5))+
        ylab("LH test")
      plotlist$LH = g_LH
      rel_heights = c(rel_heights,0.15)
    }
    
    
    # Preg tests
    if(show_all_variables | any((d$preg_test !=0)&(!is.na(d$preg_test)))){
      g_PT = g+
        geom_select+
        geom_point(data = d,aes(x = rel_date, y = preg_test, col = PT_cols))+
        scale_color_identity()+
        scale_y_continuous(breaks = c(-1,1), minor_breaks = NULL, labels = c("neg","pos"), limits = c(-1.5,1.5))+
        ylab("Preg test")
      plotlist$preg =  g_PT
      rel_heights = c(rel_heights,0.15)
    }
    
    # Sex
    if(show_all_variables | any(!is.na(d$sex))){
      g_sex = g+
        geom_select+
        geom_point(data = d,aes(x = rel_date, y = sex_y, col = factor(sex_y), shape = factor(sex_y)))+
        scale_y_continuous(breaks = 1:4, minor_breaks = NULL, labels = levels(d$sex_short),limits = c(0.8,4.2))+
        ylab("Sex")
      plotlist$sex = g_sex
      rel_heights = c(rel_heights,0.2)
    }
  }
  
  if(nrow(obs)>0){
    obs_score = obs %>%  dplyr::select(rel_date, ends_with("_score"))
    for(c in 2:ncol(obs_score)){
      col_name = colnames(obs_score)[c]
      df = data.frame(rel_date = obs_score$rel_date, y = obs_score[,c])
      df_range = range(c(0,1,df$y), na.rm = TRUE)
      df_break = floor(min(df_range)):1
      g_score = g + 
        geom_line(data = df, aes(x = rel_date, y = y), col = "gray")+
        geom_point(data = df, aes(x = rel_date, y = y, col = y), size = 0.5)+
        scale_y_continuous(breaks = df_break, limits = df_range)+
        scale_color_gradient2(low = "cornflowerblue", mid = "gray", high = "red", midpoint = 0)+
        ylab(col_name)
      
      plotlist[[col_name]] = g_score
      rel_heights = c(rel_heights, 0.2)
    }
  }
  
  # States
  if(length(states_labels_list) > 0){
    for(sl in names(states_labels_list)){
      states_labels = states_labels_list[[sl]]
      if(nrow(states_labels)>0){
        if(nrow(states_labels) == lu(states_labels$rel_date)){
          if("cs" %in% colnames(states_labels)){
            g_state_cs =  g+
              geom_select+
              geom_line(data = states_labels, aes(x = rel_date, y = cs))+
              scale_y_continuous(breaks = c(0,1), minor_breaks = NULL, limits = c(0,1))+
              ylab(paste0("CS (",gsub("_"," ",sl),")"))
            plotlist[[str_c("CS_",sl)]] = g_state_cs
            rel_heights = c(rel_heights,0.18)
          }
          if(all(c("cs_cycle_min","cs_cycle_mean") %in% colnames(states_labels))){
            g_state_LR_agg =  g+
              geom_select+
              geom_line(data = states_labels, aes(x = rel_date, y = cs_cycle_min), col = "green3")+
              geom_line(data = states_labels, aes(x = rel_date, y = cs_cycle_mean), col = "black")+
              scale_y_continuous(breaks = c(0,1), minor_breaks = NULL, limits = c(0,1))+
              ylab(paste0("CS (",gsub("_"," ",sl),")"))
            plotlist[[str_c("CS_agg_",sl)]] = g_state_LR_agg
            rel_heights = c(rel_heights,0.18)
          }
          if("color" %in% colnames(states_labels)){
            if(!("alpha" %in% colnames(states_labels))){states_labels$alpha = 1}
            if(n_days <= 200){border = "white"}else{border = NA}
            g_state_dec =  g+
              geom_select+
              geom_rect(data = states_labels, aes(xmin = rel_date-0.5,xmax = rel_date+0.5, fill = color, alpha = alpha),ymin = 0, ymax = 1, col = border)+
              scale_y_continuous(breaks = NULL, minor_breaks = NULL, limits = c(0,1))+
              scale_alpha_continuous(limits = c(0,1))+
              ylab(paste0("States\n(",gsub("_"," ",sl),")"))
            rel_heights = c(rel_heights,0.18)
            plotlist[[sl]] =  g_state_dec
          }
        }else{
          g_state_dec =  g+
            geom_select+
            geom_line(data = states_labels, aes(x = rel_date, y = prob, col = color))+
            scale_color_identity()+
            scale_y_continuous(breaks = c(0,1), minor_breaks = c(1/2), limits = c(0,1)+0.04*c(-1,1))+
            ylab(paste0("States\n(",gsub("_"," ",sl),")"))
          rel_heights = c(rel_heights,0.2)
          plotlist[[sl]] =  g_state_dec
        }
      }
    }
  }
  
  plotlist$axis =  g_axis
  rel_heights = c(rel_heights,0.15)
  g_history = plot_grid(plotlist = plotlist, ncol=1, align="v", rel_heights = rel_heights)
  
  return(g_history)
  
}





ggplot_confusion_matrix = function(true , decoded){
  
  cm = table(true = true, decoded =  decoded)
  cm
  mx = match(hsmm$states$abbr,colnames(cm)); #if(any(is.na(mx))){mx[is.na(mx)]= (max(mx,na.rm = TRUE)+1) : (length(mx))}
  my = match(hsmm$states$abbr,rownames(cm))
  cm = cm[my, mx]
  
  cm_sum = rowSums(cm, na.rm = TRUE)
  (cm/cm_sum*100) %>%  round(.,1) 
  cm = cm/cm_sum *100
  
  cm_long = data.frame(cm)
  cm_long$true = factor(cm_long$true, levels = rev(hsmm$states$abbr))
  
  confusion_matrix_viz = ggplot(cm_long, aes(x = decoded, y = true , fill = Freq))+ coord_fixed()+
    geom_tile()+
    geom_abline(slope = -1, intercept = hsmm$n_states+1, col = "black", size = 0.2, linetype = 2)+
    geom_hline(yintercept = 0:hsmm$n_states+0.5, col = "black", size = 0.1)+
    geom_vline(xintercept = c(0.5,hsmm$n_states+0.5), col = "black", size = 0.1)+
    scale_fill_gradient(name = "% of manual labels",low = "white", high = "royalblue3")+
    ylab("Manual\nlabels")+
    xlab("Decoded labels")+
    theme(legend.position = "bottom",
          legend.title = element_text(vjust = 0.75),
          axis.title.y = element_text(angle = 0, hjust = 1, vjust = 1),
          axis.title.x = element_text(hjust = 1),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  return(confusion_matrix_viz)
}



ggplot_event_accuracy = function(event_accuracy){
  event_accuracy$no_event_in_30d_window = is.na(event_accuracy$decoded_rel_date)
  
  g = ggplot(event_accuracy, aes(x = time_diff, fill = no_event_in_30d_window)) +
    geom_histogram(binwidth = 1)+
    geom_vline(xintercept = 0, linetype = 2)+
    facet_grid(event ~ ., scale = "free")+
    scale_fill_manual(values = c("turquoise","red"))+
    xlab("difference in days (decoded - manual label)")+
    xlim(c(-17.5,15.5))+
    guides(fill = FALSE)
  return(g)
}




viz_em_parm = function(p1, p2, hsmm){
  
  
  N = 1000
  
  ## 1
  r_obs = foreach(state = 1:hsmm$n_states, .combine = rbind) %do%{
    this_state_r_obs = generate_random_obs(n = N, parem = p1, state = state) %>% 
      as.data.frame() %>% set_colnames(hsmm$obs_names) %>% 
      mutate(state_name =  hsmm$states$abbr[state], color = hsmm$state$colors[state],n = 1:N) %>% 
      tidyr::pivot_longer(.,cols = hsmm$obs_names, names_to = "obs_name")
    return(this_state_r_obs)
  }
  
  r_obs$state_name = factor(r_obs$state_name, levels = hsmm$states$abbr)
  r_obs$obs_name = factor(r_obs$obs_name, levels = hsmm$obs_names)
  
  r_obs1 = r_obs
  
  ## 2
  r_obs = foreach(state = 1:hsmm$n_states, .combine = rbind) %do%{
    this_state_r_obs = generate_random_obs(n = N, parem = p2, state = state) %>% 
      as.data.frame() %>% set_colnames(hsmm$obs_names) %>% 
      mutate(state_name =  hsmm$states$abbr[state], color = hsmm$state$colors[state],n = 1:N) %>% 
      tidyr::pivot_longer(.,cols = hsmm$obs_names, names_to = "obs_name")
    return(this_state_r_obs)
  }
  
  r_obs$state_name = factor(r_obs$state_name, levels = hsmm$states$abbr)
  r_obs$obs_name = factor(r_obs$obs_name, levels = hsmm$obs_names)
  
  r_obs2 = r_obs
  
  ## combining them
  
  r_obs1$type = "1"
  r_obs2$type = "2"
  r_obs = rbind(r_obs1, r_obs2)
  
  g_par = ggplot(r_obs, aes(x = type, y = value, fill = color, col = color))+  
    geom_hline(yintercept = c(0,1),col = "gray80")+
    geom_violin(bw = 0.05, size = 0.5)+
    facet_grid( obs_name ~ state_name, scale = "free")+
    scale_fill_identity()+ scale_color_identity()+
    coord_cartesian(ylim = c(-1.5,1.5))+
    ylab("emission parameter value")
  
  return(g_par)
}




viz_em_parm_old_style = function(p1, p2, hsmm){
  
  emis1 = expand.grid(variable = names(p1$mu[[1]]),state_num = 1:hsmm$n_states, model = "initial")
  emis1$mu = unlist(p1$mu)
  emis1$sigma = lapply(1:nrow(hsmm$states),function(i){diag(p1$sigma[[i]])}) %>%  unlist()
  
  emis2 = expand.grid(variable = names(p1$mu[[1]]), state_num =  1:hsmm$n_states, model = "fitted")
  emis2$mu = unlist(p2$mu)
  emis2$sigma = lapply(1:nrow(hsmm$states),function(i){diag(p2$sigma[[i]])}) %>%  unlist()
  
  emis = rbind(emis1, emis2)
  emis$state = hsmm$states$abbr[emis$state_num]
  emis$color = hsmm$states$colors[emis$state_num]
  emis$x = emis$state_num + 0.25*(emis$model == "fitted")
  
  g = ggplot(emis, aes(x = x, col = color))+
    geom_hline(yintercept = c(0,1),col = "gray80")+
    geom_point(aes(y = mu, shape = model))+
    geom_segment(aes(xend = x, y = mu-sigma, yend = mu+sigma))+
    scale_color_identity()+
    ylab("emission parameter value")+
    scale_x_continuous(breaks = unique(emis$state_num), labels = hsmm$states$abbr)+
    scale_shape_manual(values = c(16,4))+
    coord_cartesian(ylim = c(-1.5,1.5))+
    facet_grid(variable ~., scale = "free_y")+
    theme(strip.text.y = element_text(angle = 0, hjust = 0),
          legend.position = "left")
  
  return(g)
}

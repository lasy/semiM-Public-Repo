

plot_sim_performances = function(m , 
                                 model_base, model_full,
                                 ground_truth, ground_truth_base,
                                 vit_base_init, vit_base_fitted, 
                                 vit_full_init, vit_full_fitted,
                                 base_fit_param, 
                                 full_fit_param){
  
  g_title = ggplot()+ ggtitle(str_c("Simulated data experiment: m = ",m))
  
  g_base = plot_perf(name = "Base", model = model_base, fit_param = base_fit_param, 
                     GT = ground_truth_base, 
                     vit_init = vit_base_init, vit_fitted = vit_base_fitted)
  
  g_full = plot_perf(name = "Full", model = model_full, fit_param = full_fit_param, 
                     GT = ground_truth, 
                     vit_init = vit_full_init, vit_fitted = vit_full_fitted)
  
  g_seq = plot_seq(ground_truth = ground_truth, ground_truth_base = ground_truth_base,
                   vit_base_init = vit_base_init, vit_base_fitted = vit_base_fitted, 
                   vit_full_init = vit_full_init, vit_full_fitted = vit_full_fitted)
  
  g = ggarrange(
    g_title, 
    g_base,
    g_full,
    g_seq,
    heights = c(1,4,4,4), ncol = 1
  )
  
  g
  
}




plot_perf = function(name, 
                     model,
                     fit_param, 
                     GT, 
                     vit_init, vit_fitted){
  
  df_acc = bind_rows(vit_init   %>%  select(-state_deprecated) %>% mutate(model = "init"),
                     vit_fitted %>%  select(-state_deprecated) %>% mutate(model = "fitted")) 
  df_acc = df_acc %>% full_join(GT %>% rename(state_GT = state), by = c("seq_id","t"))
  df_acc_overall = df_acc %>% group_by(model) %>% 
    summarize(value = mean(state == state_GT, na.rm = TRUE)) %>% mutate(accuracy = "overall accuracy")
  df_acc_M = df_acc %>% filter(state_GT == 1) %>% group_by(model) %>% 
    summarize(value = mean(state == state_GT, na.rm = TRUE)) %>% mutate(accuracy = "accuracy on M")
  df_acc_state_av = df_acc %>% group_by(model, state_GT) %>% 
    summarize(value = mean(state == state_GT, na.rm = TRUE)) %>%  group_by(model) %>%
    summarize(value = mean(value,na.rm = TRUE)) %>% mutate(accuracy = "average state accuracy") 
  df = bind_rows(df_acc_overall, df_acc_M, df_acc_state_av)
  df = df %>% mutate(model = model %>% factor(.,levels = c("init","fitted")),
                     accuracy = accuracy %>% factor(., levels = unique(df$accuracy) ))
  
  g_acc = ggplot(df, aes(x = model, y = value, fill = value))+ 
    geom_bar(stat = "identity") + 
    geom_text(aes(label = round(100*value,1), 
                  y = value + ifelse(value < 0.75, 0.05, -0.05),
                  vjust = ifelse(value < 0.75, 0, 1))) +
    ylim(c(0,1))+ guides(fill = FALSE)+
    scale_fill_gradient(low = "tomato",high = "steelblue", limits = c(0,1))+
    facet_grid(accuracy ~ .)+
    ggtitle(str_c("Model: ", name))
  
  g_fit = ggplot(data.frame(iter = 1:fit_param$n_iter, ll = fit_param$ll), 
                 aes(x = iter, y = ll))+ 
    geom_line() + geom_point()+
    ggtitle(fit_param$message)
  
  GT_matched = full_join(vit_init %>% select(seq_id, t) , GT, by = c("seq_id","t"))
  
  g_mat_init = ggplot_confusion_matrix(true = model$state_names[GT_matched$state], 
                                       decoded = model$state_names[vit_init$state], 
                                       states = model$state_names )
  
  g_mat_fitted = ggplot_confusion_matrix(true = model$state_names[GT_matched$state], 
                                         decoded = model$state_names[vit_fitted$state], 
                                         states = model$state_names )
  
  g = ggarrange(g_acc, g_fit, g_mat_init, g_mat_fitted, nrow = 1, widths = c(1,1,2,2))
  
  g
  
}


plot_seq = function(ground_truth, ground_truth_base,
                    vit_base_init, vit_base_fitted, 
                    vit_full_init, vit_full_fitted){
  
  df_seq = bind_rows(ground_truth_base %>% 
                       mutate(name = "GT_base", state_col = FAM_base_init$state_colors[state]),
                     vit_base_init %>% select(seq_id ,t , state) %>% 
                       mutate(name = "vit_base_I", state_col = FAM_base_init$state_colors[state]),
                     vit_base_fitted %>% select(seq_id ,t , state) %>% 
                       mutate(name = "vit_base_F", state_col = FAM_base_init$state_colors[state]),
                     ground_truth %>% 
                       mutate(name = "GT_full", state_col = FAM_full_init$state_colors[state]),
                     vit_full_init %>% select(seq_id ,t , state) %>% 
                       mutate(name = "vit_full_I", state_col = FAM_full_init$state_colors[state]),
                     vit_full_fitted %>% select(seq_id ,t , state) %>% 
                       mutate(name = "vit_full_F", state_col = FAM_full_init$state_colors[state]))
  df_seq = df_seq %>% 
    mutate(name = name %>% factor(., levels = rev(unique(df_seq$name))))
  g_seq = ggplot(df_seq, aes(x = t, y =  name, fill = state_col))+
    geom_tile()+ scale_fill_identity()+
    facet_grid(seq_id ~ . )
  
  g_seq
}






###------- 



plot_experiment_performances = function(output){
  
  # internal variables
  decodings = c("state_VI", "state_SI", "state_VF", "state_SF" )
  decoding_cols = c("steelblue2","steelblue3","pink2","pink3")
  
  
  # convergence of fit
  df = output$fitting_par
  df = df %>% 
    mutate(
      message_col = case_when(
        str_detect(message, "^Converged") ~ "steelblue",
        str_detect(message, "^Reached") ~ "orange",
        TRUE ~ "red"
  ),
  message_short = str_extract(message, "(\\w+)"),
  message_short = ifelse(message_short == "Reached", "Reached max iter",message_short),
  message_short = message_short %>% factor(.,levels = c("Converged","Reached max iter","Error")))
  
  g_fit_conv = ggplot(df, aes(x = iter, y = ll, col = message_short, group = seq_id))+
    geom_line()+
    geom_point()+
    xlab("EM Iteration")+ylab("Log Likelihood")+
    scale_x_continuous(breaks = df$iter)+
    scale_color_manual(name = "", values = c("steelblue1","orange","red"))+
    theme(legend.position = "bottom")+
    ggtitle("Fitting: EM status")
  
  # overal accuracies
  output_long = output$X_ad %>% select(seq_id, t, state_GT, all_of(decodings)) %>% 
    pivot_longer(cols = all_of(decodings), names_to = "decoding",values_to = "state_decoded") %>% 
    mutate(decoding = decoding %>% factor(., levels = decodings))
  
  df = output_long %>% 
    group_by(decoding) %>% 
    summarize(accuracy = mean(state_GT == state_decoded, na.rm = TRUE)*100) %>% 
    mutate(decoding_short = str_remove(decoding, "state_") %>% factor(.,levels = c("VI","SI","VF","SF")))
  
  g_overal_acc = ggplot(df, aes(x = decoding_short, y = accuracy, fill = decoding))+
    geom_bar(stat = "identity")+
    geom_text(aes(label = round(accuracy), col = decoding), vjust = 0, nudge_y = 1)+
    scale_y_continuous(limits = c(0,100), breaks = seq(0,100,by = 20))+
    xlab("Decoding")+
    scale_fill_manual(values = decoding_cols)+
    scale_color_manual(values = decoding_cols)+
    guides(col = FALSE, fill = FALSE)+
    ggtitle("Overal accuracies")
  
  # overal accuracies (average by state)
  df = output_long %>% 
    group_by(decoding, state_GT) %>% 
    summarize(accuracy = mean(state_GT == state_decoded, na.rm = TRUE)*100) %>% 
    group_by(decoding) %>% 
    summarize(accuracy = mean(accuracy, na.rm = TRUE)) %>% 
    mutate(decoding_short = str_remove(decoding, "state_") %>% factor(.,levels = c("VI","SI","VF","SF")))
  
  g_acc_mean_state_acc = ggplot(df, aes(x = decoding_short, y = accuracy, fill = decoding))+
    geom_bar(stat = "identity")+
    geom_text(aes(label = round(accuracy), col = decoding), vjust = 0, nudge_y = 1)+
    scale_y_continuous(limits = c(0,100), breaks = seq(0,100,by = 20))+
    xlab("Decoding")+
    scale_fill_manual(values = decoding_cols)+
    scale_color_manual(values = decoding_cols)+
    guides(col = FALSE, fill = FALSE)+
    ggtitle("Average of state accuracies")
  
  

  
  # sequence accuracies
  
  df = output_long %>% 
    group_by(seq_id, decoding) %>% 
    summarize(accuracy = mean(state_GT == state_decoded, na.rm = TRUE)*100) %>% 
    mutate(decoding_short = str_remove(decoding, "state_") %>% factor(.,levels = c("VI","SI","VF","SF")))
  
  g_seq_acc = ggplot(df, aes(x = decoding_short, y = accuracy, fill = decoding))+
    geom_bar(stat = "identity")+
    geom_text(aes(label = round(accuracy), col = decoding), vjust = 0, nudge_y = 1)+
    scale_y_continuous(limits = c(0,100), breaks = seq(0,100,by = 20))+
    facet_grid(. ~ seq_id)+
    xlab("Decoding")+
    scale_fill_manual(values = decoding_cols)+
    scale_color_manual(values = decoding_cols)+
    guides(col = FALSE, fill = FALSE)+
    ggtitle("Accuracies by sequence")
  
  # confusion matrices
  g_conf_mat = plot_exp_confusion_matrices(X_ad = output$X_ad, model = output$model, decodings = decodings)
  
  # event timing
  # g_event = plot_timing_accuracy(X_ad = output$X_ad, state_names = output$model$state_names, decodings = decodings)
  
  # state_probs
  g_state_probs = plot_state_prob(output = output, model = output$model)
  
  g = ggarrange(
    ggarrange(g_fit_conv, g_overal_acc, g_acc_mean_state_acc, g_seq_acc, 
              ncol = 4, nrow = 1, widths = c(1,1,1,2)),
    g_conf_mat,
    #g_event,
    g_state_probs,
    ncol = 1, nrow = 3
  )
  
  g
}



plot_user_history= function(d = data.frame(),
                            o = data.frame(),
                            states_labels_list = list(), 
                            current_selection = data.frame(), 
                            show_all_variables = FALSE){
  
  if((lu(d$user_id)>1)|(lu(o$user_id)>1)){stop("there is more than one user in the data frame\n")}
  user_id = unique(c(d$user_id, o$user_id))
  if((nrow(d) + nrow(o)) == 0){stop("both 'd' and 'o' are empty \n")}
  
  # AXES
  start = min(d$rel_date, o$rel_date)-0.55
  end = max(d$rel_date, o$rel_date)+0.55
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
  
  if(nrow(o)>0){
    for(c in 3:ncol(o)){
      col_name = colnames(o)[c]
      df = data.frame(rel_date = o$rel_date, y = o[,c] %>%  unlist() %>%  set_names(NULL))
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





ggplot_confusion_matrix = function(true , decoded, states){
  
  J = length(states)
  
  true = true %>% factor(., levels = states)
  decoded = decoded %>% factor(., levels = states)
  
  cm = table(true = true, decoded =  decoded)
  
  cm_sum = rowSums(cm, na.rm = TRUE)
  (cm/cm_sum*100) %>%  round(.,1) 
  cm = cm/cm_sum *100
  
  cm_long = data.frame(cm)
  cm_long$true = factor(cm_long$true, levels = rev(states))
  
  confusion_matrix_viz = ggplot(cm_long, aes(x = decoded, y = true , fill = Freq))+ coord_fixed()+
    geom_tile()+
    geom_abline(slope = -1, intercept = J+1, col = "black", size = 0.2, linetype = 2)+
    geom_hline(yintercept = 0:J+0.5, col = "black", size = 0.1)+
    geom_vline(xintercept = c(0.5,J+0.5), col = "black", size = 0.1)+
    scale_fill_gradient(name = "% of manual labels",low = "white", high = "royalblue3")+
    ylab("Ground\ntruth")+
    xlab("Decoded labels")+
    theme(legend.position = "bottom",
          legend.key.height = unit(1/70,units = 'npc'),
          legend.title = element_text(vjust = 0.75),
          axis.title.y = element_text(angle = 0, hjust = 1, vjust = 1),
          axis.title.x = element_text(hjust = 1),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  return(confusion_matrix_viz)
}


plot_exp_confusion_matrices = function(X_ad, model = FAM_init, decodings = decodings){
  
  state_names = model$state_names
  
  gVI = ggplot_confusion_matrix(true = state_names[X_ad$state_GT], decoded = state_names[X_ad$state_VI], states = state_names)+ggtitle("Viterbi Initial Model")+theme(legend.position = "none")
  gSI = ggplot_confusion_matrix(true = state_names[X_ad$state_GT], decoded = state_names[X_ad$state_SI], states = state_names)+ggtitle("Smoothed Initial Model")+theme(legend.position = "none")
  gVF = ggplot_confusion_matrix(true = state_names[X_ad$state_GT], decoded = state_names[X_ad$state_VF], states = state_names)+ggtitle("Viterbi Fitted Model")+theme(legend.position = "none")
  gSF = ggplot_confusion_matrix(true = state_names[X_ad$state_GT], decoded = state_names[X_ad$state_SF], states = state_names)+ggtitle("Smoothed Fitted Model")+theme(legend.position = "none")
  
  g = cowplot::plot_grid(plotlist = list(gVI, gSI, gVF, gSF), nrow = 1)
  g
}


plot_average_state_accuracy = function(output_long, model, state_cols = NULL){
  
  if(is.null(state_cols)) state_cols = viridis_pal(option = "C")(model$J)
  
  df_conf_mat = output_long %>% 
    group_by(decoding, state_GT) %>% 
    summarize(state_accuracy = mean(state_decoded == state_GT, na.rm = TRUE))
  
  df_conf_mat_metric = df_conf_mat %>% 
    group_by(decoding) %>% 
    summarize(average_state_accuracy = mean(state_accuracy, na.rm = TRUE)) %>% 
    mutate(decoding_with_acc = 
             str_c(decoding, "\n av. state acc: ", round(100*average_state_accuracy), "%")) %>% 
    mutate(decoding_with_acc = decoding_with_acc %>%  factor(.,levels = decoding_with_acc))
  
  df_conf_mat_metric
  
  df_conf_mat = df_conf_mat %>% 
    mutate(state_col = state_cols[state_GT]) %>% 
    full_join(.,df_conf_mat_metric, by = "decoding") 
  
  g = ggplot(df_conf_mat, aes(x = state_GT, y = state_accuracy, col = state_col))+
    geom_point()+ 
    scale_color_identity()+
    xlab("States (GT)")+ ylab("State accuracy")+
    facet_grid(. ~ decoding_with_acc)
  
  g
  
}





ggplot_event_accuracy = function(event_accuracy){
  event_accuracy$no_event_in_30d_window = is.na(event_accuracy$decoded_rel_date)
  
  g = ggplot(event_accuracy, aes(x = time_diff, fill = no_event_in_30d_window)) +
    geom_histogram(binwidth = 1)+
    geom_vline(xintercept = 0, linetype = 2)+
    facet_grid(event ~ ., scale = "free")+
    scale_fill_manual(values = c("turquoise","red"))+
    xlab("difference in days\n(decoded - manual label)")+
    xlim(c(-17.5,15.5))+
    guides(fill = FALSE)
  return(g)
}


plot_timing_accuracy = function(X_ad = X_ad, state_names = hsmm$states$names, decodings = decodings){
  if(is.null(state_names)) state_names = 1:model$J
  
  plotlist = list()
  
  for(decoding in decodings){
    
    ev_acc = event_accuracy(
      decoding = X_ad %>% select(seq_id, t, all_of(decoding)) %>% 
        set_colnames(c("user_id","rel_date","state")) %>% 
        mutate(state_name = state_names[state]), 
      manual_labels = X_ad %>% select(seq_id, t, state_GT) %>% 
        set_colnames(c("user_id","rel_date","state")) %>% 
        mutate(state_name = state_names[state]),
      events = c("O","B","L"))
    
    ev_acc_summary = ev_acc %>% 
      group_by(event) %>% 
      summarize(timing_error = mean(time_diff^2) %>% sqrt()) %>% 
      mutate(title = str_c(event,": ",round(timing_error,2)))
    
    metric = mean(ev_acc_summary$timing_error)
    
    g = ggplot_event_accuracy(event_accuracy = ev_acc)+
      ggtitle(str_c(decoding," (mean error: ",round(metric, 2),")",
                    "\n", str_c(ev_acc_summary$title,collapse = " | ")))
    plotlist[[decoding]] = g
    
  }
  
  g = cowplot::plot_grid(plotlist = plotlist, nrow = 1)
  g
  
}



plot_state_prob = function(output, model){
  
  state_names = model$state_names
  state_cols = model$state_colors
  
  state_probs = bind_rows(
    output$smoothed_prob$init %>% mutate(decoding = "I"),
    output$smoothed_prob$fitted %>% mutate(decoding = "F")
  ) %>% mutate(
    decoding = decoding %>% factor(.,levels = c("I", "F"))
  )
  
  state_probs = state_probs %>% 
    mutate(state_name = state_names[state] %>% factor(.,levels = state_names),
           state_col = state_cols[state])
  
  state_probs  = left_join(output$X_ad %>% select(seq_id, t, state_GT) %>% rename(state = state_GT) %>% filter(!is.na(state)),
                           state_probs, 
                           by = c("seq_id","t","state"))
  
  ggplot(state_probs, aes(x = state_name, y = posterior, fill = state_col, color = state_col))+
    geom_violin()+
    scale_fill_identity()+
    scale_color_identity()+
    facet_grid(. ~ state_name + decoding, scale = "free")
  
}


percentage_of_incorrect_decoding_viz = function(decoding = rule_based_decoding){
  
  decoding_agg = decoding %>% 
    arrange(user_id) %>% 
    group_by(user_id) %>% 
    dplyr::summarise(n_days = n(),
                     n_days_incorrect_decoding = sum(!correct_decoding, na.rm = TRUE),
                     percent_incorrect_decoding = n_days_incorrect_decoding/n_days *100,
                     n_days_with_manual_labels = sum(!is.na(correct_decoding)),
                     percent_with_manual_labels = n_days_with_manual_labels/n_days *100) %>% 
    arrange(percent_incorrect_decoding, percent_with_manual_labels) %>% 
    ungroup() %>% 
    dplyr::mutate(user_id_short = substr(user_id,1,5),
                  user_id_short = factor(user_id_short, levels = user_id_short))
  
  
  
  perc_of_incorrect_decoding_per_user = 
    ggplot(decoding_agg, aes(x = user_id_short, y = percent_incorrect_decoding))+
    geom_bar(stat = "identity")+
    geom_point(aes(y = -0.5, col = percent_with_manual_labels))+
    scale_color_gradient(name = "% of sequence manually labelled" ,low = "white",high = "blue4", limits = c(0,100))+
    coord_flip()+ xlab("user ID")+ ylab("% of the user's time-series with incorrect decoding")+
    guides(fill = FALSE)+
    theme(legend.position = "bottom", legend.key.height = unit(5,"pt"))
  perc_of_incorrect_decoding_per_user
  
  
  return(list(decoding_agg = decoding_agg, plot = perc_of_incorrect_decoding_per_user))
}



viz_em_parm = function(p1, p2, hsmm){
  
  N = 1000
  
  ## 1
  r_obs = foreach(state = 1:hsmm$n_states, .combine = rbind) %do%{
    this_state_r_obs = generate_random_obs(n = N, parem = p1, state = state) %>% 
      as.data.frame() %>% set_colnames(hsmm$obs_names) %>% 
      mutate(state_abbr =  hsmm$states$abbr[state], color = hsmm$state$colors[state],n = 1:N) %>% 
      tidyr::pivot_longer(.,cols = hsmm$obs_names, names_to = "obs_name")
    return(this_state_r_obs)
  }
  
  r_obs$state_abbr = factor(r_obs$state_abbr, levels = hsmm$states$abbr)
  r_obs$obs_name = factor(r_obs$obs_name, levels = hsmm$obs_names)
  
  r_obs1 = r_obs
  
  ## 2
  r_obs = foreach(state = 1:hsmm$n_states, .combine = rbind) %do%{
    this_state_r_obs = generate_random_obs(n = N, parem = p2, state = state) %>% 
      as.data.frame() %>% set_colnames(hsmm$obs_names) %>% 
      mutate(state_abbr =  hsmm$states$abbr[state], color = hsmm$state$colors[state],n = 1:N) %>% 
      tidyr::pivot_longer(.,cols = hsmm$obs_names, names_to = "obs_name")
    return(this_state_r_obs)
  }
  
  r_obs$state_abbr = factor(r_obs$state_abbr, levels = hsmm$states$abbr)
  r_obs$obs_name = factor(r_obs$obs_name, levels = hsmm$obs_names)
  
  r_obs2 = r_obs
  
  ## combining them
  
  r_obs1$type = 1
  r_obs2$type = 2
  r_obs = rbind(r_obs1, r_obs2)
  
  g_par = ggplot_viz_em_parm(r_obs)
  
  return(g_par)
}

ggplot_viz_em_parm = function(r_obs){
  
  if("type" %in% colnames(r_obs)){r_obs$type = as.numeric(r_obs$type)}else{r_obs$type = 1}
  
  precision = 1/5
  r_obs_h = r_obs %>% mutate(rounded_value = round(value/precision)*precision) %>% 
    group_by(type, state_abbr, color, obs_name, rounded_value) %>% 
    dplyr::summarise(n = n()) %>% 
    group_by(obs_name) %>%  
    mutate(max_n = max(n), f = n/max_n)
  
  g_par = ggplot(r_obs_h, aes(xmin = type-1, xmax = type-1+f, ymin = rounded_value-precision/2, ymax = rounded_value+precision/2, fill = color))
  g_par = g_par + 
    geom_vline(xintercept = seq(0,max(r_obs$type)-1, by = 1), col = "gray85")+
    geom_rect()+
    scale_color_identity()+scale_fill_identity()+
    facet_grid( obs_name ~ state_abbr, scales = "free")+
    xlab("pdf (scaled by its max. value per variable for clarity of the plot)")+
    ylab("")+
    scale_y_continuous(breaks = seq(-5,5,by = 1), limits = function(x){x[1] = min(x[1],-0.75); x[2] = max(x[2],1.75); return(x)})+
    scale_x_continuous(breaks = 0)+
    theme(strip.text.x = element_text(angle = 90, hjust = 0)) #,strip.text.y = element_text(angle = 0, hjust = 0)
  
  g_par
  
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

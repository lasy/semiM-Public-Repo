

A = matrix(sample(1:10,10*3, replace = TRUE),10,3)
c = matrix(sample(c(1,1,1,1,0),10*3, replace = TRUE),10,3)
B = A
B[c==0] = NA

not_censored = cov.wt.censored(A)
not_censored$center
not_censored$cov

censored = cov.wt.censored(B)
censored$center
censored$cov


###########################





most_probable_value = function(par, state){
  if(par$type == "norm"){
    val = par$param$mean[state]
  }else if(par$type == "binom"){
    val = par$param$prob[state] %>%  round()
  }else{stop("This distribution has not been implemented yet")}
}

compute_prob_obs = function(x, state, par){
  ix = which(is.na(x))
  if(length(ix)>0){x[ix] = most_probable_value(par,state)}
  if(!(par$type %in% c("norm","binom"))){stop("This distribution has not been implemented yet")}
  cmd = str_c("p = d",par$type,"(x,",
              str_c(names(par$param)," = par$param$",names(par$param),"[state]") %>% str_c(collapse = ",")
              ,")")
  eval(parse(text = cmd))
  return(p)
}


compute_prob = function(obs, state, parem){
  p = sapply(1:length(parem), function(o) compute_prob_obs(obs[,o], state, parem[[o]])) 
  pp = apply(p, 1, FUN = prod)
  return(pp)
}



parem = list(bleeding = list(type = "norm", param = list(mean = 1:2, sd = rep(1,2))),
             LH = list(type = "binom", param = list(size = rep(1,2),prob = seq(0.2,0.8,len = 2)))
                             )

obs = data.frame(bleeding = c(1,2,1,2,NA),LH = c(0,0,1,1,1)) %>%  as.matrix()


matplot(obs, type = "b", lty = 1, pch = 16)


p = sapply(1:2, function(state) compute_prob(obs, state, parem))

matplot(p, type = "l", col = hsmm$states$colors, lty = 1)


w = p/rowSums(p)


compute_new_em_parms = function(obs, w, parem){
  new.parem = parem
  for(i in 1:length(parem)){
    cat(i, "\n")
    if(new.parem[[i]]$type == "norm"){
      new.parem[[i]]$param$mean = sapply(1:ncol(w),function(state) weighted.mean(obs[,i],w = w[,state], na.rm = TRUE))
      new.parem[[i]]$param$sd = sapply(1:ncol(w),function(state) wt.sd(obs[,i],w = w[,state]))
    }else if(new.parem[[i]]$type == "binom"){
      new.parem[[i]]$param$prob = sapply(1:ncol(w),function(state) weighted.mean(obs[,i],w = w[,state], na.rm = TRUE))
    }else{ stop("This type of distribution has not been handled yet")}
  }
  return(new.parem)
}

new.parem = compute_new_em_parms(obs, w, parem)


p = sapply(1:2, function(state) compute_prob(obs, state, new.parem))

matplot(p, type = "l", col = hsmm$states$colors, lty = 1)

w = p/rowSums(p)


new.parem = compute_new_em_parms(obs, w, new.parem)


p = sapply(1:2, function(state) compute_prob(obs, state, new.parem))

matplot(p, type = "l", col = hsmm$states$colors, lty = 1)



r_obs = generate_random_obs(100, parem = parem, state = 1)






############


oo = obs %>% mutate(state = labels$state) %>% dplyr::filter(!is.na(state)) %>% dplyr::select(str_c(hsmm$obs_names,"_score"), state)
ool = tidyr::pivot_longer(oo, cols = contains("_score"))

ggplot(ool, aes(x = value, fill = name))+
  geom_histogram(binwidth = 0.1)+
  facet_grid(state ~ name, scale = "free")



# DENSITY FUNCTIONS

compute_prob = function(obs, state, parem, w = NULL){
  if(length(obs) == length(parem)){obs = matrix(obs, nrow = 1, ncol = length(parem))}
  p = sapply(1:length(parem), function(o) compute_prob_obs(obs[,o], state, parem[[o]]))
  if(length(p)>length(parem)){pp = apply(p, 1, FUN = prod)}else{pp = prod(p)}
  if(is.null(w)){w = rep(1, length(parem)); w = w/sum(w)}
  ww = colSums(w * t(!is.na(obs)))
  pp = ww*pp
  return(pp)
}


compute_prob_obs = function(x, state, par){
  ix = which(is.na(x))
  if(length(ix)>0){x[ix] = most_probable_value(par,state)}
  if(!(par$type %in% c("norm","binom","non-par"))){stop("This distribution has not been implemented yet")}
  if(!(par$type == "non-par")){
  cmd = str_c("p = d",par$type,"(x,",
              str_c(names(par$param)," = par$param$",names(par$param),"[state]") %>% str_c(collapse = ",")
              ,")")
  eval(parse(text = cmd))
  }else{
    j = match(as.character(x), as.character(par$param$values))
    p = par$param$probs[j, state]
  }
  return(p)
}


most_probable_value = function(par, state){
  if(par$type == "norm"){
    val = par$param$mean[state]
  }else if(par$type == "binom"){
    val = par$param$prob[state] %>%  round()
  }else if(par$type == "non-par"){
    val = par$param$values[max.col(t(par$param$probs[,state]))]
  }else{stop("This distribution has not been implemented yet")}
  return(val)
}


# RANDOM GENERATION FUNCTION

generate_random_obs = function(n = 10, parem, state){
  r_obs = sapply(1:length(parem), function(i) generate_random_obs_par(n, parem[[i]],state))
  #colnames(r_obs) = names(parem)
  return(r_obs)
}


generate_random_obs_par = function(n, par, state){
  if(par$type == "norm"){
    x = rnorm(n, mean = par$param$mean[state], sd = par$param$sd[state])
  }else if(par$type == "binom"){
    x = rbinom(n, size = par$param$size[state], prob = par$param$prob[state])
  }else if(par$type == "non-par"){
    x = sample(rep(par$param$values, round(par$param$probs[,state]*1000)) , n)
  }else{stop("This distributions hasn't been implemented yet")}
  return(x)
}



# MSTEP (estimate parameters functions)

library(SDMTools)

compute_new_em_parms = function(obs, w, parem){
  new.parem = parem
  for(i in 1:length(parem)){
    #cat(i, "\n")
    if(new.parem[[i]]$type == "norm"){
      new.parem[[i]]$param$mean = sapply(1:ncol(w),function(state) weighted.mean(obs[,i],w = w[,state], na.rm = TRUE))
      #new.parem[[i]]$param$sd = sapply(1:ncol(w),function(state) sqrt(sum(w[,state] * (obs[,i] - new.parem[[i]]$param$mean[state])^2)))
      new.parem[[i]]$param$sd = sapply(1:ncol(w),function(state) wt.sd(obs[,i],w = w[,state])) %>% pmax(.,0.1)
    }else if(new.parem[[i]]$type == "binom"){
      new.parem[[i]]$param$prob = sapply(1:ncol(w),function(state) weighted.mean(obs[,i],w = w[,state], na.rm = TRUE)) %>% pmin(.,0.99) %>% pmax(.,0.01)
    }else if(new.parem[[i]]$type == "non-par"){
      new.parem[[i]]$param$probs = sapply(1:ncol(w), 
                                          function(state) {
                                            x = rep(obs[,i], round(w[,state]*5)) %>% factor(.,levels = parem[[i]]$param$values); 
                                            tt = table(x);return(t(tt/max(tt)))})
    }else{ stop("This type of distribution has not been handled yet")}
  }
  return(new.parem)
}





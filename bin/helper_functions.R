gen_corr_mats <- function(n, d, alpha_func = NULL){
  
  corr_mats = list()
  pb = progress_bar$new(total = n, format = "[:bar] :current/:total :eta")
  for(i in 1:n){
    pb$tick()
    corr_mats[[i]] = rand_corr_wrapper(d, alpha_func)
    perm = sample(1:d)
    corr_mats[[i]] = corr_mats[[i]][perm, perm]
  }
  
  return(corr_mats)
  
}
.rbo.ext <- function(x, y, p, k, uneven.lengths = TRUE) {
  if (length(x) <= length(y)) {
    S <- x
    L <- y
  } else {
    S <- y
    L <- x
  }
  l <- min(k, length(L))
  s <- min(k, length(S))
  
  if (uneven.lengths) {
    Xd <- sapply(1:l, function(i) length(intersect(S[1:i], L[1:i])))
    ((1-p) / p) *
      ((sum(Xd[seq(1, l)] / seq(1, l) * p^seq(1, l))) +
         (sum(Xd[s] * (seq(s+1, l) - s) / (s * seq(s+1, l)) * p^seq(s+1, l)))) +
      ((Xd[l] - Xd[s]) / l + (Xd[s] / s)) * p^l  
  } else {
    #stopifnot(l == s)
    k <- min(s, k)
    Xd <- sapply(1:k, function(i) length(intersect(x[1:i], y[1:i])))
    Xk <- Xd[k]
    (Xk / k) * p^k + (((1-p)/p) * sum((Xd / seq(1,k)) * p^seq(1,k)))
  }
}
rbo <- function(s, t, p, k=floor(max(length(s), length(t))/2), side=c("top", "bottom"), mid=NULL, uneven.lengths = TRUE) {
  side <- match.arg(side)
  if (!is.numeric(s) | !is.numeric(t))
    stop("Input vectors are not numeric.")
  if (is.null(names(s)) | is.null(names(t)))
    stop("Input vectors are not named.")
  ids <- switch(side,
                "top"=list(s=.select.ids(s, "top", mid), t=.select.ids(t, "top", mid)),
                "bottom"=list(s=.select.ids(s, "bottom", mid), t=.select.ids(t, "bottom", mid))
  )
  min(1, .rbo.ext(ids$s, ids$t, p, k, uneven.lengths = uneven.lengths))
}
thres_abs <- function(net, thres, weighted){
  
  if(!weighted){
    return(1*(net > thres))
  }else{
    net[which(net < thres, arr.ind = T)] = 0
    return(net)
  }
}
thres_prop <- function(net, prop){
  target = quantile(net[upper.tri(net)], probs = 1-prop)
  return(1*(net > target))
}
thres_perc <- function(net, weighted){
  net[which(net < 0, arr.ind = T)] = 0
  edge_weights = sort(net[upper.tri(net)],decreasing =  F)
  toggle = T
  counter = 1
  temp = graph.adjacency(net, mode = "undirected", weighted = T, diag = F)
  nconn = components(temp, mode = "weak")$csize[1]
  for(i in 1:length(edge_weights)){
    temp_net = net
    temp_net[which(temp_net < edge_weights[i])] = 0
    temp = graph.adjacency(temp_net, mode = "undirected", weighted = T, diag = F)
    nconn2 = components(temp, mode = "weak")$csize[1]
    
    if(nconn2 < nconn){
      temp_net = net
      temp_net[which(net < edge_weights[i-1])] = 0
      break
      
    }
    counter = counter + 1
  }
  if(weighted){
    temp_net[which(temp_net < 0, arr.ind = T)] = 0
  }else{
    temp_net = 1*(temp_net > 0)
  }
  return(temp_net)
}
thres_pos <- function(net, weighted){
  temp_net = net
  if(weighted){
    temp_net[which(temp_net < 0, arr.ind = T)] = 0
  }else{
    temp_net = 1*(net > 0)
    print(temp_net)
  }
  return(temp_net)
  
}
extract_RBO = function(res_df,weighted_var_name ,var_name, param_name){
  
  group_levels = unique(res_df[,param_name])
  res_vec = vector()
  pb = progress_bar$new(total = length(group_levels), format = "[:bar] :current/:total :eta")
  for(i in 1:length(group_levels)){
    subset = res_df[which(res_df[,param_name] == group_levels[i]),]
    weigh =  subset[,weighted_var_name]
    names(weigh) = subset[, "ID"]
    comp = subset[, var_name]
    names(comp) = subset[, "ID"]
    res_vec[i] = rbo(weigh, comp,p = .50)
    pb$tick()
  }
  
  return(data.frame(param = group_levels, rbo = res_vec)) 
  
}
apply_modularity = function(nets, abs_thres = seq(.01, .7, .01), prop_thres = seq(.01, .5, .01)){
  toReturn_list = list() 
  pb = progress_bar$new(total = length(nets), format = "[:bar] :current/:total :eta")
  pb$tick(0)
  for(i in 1:length(nets)){
    pos_thres = thres_pos(nets[[i]], weighted = T)
    
    weighted_net = graph.adjacency(pos_thres, mode = "undirected", weighted = T, diag = F)
    
    
    abs_thres_nets = list()
    for(j in 1:length(abs_thres)){
      
      temp = thres_abs(net = nets[[i]], thres = abs_thres[j], weighted = F)
      abs_thres_nets[[j]] = graph.adjacency(temp, mode = "undirected",  diag = F)
      
    } 
    
    
    prop_thres_nets = list()
    for(j in 1:length(prop_thres)){
      
      temp = thres_prop(net = nets[[i]], prop = prop_thres[j])
      prop_thres_nets[[j]] = graph.adjacency(temp, mode = "undirected",  diag = F)
      
    } 
    
    perc_thres_net = thres_perc(net = nets[[i]], weighted = F)
    perc_thres_net =  graph.adjacency(perc_thres_net, mode = "undirected", diag = F)
    
    weighted_comms = walktrap.community(weighted_net)
    
    weighted_comms = walktrap.community(weighted_net)
    perc_thres_comms = walktrap.community(perc_thres_net)
    abs_thres_comms = lapply(abs_thres_nets, FUN = function(x){return(walktrap.community(x))})
    prop_thres_comms = lapply(prop_thres_nets, FUN = function(x){return(walktrap.community(x))})
    
    weighted_mod = modularity(weighted_comms)
    perc_thres_mod = modularity(perc_thres_comms)
    abs_thres_mod = lapply(abs_thres_comms, FUN = function(x){return(modularity(x))})
    prop_thres_mod = lapply(prop_thres_comms, FUN = function(x){return(modularity(x))})
    
    perc_thres_ARI = aricode::ARI(c1 = as.numeric(membership(weighted_comms)), c2=as.numeric(membership(perc_thres_comms)))
    abs_thres_ARI = lapply(abs_thres_comms, FUN = function(x){aricode::ARI(c1 = as.numeric(membership(weighted_comms)), c2=as.numeric(membership(x)))})
    prop_thres_ARI = lapply(prop_thres_comms, FUN = function(x){aricode::ARI(c1 = as.numeric(membership(weighted_comms)), c2=as.numeric(membership(x)))})
    
    
    abs_labs = paste("abs_thres_", abs_thres,"_mod",  sep = "")
    prop_labs = paste("prop_thres_", prop_thres, "_mod",  sep = "")
    abs_ARIlabs = paste("abs_thres_", abs_thres,"_ARI",  sep = "")
    prop_ARIlabs = paste("prop_thres_", prop_thres, "_ARI",  sep = "")
    
    
    abs_df = data.frame(thres = abs_thres,thres_mod = as.numeric(abs_thres_mod),thres_ARI = as.numeric(abs_thres_ARI),thres_type = "Abs")
    
    prop_df = data.frame(thres = prop_thres, thres_mod = as.numeric(prop_thres_mod),thres_ARI = as.numeric(prop_thres_ARI), thres_type = "Prop")
    
    
    toReturn = data.frame(ID = i, weighted_mod = weighted_mod, perc_thres_mod = perc_thres_mod, perc_thres_ARI = perc_thres_ARI)
    toReturn = rbind(cbind(toReturn, abs_df),cbind(toReturn, prop_df))
    toReturn_list[[i]] = toReturn
    pb$tick()
  } 
  toReturn = do.call("rbind", toReturn_list)
  return(toReturn)   
}
apply_GE = function(nets, abs_thres = seq(.01, .7, .01), prop_thres = seq(.01, .5, .01)){
  toReturn_list = list() 
  pb = progress_bar$new(total = length(nets), format = "[:bar] :current/:total :eta")
  pb$tick(0)
  for(i in 1:length(nets)){
    pos_thres = thres_pos(nets[[i]], weighted = T)
    temp = 1/pos_thres
    temp[which(is.infinite(temp))] = 0
    weighted_net = graph.adjacency(temp, mode = "undirected", weighted = T, diag = F)
    
    
    abs_thres_nets = list()
    for(j in 1:length(abs_thres)){
      
      temp = thres_abs(net = nets[[i]], thres = abs_thres[j], weighted = F)
      abs_thres_nets[[j]] = graph.adjacency(temp, mode = "undirected",  diag = F)
      
    } 
    
    
    prop_thres_nets = list()
    for(j in 1:length(prop_thres)){
      
      temp = thres_prop(net = nets[[i]], prop = prop_thres[j])
      prop_thres_nets[[j]] = graph.adjacency(temp, mode = "undirected",  diag = F)
      
    } 
    
    perc_thres_net = thres_perc(net = nets[[i]], weighted = F)
    perc_thres_net =  graph.adjacency(perc_thres_net, mode = "undirected", diag = F)
    
    weighted_mod = efficiency(weighted_net, type = "global")
    perc_thres_mod = efficiency(perc_thres_net, type = "global")
    abs_thres_mod = lapply(abs_thres_nets, FUN = function(x){return(efficiency(x, type = "global"))})
    prop_thres_mod = lapply(prop_thres_nets, FUN = function(x){return(efficiency(x, type = "global"))})
    
    abs_labs = paste("abs_thres_", abs_thres,"_mod",  sep = "")
    prop_labs = paste("prop_thres_", prop_thres, "_mod",  sep = "")
    
    abs_df = data.frame(thres = abs_thres,thres_ge = as.numeric(abs_thres_mod),thres_type = "Abs")
    
    prop_df = data.frame(thres = prop_thres, thres_ge = as.numeric(prop_thres_mod), thres_type = "Prop")
    
    
    toReturn = data.frame(ID = i, weighted_ge = weighted_mod, perc_thres_ge = perc_thres_mod)
    toReturn = rbind(cbind(toReturn, abs_df),cbind(toReturn, prop_df))
    toReturn_list[[i]] = toReturn
    pb$tick()
  } 
  toReturn = do.call("rbind", toReturn_list)
  return(toReturn)   
}
bassett_lattice = function(net){
  vals = net[upper.tri(net)]
  
  vals = sort(vals, decreasing = T)
  indices = list()
  counter = 1
  dist = vector()
  mat = matrix(0, dim(net)[1], dim(net)[2])
  indices = list()
  
  
  for(d in 1:(ceiling((dim(net)[1]-1)/2))){
    
    for(i in 1:(dim(net)[1]-1)){ 
      dist = c(dist,rep(d,2))
      target_indices =  ((i +c( -d,d)) %% dim(net)[1])
      indices[[counter]] = cbind(i, target_indices)
      counter = counter + 1
    }
  }
  
  indices = do.call("rbind", indices)+1
  indices = cbind(indices,dist)
  indices = indices[which(indices[,1]> indices[,2]),]
  indices = unique(indices)
  
  
  for(d in 1:(ceiling((dim(net)[1]-1)/2))){
    indices_c = which(indices[,3] == d)
    pindices_c = sample(indices_c)
    indices[indices_c,] = indices[pindices_c,]
  }
  
  
  counter = 1
  while(length(vals) >0){
    mat[indices[counter,1], indices[counter,2]] = vals[1]
    vals = vals[-1]
    counter = counter+1
    
    
  }
  return(mat+t(mat))
}
unweighted_cc = function(net){
  
  return(mean(transitivity(net, type = "local", isolates = "zero")))
  
}
swp_unweighted = function(net){
  
  t_cc = mean(transitivity(net, type = "local", isolates = "zero"))
  t_cpl = mean_distance(net, directed = F, unconnected = F)
  
  lattices = lapply(1:10,FUN = function(x){
    
    return(graph.adjacency(bassett_lattice(as.matrix(as_adjacency_matrix(net))), mode = "undirected", diag = F))
    
  } )
  
  rands = lapply(1:10,FUN = function(x){
    
    return(rewire(net,each_edge(p = 1, loops = FALSE)))
    
  } )
  
  lat_cc = mean(sapply(lattices, unweighted_cc))
  lat_cpl = mean(sapply(lattices, mean_distance, directed = F, unconnected = F))
  
  rand_cc = mean(sapply(rands, unweighted_cc))
  rand_cpl = mean(sapply(rands, mean_distance, directed = F, unconnected = F))
  
  dC_true  =(lat_cc - t_cc)/(lat_cc - rand_cc)
  dL_true =(t_cpl - rand_cpl)/(lat_cpl-rand_cpl)
  
  dC = pmax(pmin(dC_true, 1), 0)
  dL = pmax(pmin(dL_true, 1), 0)
  swp = 1-sqrt((dC^2 + dL^2)/2)
  if(all(c(!is.na(dC), !is.na(dL)))){
  if(dC == 1 & dL == 1){
    swp = 2
  }
  
  if(dC == 1 & dL == 0){
    swp = 1
  }
  
  if(dC == 0 & dL == 1){
    swp = -1
  }
  
  if(dC == 0 & dL == 0){
    swp = 0
  }
  }
  return(swp)
}
w_clust= function(net){
  net_adj = as_adj(net, type = "both", attr = "weight", sparse = F)
  return(ClustF(abs(net_adj))$GlobalCC)
}
swp_weighted = function(net){
  net_mat = as.matrix(as_adjacency_matrix(graph = net, attr = "weight"))
  net_mat[which(abs(net_mat) < .001)] = 0
  net = graph.adjacency(abs(net_mat), mode = "undirected", weighted = T, diag = F)
  t_cc = w_clust(net)
  inv_net = abs(1/net_mat)
  inv_net[is.infinite(inv_net)] = 0
  
  inv_net = graph.adjacency(inv_net, mode = "undirected", weighted = T, diag = F)
  
  t_cpl = mean_distance_wt(inv_net)
  
  lattices = lapply(1:10,FUN = function(x){
    
    return(graph.adjacency(bassett_lattice(as.matrix(as_adjacency_matrix(net,attr = "weight"))), mode = "undirected", diag = F, weighted = T))
    
  } )
  
  rands = lapply(1:10,FUN = function(x){
    
    return(rewire(net,each_edge(p = 1, loops = FALSE)))
    
  } )
  
  lat_cc = mean(sapply(lattices, w_clust))
  lat_cpl = mean(sapply(lattices, mean_distance_wt, xfm = T, xfm.type = "1/w"))
  
  rand_cc = mean(sapply(rands, w_clust))
  rand_cpl = mean(sapply(rands, mean_distance_wt, xfm = T, xfm.type = "1/w"))
  dC  =(lat_cc - t_cc)/(lat_cc - rand_cc)
  dL =(t_cpl - rand_cpl)/(lat_cpl-rand_cpl)
  dC = pmax(pmin(dC, 1), 0)
  dL = pmax(pmin(dL, 1), 0)
  
  swp = 1-sqrt((dC^2 + dL^2)/2)
  if(all(c(!is.na(dC), !is.na(dL)))){
    if(dC == 1 & dL == 1){
      swp = 2
    }
    
    if(dC == 1 & dL == 0){
      swp = 1
    }
    
    if(dC == 0 & dL == 1){
      swp = -1
    }
    
    if(dC == 0 & dL == 0){
      swp = 0
    }
  }
  return(swp)
}
apply_SWP = function(nets, abs_thres = seq(.1, .7, .05), prop_thres = seq(.01, .5, .05)){
  toReturn_list = list() 
  pb = progress_bar$new(total = length(nets), format = "[:bar] :current/:total :eta")
  pb$tick(0)
  
  for(i in 1:length(nets)){
    
    pos_thres = thres_pos(nets[[i]], weighted = T)
    temp = 1/pos_thres
    temp[which(is.infinite(temp))] = 0
    weighted_net = graph.adjacency(temp, mode = "undirected", weighted = T, diag = F)
    
    
    abs_thres_nets = list()
    for(j in 1:length(abs_thres)){
      
      temp = thres_abs(net = nets[[i]], thres = abs_thres[j], weighted = F)
      abs_thres_nets[[j]] = graph.adjacency(temp, mode = "undirected",  diag = F)
      
    } 
    
    
    prop_thres_nets = list()
    for(j in 1:length(prop_thres)){
      
      temp = thres_prop(net = nets[[i]], prop = prop_thres[j])
      prop_thres_nets[[j]] = graph.adjacency(temp, mode = "undirected",  diag = F)
      
    } 
    
    perc_thres_net = thres_perc(net = nets[[i]], weighted = F)
    perc_thres_net =  graph.adjacency(perc_thres_net, mode = "undirected", diag = F)
    
    weighted_mod = swp_weighted(weighted_net)
    perc_thres_mod = swp_unweighted(perc_thres_net)
    
    abs_thres_mod = lapply(abs_thres_nets, FUN = function(x){return(swp_unweighted(x))})
    prop_thres_mod = lapply(prop_thres_nets, FUN = function(x){return(swp_unweighted(x))})
    
    abs_labs = paste("abs_thres_", abs_thres,"_mod",  sep = "")
    prop_labs = paste("prop_thres_", prop_thres, "_mod",  sep = "")
    
    abs_df = data.frame(thres = abs_thres,thres_swp = as.numeric(abs_thres_mod),thres_type = "Abs")
    
    prop_df = data.frame(thres = prop_thres, thres_swp = as.numeric(prop_thres_mod), thres_type = "Prop")
    
    
    toReturn = data.frame(ID = i, weighted_swp = weighted_mod, perc_thres_swp = perc_thres_mod)
    toReturn = rbind(cbind(toReturn, abs_df),cbind(toReturn, prop_df))
    toReturn_list[[i]] = toReturn
    pb$tick()
  }
  return(do.call("rbind", toReturn_list))
}

simulate_mean = function(real_mean){
  mu_mean = mean(real_mean,na.rm = TRUE)
  mu_s2 = var(real_mean, na.rm = TRUE)
  scale = mu_s2/mu_mean
  shape = mu_mean^2/mu_s2
  simu_mean = rgamma(length(real_mean), scale = scale, shape = shape)
  simu_mean[simu_mean < log10(1.01)] = log10(1.01) 
  return(simu_mean)
}


simulate_sd = function(real_mean, real_sd, simu_mean, nbin = 15){
  mean_bin = seq(min(real_mean), max(real_mean), length.out = nbin)
  mean_bin[1] = log10(1.01)
  mean_bin[nbin] = Inf
  real_bin_assign = findInterval(real_mean, mean_bin)
  real_bin_unq = unique(real_bin_assign)
  
  simu_bin_assign = findInterval(simu_mean, mean_bin)
  simu_sd = rep(0, length(simu_mean))
  for(b in 1:nbin){
    if(sum(simu_bin_assign == b) == 0) next
    if(b %in% real_bin_unq){
      pop = real_sd[real_bin_assign == b]
    }else{
      ind = which.min(abs(real_bin_unq - b))
      pop = real_sd[real_bin_assign == real_bin_unq[ind]]
    }
    samp = sample(pop, size = sum(simu_bin_assign == b), replace = TRUE)
    simu_sd[simu_bin_assign == b] = samp
  }
  return(simu_sd)
}


simulate_sj = function(real_sj, Js){
  mean = mean(real_sj)
  sd = sd(real_sj)
  simu_sj = rnorm(Js, mean = mean, sd = sd)
  simu_sj[simu_sj < min(real_sj)] = min(real_sj)
  simu_sj[simu_sj > max(simu_sj)] = max(simu_sj)
  return(simu_sj)
}


simulate_dropout = function(gene_drate, cell_drate){
  J = length(cell_drate)

  cell_dnum = rbinom(length(gene_drate), size = J, prob = gene_drate)
  cell_drate_norm = cell_drate/sum(cell_drate)
  
  drop_inds = t(sapply(1:length(gene_drate), function(ii){
    x = rep(0, J)
    x[sample(1:J, cell_dnum[ii], prob = cell_drate_norm)] = 1
    return(x)
  }))
  return(drop_inds)
}


simulate_newdropout = function(real_mean, gene_drate, simu_mean, cell_drate, nbin = 15){
  J = length(cell_drate)
  mean_bin = seq(min(real_mean), max(real_mean), length.out = nbin)
  mean_bin[1] = log10(1.01)
  mean_bin[nbin] = Inf
  real_bin_assign = findInterval(real_mean, mean_bin)
  real_bin_unq = unique(real_bin_assign)
  
  simu_bin_assign = findInterval(simu_mean, mean_bin)
  simu_drate = rep(0, length(simu_mean))
  for(b in 1:nbin){
    if(sum(simu_bin_assign == b) == 0) next
    if(b %in% real_bin_unq){
      pop = gene_drate[real_bin_assign == b]
    }else{
      ind = which.min(abs(real_bin_unq - b))
      pop = gene_drate[real_bin_assign == real_bin_unq[ind]]
    }
    samp = sample(pop, size = sum(simu_bin_assign == b), replace = TRUE)
    simu_drate[simu_bin_assign == b] = samp
  }
  cell_dnum = rbinom(length(simu_drate), size = J, prob = simu_drate)
  cell_drate_norm = cell_drate/sum(cell_drate)
  
  drop_inds = t(sapply(1:length(simu_mean), function(ii){
    x = rep(0, J)
    x[sample(1:J, cell_dnum[ii], prob = cell_drate_norm)] = 1
    return(x)
  }))
  return(drop_inds)
}


simulate_gene_drate = function(real_mean, gene_drate, simu_mean_list, nbin = 15){
  simu_drate_list = list()
  
  mean_bin = seq(min(real_mean), max(real_mean), length.out = nbin)
  mean_bin[1] = log10(1.01)
  mean_bin[nbin] = Inf
  real_bin_assign = findInterval(real_mean, mean_bin)
  real_bin_unq = unique(real_bin_assign)
  
  ### simulate gene_drate1
  csimu_mean = simu_mean_list[[1]]
  simu_drate1 = rep(0, length(csimu_mean))
  
  simu_bin_assign = findInterval(csimu_mean, mean_bin)
  for(b in 1:nbin){
    if(sum(simu_bin_assign == b) == 0) next
    if(b %in% real_bin_unq){
      pop = gene_drate[real_bin_assign == b]
    }else{
      ind = which.min(abs(real_bin_unq - b))
      pop = gene_drate[real_bin_assign == real_bin_unq[ind]]
    }
    samp = sample(pop, size = sum(simu_bin_assign == b), replace = TRUE)
    simu_drate1[simu_bin_assign == b] = samp
  }
  simu_drate_list[[1]] = simu_drate1
  csimu_drate = simu_drate1
  
  for(g in 2:length(simu_mean_list)){
    ### simulate gene_drateg
    simu_mean2 = simu_mean_list[[g]]
    simu_drate2 = rep(0, length(simu_mean2))
    ind_same = which(csimu_mean == simu_mean2)
    simu_drate2[ind_same] = csimu_drate[ind_same]
    ind_diff = which(csimu_mean != simu_mean2)
    simu_tp = rep(0, length(ind_diff))
    if(length(ind_diff) > 0){
      simu_bin_assign = findInterval(simu_mean2[ind_diff], mean_bin)
      for(b in 1:nbin){
        if(sum(simu_bin_assign == b) == 0) next
        if(b %in% real_bin_unq){
          pop = gene_drate[real_bin_assign == b]
        }else{
          ind = which.min(abs(real_bin_unq - b))
          pop = gene_drate[real_bin_assign == real_bin_unq[ind]]
        }
        samp = sample(pop, size = sum(simu_bin_assign == b), replace = TRUE)
        simu_tp[simu_bin_assign == b] = samp
      }
      simu_drate2[ind_diff] = simu_tp
    }
    simu_drate_list[[g]] = simu_drate2
    csimu_drate = simu_drate2
    csimu_mean = simu_mean2
  }

  return(simu_drate_list)
}
simulate_dedropout = function(real_mean, gene_drate, simu_drate, cell_drate){
  J = length(cell_drate)
 
  cell_dnum = rbinom(length(simu_drate), size = J, prob = simu_drate)
  cell_drate_norm = cell_drate/sum(cell_drate)
  
  drop_inds = t(sapply(1:length(simu_drate), function(ii){
    x = rep(0, J)
    x[sample(1:J, cell_dnum[ii], prob = cell_drate_norm)] = 1
    return(x)
  }))
  return(drop_inds)
}



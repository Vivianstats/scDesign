estimate_pa = function(rcount, ncores = 10){
  rcount = as.matrix(rcount)
  colsums = colSums(rcount)
  rcount = sweep(rcount, 2, median(colsums)/colsums, FUN = "*")

  print("estimate expression parameters")
  lcount = log10(rcount+1.01)
  pa = get_mix_parameters(count = lcount, ncores = ncores)
  rownames(pa) = rownames(rcount)

  ngenes = nrow(rcount)
  ncells = ncol(rcount)

  ### calculate cell-wise dropout rate
  droprate = t(sapply(1:ngenes, function(i) {
    wt = calculate_weight(lcount[i, ], pa[i, ])
    return(wt[, 1])
  }))
  cell_drate = colSums(droprate > 0.99, na.rm = TRUE)/ngenes
  maxlog = max(lcount)
  return(list(pa = pa, ngenes = ngenes, ncells = ncells,
              cell_drate = cell_drate, maxlog = maxlog, colsums = colsums))
}

# single count matrix; keep real genes
simulate_ofo = function(estpa, Js, S = NULL){
  ngenes = estpa$ngenes
  cell_drate = estpa$cell_drate
  maxlog = estpa$maxlog
  colsums = estpa$colsums

  pa = estpa$pa
  genenames = rownames(pa)
  pa_rowind = which(complete.cases(pa))
  pa = pa[pa_rowind, , drop = FALSE]
  real_mean = pa[,"mu"]
  real_sd = pa[,"sigma"]
  gene_drate = pa[,"rate"]

  simucount = matrix(0, ncol = Js, nrow = ngenes)
  simu_mean = real_mean
  simu_sd = real_sd

  ### simulate ideal counts
  simu_ideal = t(sapply(1:nrow(pa), function(ii){
    x = rnorm(Js, mean = simu_mean[ii], sd = simu_sd[ii])
    v = min(x[x>log10(1.01)])
    x[x<v] = v
    return(x)
  }))
  simu_ideal[simu_ideal > maxlog] = maxlog

  ### print("introduce dropouts")
  cell_drate = sample(cell_drate, Js, replace = TRUE)
  drop_inds = simulate_dropout(gene_drate = gene_drate, cell_drate = cell_drate)
  simu_ideal[drop_inds == 1] = log10(1.01)

  ### translate back to count data
  simu_trans = 10^simu_ideal - 1.01
  ### adjust for library size
  simu_sj = simulate_sj(real_sj = colsums, Js)

  simu_trans = sweep(simu_trans, 2, simu_sj/median(colsums), FUN = "*")
  ### convert count to proportions
  simu_prop = simu_trans/sum(simu_trans)
  if(is.null(S)){
    simu_multi = rmultinom(1, size = sum(colsums), prob = simu_prop)
  }else{
    simu_multi = rmultinom(1, size = S, prob = simu_prop)
  }
  simu_multi = matrix(simu_multi, ncol = Js)
  simucount[pa_rowind, ] = simu_multi
  rownames(simucount) = genenames
  colnames(simucount) = paste0("cell", 1:Js)
  return(simucount)
}

# single count matrix; new genes
simulate_new_ofo = function(rcount, estpa, Js, S = NULL){
  ngenes = estpa$ngenes
  cell_drate = estpa$cell_drate
  maxlog = estpa$maxlog
  colsums = estpa$colsums

  pa = estpa$pa
  genenames = rownames(pa)
  pa_rowind = which(complete.cases(pa))
  rcount = rcount[pa_rowind, , drop = FALSE]
  pa = pa[pa_rowind, , drop = FALSE]
  real_mean = pa[,"mu"]
  real_sd = pa[,"sigma"]
  gene_drate = pa[,"rate"]

  simucount = matrix(0, ncol = Js, nrow = ngenes)
  ### draw gene means from Gamma distribution
  simu_mean = simulate_mean(real_mean)
  ### draw gene sds from empirical bins
  simu_sd = simulate_sd(real_mean, real_sd, simu_mean)

  ### simulate ideal counts
  simu_ideal = t(sapply(1:nrow(pa), function(ii){
    x = rnorm(Js, mean = simu_mean[ii], sd = simu_sd[ii])
    v = min(x[x>log10(1.01)])
    x[x<v] = v
    return(x)
  }))
  # simu_ideal[simu_ideal < log10(1.01)] = log10(1.01)
  simu_ideal[simu_ideal > maxlog] = maxlog

  ### print("introduce dropouts")
  cell_drate = sample(cell_drate, Js, replace = TRUE)
  drop_inds = simulate_newdropout(real_mean, gene_drate, simu_mean, cell_drate)
  simu_ideal[drop_inds == 1] = log10(1.01)

  ### translate back to count data
  simu_trans = 10^simu_ideal - 1.01
  ### adjust for library size
  simu_sj = simulate_sj(real_sj = colsums, Js)

  simu_trans = sweep(simu_trans, 2, simu_sj/median(colsums), FUN = "*")
  ### convert count to proportions
  simu_prop = simu_trans/sum(simu_trans)
  ### print("draw multinom reads")
  if(is.null(S)){
    simu_multi = rmultinom(1, size = sum(colsums), prob = simu_prop)
  }else{
    simu_multi = rmultinom(1, size = S, prob = simu_prop)
  }
  simu_multi = matrix(simu_multi, ncol = Js)
  simucount[pa_rowind, ] = simu_multi
  rownames(simucount) = paste0("gene", 1:ngenes)
  colnames(simucount) = paste0("cell", 1:Js)
  return(simucount)
}

# multiple count matrices; new de genes
simulate_de_mfo = function(rcount, estpa, Js, ngroup = 2, S,
                           pUp = 0.05, pDown = 0.05,
                           fU = 5, fL = 3){
  if(ngroup < 2){stop("number of cell groups < 2!")}
  ngenes = estpa$ngenes
  cell_drate = estpa$cell_drate
  maxlog = estpa$maxlog
  colsums = estpa$colsums

  pa = estpa$pa
  pa_rowind = which(complete.cases(pa))
  rcount = rcount[pa_rowind, , drop = FALSE]
  pa = pa[pa_rowind, , drop = FALSE]
  real_mean = pa[,"mu"]
  real_sd = pa[,"sigma"]
  gene_drate = pa[,"rate"]

  ### cell condition 1
  simu_mean1 = simulate_mean(real_mean)

  parlist = list()
  cmean = simu_mean1
  for(g in 1:(ngroup-1)){
    ### cell condition g+1
    simu_mean = cmean
    nup = round(ngenes * pUp)
    ndown = round(ngenes * pDown)
    degenes = sample(1:length(cmean), nup+ndown, replace = FALSE)
    upgenes = sample(degenes, nup, replace = FALSE)
    downgenes = setdiff(degenes, upgenes)
    upfds = log10(runif(nup, fL, fU))
    downfds = -log10(runif(ndown, fL, fU))
    simu_mean[upgenes] = simu_mean[upgenes] + upfds
    simu_mean[downgenes] = simu_mean[downgenes] + downfds
    simu_mean[simu_mean >= max(simu_mean1)*1.1] = max(simu_mean1)
    simu_mean[simu_mean <= log10(1.01)] = 0.1
    plist = list(simu_mean = simu_mean, upgenes = upgenes, downgenes = downgenes)
    parlist[[g]] = plist
    cmean = simu_mean
  }

  simu_mean_list = lapply(1:ngroup, function(g){
    if(g == 1) return(simu_mean1)
    return(parlist[[g-1]]$simu_mean)
  })
  gene_drate_list = simulate_gene_drate(real_mean, gene_drate, simu_mean_list)

  res = lapply(1:ngroup, function(cc){
    if(cc == 1){simu_mean = simu_mean1}
    if(cc > 1){simu_mean = parlist[[cc-1]]$simu_mean}
    simu_sd = simulate_sd(real_mean, real_sd, simu_mean) #, bw = bw_sd)
    simucount = matrix(0, ncol = Js[cc], nrow = ngenes)
    ### simulate ideal counts
    simu_ideal = t(sapply(1:nrow(pa), function(ii){
      x = rnorm(Js[cc], mean = simu_mean[ii], sd = simu_sd[ii])
      v = min(x[x>log10(1.01)])
      x[x<v] = v
      return(x)
    }))
    # simu_ideal[simu_ideal < log10(1.01)] = log10(1.01)
    simu_ideal[simu_ideal > maxlog] = maxlog

    ### print("introduce dropouts")
    sim_cell_drate = sample(cell_drate, Js[cc], replace = TRUE)

    drop_inds = simulate_dedropout(real_mean, gene_drate,
                                   simu_drate = gene_drate_list[[cc]],
                                   cell_drate = sim_cell_drate)
    simu_ideal[drop_inds == 1] = log10(1.01)

    ### translate back to count data
    simu_trans = 10^simu_ideal - 1.01
    ### adjust for library size
    simu_sj = simulate_sj(real_sj = colsums, Js[cc])

    simu_trans = sweep(simu_trans, 2, simu_sj/median(colsums), FUN = "*")
    ### convert count to proportions
    simu_prop = simu_trans/sum(simu_trans)
    ### print("draw multinom reads")
    if(is.null(S)){
      simu_multi = rmultinom(1, size = sum(colsums), prob = simu_prop)
    }else{
      simu_multi = rmultinom(1, size = S[cc], prob = simu_prop)
    }
    simu_multi = matrix(simu_multi, ncol = Js[cc])
    simucount[pa_rowind, ] = simu_multi
    rownames(simucount) = paste0("gene", 1:ngenes)
    colnames(simucount) = paste0("C", cc, "_", 1:(Js[cc]))
    return(simucount)
  })
  names(res) = paste0("count", 1:ngroup)
  genesUp = lapply(1:ngroup, function(g){
    if(g==1) return(NULL)
    paste0("gene", pa_rowind[parlist[[g-1]]$upgenes])
  })

  genesDown = lapply(1:ngroup, function(g){
    if(g==1) return(NULL)
    paste0("gene", pa_rowind[parlist[[g-1]]$downgenes])
  })

  return(list(count = res, genesUp = genesUp, genesDown = genesDown))
}

# two count matrices; real de genes
simulate_mfm = function(estpa1, estpa2, Js, S = NULL, p1 = 0.3, p2 = 0.3){
  pa1 = estpa1$pa
  pa2 = estpa2$pa
  if(nrow(pa1) != nrow(pa2)) stop("two real datasets have different gene numbers!")
  if(sum(rownames(pa1) == rownames(pa2)) != nrow(pa2)) stop("gene names in the two datasets don't match!")

  cn = rmultinom(1, Js, c(p1, p2, 1-p1-p2))
  J1 = cn[1]; J2 = cn[2]
  Ssub = round(S*(J1+J2)/Js)
  res = lapply(1:2, function(cc){
    if(cc == 1){
      estpa = estpa1; JJ = J1
    }else{estpa = estpa2; JJ = J2}

    ngenes = estpa$ngenes
    cell_drate = estpa$cell_drate
    maxlog = estpa$maxlog
    colsums = estpa$colsums
    pa = estpa$pa
    pa_rowind = which(complete.cases(pa))
    pa = pa[pa_rowind, , drop = FALSE]
    real_mean = pa[,"mu"]
    real_sd = pa[,"sigma"]
    gene_drate = pa[,"rate"]

    simu_mean = real_mean
    simu_sd = real_sd
    simucount = matrix(0, ncol = JJ, nrow = ngenes)
    ### simulate ideal counts
    simu_ideal = t(sapply(1:nrow(pa), function(ii){
      x = rnorm(JJ, mean = simu_mean[ii], sd = simu_sd[ii])
      v = min(x[x>log10(1.01)])
      x[x<v] = v
      return(x)
    }))
    simu_ideal[simu_ideal > maxlog] = maxlog

    ### print("introduce dropouts")
    cell_drate = sample(cell_drate, JJ, replace = TRUE)
    drop_inds = simulate_newdropout(real_mean, gene_drate,
                                    simu_mean, cell_drate)
    simu_ideal[drop_inds == 1] = log10(1.01)

    ### translate back to count data
    simu_trans = 10^simu_ideal - 1.01
    ### adjust for library size
    simu_sj = simulate_sj(real_sj = colsums, JJ)

    simu_trans = sweep(simu_trans, 2, simu_sj/median(colsums), FUN = "*")
    simures = matrix(0, ncol = ncol(simu_trans), nrow = ngenes)
    simures[pa_rowind, ] = simu_trans
    if(cc == 1){colnames(simures) = paste0("C1_", 1:JJ)}
    if(cc == 2){colnames(simures) = paste0("C2_", 1:JJ)}
    return(simures)
  })

  simu_trans = matrix(0, nrow = nrow(pa1), ncol = J1+J2)
  simu_trans[, 1:J1] = res[[1]]
  simu_trans[, (J1+1):(J1+J2)] = res[[2]]
  colnames(simu_trans) = c(colnames(res[[1]]), colnames(res[[2]]))
  ### convert count to proportions
  simu_prop = simu_trans/sum(simu_trans)
  ### print("draw multinom reads")
  simu_multi = rmultinom(1, size = Ssub, prob = simu_prop)
  simu_multi = matrix(simu_multi, ncol = J1+J2)
  rownames(simu_multi) = rownames(pa1)
  colnames(simu_multi) = colnames(simu_trans)

  return(simu_multi)
}





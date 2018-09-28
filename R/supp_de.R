de_scimpute = function(pa1, pa2, ncores, by = "factor", factor = 1, rank = 500){
  pa1 = pa1[complete.cases(pa1), , drop = FALSE]
  pa2 = pa2[complete.cases(pa2), , drop = FALSE]
  genesint = intersect(rownames(pa1), rownames(pa2))
  pa1 = pa1[genesint, , drop = FALSE]
  pa2 = pa2[genesint, , drop = FALSE]

  mean1 = pa1[,"mu"]
  sd1 = pa1[,"sigma"]
  mean2 = pa2[,"mu"]
  sd2 = pa2[,"sigma"]
  facs = abs(mean1 - mean2)/(sd1 + sd2)
  if(by == "factor"){degenes = genesint[facs >= factor]}
  if(by == "rank"){degenes = genesint[order(facs, decreasing = TRUE)[1:rank]]}
  return(degenes)
}



# MAST
# multiple fdr
# return a  matrix,
# rows:fdr values, 1st column:precision, 2nd column:recall
de_mast_mfdr = function(count1, count2, by = "fdr", fdr_thres = 10^seq(-2,-6,-1),
                   logfc_thre = log2(1.2), truedegenes){
  freq_expressed = 0.1
  genesint = intersect(rownames(count1), rownames(count2))
  count1 = count1[genesint, , drop = FALSE]
  count2 = count2[genesint, , drop = FALSE]
  count = cbind(count1, count2)
  labels = c(rep("C1", ncol(count1)),rep("C2", ncol(count2)))

  cdat = data.frame(wellKey = colnames(count), condition = labels)
  fdat = data.frame(primerid = rownames(count), row.names = rownames(count))
  scRaw = FromMatrix(count, cdat, fdat)
  sca = scRaw[which(freq(scRaw) > 0), ]
  cdr2 = colSums(assay(sca) > 0)
  colData(sca)$cngeneson = scale(cdr2)
  thres = thresholdSCRNACountMatrix(assay(sca), nbins = 20, min_per_bin = 30, data_log=FALSE)
  assays(sca) = list(thresh=thres$counts_threshold, cpm=assay(sca))
  expressed_genes = freq(sca) > freq_expressed
  sca = sca[expressed_genes, ]

  zlmCond = zlm(~condition + cngeneson, sca)
  summaryCond = summary(zlmCond, doLRT='conditionC2')
  summaryDt = summaryCond$datatable
  fcHurdle = merge(summaryDt[contrast=='conditionC2' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                    summaryDt[contrast=='conditionC2' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  if(by == "fdr"){
    mat = t(sapply(fdr_thres, function(fdr_thre){
      fcHurdleSig = merge(fcHurdle[fdr<fdr_thre & abs(coef)>logfc_thre],
                          data.table::as.data.table(mcols(sca)), by='primerid')
      data.table::setorder(fcHurdleSig, fdr)
      myde = fcHurdleSig$primerid
      tp = intersect(truedegenes, myde)
      prec = length(tp)/length(myde)
      if(length(myde) == 0) prec = 1
      reca = length(tp)/length(truedegenes)
      # 1 - false positive rate
      truenegative = setdiff(genesint, truedegenes)
      fp = 1 - length(intersect(myde, truenegative))/length(truenegative)
      # if(fdr_thre == 1e-2){print(fp)}
      F_pr = 2*prec*reca/(prec + reca)
      F_roc = 2*fp*reca/(fp + reca)
      return(c(prec, reca, fp, F_pr, F_roc))
    }))
  }
  return(mat)
}


# t test
# multiple fdr
# return a  matrix,
# rows:fdr values, 1st column:precision, 2nd column:recall
de_ttest_mfdr = function(count1, count2, by = "fdr", fdr_thres = 10^seq(-2,-6,-1),
                        truedegenes){
  colsums1 = colSums(count1)
  count1 = sweep(count1, 2, 1e6/colsums1, FUN = "*")
  colsums2 = colSums(count2)
  count2 = sweep(count2, 2, 1e6/colsums2, FUN = "*")

  nexpressed1 = max(5, round(0.05*ncol(count1)))
  nexpressed2 = max(5, round(0.05*ncol(count2)))
  expressed1 = rowSums(count1 > 0)
  count1 = count1[expressed1 > nexpressed1, ]
  expressed2 = rowSums(count2 > 0)
  count2 = count2[expressed2 > nexpressed2, ]

  genesint = intersect(rownames(count1), rownames(count2))
  count1 = count1[genesint, , drop = FALSE]
  count2 = count2[genesint, , drop = FALSE]

  pvals = sapply(1:length(genesint), function(i){
    x = count1[i, ]; x = log10(x[x>0]+1)
    y = count2[i, ]; y = log10(y[y>0]+1)
    tobs = try(t.test(x = x, y = y), silent = TRUE)
    if(class(tobs) == "try-error"){pval = 1
    }else{ pval = tobs$p.value}
    return(pval)
  })
  pvals.fdr = stats::p.adjust(pvals, "fdr")
  if(by == "fdr"){
    mat = t(sapply(fdr_thres, function(fdr_thre){
      # print(fdr_thre)
      myde = genesint[pvals.fdr < fdr_thre]
      tp = intersect(truedegenes, myde)
      prec = length(tp)/length(myde)
      if(length(myde) == 0) prec = 1
      reca = length(tp)/length(truedegenes)
      # 1 - false positive rate
      truenegative = setdiff(genesint, truedegenes)
      fp = 1 - length(intersect(myde, truenegative))/length(truenegative)
      # if(fdr_thre == 1e-2){print(fp)}
      F_pr = 2*prec*reca/(prec + reca)
      F_roc = 2*fp*reca/(fp + reca)
      return(c(prec, reca, fp, F_pr, F_roc))
    }))
  }
  return(mat)
}


get_pr_mfdr = function(estpa1, estpa2, B, ncell, S1, S2,
                       de_method = "mast", degenes, by = "fdr",
                       fdr_thres = 10^seq(-2,-6,-1), ncores = 20){

  prls = lapply(1:nrow(ncell), function(ii){
    J1 = ncell[ii,1]; J2 = ncell[ii,2]
    message(paste("state 1 cell number:", J1))
    message(paste("state 2 cell number:", J2))

    prreps = mclapply(1:B, function(bb){
      set.seed(bb)
      if(bb %% 10 == 0) gc()
      count1 = simulate_ofo(estpa1, J1, S1)
      count2 = simulate_ofo(estpa2, J2, S2)
      colnames(count1) = paste0(colnames(count1), "_c1")
      colnames(count2) = paste0(colnames(count2), "_c2")
      if(de_method == "mast"){
        # print("check")
        de_simu = try(de_mast_mfdr(count1, count2, by = by,
                                   fdr_thres = fdr_thres, truedegenes = degenes),
                      silent = TRUE)
        if(class(de_simu) == "try-error") {
          de_simu = matrix(NA, nrow = length(fdr_thres), ncol = 5)
        }
      }
      if(de_method == "ttest"){
        # print("check")
        de_simu = try(de_ttest_mfdr(count1, count2, by = by, fdr_thres = fdr_thres,
                                    truedegenes = degenes),
                      silent = TRUE)
        if(class(de_simu) == "try-error") {
          de_simu = matrix(NA, nrow = length(fdr_thres), ncol = 5)
        }
      }
      return(de_simu)
    }, mc.cores = ncores)
    prreps = prreps[!sapply(prreps, is.null)]
    pr = try(apply(simplify2array(prreps), 1:2, mean, na.rm = TRUE), silent = TRUE)
    if(class(pr) == "try-error"){
      message(pr)
      pr = matrix(NA, nrow = length(fdr_thres), ncol = 5)
    }
    return(pr)
  })

  res = lapply(1:5, function(k){
    # rows for pval, columns for cell number
    mat = sapply(1:nrow(ncell), function(s1){ prls[[s1]][,k]})
    colnames(mat) = paste0(ncell[,1], "vs", ncell[,2])
    rownames(mat) = paste0(fdr_thres)
    return(mat)
  })
  names(res) = c("precision", "recall", "TN", "F1", "F2")
  return(res)
}




# two cell types sequenced together
# multiple fdr thresholds
# fixed proportions
# pmat(rmat), rows: fdr values, columns: cell numbers
get_pr_mfdr_joint = function(estpa1, estpa2, B, Js_vec, prop1, prop2, S,
                              de_method = "mast", degenes, by = "fdr",
                              fdr_thres = 10^seq(-2,-6,-1), ncores = 20){
  prls = lapply(Js_vec, function(Js){
    print(paste("cell number", Js))
    prreps = mclapply(1:B, function(bb){
      set.seed(bb)
      if(bb %% 50 == 0) print(paste("Iteration: ", bb)); gc()
      simucount = simulate_mfm(estpa1, estpa2, Js, S = S, p1 = prop1, p2 = prop2)
      labels = colnames(simucount)
      labels = sapply(strsplit(labels, split = "_"), function(x) x[1])
      count1 = simucount[,labels == "C1"]
      count2 = simucount[,labels == "C2"]

      if(de_method == "mast"){
        # print("check")
        de_simu = try(de_mast_mfdr(count1, count2, by = by,
                                   fdr_thres = fdr_thres, truedegenes = degenes),
                      silent = TRUE)
        if(class(de_simu) == "try-error") {
          de_simu = matrix(NA, nrow = length(fdr_thres), ncol = 5)
        }
      }
      if(de_method == "ttest"){
        # print("check")
        de_simu = try(de_ttest_mfdr(count1, count2, by = by,
                                    truedegenes = degenes),
                      silent = TRUE)
        if(class(de_simu) == "try-error") {
          de_simu = matrix(NA, nrow = length(fdr_thres), ncol = 5)
        }
      }
      return(de_simu)
    }, mc.cores = ncores)
    prreps = prreps[!sapply(prreps, is.null)]
    pr = try(apply(simplify2array(prreps), 1:2, mean, na.rm = TRUE), silent = TRUE)
    if(class(pr) == "try-error"){
      message(pr)
      pr = matrix(NA, nrow = length(fdr_thres), ncol = 5)
    }
    return(pr)
  })
  res = lapply(1:5, function(k){
    # rows for pval, columns for cell number
    mat = sapply(1:length(Js_vec), function(s1){ prls[[s1]][,k]})
    colnames(mat) = as.character(ncell)
    rownames(mat) = as.character(fdr_thres)
    return(mat)
  })
  names(res) = c("precision", "recall", "TN", "F1", "F2")

  return(res)
}





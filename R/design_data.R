#' use scDesign to simulate scRNA-seq data
#' @param realcount A numeric matrix with rows representing genes and columns representing cells. 
#' Gene names are given as row names.
#' @param S A number specifying the total number of RNA-seq reads. Default to 1e8.
#' @param ncell An integer specifying the number of cells. When \code{ngroup > 1},
#' \code{ncell} is the number of cells in each cell state.
#' @param ngroup An integer giving the number of cell states to simulate. Defaults to 1.
#' @param pUp A value between 0 and 1 specifying the proportion of up regulated genes
#' between two adjacent cell states. Defaults to 0.05 and only used when \code{ngroup > 1}.
#' @param pDown A value between 0 and 1 specifying the proportion of down regulated genes
#' between two adjacent cell states. Defaults to 0.05 and only used when \code{ngroup > 1}.
#' @param fU A value specifying the upper bound of fold changes of differentially expressed genes.
#' Deaults to 5.
#' @param fL A value specifying the lower bound of fold changes of differentially expressed genes.
#' Deaults to 1.5.
#' @param ncores An integer specifying the number of cores used for parallel computation.
#' Defaults to 1.
#' @return When \code{ngroup = 1}, \code{design_data} returns a simulated count matrix with
#' rows representing genes and columns representing cells.
#' When \code{ngroup > 1}, \code{design_data} returns a list of \code{ngroup} elements. 
#' The g-th element corresponds to the g-th cell state, and is a list containing three elements:
#' \describe{
#'   \item{count:}{a count matrix with rows representing genes and columns representing cells;}
#'   \item{genesUp:}{a character vector giving the names of up-regulated genes from state g-1 to g;}
#'   \item{genesDown:}{a character vector giving the names of down-regulated genes from state g-1 to g.}
#' }
#' @export
#' @import parallel
#' @importFrom stats complete.cases dgamma dnorm prcomp quantile rgamma rnorm sd uniroot
#' @author Wei Vivian Li, \email{liw@ucla.edu}
#' @author Jingyi Jessica Li, \email{jli@stat.ucla.edu}
design_data = function(realcount, S = 1e8, ncell, ngroup = 1, 
                       pUp = 0.05, pDown = 0.05, fU = 5, fL = 1.5, ncores = 1){
  estpa = estimate_pa(realcount, ncores = ncores)
  if(ngroup == 1){
    simu_res = simulate_new_ofo(realcount, Js = ncell, estpa = estpa, S = S)
  }else{
    simu_res = simulate_de_mfo(realcount, estpa, Js = ncell, 
                               ngroup = ngroup, S = S,
                               pUp = pUp, pDown = pDown, fU = fU, fL = fL)
  }
  return(simu_res)
}


design_plot_sep = function(prlist, ncell, plot_dir, pvals){
  theme_set(theme_bw())
  measures = c("precision", "recall", "TN", "F1", "F2")
  bestnum = sapply(prlist, function(mat){
    nums = sapply(1:nrow(mat), function(i){
      paste0(ncell[which.max(mat[i, ]), ], collapse = "vs")
    })
    num = names(sort(table(nums), decreasing=TRUE)[1])
    return(num)
  })
  measures = c("precision", "recall (TP)", "TN", 
               "F1 (precision vs. recall)", "F2 (TP vs. TN)")
  numdata = data.frame(metric = measures, "cell number" = bestnum)
  
  prmat = lapply(1:5, function(k){
    mat = prlist[[k]]
    data = as.data.frame(mat)
    colnames(data) = paste0(ncell[,1], "vs", ncell[,2])
    data$pval = -log10(pvals)
    data$measure = measures[k]
    data = data %>% gather(ncell, value, -(measure:pval))
    data$ncell = factor(data$ncell, levels = paste0(ncell[,1], "vs", ncell[,2]))
    return(data)
  })
  red = "#CC0C00"; lred = "#EEAEAA"
  blue = "#5C88DA"; lblue = "#C8D7F2"
  yellow = "#FFB200"; lyellow = "#FFE5AA"
  
  ylabs = c("precision", "recall (TP)", "TN", "F1", "F2")
  plots = lapply(1:5, function(k){
    if(k %in% c(1,3)){low = lblue; high = blue}
    if(k == 2){low = lred; high = red}
    if(k %in% c(4,5)){low = lyellow; high = yellow}
    
    g1 = ggplot(prmat[[k]], aes(x = ncell, y = value, group = pval, color = pval)) +
      geom_line()+ scale_color_gradient(low = low, high = high, name = "-log(p)") +
      xlab("") + ylab(ylabs[k]) +
      # labs(color = "-log(threshold)") +
      theme(text = element_text(size=16), axis.text = element_text(size=16),
            legend.title = element_text(size=16),
            legend.text = element_text(size=16),
            legend.key.size = unit(2, 'lines'),
            axis.text.x = element_text(angle = 45, hjust = 1))
    return(g1)
  })
  
  mytheme = gridExtra::ttheme_default(base_size = 20)
  gtable = tableGrob(numdata, rows = NULL, theme = mytheme)
  width = nrow(ncell) * 1.5
  if(nrow(ncell) <= 3) width = 5
  height = 5.5
  if(max(ncell)>=1e4) height = height + 0.5*(log10(max(ncell))-4)
  pdf(paste0(plot_dir, "design_summary.pdf"),
      width = width, height = height)
  grid.draw(plots[[1]])
  grid.draw(plots[[2]])
  grid.draw(plots[[3]])
  grid.draw(plots[[4]])
  grid.draw(plots[[5]])
  grid.arrange(gtable)
  dev.off()
  return(0)
}


#' use scDesign to make experimental design assuming two cell states are sequenced independently
#' @param realcount1 A numeric matrix with rows representing genes and columns representing cells (cell state 1). 
#' Gene names are given as row names.
#' @param realcount2 A numeric matrix with rows representing genes and columns representing cells (cell state 2). 
#' Gene names are given as row names.
#' @param S1 A number specifying the total number of RNA-seq reads for cell state 1. Default to 1e8.
#' @param S2 A number specifying the total number of RNA-seq reads for cell state 2. Default to 1e8.
#' @param ncell A two-column matrix specifying the numbers of cells. Column 1 is for cell state 1 and
#' column 2 is for cell state 2. By default, the following cell number matrix is used:
#' \tabular{rr}{
#' 64 \tab 64 \cr
#' 128 \tab 128 \cr
#' 256 \tab 256 \cr
#' 512 \tab 512 \cr
#' 1024 \tab 1024 \cr
#' 2048 \tab 2048 \cr
#' 4096 \tab 4096 \cr
#' }
#' @param B An integer giving the number of experiments to repeat in order the calculate 
#' the average DE analysis accuracy. Defaults to 100.
#' @param de_method A character specifying the differential expression analysis method to use.
#' Currently supports "ttest" (default) or "mast".
#' @param p_thre A numeric vector specifying the FDR thresholds used to identify 
#' differentially expressed genes. Defaults to \code{10^seq(-2,-6,-1)}.
#' @param plot_dir A character giving the directory to save experimental design results
#' Defaults to "./".
#' @param ncores An integer specifying the number of cores used for parallel computation.
#' Defaults to 1.
#' @param rank An integer specifying the number of top DE genes to identify 
#' from scImpute's results as the standard in DE analysis. Defaults to 1000.
#' @return A list of five elments:
#' \describe{
#'   \item{precision:}{a matrix of precision. }
#'   \item{recall:}{a matrix of recall (true positive rate).}
#'   \item{TN:}{a matrix of TN (true negative rate).}
#'   \item{F1:}{a matrix of F1 (precision vs. recall).}
#'   \item{F2:}{a matrix of F2 (TN vs. recall).}
#' }
#' In all the matrices, rows correspond to different FDR thresholds and 
#' columns correspond to the cell numbers  specified in \code{ncell}.
#' \code{design_sep} also writes the list to design_summary.txt and
#' saves it to \code{plot_dir}. 
#' The corresponding plots are also saved to \code{plot_dir}.
#' @export
#' @import parallel
#' @import grid
#' @import MAST
#' @importFrom SummarizedExperiment assay colData
#' @importFrom GenomicRanges mcols
#' @importFrom gridExtra ttheme_default tableGrob grid.arrange
#' @import ggplot2
#' @importFrom dplyr %>%
#' @importFrom tidyr gather
#' @importFrom utils write.table
#' @importFrom grDevices dev.off pdf
#' @importFrom data.table as.data.table setorder
#' @importFrom stats coef p.adjust rbinom rmultinom runif t.test var
#' @author Wei Vivian Li, \email{liw@ucla.edu}
#' @author Jingyi Jessica Li, \email{jli@stat.ucla.edu}
design_sep = function(realcount1, realcount2, 
                        S1 = 1e8, S2 = 1e8, ncell = NULL, B = 100,
                        de_method = "ttest", p_thre = 10^seq(-2,-6,-1),
                        plot_dir = "./",
                        ncores = 1, rank = 1000){
  if(is.null(ncell)){
    message("Using default cell numbers ... ")
    ncell = cbind(round(2^seq(6,12,1)), round(2^seq(6,12,1)))
  }
  if(class(ncell) != "matrix") stop("ncell must be a matrix!")
  if(ncol(ncell) != 2) stop("ncell must have two columns!")
  if(is.null(rownames(realcount1))) stop("realcount1 must have rownames!")
  if(is.null(rownames(realcount1))) stop("realcount1 must have rownames!")
    
  ncell = data.frame(ncell)
  colnames(ncell) = c("n1", "n2")
  ncell = ncell[with(ncell, order(n1, n2)), , drop = FALSE]
  
  estpa1 = estimate_pa(realcount1, ncores = ncores)
  estpa2 = estimate_pa(realcount2, ncores = ncores)
  
  truede = de_scimpute(estpa1$pa, estpa2$pa, ncores = ncores, by = "rank", rank = rank)
  
  prlist = get_pr_mfdr(estpa1, estpa2, B=B, ncell=ncell, S1=S1, S2=S2, 
                       de_method = de_method, degenes = truede, 
                       fdr_thres = p_thre, ncores = ncores)
  dir.create(plot_dir, recursive = TRUE)
  design_plot_sep(prlist, ncell, plot_dir = plot_dir, pvals = p_thre)
  
  for(i in 1:5){
    cat(names(prlist)[i], file=paste0(plot_dir,"design_summary.txt"), 
        sep="\n", append = TRUE)
    mat = round(prlist[[i]], digits = 3)
    mat1 = data.frame(rownames(mat), mat)
    colnames(mat1) = c("p_thre", colnames(mat))
    write.table(mat1, file=paste0(plot_dir,"design_summary.txt"), 
                row.names=FALSE, col.names=TRUE, append = TRUE, quote = FALSE)
    cat("\n", file=paste0(plot_dir,"design_summary.txt"), append = TRUE)
  }
  return(prlist)
}


design_plot_joint = function(prlist, ncell, plot_dir, pvals){
  theme_set(theme_bw())
  measures = c("precision", "recall", "TN", "F1", "F2")
  bestnum = sapply(prlist, function(mat){
    nums = sapply(1:nrow(mat), function(i){
      as.character(ncell[which.max(mat[i, ])])
    })
    num = names(sort(table(nums), decreasing=TRUE)[1])
    return(num)
  })
  measures = c("precision", "recall (TP)", "TN", 
               "F1 (precision vs. recall)", "F2 (TP vs. TN)")
  numdata = data.frame(metric = measures, "cell number" = bestnum)
  
  prmat = lapply(1:5, function(k){
    mat = prlist[[k]]
    data = as.data.frame(mat)
    colnames(data) = ncell
    data$pval = -log10(pvals)
    data$measure = measures[k]
    data = data %>% gather(ncell, value, -(measure:pval))
    data$ncell = factor(data$ncell, levels = ncell)
    return(data)
  })
  red = "#CC0C00"; lred = "#EEAEAA"
  blue = "#5C88DA"; lblue = "#C8D7F2"
  yellow = "#FFB200"; lyellow = "#FFE5AA"
  
  ylabs = c("precision", "recall (TP)", "TN", "F1", "F2")
  plots = lapply(1:5, function(k){
    if(k %in% c(1,3)){low = lblue; high = blue}
    if(k == 2){low = lred; high = red}
    if(k %in% c(4,5)){low = lyellow; high = yellow}
    
    g1 = ggplot(prmat[[k]], aes(x = ncell, y = value, group = pval, color = pval)) +
      geom_line()+ scale_color_gradient(low = low, high = high, name = "-log(p)") +
      xlab("") + ylab(ylabs[k]) +
      # labs(color = "-log(threshold)") +
      theme(text = element_text(size=16), axis.text = element_text(size=16),
            legend.title = element_text(size=16),
            legend.text = element_text(size=16),
            legend.key.size = unit(2, 'lines'),
            axis.text.x = element_text(angle = 45, hjust = 1))
    return(g1)
  })
  
  mytheme = gridExtra::ttheme_default(base_size = 20)
  gtable = tableGrob(numdata, rows = NULL, theme = mytheme)
  width = length(ncell) * 1.5
  if(length(ncell) <= 3) width = 5
  height = 5.5
  if(max(ncell)>=1e4) height = height + 0.5*(log10(max(ncell))-4)
  pdf(paste0(plot_dir, "design_summary.pdf"),
      width = width, height = height)
  grid.draw(plots[[1]])
  grid.draw(plots[[2]])
  grid.draw(plots[[3]])
  grid.draw(plots[[4]])
  grid.draw(plots[[5]])
  grid.arrange(gtable)
  dev.off()
  return(0)
}



#' use scDesign to make experimental design assuming two cell states are sequenced together
#' @param realcount1 A numeric matrix with rows representing genes and columns representing cells (cell state 1). 
#' Gene names are given as row names.
#' @param realcount2 A numeric matrix with rows representing genes and columns representing cells (cell state 2). 
#' Gene names are given as row names.
#' @param prop1 A number giving the proportion of state 1 cells in the cell population.
#' @param prop2 A number giving the proportion of state 2 cells in the cell population.
#' @param S A number specifying the total number of RNA-seq reads for the cell population. Default to 1e8.
#' @param ncell An integer vector specifying the total number of cells to sequence. 
#' Defaults to \code{2^seq(6,13,1)}.
#' @param B An integer giving the number of experiments to repeat in order the calculate 
#' the average DE analysis accuracy. Defaults to 100.
#' @param de_method A character specifying the differential expression analysis method to use.
#' Currently supports "ttest" (default) or "mast".
#' @param p_thre A numeric vector specifying the FDR thresholds used to identify 
#' differentially expressed genes. Defaults to \code{10^seq(-2,-6,-1)}.
#' @param plot_dir A character giving the directory to save experimental design results
#' Defaults to "./".
#' @param ncores An integer specifying the number of cores used for parallel computation.
#' Defaults to 1.
#' @param rank An integer specifying the number of top DE genes to identify 
#' from scImpute's results as the standard in DE analysis. Defaults to 1000.
#' @return A list of five elments:
#' \describe{
#'   \item{precision:}{a matrix of precision. }
#'   \item{recall:}{a matrix of recall (true positive rate).}
#'   \item{TN:}{a matrix of TN (true negative rate).}
#'   \item{F1:}{a matrix of F1 (precision vs. recall).}
#'   \item{F2:}{a matrix of F2 (TN vs. recall).}
#' }
#' In all the matrices, rows correspond to different FDR thresholds and 
#' columns correspond to the cell numbers  specified in \code{ncell}.
#' \code{design_joint} also writes the list to design_summary.txt and
#' saves it to \code{plot_dir}. 
#' The corresponding plots are also saved to \code{plot_dir}.
#' @export
#' @import parallel
#' @import grid
#' @import MAST
#' @importFrom SummarizedExperiment assay colData
#' @importFrom GenomicRanges mcols
#' @importFrom gridExtra ttheme_default tableGrob grid.arrange
#' @import ggplot2
#' @importFrom dplyr %>%
#' @importFrom tidyr gather
#' @importFrom utils write.table
#' @importFrom grDevices dev.off pdf
#' @importFrom data.table as.data.table setorder
#' @importFrom stats coef p.adjust rbinom rmultinom runif t.test var
#' @author Wei Vivian Li, \email{liw@ucla.edu}
#' @author Jingyi Jessica Li, \email{jli@stat.ucla.edu}
design_joint = function(realcount1, realcount2, prop1, prop2,
                        S = 1e8, ncell = round(2^seq(6,13,1)), B = 100,
                        de_method = "ttest", p_thre = 10^seq(-2,-6,-1),
                        plot_dir = "./",
                        ncores = 1, rank = 1000){

  if(!is.numeric(ncell)) stop("ncell must be an integer vector!")
  if(length(ncell) == 1) stop("ncell must have at least two elements!")
  if(prop1 + prop2 >1) stop("prop1 + prop2 > 1!")
  genes_int = intersect(rownames(realcount1), rownames(realcount2))
  if(length(genes_int) < 1000) stop("two real datasets have less than 1000 common genes!")
  realcount1 = realcount1[genes_int, ]
  realcount2 = realcount2[genes_int, ]
  
  estpa1 = estimate_pa(realcount1, ncores = ncores)
  estpa2 = estimate_pa(realcount2, ncores = ncores)
  
  truede = de_scimpute(estpa1$pa, estpa2$pa, ncores = ncores, by = "rank", rank = rank)
  
  prlist = get_pr_mfdr_joint(estpa1, estpa2, B=B, Js_vec = ncell, 
                             prop1 = prop1, prop2 = prop2, S=S, 
                             de_method = de_method, degenes = truede, by = "fdr", 
                             fdr_thres = p_thre, ncores = ncores)
  dir.create(plot_dir, recursive = TRUE)
  design_plot_joint(prlist, ncell, plot_dir = plot_dir, pvals = p_thre)
  
  for(i in 1:5){
    cat(names(prlist)[i], file=paste0(plot_dir,"design_summary.txt"), 
        sep="\n", append = TRUE)
    mat = round(prlist[[i]], digits = 3)
    mat1 = data.frame(rownames(mat), mat)
    colnames(mat1) = c("p_thre", colnames(mat))
    write.table(mat1, file=paste0(plot_dir,"design_summary.txt"), 
                row.names=FALSE, col.names=TRUE, append = TRUE, quote = FALSE)
    cat("\n", file=paste0(plot_dir,"design_summary.txt"), append = TRUE)
  }
  
  return(prlist)
}
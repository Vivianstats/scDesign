### install scDesign -------------------------------------------------
# library(devtools)
# ## install MAST package
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("MAST")

# install_github("Vivianstats/scDesign")
library(scDesign)


set.seed(2020)
### simulate one scRNA-seq dataset -----------------------------------
realcount1 = readRDS(system.file("extdata", "astrocytes.rds", package = "scDesign"))
dim(realcount1)
realcount1[1:10, 1:2]

simcount1 = design_data(realcount1, S = sum(realcount1), ncell = 49, ncores = 1)
dim(simcount1)


### compare real and simulated data ----------------------------------
get_stats = function(mat, data){
  libsize = colSums(mat)
  mat = sweep(mat, 2, 1e6/libsize, FUN = "*")
  mat[, libsize == 0] = 0
  mat = log10(mat + 1)
  mean = rowMeans(mat)
  var = apply(mat,1,var)
  cv = sqrt(var)/mean

  zero_gene = rowSums(mat == 0)/ncol(mat)
  zero_cell = colSums(mat == 0)/nrow(mat)

  summs = list(libsize = libsize, mean = mean, var = var,
               cv = cv, drop_gene = zero_gene, drop_cell = zero_cell)
  summs = lapply(1:length(summs), function(i){
    data.frame(value = summs[[i]], measure = names(summs)[i], data = data,
               stringsAsFactors = FALSE)
  })
  summs = Reduce(rbind, summs)
  return(summs)
}
stats0 = get_stats(realcount1, data = "real")
stats1 = get_stats(simcount1, data = "simulated1")

dt = rbind(stats0, stats1)
ggplot(dt, aes(x = data, y = value, color = data)) +
  geom_boxplot(outlier.size = .5) + facet_wrap(~measure, scales = "free", ncol = 3) +
  theme_bw() + theme(legend.position = "none")


simcount2 = design_data(realcount1, S = sum(realcount1), ncell = 100, ncores = 1)
stats2 = get_stats(simcount2, data = "simulated2")
dt = rbind(dt, stats2)
ggplot(dt, aes(x = data, y = value, color = data)) +
  geom_boxplot(outlier.size = .5) + facet_wrap(~measure, scales = "free", ncol = 3) +
  theme_bw() + theme(legend.position = "none")


### simulate multiple scRNA-seq datasets -----------------------------------
simdata = design_data(realcount1, ngroup = 3, S = rep(1e7,3), ncell = rep(50,3), ncores = 1)

names(simdata)
names(simdata$count)
sapply(simdata$count, dim)

# up-regulated genes from state 1 to state 2
length(simdata$genesUp[[2]])
head(simdata$genesUp[[2]])


# visualize the result
count = cbind(simdata$count[[1]], simdata$count[[2]], simdata$count[[3]])
pca = prcomp(t( log10(count+1) ))
dt = data.frame(pc1 = pca$x[,1], pc2 = pca$x[,2], group = rep(1:3, each = 50))
ggplot(dt, aes(x = pc1, y = pc2, color = as.factor(group))) +
  geom_point() + theme_bw()

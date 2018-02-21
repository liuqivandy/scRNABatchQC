library(ggplot2)

checkClusterSeparateness <- function(object) {
  ####check the separateness of clusters using the silhouette width
  
  clust.col <- scater:::.get_palette("tableau10medium") # hidden scater colours
  sil <- silhouette(object$cluster, dist = object$dist)
  
  sil.cols <- clust.col[ifelse(sil[, 3] > 0, sil[, 1], sil[, 2])]
  
  sil.cols <- sil.cols[order(-sil[, 1], sil[, 3])]
  
  plot(sil, main = paste(length(unique(object$cluster)), "clusters "), border = sil.cols, col = sil.cols, do.col.sort = FALSE)
}

plot.multi.dens <- function(s, feature = "") {
  junk.x = NULL
  junk.y = NULL
  
  fNo <- which(colnames(colData(s[[1]]$sce)) == feature)
  
  stopifnot(fNo > 0)
  
  for (i in 1:length(s)) {
    junk.x = c(junk.x, density(colData(s[[i]]$sce)[, fNo])$x)
    junk.y = c(junk.y, density(colData(s[[i]]$sce)[, fNo])$y)
  }
  xr <- range(junk.x)
  yr <- range(junk.y)
  plot( density(colData(s[[1]]$sce)[, fNo]), xlim = xr, ylim = yr, main = "", xlab = feature)
  for (i in 1:length(s)) {
    lines(density(colData(s[[i]]$sce)[, fNo]), xlim = xr, ylim = yr, col = i)
  }
}

plot_multi_dens <- function(s, feature = "", cex = 2) {

  fNo <- which(colnames(colData(s[[1]]$sce)) == feature)
  
  stopifnot(fNo > 0)
  
  pdata <- data.frame()
  
  for (i in 1:length(s)) {
    pdata <- rbind(pdata, data.frame(data = colData(s[[i]]$sce)[, fNo], 
                                     Sample = rep(names(s)[i], length((colData(s[[i]]$sce)[, fNo])))))
  }
  
  p <- ggplot(pdata, aes(x = data, color = Sample)) + 
    geom_density(cex = cex) + xlab(feature) + scale_color_manual(values = 1:length(s))
  
  return(p)
}

plotMultiSamplesOneExplanatoryVariables <- function(s, var = "", size = 2) {
  d <- list()
  pct_var_explained <- c()
  
  for (i in 1:length(s)) {
    d[[i]] <- plotExplanatoryVariables(s[[i]]$sce, variables = c(var)) ##scater
    pct_var_explained <- c(pct_var_explained, d[[i]]$data$Pct_Var_Explained)
  }

  dat <- data.frame(Pct_Var_Explained = pct_var_explained,
                    Sample = rep(names(s), each = length(d[[1]]$data$Pct_Var_Explained)))
  
  p <- ggplot(dat, aes(x = Pct_Var_Explained, colour = Sample)) + 
    geom_line(stat = "density", alpha = 0.7, size = size, trim = T) + 
    geom_vline(xintercept = 1, linetype = 2) + 
    scale_x_log10(breaks = 10 ^ (-3:2), labels = c(0.001, 0.01, 0.1, 1, 10, 100)) + 
    xlab(paste0("% variance explained (log10-scale)")) + 
    ylab("Density") + ggtitle(var) + scale_color_manual(values = 1:length(s)) +
    coord_cartesian(xlim = c(10 ^ (-3), 100)) 
  return(p)
}

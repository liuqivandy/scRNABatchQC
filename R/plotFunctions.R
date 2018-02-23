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

### top 500 genes count distribution
plotGeneCountDistribution <- function(dat, scolors, nfeatures = 500) {
  prop_mat <- c()
  
  for (i in 1:length(dat)) {
    aveCountSum <- sum(rowData(dat[[i]]$sce)$ave.count)
    aveCountProp <- cumsum(sort(rowData(dat[[i]]$sce)$ave.count, decreasing = TRUE)) / aveCountSum
    prop_mat <- cbind(prop_mat, aveCountProp)
    colnames(prop_mat)[i] <- names(dat)[i]
  }
  rownames(prop_mat) <- seq_len(nrow(prop_mat))
  
  prop_to_plot <- reshape2::melt(prop_mat[seq_len(nfeatures), ])
  colnames(prop_to_plot) <- c("Feature", "Sample", "Proportion_Library")
  
  p <- ggplot(prop_to_plot, 
              aes_string(x = "Feature", y = "Proportion_Library", 
                         group = "Sample", colour = "Sample")) +
    geom_line() + 
    xlab("Number of features") + ylab("Cumulative proportion of library") +
    scale_color_manual(values = scolors) +
    theme_classic()
  
  print(p)
}

####averge count vs. detection rate

plotAveCountVSdetectRate <- function(dat, scolors) {
  avedetect <- data.frame()
  for (i in 1:length(dat)) {
    tmpavedec <- data.frame(avecount = log10(rowData(dat[[i]]$sce)$ave.count), 
                            detectrate = rowData(dat[[i]]$sce)$num.cells / dim(dat[[i]]$sce)[2], 
                            Sample = rep(names(dat)[i], length(log10(rowData(dat[[i]]$sce)$ave.count))))
    avedetect <- rbind(avedetect, tmpavedec)
  }
  
  p <- ggplot(avedetect, aes_string(x = "avecount", y = "detectrate", 
                                    group = "Sample", colour = "Sample")) + 
    geom_point() + xlab("Average count of genes") + ylab("Detecting rate") +
    scale_color_manual(values = scolors) +
    theme_classic()
  return(p)
}

##variance trend

plotVarianceTrend <- function(dat, scolors) {
  vartrend_dat <- data.frame()
  for (i in 1:length(dat)) {
    tmpvartrend <- data.frame(mean = dat[[i]]$hvg$mean, 
                              total = dat[[i]]$hvg$total, 
                              trend = dat[[i]]$var.fit$trend(dat[[i]]$hvg$mean),
                              Sample = rep(names(dat)[i], length(dat[[i]]$hvg$mean)))
    vartrend_dat <- rbind(vartrend_dat, tmpvartrend)
  }
  
  pp <- ggplot(vartrend_dat, aes_string(x = "mean", y = "total", group = "Sample", colour = "Sample")) + geom_point()
  pl <- ggplot(vartrend_dat, aes_string(x = "mean", y = "trend", group = "Sample", colour = "Sample")) + geom_line(alpha = 0.3, size = 1.5)
  p <- ggplot(vartrend_dat) + 
    geom_point(pp$mapping) + 
    geom_line(pl$mapping) + 
    scale_color_manual(values = scolors) + 
    xlab("Mean log-expression") + 
    ylab("Variance of log-expression") +
    theme_classic()
  
  return(p)
}


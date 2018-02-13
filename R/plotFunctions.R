

checkClusterSeparateness <- function(object) {
  ####check the separateness of clusters using the silhouette width
  
  clust.col <- scater:::.get_palette("tableau10medium") # hidden scater colours
  sil <- silhouette(object$cluster, dist = object$dist)
  
  sil.cols <- clust.col[ifelse(sil[, 3] > 0, sil[, 1], sil[, 2])]
  
  sil.cols <- sil.cols[order(-sil[, 1], sil[, 3])]
  
  plot(sil, main = paste(length(unique(object$cluster)), "clusters "), border = sil.cols, col = sil.cols, do.col.sort = FALSE)
}


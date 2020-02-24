#' An evolClustR Function
#'
#' Function to count number of neighbors for any given radius
#' @param Cdata PDB file parsed to contain only central Carbon atoms
#' @param radius Radius size of cluster of interest
#' @keywords cluster
#' @export
#' @examples
#' countneighbors()

count.neighbors<-function(Cdata,radius){
  all.pos <- Cdata[,c(7,8,9)]
  all.dists <- as.matrix(dist(all.pos))   
  all.insphere <- all.dists<radius

  # minus 1 to remove itself; just counting neighbors
  count <- apply(all.insphere, 2, sum) - 1

  return(as.vector(count))
}



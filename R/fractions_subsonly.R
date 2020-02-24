#' An evolClustR Function
#'
#' Function to return fractions in each cluster for any given radius
#' @param Cdata PDB file parsed to contain only central Carbon atoms
#' @param subs substited sites
#' @param radius Radius size of cluster of interest
#' @keywords cluster
#' @export
#' @examples
#' fractions()

# updated 5/20/19 to count central residue
# still in progress, vectorized distance calculation
# possibly done but check it

fractions_subsonly<-function(Cdata,subs,radius,fraction=TRUE){


  all.pos <- Cdata[,c(7,8,9)]
  all.dists <- as.matrix(dist(all.pos))
  all.insphere <- all.dists<radius
  sub.i <- which(is.element(Cdata$pos, subs))

  count <- apply(all.insphere, 2, sum)

  subbed <- all.insphere[sub.i,]
  clust.count <- apply(subbed, 2, sum)

  if(fraction){
    return((clust.count/count)[which(is.element(Cdata$pos, subs))])
  } else {
    return((clust.count)[which(is.element(Cdata$pos, subs))])
  }
}


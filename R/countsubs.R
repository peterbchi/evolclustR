#' An evolClustR Function
#'
#' Function to return the count of substitutions in each sphere for any given radius
#' @param Cdata PDB file parsed to contain only central Carbon atoms
#' @param subs all substited sites
#' @param radius Radius size of cluster of interest
#' @keywords cluster
#' @export
#' @examples
#' count.subs()

# Yu/Thorne excludes the centers
# pg. 683: Imagine that the numbers of amino acid replacements
# within the four balls (excluding the centers) are 2, 2, 2, and 0. In
# this case, the average number of replacements within a ball would
# be 1.5.

count.subs<-function(Cdata,subs,radius){

  # this is going to re-order them, not sure if I care or not (I don't think I do since N_i is just the average)
  subs.pos <- Cdata[which(is.element(Cdata$pos, subs)),c(7,8,9)]
  subs.dists <- as.matrix(dist(subs.pos))
  subs.insphere <- subs.dists<radius

  subs.count <- apply(subs.insphere, 2, sum) - 1 # to exclude centers
  return(as.vector(subs.count))  
     
}


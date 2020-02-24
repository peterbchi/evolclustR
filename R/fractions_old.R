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

fractions<-function(Cdata,subs,radius){
   count<-rep(0,dim(Cdata)[1])
   clust.count<-rep(0,dim(Cdata)[1])

   for(i in Cdata$pos){
      center<-Cdata[which(Cdata$pos==i),c(7,8,9)]		# i is pos, not row number
      dists <- apply(Calphas[,c(7,8,9)], 1, function(x,y) return(sqrt(sum((x-y)^2))), center)

      count[which(Cdata$pos==i)] <- sum(dists<radius)
      clust.count[which(Cdata$pos==i)] <- sum(is.element(subs, Cdata$pos[dists<radius]))
   }
   return(clust.count/count)
}


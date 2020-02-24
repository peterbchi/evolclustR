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

fractions<-function(Cdata,subs,radius){
   count<-rep(0,length(subs))
   clust.count<-rep(0,length(subs))

   for(i in subs){
      center<-Cdata[which(Cdata[,6]==i),c(7,8,9)]
   #   dists<-rep(NA,dim(Cdata)[1]-1)
      for(j in Cdata[,6][-which(Cdata[,6]==i)]){
         # print(j) sanity check
         other<-Cdata[which(Cdata[,6]==j),c(7,8,9)]
         dist.temp<-dist(rbind(center,other))[1]
         if(dist.temp<radius){
            count[which(subs==i)]<-count[which(subs==i)]+1
         }
         if(dist.temp<radius & is.element(j, subs)){
            clust.count[which(subs==i)]<-clust.count[which(subs==i)]+1
         }
      }
   }
   return(clust.count/count)
}


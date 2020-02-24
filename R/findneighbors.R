#' An evolClustR Function
#'
#' Function to return fractions in each cluster for any given radius
#' @param Cdata PDB file parsed to contain only central Carbon atoms
#' @param radius Radius of interest
#' @param site location of central site
#' @keywords cluster
#' @export
#' @examples
#' find.neighbors()


find.neighbors<-function(Cdata,radius,site){
  neighbors<-NA
  row.ind<-which(Cdata[,6]==site)
  center<-Cdata[row.ind,c(7,8,9)]
  for(j in Cdata[,6][-row.ind]){	# 
    other<-Cdata[which(Cdata[,6]==j),c(7,8,9)]
    dist.temp<-dist(rbind(center,other))[1]
    if(dist.temp<radius){
      neighbors<-c(neighbors,j)
    }         
  }
  return(neighbors[-1])
}


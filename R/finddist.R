#' An evolClustR Function
#'
#' Function to find the distance from the center of any given residue
#' @param Cdata PDB file parsed to contain only central Carbon atoms
#' @keywords cluster
#' @export
#' @examples
#' finddist()

# Note: indexing was changed on 4/6/18; make sure this doesn't fuck anything else up
# Index i below doesn't make any sense. Not sure why it doesn't throw an error

find.dist<-function(Cdata, sub){		# sub represents row number of subbed site
   all.dist<-rep(NA,dim(Cdata)[1])
   n.rows<-dim(Cdata)[1]
   count<-rep(0,n.rows)
   center<-Cdata[sub,c(7,8,9)]
   for(j in 1:n.rows){			# index represents row number (NOT residue number)
      other<-Cdata[j,c(7,8,9)]			# 7/26/18: took out [-sub] -- i think this was wrong.
      all.dist[j]<-dist(rbind(center,other))[1]
   }
   return(all.dist)
}



